import abc
import collections
import itertools
import json
import re
import sqlite3
from collections.abc import Mapping
from typing import Final

import numpy as np

from . import compression
from . import gemini_utils as util
from .gemini_constants import (
    AFFECTED,
    HET,
    HOM_ALT,
    HOM_REF,
    MISSING,
    UNAFFECTED,
    UNKNOWN,
)
from .gemini_region import add_region_to_query
from .gemini_subjects import Subject, get_subjects
from .pdict import PDict
from .sql_utils import ensure_columns, get_select_cols_and_rest

flatten = itertools.chain.from_iterable


class GeminiError(Exception):
    pass


class RowFormat(abc.ABC):
    """Abstract base class for custom query output formatters."""

    @property
    @abc.abstractmethod
    def name(self):
        pass

    @abc.abstractmethod
    def format(self, row):
        return str(row.print_fields)

    def format_query(self, query):
        return query

    def predicate(self, row):
        return True

    def header(self, fields):
        return "\t".join(fields)


class DefaultRowFormat(RowFormat):
    name = "default"

    def __init__(self, args=None):
        self.args = args

    def format(self, row):
        return str(row.print_fields)


class JSONRowFormat(RowFormat):
    name = "json"

    def __init__(self, args=None):
        self.args = args

    def format(self, row):
        from .pdict import to_json

        return json.dumps(row.print_fields, default=to_json)

    def header(self, fields):
        return None


class CarrierSummary(RowFormat):
    """Generates counts of carrier/noncarrier status for a given column of the samples table."""

    name = "carrier_summary"

    def __init__(self, gq_or_args=None, subjects=None, carrier_summary=None):
        self.carrier_summary = carrier_summary or getattr(gq_or_args, "carrier_summary", "affected")
        if subjects is not None:
            self.subjects = subjects
        elif hasattr(gq_or_args, "sample_to_sample_object") and gq_or_args.sample_to_sample_object:
            self.subjects = gq_or_args.sample_to_sample_object
        elif gq_or_args is not None:
            self.subjects = get_subjects(gq_or_args)
        else:
            self.subjects = {}

        self.column_types = sorted(
            {
                getattr(x, self.carrier_summary)
                for x in self.subjects.values()
                if getattr(x, self.carrier_summary, None) is not None
            }
        )
        self.column_counters = {None: set()}
        for ct in self.column_types:
            self.column_counters[ct] = {
                k for k, v in self.subjects.items() if getattr(v, self.carrier_summary, None) == ct
            }

    def format(self, row):
        have_variant = set(row["variant_samples"])
        have_reference = set(row["hom_ref_samples"])
        unknown = len(set(row["unknown_samples"]).union(self.column_counters[None]))
        carrier_counts = []
        for ct in self.column_types:
            carrier_counts.append(len(self.column_counters[ct].intersection(have_variant)))
        for ct in self.column_types:
            carrier_counts.append(len(self.column_counters[ct].intersection(have_reference)))
        carrier_counts.append(unknown)

        r = row.print_fields
        prefix = "\t".join(
            ",".join(str(s) for s in r[c]) if isinstance(r[c], np.ndarray) else str(r[c]) for c in r
        )
        suffix = "\t".join(str(c) for c in carrier_counts)
        return f"{prefix}\t{suffix}" if prefix else suffix

    def header(self, fields):
        header_columns = self.column_types
        if self.carrier_summary == "affected":
            header_columns = [
                "affected" if ct is True else "unaffected" if ct is False else str(ct)
                for ct in self.column_types
            ]
        carriers = [f"{x}_carrier" for x in header_columns]
        noncarriers = [f"{x}_noncarrier" for x in header_columns]
        res = list(fields) + carriers + noncarriers + ["unknown"]
        return "\t".join(res)


class TPEDRowFormat(RowFormat):
    name = "tped"
    NULL_GENOTYPES: Final[list[str]] = ["."]
    PED_MISSING: Final[list[str]] = ["0", "0"]
    VALID_CHROMOSOMES: Final[list[str]] = [str(x) for x in range(1, 23)] + ["X", "Y", "XY", "MT"]
    POSSIBLE_HAPLOID: Final[list[str]] = ["X", "Y"]
    X_PAR_REGIONS: Final[list[tuple[int, int]]] = [(60001, 2699520), (154931044, 155260560)]
    Y_PAR_REGIONS: Final[list[tuple[int, int]]] = [(10001, 2649520), (59034050, 59363566)]

    def __init__(self, gq_or_args):
        if isinstance(gq_or_args, GeminiQuery):
            gq = gq_or_args
        else:
            gq = GeminiQuery(gq_or_args.db if hasattr(gq_or_args, "db") else gq_or_args)
        self.gq = gq
        self.samples = [gq.idx_to_sample_object[x] for x in range(len(gq.sample_to_idx))]

    def format_query(self, query):
        needed = ["chrom", "start", "ref", "alt", "gts", "type", "variant_id"]
        if hasattr(self.gq, "variant_column_names") and "rs_ids" in self.gq.variant_column_names:
            needed.insert(1, "rs_ids")
        return ensure_columns(query, needed)

    def predicate(self, row, _splitter=re.compile(r"[|/]")):
        geno = list(flatten(_splitter.split(util.to_str(x)) for x in row["gts"]))
        num_alleles = len(set(geno).difference(self.NULL_GENOTYPES))
        return 0 < num_alleles <= 2 and row.row["type"] != "sv"

    def format(self, row):
        chrom = row["chrom"].removeprefix("chr")
        chrom = chrom if chrom in self.VALID_CHROMOSOMES else "0"
        start = str(row.row["start"])
        geno = [re.split(r"[|/]", util.to_str(x)) for x in row["gts"]]
        geno = [
            self._fix_genotype(chrom, int(start), g, getattr(self.samples[i], "sex", None))
            for i, g in enumerate(geno)
        ]
        genotypes = " ".join(flatten(geno))
        name = str(row.get("rs_ids") or row["variant_id"])
        return f"{chrom} {name} 0 {start} {genotypes}"

    def _fix_genotype(self, chrom, start, genotype, sex):
        if sex == "2":
            if len(genotype) < 2 or chrom == "Y" or any(a in self.NULL_GENOTYPES for a in genotype):
                return self.PED_MISSING
            return genotype
        if chrom in self.POSSIBLE_HAPLOID and sex == "1":
            genotype = [x for x in genotype if x not in self.NULL_GENOTYPES]
            if not genotype:
                return self.PED_MISSING
            if len(genotype) == 2 and genotype[0] != genotype[1] and not self._in_par(chrom, start):
                return self.PED_MISSING
            if len(genotype) < 2:
                return [genotype[0], genotype[0]]
        if any(a in self.NULL_GENOTYPES for a in genotype) or len(genotype) < 2:
            return self.PED_MISSING
        return genotype

    def _in_par(self, chrom, start):
        regions = self.X_PAR_REGIONS if chrom == "X" else self.Y_PAR_REGIONS if chrom == "Y" else []
        return any(r[0] < start < r[1] for r in regions)

    def header(self, fields):
        return None


class VCFRowFormat(RowFormat):
    name = "vcf"

    def __init__(self, gq_or_args):
        if isinstance(gq_or_args, GeminiQuery):
            self.gq = gq_or_args
        else:
            self.gq = GeminiQuery(gq_or_args.db if hasattr(gq_or_args, "db") else gq_or_args)

    def format_query(self, query):
        needed = [
            "chrom",
            "start",
            "vcf_id",
            "ref",
            "alt",
            "qual",
            "filter",
            "info",
            "gts",
            "gt_types",
            "gt_phases",
        ]
        if hasattr(self.gq, "variant_column_names") and self.gq.variant_column_names:
            needed = [c for c in needed if c in self.gq.variant_column_names]
        return ensure_columns(query, needed)

    def format(self, row):
        istr = self.gq._info_dict_to_string
        vcf_rec = [row.row["chrom"], row.row["start"] + 1]
        vcf_rec.append("." if row.row["vcf_id"] is None else row.row["vcf_id"])
        vcf_rec.extend([row.row["ref"], row.row["alt"], row.row["qual"]])
        vcf_rec.append("PASS" if row.row["filter"] is None else row.row["filter"])
        info_val = istr(row["info"]) if "info" in row and row["info"] is not None else "."
        vcf_rec.extend([info_val, "GT"])

        gts = [util.to_str(r) for r in row["gts"]]
        gt_types = list(row["gt_types"])
        gt_phases = list(row["gt_phases"])
        for idx, gt_type in enumerate(gt_types):
            phase_char = "|" if gt_phases[idx] else "/"
            gt = gts[idx]
            alleles = gt.split(phase_char)
            if gt_type == HOM_REF:
                vcf_rec.append(f"0{phase_char}0")
            elif gt_type == HET:
                if gt_phases[idx] and len(alleles) > 0 and alleles[0] != row.row["ref"]:
                    vcf_rec.append(f"1{phase_char}0")
                else:
                    vcf_rec.append(f"0{phase_char}1")
            elif gt_type == HOM_ALT:
                vcf_rec.append(f"1{phase_char}1")
            else:
                vcf_rec.append(f".{phase_char}.")

        return "\t".join(str(c) if c is not None else "." for c in vcf_rec)

    def header(self, fields):
        cursor = self.gq.conn.cursor()
        try:
            cursor.execute("SELECT vcf_header FROM vcf_header")
            row = cursor.fetchone()
            if row:
                return str(row[0]).strip()
        except sqlite3.OperationalError:
            pass
        return None


class SampleDetailRowFormat(RowFormat):
    name = "sampledetail"

    def __init__(self, gq_or_args, sample_delim=","):
        if isinstance(gq_or_args, GeminiQuery):
            self.gq = gq_or_args
        else:
            self.gq = GeminiQuery(gq_or_args.db if hasattr(gq_or_args, "db") else gq_or_args)
        self.sample_delim = sample_delim
        cursor = self.gq.conn.cursor()
        cursor.execute("SELECT * FROM samples")
        self.cols = [col[0] for col in cursor.description if col[0] != "name"]
        self.samples = {}
        for row in cursor.fetchall():
            self.samples[row["name"]] = [str(row[x]) for x in self.cols]

    def format(self, row):
        r = row.print_fields
        variant_samples = str(r.get("variant_samples", "")).split(self.sample_delim)
        samples = [s for s in variant_samples if s]
        out = []
        for sample in samples:
            out.append("\t".join([str(r[c]) for c in r] + self.samples.get(sample, [])))
        return "\n".join(out)

    def header(self, fields):
        return "\t".join(fields + self.cols)


class GeminiRow:
    """Represents a single query result row from a GEMINI database."""

    __slots__ = (
        "cache",
        "formatter",
        "genotype_dict",
        "print_fields",
        "query",
        "row",
        "unpacker",
    )

    def __init__(
        self,
        row,
        query,
        formatter=None,
        print_fields=None,
        unpacker=None,
    ):
        self.cache = getattr(row, "cache", {}).copy()
        self.genotype_dict = getattr(row, "genotype_dict", {}).copy()
        self.row = getattr(row, "row", row)
        self.unpacker = unpacker or compression.unpack_genotype_blob
        self.formatter = formatter or DefaultRowFormat()
        self.query = query
        self.print_fields = print_fields or PDict()

    @property
    def info(self):
        return self["info"]

    def __getitem__(self, key):
        if isinstance(key, (int, slice)):
            try:
                return self.row[key]
            except IndexError:
                raise IndexError(key) from None

        key = str(key)
        if key in (
            "het_samples",
            "hom_alt_samples",
            "unknown_samples",
            "variant_samples",
            "hom_ref_samples",
        ):
            if not self.genotype_dict:
                self.genotype_dict = self.query._group_samples_by_genotype(self["gt_types"])
            if key == "het_samples":
                return self.genotype_dict[HET]
            if key == "hom_alt_samples":
                return self.genotype_dict[HOM_ALT]
            if key == "hom_ref_samples":
                return self.genotype_dict[HOM_REF]
            if key == "unknown_samples":
                return self.genotype_dict[UNKNOWN]
            if key == "variant_samples":
                return self.genotype_dict[HET] + self.genotype_dict[HOM_ALT]

        if key in self.cache:
            return self.cache[key]

        if key == "info":
            if "info" not in self.cache:
                self.cache["info"] = compression.unpack_ordereddict_blob(self.row["info"])
            return self.cache["info"]

        # Named genotype columns like gts.NA12878
        if "." in key:
            col, sample = key.split(".", 1)
            if col in self.query.gt_cols and sample in self.query.sample_to_idx:
                idx = self.query.sample_to_idx[sample]
                arr = self[col]
                return arr[idx] if arr is not None else None

        if key in self.query.gt_cols:
            if key not in self.cache:
                blob = self.row[key]
                arr = self.unpacker(blob)
                if arr is not None and key == "gts":
                    arr = arr.astype(str)
                self.cache[key] = arr
            return self.cache[key]

        try:
            return self.row[key]
        except IndexError, KeyError:
            raise KeyError(key) from None

    def __iter__(self):
        return iter(self.keys())

    def __len__(self):
        return len(self.keys())

    def __contains__(self, key):
        return key in self.keys()

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def keys(self):
        row_keys = list(self.row.keys()) if hasattr(self.row, "keys") else []
        return row_keys + list(self.cache.keys())

    def values(self):
        return [self[k] for k in self.keys()]

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def __repr__(self):
        return self.formatter.format(self)

    __str__ = __repr__


Mapping.register(GeminiRow)


class GeminiQuery:
    """Submit queries to an existing GEMINI database and iterate over the results."""

    def __init__(
        self,
        db,
        include_gt_cols=False,
        out_format=None,
    ):
        self.db = db
        self.include_gt_cols = include_gt_cols
        self._connect_to_database()

        cursor = self.conn.cursor()
        self.gt_cols = util.get_gt_cols(cursor)
        self._collect_sample_table_columns(cursor)
        self._collect_variant_table_columns(cursor)

        self.sample_info = collections.defaultdict(list)
        self.sample_to_idx = {}
        self.idx_to_sample = {}
        self.sample_to_sample_object = {}
        self.idx_to_sample_object = {}

        try:
            cursor.execute("SELECT * FROM samples")
            samples = cursor.fetchall()
            for s in samples:
                idx = int(s["sample_id"]) - 1
                name = str(s["name"])
                self.sample_to_idx[name] = idx
                self.idx_to_sample[idx] = name
                subj = Subject(s)
                self.sample_to_sample_object[name] = subj
                self.idx_to_sample_object[idx] = subj
        except sqlite3.OperationalError:
            pass

        self.formatter = self._resolve_formatter(out_format)

        self.unpacker = compression.unpack_genotype_blob
        try:
            cursor.execute("SELECT feature FROM features")
            features = [r["feature"] for r in cursor.fetchall()]
            if "snappy_compression" in features:
                self.unpacker = compression.snappy_unpack_blob
        except sqlite3.OperationalError:
            pass

        self.predicates = [self.formatter.predicate]
        self.sample_show_fields = ["variant_samples", "het_samples", "hom_alt_samples"]
        self.result_proxy = None

    def _resolve_formatter(self, out_format):
        if out_format is None:
            return DefaultRowFormat()
        if isinstance(out_format, RowFormat):
            return out_format
        if isinstance(out_format, str):
            fmt = out_format.lower()
            if fmt == "default":
                return DefaultRowFormat()
            if fmt == "json":
                return JSONRowFormat()
            if fmt == "vcf":
                return VCFRowFormat(self)
            if fmt == "tped":
                return TPEDRowFormat(self)
            if fmt == "carrier_summary":
                return CarrierSummary(self)
            if fmt == "sampledetail":
                return SampleDetailRowFormat(self)
        raise ValueError(f"Unknown row format: {out_format}")

    def _connect_to_database(self):
        if isinstance(self.db, sqlite3.Connection):
            self.conn = self.db
        else:
            self.conn = sqlite3.connect(self.db)
        self.conn.row_factory = sqlite3.Row

    def _collect_sample_table_columns(self, cursor):
        try:
            cursor.execute("PRAGMA table_info('samples')")
            self.sample_column_names = [r["name"] for r in cursor.fetchall()]
        except sqlite3.OperationalError:
            self.sample_column_names = []

    def _collect_variant_table_columns(self, cursor):
        try:
            cursor.execute("PRAGMA table_info('variants')")
            self.variant_column_names = [r["name"] for r in cursor.fetchall()]
        except sqlite3.OperationalError:
            self.variant_column_names = []

    def run(
        self,
        query,
        gt_filter=None,
        show_variant_samples=False,
        variant_samples_delim=",",
        predicates=None,
        needs_genotypes=False,
        needs_genes=False,
        show_families=False,
        subjects=None,
        region=None,
    ):
        """Execute a query against the GEMINI database."""
        if region:
            query = add_region_to_query(query, region)

        self.query = self.formatter.format_query(query)
        self.gt_filter = gt_filter
        if not self._is_gt_filter_safe():
            raise ValueError("ERROR: invalid --gt-filter expression.")

        self.show_variant_samples = show_variant_samples
        self.variant_samples_delim = variant_samples_delim
        self.needs_genotypes = needs_genotypes
        self.needs_vcf_columns = isinstance(self.formatter, VCFRowFormat)
        self.needs_genes = needs_genes
        self.show_families = show_families
        self.subjects = subjects
        if predicates:
            self.predicates = [self.formatter.predicate] + list(predicates)

        self.query = self.query.replace(",", ", ")
        self.query_pieces = self.query.split()

        has_gt_in_query = any(s.startswith(("gt", "(gt")) or ".gt" in s for s in self.query_pieces)

        if not has_gt_in_query:
            if self.gt_filter is None:
                self.query_type = "no-genotypes"
            else:
                self.gt_filter = self._correct_genotype_filter()
                self.query_type = "filter-genotypes"
        else:
            if self.gt_filter is not None:
                self.gt_filter = self._correct_genotype_filter()
            self.query_type = "select-genotypes"

        if self.gt_filter:
            self.gt_filter_compiled = compile(self.gt_filter, "<gt_filter>", "eval")

        self.result_proxy = iter(self._apply_query())
        return self

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            row_raw = next(self.result_proxy)
            row = GeminiRow(row_raw, self, unpacker=self.unpacker)

            if self.gt_filter:
                try:
                    eval_context = {
                        "sample_info": self.sample_info,
                        "HET": HET,
                        "HOM_REF": HOM_REF,
                        "HOM_ALT": HOM_ALT,
                        "UNKNOWN": UNKNOWN,
                        "MISSING": MISSING,
                        "AFFECTED": AFFECTED,
                        "UNAFFECTED": UNAFFECTED,
                        "all": all,
                        "any": any,
                        "sum": sum,
                    }
                    for col in self.gt_cols:
                        if col in self.gt_filter:
                            eval_context[col] = row[col]

                    if not eval(self.gt_filter_compiled, eval_context):
                        continue
                except TypeError:
                    continue

            fields = PDict()
            for col in self.report_cols:
                if col == "*":
                    continue
                if not col.lower().startswith("gt"):
                    fields.add(col, row[col])
                else:
                    if "[" in col:
                        orig_col = self.gt_idx_to_name_map[col]
                        source, extra = col.split("[", 1)
                        idx = int(extra[:-1])
                        val = row[source][idx] if row[source] is not None else None
                        if isinstance(val, (np.integer, np.bool_)):
                            fields.add(orig_col, int(val))
                        elif isinstance(val, np.floating):
                            fields.add(orig_col, float(val))
                        else:
                            fields.add(orig_col, val)
                    else:
                        fields[col] = row[col] if row[col] is not None else str(None)

            if self.show_variant_samples:
                fields["variant_samples"] = self.variant_samples_delim.join(
                    self._filter_samples(row["variant_samples"])
                )
                fields["het_samples"] = self.variant_samples_delim.join(
                    self._filter_samples(row["het_samples"])
                )
                fields["hom_alt_samples"] = self.variant_samples_delim.join(
                    self._filter_samples(row["hom_alt_samples"])
                )

            if self.show_families:
                families = sorted(
                    {
                        str(self.sample_to_sample_object[x].family_id)
                        for x in row["variant_samples"]
                        if x in self.sample_to_sample_object
                        and getattr(self.sample_to_sample_object[x], "family_id", None)
                    }
                )
                fields["families"] = self.variant_samples_delim.join(families)

            if not all(predicate(row) for predicate in self.predicates):
                continue

            return GeminiRow(
                row,
                self,
                formatter=self.formatter,
                print_fields=fields,
                unpacker=self.unpacker,
            )

    @property
    def header(self):
        if self.query_type == "no-genotypes":
            h = list(self.all_query_cols)
        else:
            h = list(self.all_query_cols) + [
                col
                for col in util.OrderedSet(self.all_columns_orig)
                - util.OrderedSet(self.select_columns)
            ]
        if self.show_variant_samples:
            h += self.sample_show_fields
        if self.show_families:
            h += ["families"]
        return self.formatter.header(h)

    @property
    def sample2index(self):
        return self.sample_to_idx

    @property
    def index2sample(self):
        return self.idx_to_sample

    def _filter_samples(self, samples):
        if self.subjects is not None:
            return [x for x in samples if x in self.subjects]
        return samples

    def _group_samples_by_genotype(self, gt_types):
        d = {HOM_REF: [], HET: [], UNKNOWN: [], HOM_ALT: []}
        lookup = self.idx_to_sample
        for i, val in enumerate(gt_types):
            sample_name = lookup.get(i)
            if sample_name is not None and val in d:
                d[val].append(sample_name)
        return d

    def _is_gt_filter_safe(self, gt_filter=None):
        gt_filter = gt_filter or self.gt_filter
        if not gt_filter or not gt_filter.strip():
            return True
        if "__" in gt_filter or "import" in gt_filter:
            return False
        evil = ["rm ", "os.", "sys.", "exec", "eval", "shutil"]
        if any(s in gt_filter for s in evil):
            return False
        valid_cols = list(flatten((f"{gtc}.", f"({gtc}).") for gtc in self.gt_cols))
        return any(s in gt_filter for s in valid_cols)

    def _execute_query(self):
        cursor = self.conn.cursor()
        try:
            cursor.execute(self.query)
        except sqlite3.OperationalError as e:
            raise ValueError(f"SQL error in query '{self.query}': {e}") from e
        return cursor

    def _apply_query(self):
        if self.needs_genes:
            self.query = self._add_gene_col_to_query()
        if self.needs_vcf_columns:
            self.query = self._add_vcf_cols_to_query()

        if self._query_needs_genotype_info():
            self._split_select()
            self.query = self._add_gt_cols_to_query()
            cursor = self._execute_query()
            self.all_query_cols = [
                str(t[0])
                for t in cursor.description
                if not t[0].lower().startswith("gt") and ".gt" not in t[0].lower()
            ]

            if "*" in self.select_columns:
                self.select_columns.remove("*")
                self.all_columns_orig.remove("*")
                self.all_columns_new.remove("*")
                self.select_columns += self.all_query_cols

            self.report_cols = self.all_query_cols + list(
                util.OrderedSet(self.all_columns_new) - util.OrderedSet(self.select_columns)
            )
        else:
            cursor = self._execute_query()
            self.all_query_cols = [
                str(t[0]) for t in cursor.description if not t[0].lower().startswith("gt")
            ]
            self.report_cols = self.all_query_cols

        return cursor

    def _correct_genotype_col(self, raw_col):
        if raw_col == "*":
            return raw_col
        if "." in raw_col:
            column, sample = raw_col.split(".", 1)
            sample_idx = self.sample_to_idx.get(sample)
            if sample_idx is None:
                raise GeminiError(f"Unknown sample name '{sample}' in genotype column '{raw_col}'")
            return f"{column.lower()}[{sample_idx}]"
        return raw_col

    def _get_matching_sample_ids(self, wildcard):
        cursor = self.conn.cursor()
        if wildcard.strip() == "*":
            q = "SELECT sample_id, name FROM samples"
        else:
            q = f"SELECT sample_id, name FROM samples WHERE {wildcard}"
        cursor.execute(q)
        return [(int(row["sample_id"]) - 1, str(row["name"])) for row in cursor.fetchall()]

    def _correct_genotype_filter(self, gt_filter=None):
        def _swap_genotype_for_number(token):
            token = token.replace("HOM_REF", str(HOM_REF))
            token = token.replace("HET", str(HET))
            token = token.replace("HOM_ALT", str(HOM_ALT))
            token = token.replace("UNKNOWN", str(UNKNOWN))
            return token

        corrected_gt_filter = []
        raw_filter = str(gt_filter or self.gt_filter)
        wildcard_tokens = re.split(r"(\(\s*gt\w+\s*\)\.\(.+?\)\.\(.+?\)\.\(.+?\))", raw_filter)

        for token_idx, token in enumerate(wildcard_tokens):
            if "gt" in token.lower() and ".(" not in token and ")self." not in token:
                tokens = re.split(r"[\s+]+", token)
                for t in tokens:
                    if not t:
                        continue
                    if "gt" in t.lower():
                        corrected_gt_filter.append(self._correct_genotype_col(t))
                    else:
                        t = _swap_genotype_for_number(t)
                        if t.strip() in ("AND", "OR"):
                            t = t.lower()
                        corrected_gt_filter.append(t)
            elif "gt" in token.lower() and ".(" in token and ")." in token:
                if token.count(".") != 3 or token.count("(") != 4 or token.count(")") != 4:
                    raise ValueError(f"Wildcard filter should consist of 4 elements: {token}")
                column, wildcard, wildcard_rule, wildcard_op = token.split(".")
                column = column.strip("() ")
                wildcard = wildcard.strip("() ")
                is_gt_expression = wildcard in ("=HET", "=HOM_REF", "=HOM_ALT")
                wildcard_rule = wildcard_rule.strip("() ")
                wildcard_op = wildcard_op.strip("() ")
                extra_filter = None

                if is_gt_expression:
                    self.sample_info[token_idx] = self._get_matching_sample_ids("*")
                    extra_filter = f" if gt_types[sample[0]] == {wildcard[1:]}"
                else:
                    self.sample_info[token_idx] = self._get_matching_sample_ids(wildcard)

                wildcard_rule = _swap_genotype_for_number(wildcard_rule)

                if wildcard_op in ("all", "any"):
                    rule = f"{wildcard_op}({column}[sample[0]]{wildcard_rule} for sample in sample_info[{token_idx}]{extra_filter or ''})"
                elif wildcard_op == "none":
                    rule = f"not any({column}[sample[0]]{wildcard_rule} for sample in sample_info[{token_idx}]{extra_filter or ''})"
                elif "count" in wildcard_op:
                    count_comp = wildcard_op.split("count")[-1]
                    rule = f"sum({column}[sample[0]]{wildcard_rule} for sample in sample_info[{token_idx}]{extra_filter or ''}){count_comp}"
                else:
                    raise ValueError(f"Unsupported wildcard operation: ({wildcard_op})")

                corrected_gt_filter.append(rule)
            else:
                if token:
                    if token.strip() in ("AND", "OR"):
                        token = token.lower()
                    corrected_gt_filter.append(token)

        return " ".join(corrected_gt_filter)

    def _add_gt_cols_to_query(self):
        select_tokens, rest_of_query = get_select_cols_and_rest(self.query)
        select_clause_list = [
            token
            for token in select_tokens
            if not token.lower().startswith("gt")
            and not token.lower().startswith("(gt")
            and ".gt" not in token.lower()
        ]
        select_clause = ", ".join(select_clause_list).strip()
        if select_clause:
            select_clause += ", "
        select_clause += ", ".join(self.gt_cols)
        return f"SELECT {select_clause} {rest_of_query}"

    def _add_gene_col_to_query(self):
        select_tokens, rest_of_query = get_select_cols_and_rest(self.query)
        if not any("gene" in s for s in select_tokens):
            select_tokens.append("gene")
            return f"SELECT {', '.join(select_tokens)} {rest_of_query}"
        return self.query

    def _add_vcf_cols_to_query(self):
        select_tokens, rest_of_query = get_select_cols_and_rest(self.query)
        vcf_cols = [
            "chrom",
            "start",
            "vcf_id",
            "ref",
            "alt",
            "qual",
            "filter",
            "info",
            "gts",
            "gt_types",
            "gt_phases",
        ]
        if self.variant_column_names:
            vcf_cols = [c for c in vcf_cols if c in self.variant_column_names]
        cols_to_add = [col for col in vcf_cols if not any(col in s for s in select_tokens)]
        select_clause = ", ".join(select_tokens + cols_to_add)
        return f"SELECT {select_clause} {rest_of_query}"

    def _split_select(self):
        self.select_columns = []
        self.all_columns_new = []
        self.all_columns_orig = []
        self.gt_name_to_idx_map = {}
        self.gt_idx_to_name_map = {}

        select_tokens, _ = get_select_cols_and_rest(self.query)
        for token in select_tokens:
            if "gt" in token.lower() and ".(" in token and ")." in token:
                column, wildcard = token.split(".")
                wildcard = wildcard.strip("()")
                column = column.strip("()")
                sample_info = self._get_matching_sample_ids(wildcard)
                for sample in sample_info:
                    wildcard_display_col = f"{column}.{sample[1]}"
                    wildcard_mask_col = f"{column}[{sample[0]}]"
                    self.all_columns_new.append(wildcard_mask_col)
                    self.all_columns_orig.append(wildcard_display_col)
                    self.gt_name_to_idx_map[wildcard_display_col] = wildcard_mask_col
                    self.gt_idx_to_name_map[wildcard_mask_col] = wildcard_display_col
            elif "gt" in token.lower() and ".(" not in token and "length" not in token:
                new_col = self._correct_genotype_col(token)
                self.all_columns_new.append(new_col)
                self.all_columns_orig.append(token)
                self.gt_name_to_idx_map[token] = new_col
                self.gt_idx_to_name_map[new_col] = token
            else:
                self.select_columns.append(token)
                self.all_columns_new.append(token)
                self.all_columns_orig.append(token)

    def _info_dict_to_string(self, info):
        if not info:
            return "."
        parts = []
        for key, value in info.items():
            if isinstance(value, list):
                parts.append(f"{key}={','.join(str(v) for v in value)}")
            else:
                parts.append(f"{key}={value}")
        return ";".join(parts)

    def _query_needs_genotype_info(self):
        if (
            self.include_gt_cols
            or self.show_variant_samples
            or self.needs_genotypes
            or bool(self.gt_filter)
            or isinstance(
                self.formatter,
                (CarrierSummary, TPEDRowFormat, VCFRowFormat, SampleDetailRowFormat),
            )
        ):
            return True
        tokens = [t.strip(",() ") for t in self.query.split()]
        return any(t.lower().startswith("gt") or ".gt" in t.lower() for t in tokens)


def select_formatter(args):
    supported = {x.name.lower(): x for x in util.itersubclasses(RowFormat)}
    fmt_name = getattr(args, "format", "default").lower()
    if fmt_name not in supported:
        raise NotImplementedError(f"Unsupported format: {fmt_name}")
    return supported[fmt_name](args)


def add_variant_ids_to_query(query, vids):
    if not vids:
        return query
    vid_str = ",".join(str(int(x)) for x in vids)
    extra = f"variant_id IN ({vid_str})"
    where_idx = query.lower().find("where")
    if where_idx == -1:
        tail_match = re.search(r"\b(order\s+by|limit)\b", query, re.IGNORECASE)
        if tail_match:
            idx = tail_match.start()
            return f"{query[:idx]} WHERE {extra} {query[idx:]}"
        return f"{query} WHERE {extra}"
    return f"{query[:where_idx]} WHERE {extra} AND {query[where_idx + len('where') :]}"
