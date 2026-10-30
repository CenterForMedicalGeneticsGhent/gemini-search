import sqlite3

import numpy as np
import pytest

from gemini import (
    HET,
    HOM_ALT,
    HOM_REF,
    JSONRowFormat,
)
from gemini.compression import pack_blob, zdumps
from gemini.gemini_query import GeminiQuery


@pytest.fixture
def test_db(tmp_path):
    """Create a minimal SQLite GEMINI database with samples and variants."""
    db_path = str(tmp_path / "test.db")
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # Create samples table
    cur.execute("""
        CREATE TABLE samples (
            sample_id INTEGER PRIMARY KEY,
            family_id TEXT,
            name TEXT,
            paternal_id TEXT,
            maternal_id TEXT,
            sex TEXT,
            phenotype INTEGER
        )
    """)
    samples = [
        (1, "fam1", "s1", "0", "0", "1", 1),  # male unaffected
        (2, "fam1", "s2", "0", "0", "2", 2),  # female affected
        (3, "fam2", "s3", "0", "0", "1", 2),  # male affected
    ]
    cur.executemany("INSERT INTO samples VALUES (?, ?, ?, ?, ?, ?, ?)", samples)

    # Create variants table with BLOB genotype columns
    cur.execute("""
        CREATE TABLE variants (
            variant_id INTEGER PRIMARY KEY,
            chrom TEXT,
            start INTEGER,
            end INTEGER,
            vcf_id TEXT,
            ref TEXT,
            alt TEXT,
            qual REAL,
            filter TEXT,
            type TEXT,
            gene TEXT,
            info BLOB,
            gts BLOB,
            gt_types BLOB,
            gt_phases BLOB,
            gt_depths BLOB
        )
    """)

    # Variant 1: s1 is HOM_REF, s2 is HET, s3 is HOM_ALT
    v1_gts = np.array(["A/A", "A/G", "G/G"])
    v1_types = np.array([HOM_REF, HET, HOM_ALT], dtype=np.int32)
    v1_phases = np.array([False, False, False], dtype=bool)
    v1_depths = np.array([30, 25, 40], dtype=np.int32)
    v1_info = {"AF": 0.5, "DP": 95}

    # Variant 2: s1 is HET, s2 is HET, s3 is HOM_REF
    v2_gts = np.array(["C/T", "C/T", "C/C"])
    v2_types = np.array([HET, HET, HOM_REF], dtype=np.int32)
    v2_phases = np.array([False, True, False], dtype=bool)
    v2_depths = np.array([20, 22, 18], dtype=np.int32)
    v2_info = {"AF": 0.33, "DP": 60}

    variants = [
        (
            1,
            "chr1",
            100,
            101,
            "rs001",
            "A",
            "G",
            99.0,
            "PASS",
            "snp",
            "GENE1",
            sqlite3.Binary(zdumps(v1_info)),
            sqlite3.Binary(pack_blob(v1_gts)),
            sqlite3.Binary(pack_blob(v1_types)),
            sqlite3.Binary(pack_blob(v1_phases)),
            sqlite3.Binary(pack_blob(v1_depths)),
        ),
        (
            2,
            "chr1",
            500,
            501,
            "rs002",
            "C",
            "T",
            85.5,
            "PASS",
            "snp",
            "GENE2",
            sqlite3.Binary(zdumps(v2_info)),
            sqlite3.Binary(pack_blob(v2_gts)),
            sqlite3.Binary(pack_blob(v2_types)),
            sqlite3.Binary(pack_blob(v2_phases)),
            sqlite3.Binary(pack_blob(v2_depths)),
        ),
    ]

    cur.executemany(
        """
        INSERT INTO variants VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """,
        variants,
    )

    conn.commit()
    conn.close()
    return db_path


def test_basic_query(test_db):
    gq = GeminiQuery(test_db)
    gq.run("SELECT chrom, start, end, gene FROM variants ORDER BY start")
    rows = list(gq)
    assert len(rows) == 2
    assert rows[0]["chrom"] == "chr1"
    assert rows[0]["start"] == 100
    assert rows[0]["gene"] == "GENE1"
    assert rows[1]["start"] == 500
    assert rows[1]["gene"] == "GENE2"


def test_info_tag_unpacking(test_db):
    gq = GeminiQuery(test_db)
    gq.run("SELECT chrom, start, info FROM variants")
    rows = list(gq)
    assert len(rows) == 2
    assert rows[0].info["AF"] == 0.5
    assert rows[0].info["DP"] == 95
    assert rows[1].info["AF"] == 0.33


def test_genotype_column_selection(test_db):
    gq = GeminiQuery(test_db)
    gq.run("SELECT chrom, gts.s1, gts.s2, gt_types.s3 FROM variants ORDER BY start")
    rows = list(gq)
    assert len(rows) == 2

    # Variant 1
    assert rows[0]["gts.s1"] == "A/A"
    assert rows[0]["gts.s2"] == "A/G"
    assert rows[0]["gt_types.s3"] == HOM_ALT

    # Variant 2
    assert rows[1]["gts.s1"] == "C/T"
    assert rows[1]["gts.s2"] == "C/T"
    assert rows[1]["gt_types.s3"] == HOM_REF


def test_genotype_filter(test_db):
    gq = GeminiQuery(test_db)
    gq.run(
        "SELECT chrom, start, gene FROM variants",
        gt_filter="gt_types.s1 == HOM_REF and gt_types.s2 == HET",
    )
    rows = list(gq)
    assert len(rows) == 1
    assert rows[0]["gene"] == "GENE1"


def test_sample_status_grouping(test_db):
    gq = GeminiQuery(test_db)
    gq.run("SELECT chrom, start FROM variants", show_variant_samples=True)
    rows = list(gq)
    assert len(rows) == 2

    # Variant 1: s1=HOM_REF, s2=HET, s3=HOM_ALT
    assert set(rows[0]["variant_samples"]) == {"s2", "s3"}
    assert rows[0]["het_samples"] == ["s2"]
    assert rows[0]["hom_alt_samples"] == ["s3"]
    assert rows[0]["hom_ref_samples"] == ["s1"]


def test_wildcard_filter(test_db):
    gq = GeminiQuery(test_db)
    # phenotype==2 means affected: s2 and s3
    # Variant 1: s2 is HET, s3 is HOM_ALT -> both have variant (!= HOM_REF)
    gq.run(
        "SELECT chrom, start, gene FROM variants",
        gt_filter="(gt_types).(phenotype == 2).(!=HOM_REF).(all)",
    )
    rows = list(gq)
    assert len(rows) == 1
    assert rows[0]["gene"] == "GENE1"


def test_region_filtering(test_db):
    gq = GeminiQuery(test_db)
    gq.run("SELECT chrom, start, end FROM variants", region="chr1:50-150")
    rows = list(gq)
    assert len(rows) == 1
    assert rows[0]["start"] == 100


def test_json_formatter(test_db):
    gq = GeminiQuery(test_db, out_format=JSONRowFormat())
    gq.run("SELECT chrom, start, gene FROM variants WHERE start = 100")
    row = next(gq)
    row_str = str(row)
    assert '"chrom": "chr1"' in row_str
    assert '"gene": "GENE1"' in row_str


def test_vcf_formatter(test_db):
    gq = GeminiQuery(test_db, out_format="vcf")
    gq.run("SELECT chrom, start, ref, alt, gene FROM variants WHERE start = 100")
    row = next(gq)
    vcf_line = str(row)
    parts = vcf_line.split("\t")
    assert parts[0] == "chr1"
    assert parts[1] == "101"  # 1-based start
    assert parts[3] == "A"  # ref
    assert parts[4] == "G"  # alt


def test_carrier_summary(test_db):
    gq = GeminiQuery(test_db, out_format="carrier_summary")
    gq.run("SELECT chrom, start FROM variants WHERE start = 100")
    row = next(gq)
    assert "chr1" in str(row)


def test_sample_mapping_properties(test_db):
    gq = GeminiQuery(test_db)
    assert gq.sample2index["s1"] == 0
    assert gq.sample2index["s2"] == 1
    assert gq.sample2index["s3"] == 2
    assert gq.index2sample[0] == "s1"
    assert gq.index2sample[1] == "s2"
    assert gq.index2sample[2] == "s3"


def test_unsafe_filter_raises(test_db):
    gq = GeminiQuery(test_db)
    with pytest.raises(ValueError, match="invalid --gt-filter"):
        gq.run("SELECT chrom FROM variants", gt_filter="__import__('os').system('ls')")


def test_tped_formatter(test_db):
    gq = GeminiQuery(test_db, out_format="tped")
    gq.run("SELECT chrom, start, end, ref, alt FROM variants LIMIT 1")
    row = next(gq)
    parts = str(row).split(" ")
    assert len(parts) >= 8
    assert parts[0] == "1"


def test_gemini_row_get(test_db):
    gq = GeminiQuery(test_db)
    gq.run("SELECT chrom, start FROM variants LIMIT 1")
    row = next(gq)
    assert row.get("chrom") == "chr1"
    assert row.get("nonexistent", "default_val") == "default_val"


def test_gemini_query_import():
    from gemini import GeminiQuery as GQ

    assert GQ is not None
