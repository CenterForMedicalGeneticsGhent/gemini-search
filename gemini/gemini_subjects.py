import sqlite3
from collections import defaultdict


class Subject:
    """Describe a single subject in the samples table."""

    def __init__(self, row):
        if hasattr(row, "keys"):
            self.__dict__.update(dict(row))
        elif isinstance(row, dict):
            self.__dict__.update(row)
        else:
            self.__dict__.update(getattr(row, "__dict__", {}))

        self.phenotype = int(self.phenotype) if self._has_phenotype() else None
        self._set_affected_status()

    def __repr__(self):
        return "\t".join(
            map(
                str,
                [
                    getattr(self, "name", ""),
                    getattr(self, "paternal_id", ""),
                    getattr(self, "maternal_id", ""),
                    self.phenotype,
                ],
            )
        )

    def _has_phenotype(self):
        return hasattr(self, "phenotype") and self.phenotype is not None

    def _set_affected_status(self):
        pheno = str(self.phenotype) if self.phenotype is not None else None
        if pheno == "2":
            self.affected = True
        elif pheno == "1":
            self.affected = False
        else:
            self.affected = None


def get_subjects(db_or_conn, sample_filter=None):
    """Return dictionary mapping subject name to Subject object."""
    should_close = False
    if isinstance(db_or_conn, str):
        conn = sqlite3.connect(db_or_conn)
        conn.row_factory = sqlite3.Row
        should_close = True
    elif hasattr(db_or_conn, "conn"):
        conn = db_or_conn.conn
    else:
        conn = db_or_conn
        if conn.row_factory is None:
            conn.row_factory = sqlite3.Row

    try:
        cursor = conn.cursor()
        query = "SELECT * FROM samples"
        if sample_filter:
            query += f" WHERE {sample_filter}"

        cursor.execute(query)
        samples_dict = {}
        for row in cursor.fetchall():
            sub = Subject(row)
            samples_dict[sub.name] = sub
        return samples_dict
    finally:
        if should_close:
            conn.close()


def get_family_dict(db_or_conn, sample_filter=None):
    """Return dictionary mapping family_id to list of Subject objects."""
    families = defaultdict(list)
    subjects = get_subjects(db_or_conn, sample_filter=sample_filter)
    for subject in subjects.values():
        fam_id = getattr(subject, "family_id", None)
        if fam_id is not None:
            families[fam_id].append(subject)
    return families
