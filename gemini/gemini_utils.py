import collections.abc
import sqlite3

import numpy as np


def to_str(s):
    if hasattr(s, "decode"):
        return s.decode("utf-8", errors="replace")
    if isinstance(s, np.str_):
        return str(s)
    return s


def get_gt_cols(cursor):
    """Return BLOB genotype column names starting with 'gt' from variants table."""
    try:
        cursor.execute("PRAGMA table_info('variants')")
        cols = cursor.fetchall()
        return [
            row["name"]
            for row in cols
            if row["name"].startswith("gt") and row["type"].upper() == "BLOB"
        ]
    except sqlite3.OperationalError:
        return []


class OrderedSet(collections.abc.MutableSet):
    """Ordered set implementation based on dict."""

    def __init__(self, iterable=None):
        self._d = dict.fromkeys(iterable) if iterable is not None else {}

    def __contains__(self, key):
        return key in self._d

    def __iter__(self):
        return iter(self._d)

    def __len__(self):
        return len(self._d)

    def add(self, key):
        self._d[key] = None

    def discard(self, key):
        self._d.pop(key, None)

    def __repr__(self):
        return f"{self.__class__.__name__}({list(self._d)})"


def itersubclasses(cls, _seen=None):
    """Generator over all subclasses of a given class, in depth first order."""
    if _seen is None:
        _seen = set()
    for sub in cls.__subclasses__():
        if sub not in _seen:
            _seen.add(sub)
            yield sub
            yield from itersubclasses(sub, _seen)
