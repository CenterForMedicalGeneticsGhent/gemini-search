import numpy as np


def _to_str(val):
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="replace")
    if isinstance(val, np.ndarray):
        return ",".join(_to_str(x) for x in val)
    return str(val)


def to_json(obj):
    return dict(obj.items())


class PDict:
    """Ordered dictionary-like object that supports repeated keys."""

    __slots__ = ("_keys", "_vals")

    def __init__(self, items=None):
        if items is None:
            self._keys = []
            self._vals = []
        elif isinstance(items, dict):
            self._keys = list(items.keys())
            self._vals = list(items.values())
        elif isinstance(items, PDict):
            self._keys = items._keys[:]
            self._vals = items._vals[:]
        elif hasattr(items, "__iter__"):
            self._keys = [i[0] for i in items]
            self._vals = [i[1] for i in items]
        else:
            raise TypeError(f"Cannot initialize PDict from {type(items)}")

    def __repr__(self):
        return "{" + ", ".join(f"{k!r}: {v!r}" for k, v in zip(self._keys, self._vals)) + "}"

    def __str__(self):
        return "\t".join(_to_str(v) for v in self._vals)

    def add(self, key, val):
        self._keys.append(key)
        self._vals.append(val)

    def get(self, key, default=None):
        try:
            idx = self._keys.index(key)
        except ValueError:
            return default
        val = self._vals[idx]
        return val.decode("utf-8", errors="replace") if isinstance(val, bytes) else val

    def __getitem__(self, key):
        try:
            idx = self._keys.index(key)
        except ValueError:
            raise KeyError(key) from None
        val = self._vals[idx]
        return val.decode("utf-8", errors="replace") if isinstance(val, bytes) else val

    def __setitem__(self, key, val):
        try:
            idx = self._keys.index(key)
        except ValueError:
            self._keys.append(key)
            self._vals.append(val)
            return
        self._vals[idx] = val

    def keys(self):
        return list(self._keys)

    def values(self):
        return list(self._vals)

    def copy(self):
        return self.__class__(self)

    def __iter__(self):
        yield from self._keys

    def __delitem__(self, key):
        idx = self._keys.index(key)
        del self._keys[idx]
        del self._vals[idx]

    def items(self):
        return list(zip(self._keys, self._vals))
