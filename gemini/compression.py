import pickle
import zlib
from collections import OrderedDict

import numpy as np

try:
    import snappy
except ImportError:
    snappy = None


def unpack_genotype_blob(blob):
    if blob is None:
        return None
    return pickle.loads(zlib.decompress(blob), encoding="latin1")


def zdumps(obj):
    return zlib.compress(pickle.dumps(obj, pickle.HIGHEST_PROTOCOL), 9)


def pack_blob(obj):
    return memoryview(zdumps(obj))


def unpack_ordereddict_blob(blob):
    if blob is None:
        return None
    blob_val = pickle.loads(zlib.decompress(blob), encoding="latin1")
    if blob_val is not None:
        return OrderedDict(blob_val)
    return None


lookup = {
    dt(1).dtype.char: dt
    for dt in (
        np.uint8,
        np.uint32,
        np.int32,
        np.float32,
        np.int64,
        np.float64,
        np.bool_,
    )
}

SEP = "\0"


def snappy_pack_blob(obj, sep=SEP):
    if snappy is None:
        raise RuntimeError("snappy is not installed")
    if obj is None:
        return b""
    c = obj.dtype.char
    if c in ("S", "U"):
        return b"U" + snappy.compress(sep.join(obj).encode("utf-8"))
    return memoryview(c.encode("utf-8") + snappy.compress(obj.tobytes(), "utf8"))


def snappy_unpack_blob(blob, sep=SEP):
    if snappy is None:
        raise RuntimeError("snappy is not installed")
    if not blob:
        return None
    if blob[0] in (ord("S"), ord("U")):
        return np.array(snappy.decompress(blob[1:], "utf8").split(sep))
    dt = lookup[chr(blob[0])]
    arr = np.frombuffer(snappy.decompress(blob[1:]), dtype=dt)
    if blob[0] == ord("?"):
        arr = np.array(arr)
        arr.setflags(write=True)
    return arr
