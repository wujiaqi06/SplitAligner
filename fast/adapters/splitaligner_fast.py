"""Python 3.8 ctypes interface to FAST-002; contains no Python matcher.

Node lengths use the original node IDs supplied to Mapping. Rebuild the
mapping whenever parent/taxon arrays or node-to-edge bindings change.
"""
import ctypes as C
import os
from pathlib import Path

NONE = C.c_size_t(-1).value


class FastError(RuntimeError):
    def __init__(self, status, message):
        self.status = status
        super().__init__("SAF status {}: {}".format(status, message))


class InvalidArgument(FastError):
    pass


class CapacityError(FastError):
    pass


class NumericOverflow(FastError):
    pass


class CountOverflow(FastError):
    pass


class _ReferenceInfo(C.Structure):
    _fields_ = [(x, C.c_size_t) for x in ("n_taxa", "n_words", "n_axes", "n_nodes")]


class _CoverageInfo(C.Structure):
    _fields_ = [(x, C.c_size_t) for x in ("n_taxa", "n_words", "n_axes", "n_projections", "n_observed")]


class _MappingInfo(C.Structure):
    _fields_ = [(x, C.c_size_t) for x in ("n_nodes", "n_axes", "n_projections", "n_gene_splits", "n_reference_absent", "n_gene_unmatched")]


class _SplitInfo(C.Structure):
    _fields_ = [("side_count", C.c_size_t), ("other_count", C.c_size_t), ("terminal", C.c_uint8)]


class _ValueBuffer(C.Structure):
    _fields_ = [("values", C.POINTER(C.c_double)), ("valid", C.POINTER(C.c_uint8)), ("capacity", C.c_size_t)]


class _SupportCount(C.Structure):
    _fields_ = [(x, C.c_uint64) for x in ("n_valid_topologies", "n_struct", "n_fused", "n_dec", "n_rec", "n_not_recovered")]


def _dict(struct):
    return {name: getattr(struct, name) for name, _ in struct._fields_}


def _integers(values, low, high, name):
    result = tuple(values)
    if any(not isinstance(v, int) or isinstance(v, bool) or v < low or v > high for v in result):
        raise ValueError("{} must contain integers in [{}, {}]".format(name, low, high))
    return result


def _text(value, name):
    if not isinstance(value, str) or not value or "\x00" in value:
        raise ValueError("{} must be a nonempty string without NUL".format(name))
    return value.encode("utf-8")


def validate_gene_id(value):
    """Preserve whitespace: 'g' and ' g ' remain different IDs."""
    if not isinstance(value, str) or not value.strip() or any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError("gene ID must be nonblank and contain no ASCII control characters or DEL")
    return value


def _require_type(value, expected, name):
    if not isinstance(value, expected):
        raise TypeError("{} must be a {}".format(name, expected.__name__))


class Library:
    """Load an explicit shared library, or SAF_LIBRARY / build/release/libsplitaligner_fast.so."""
    def __init__(self, path=None):
        path = path or os.environ.get("SAF_LIBRARY") or str(Path(__file__).resolve().parents[1] / "build" / "release" / "libsplitaligner_fast.so")
        self.path = str(Path(path).resolve())
        self._c = C.CDLL(self.path)
        p, z = C.c_void_p, C.c_size_t
        pz, pu64, pi64 = C.POINTER(z), C.POINTER(C.c_uint64), C.POINTER(C.c_int64)
        pu8, pd = C.POINTER(C.c_uint8), C.POINTER(C.c_double)
        signatures = {
            "reference_create": [z, C.c_char_p, pi64, pi64, z, C.POINTER(p)],
            "reference_get_info": [p, C.POINTER(_ReferenceInfo)],
            "reference_axis": [p, z, pu64, z, pz, C.POINTER(_SplitInfo)],
            "reference_sources": [p, z, pz, z, pz],
            "reference_universe": [p, C.POINTER(C.c_char), z, pz],
            "coverage_create": [p, pu64, z, C.POINTER(p)],
            "coverage_get_info": [p, C.POINTER(_CoverageInfo)],
            "coverage_mask": [p, pu64, z, pz],
            "coverage_projection": [p, z, pu64, z, pz, C.POINTER(_SplitInfo)],
            "coverage_members": [p, z, pz, z, pz],
            "coverage_axes": [p, pz, pu8, pu8, z],
            "mapping_create": [p, pi64, pi64, z, C.POINTER(p)],
            "mapping_get_info": [p, C.POINTER(_MappingInfo)],
            "mapping_gene_split": [p, z, pu64, z, pz, C.POINTER(_SplitInfo)],
            "mapping_gene_sources": [p, z, pz, z, pz],
            "mapping_projection_to_gene": [p, pz, z, pz],
            "mapping_gene_to_projection": [p, pz, z, pz],
            "mapping_reference_absent": [p, pz, z, pz],
            "mapping_gene_unmatched": [p, pz, z, pz],
            "mapping_states": [p, C.POINTER(C.c_char), z],
            "mapping_recovered": [p, pu8, z],
            "mapping_apply": [p, pd, pu8, z, _ValueBuffer, _ValueBuffer, _ValueBuffer],
            "support_create": [p, C.c_char_p, C.c_char_p, C.POINTER(p)],
            "support_update": [p, p], "support_merge": [p, p],
            "support_read": [p, C.POINTER(_SupportCount), z, pz],
        }
        for name, args in signatures.items():
            fn = getattr(self._c, "saf_" + name)
            fn.argtypes, fn.restype = args, C.c_int
        for name in ("reference", "coverage", "mapping", "support"):
            fn = getattr(self._c, "saf_" + name + "_destroy")
            fn.argtypes, fn.restype = [p], None
        for name in ("version", "last_error"):
            fn = getattr(self._c, "saf_" + name)
            fn.argtypes, fn.restype = [], C.c_char_p

    @property
    def version(self):
        return self._c.saf_version().decode("utf-8")

    def _check(self, status):
        if status:
            raw = self._c.saf_last_error()
            message = raw.decode("utf-8", "replace") if raw else "native error without text"
            error = {1: InvalidArgument, 2: CapacityError, 4: NumericOverflow, 5: CountOverflow}.get(status, FastError)
            raise error(status, message)

    def reference(self, n_taxa, universe_id, parent, taxon):
        return Reference(self, n_taxa, universe_id, parent, taxon)


class _Handle:
    def _init(self, lib, kind):
        _require_type(lib, Library, "library")
        self._lib, self._kind, self._ptr = lib, kind, C.c_void_p()

    def __copy__(self):
        raise TypeError("owning native handles cannot be copied; create a new plan")

    def __deepcopy__(self, memo):
        raise TypeError("owning native handles cannot be deep-copied; create a new plan")

    def _handle(self):
        if not self._ptr.value:
            raise RuntimeError("{} handle is closed".format(self._kind))
        return self._ptr

    @property
    def closed(self):
        return not self._ptr.value

    def close(self):
        if getattr(self, "_ptr", None) and self._ptr.value:
            getattr(self._lib._c, "saf_" + self._kind + "_destroy")(self._ptr)
            self._ptr = C.c_void_p()

    def __enter__(self):
        self._handle()
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        self.close()

    def _call(self, suffix, *args):
        handle = self._handle()
        status = getattr(self._lib._c, "saf_" + self._kind + "_" + suffix)(handle, *args)
        self._lib._check(status)

    def _info(self, typ):
        result = typ()
        self._call("get_info", C.byref(result))
        return _dict(result)

    def _vector(self, suffix, typ=C.c_size_t, index=None, info=False):
        count = C.c_size_t()
        prefix = () if index is None else (index,)
        extra = _SplitInfo() if info else None
        tail = (C.byref(extra),) if info else ()
        self._call(suffix, *(prefix + (None, 0, C.byref(count)) + tail))
        array = (typ * count.value)()
        self._call(suffix, *(prefix + (array, count.value, C.byref(count)) + tail))
        if info:
            result = _dict(extra)
            result["terminal"] = bool(result["terminal"])
            result["words"] = tuple(array)
            return result
        return tuple(array)


def _topology(parent, taxon):
    parent = _integers(parent, -1, (1 << 63) - 1, "parent")
    taxon = _integers(taxon, -1, (1 << 63) - 1, "taxon")
    if len(parent) != len(taxon):
        raise ValueError("parent and taxon lengths differ")
    return (C.c_int64 * len(parent))(*parent), (C.c_int64 * len(taxon))(*taxon)


class Reference(_Handle):
    def __init__(self, library, n_taxa, universe_id, parent, taxon):
        self._init(library, "reference")
        _integers([n_taxa], 0, NONE, "n_taxa")
        p, t = _topology(parent, taxon)
        library._check(library._c.saf_reference_create(n_taxa, _text(universe_id, "universe_id"), p, t, len(p), C.byref(self._ptr)))
        self.info = self._info(_ReferenceInfo)

    def axes(self):
        result = []
        for axis in range(self.info["n_axes"]):
            entry = self._vector("axis", C.c_uint64, axis, True)
            entry.update(axis_id=axis, sources=self._vector("sources", index=axis))
            result.append(entry)
        return result

    def universe(self):
        return b"".join(self._vector("universe", C.c_char)).rstrip(b"\x00").decode("utf-8")

    def coverage(self, observed_taxa=None, mask_words=None):
        if (observed_taxa is None) == (mask_words is None):
            raise ValueError("supply exactly one of observed_taxa and mask_words")
        if observed_taxa is not None:
            taxa = _integers(observed_taxa, 0, self.info["n_taxa"] - 1, "observed_taxa")
            if len(set(taxa)) != len(taxa):
                raise ValueError("duplicate observed taxon")
            mask_words = [0] * self.info["n_words"]
            for taxon in taxa:
                mask_words[taxon // 64] |= 1 << (taxon % 64)
        return Coverage(self, mask_words)

    def support(self, unit_id="gene", cohort_id="default"):
        return Support(self, unit_id, cohort_id)


class Coverage(_Handle):
    def __init__(self, reference, mask_words):
        _require_type(reference, Reference, "reference")
        reference_handle = reference._handle()
        self._init(reference._lib, "coverage")
        words = _integers(mask_words, 0, (1 << 64) - 1, "mask_words")
        array = (C.c_uint64 * len(words))(*words)
        self._lib._check(self._lib._c.saf_coverage_create(reference_handle, array, len(words), C.byref(self._ptr)))
        self.info = self._info(_CoverageInfo)

    def truth(self):
        size = self.info["n_axes"]
        group, structural, eligible = (C.c_size_t * size)(), (C.c_uint8 * size)(), (C.c_uint8 * size)()
        self._call("axes", group, structural, eligible, size)
        projections = []
        for projection in range(self.info["n_projections"]):
            entry = self._vector("projection", C.c_uint64, projection, True)
            entry.update(projection_id=projection, members=self._vector("members", index=projection))
            projections.append(entry)
        return {"mask_words": self._vector("mask", C.c_uint64),
                "axis_to_projection": tuple(None if x == NONE else x for x in group),
                "structural": tuple(bool(x) for x in structural),
                "eligible": tuple(bool(x) for x in eligible), "projections": projections}

    def mapping(self, parent, taxon, gene_id=None):
        return Mapping(self, parent, taxon, gene_id)


class Mapping(_Handle):
    def __init__(self, coverage, parent, taxon, gene_id=None):
        _require_type(coverage, Coverage, "coverage")
        coverage_handle = coverage._handle()
        self._init(coverage._lib, "mapping")
        self.gene_id = None if gene_id is None else validate_gene_id(gene_id)
        p, t = _topology(parent, taxon)
        self._lib._check(self._lib._c.saf_mapping_create(coverage_handle, p, t, len(p), C.byref(self._ptr)))
        self.info = self._info(_MappingInfo)

    def ledger(self):
        splits = []
        for edge in range(self.info["n_gene_splits"]):
            entry = self._vector("gene_split", C.c_uint64, edge, True)
            entry.update(gene_edge_id=edge, sources=self._vector("gene_sources", index=edge))
            splits.append(entry)
        states = (C.c_char * self.info["n_axes"])()
        recovered = (C.c_uint8 * self.info["n_projections"])()
        self._call("states", states, len(states))
        self._call("recovered", recovered, len(recovered))
        result = {"gene_id": self.gene_id, "gene_splits": splits,
                  "states": bytes(states).decode("ascii"), "recovered": tuple(bool(x) for x in recovered)}
        for name in ("projection_to_gene", "gene_to_projection", "reference_absent", "gene_unmatched"):
            result[name] = tuple(None if x == NONE else x for x in self._vector(name))
        return result

    def apply(self, node_lengths=None, length_valid=None, primitive=True, projection=True, gene_edge=True):
        """Requested families map to (values, valid) tuples; skipped keys absent.

        Buffers are disjoint and Python-owned. Errors return no output. Original
        node IDs are mandatory; count equality cannot detect a wrong permutation.
        """
        self._handle()
        requested = (bool(primitive), bool(projection), bool(gene_edge))
        if any(requested):
            if node_lengths is None or len(node_lengths) != self.info["n_nodes"]:
                raise ValueError("node_lengths must use the original layout and exact node count")
            lengths = (C.c_double * len(node_lengths))(*node_lengths)
            if length_valid is None:
                valid = None
            else:
                if len(length_valid) != len(node_lengths):
                    raise ValueError("length_valid length differs from node_lengths")
                valid = (C.c_uint8 * len(length_valid))(*(1 if x else 0 for x in length_valid))
            count = len(node_lengths)
        else:
            lengths, valid, count = None, None, 0
        buffers, owners = [], []
        for need, size in zip(requested, (self.info["n_axes"], self.info["n_projections"], self.info["n_gene_splits"])):
            if need:
                values, available = (C.c_double * size)(), (C.c_uint8 * size)()
                buffers.append(_ValueBuffer(values, available, size))
                owners.append((values, available))
            else:
                buffers.append(_ValueBuffer(None, None, 0))
                owners.append(None)
        self._call("apply", lengths, valid, count, *buffers)
        return {name: (tuple(owner[0]), tuple(bool(x) for x in owner[1]))
                for name, owner in zip(("primitive", "projection", "gene_edge"), owners) if owner is not None}


class Support(_Handle):
    def __init__(self, reference, unit_id, cohort_id):
        _require_type(reference, Reference, "reference")
        reference_handle = reference._handle()
        self._init(reference._lib, "support")
        self._lib._check(self._lib._c.saf_support_create(reference_handle, _text(unit_id, "unit_id"), _text(cohort_id, "cohort_id"), C.byref(self._ptr)))

    def update(self, observation):
        _require_type(observation, Mapping, "observation")
        self._handle()
        observation_handle = observation._handle()
        if observation._lib is not self._lib:
            raise ValueError("mapping and support must use the same Library instance")
        self._call("update", observation_handle)

    def merge(self, disjoint_shard):
        _require_type(disjoint_shard, Support, "disjoint_shard")
        self._handle()
        shard_handle = disjoint_shard._handle()
        if disjoint_shard._lib is not self._lib:
            raise ValueError("support shards must use the same Library instance")
        self._call("merge", shard_handle)

    def counts(self):
        result = []
        for axis, record in enumerate(self._vector("read", _SupportCount)):
            entry = _dict(record)
            fraction = entry["n_rec"] / entry["n_dec"] if entry["n_dec"] else None
            entry.update(axis_id=axis, support_fraction=fraction,
                         support_percent=None if fraction is None else 100 * fraction)
            result.append(entry)
        return result
