#!/usr/bin/env python3
"""Bounded Python 3.8 ABI tests; Linux only, SAF_LIBRARY required."""
import ctypes as C
import copy
import math
from pathlib import Path
import sys
import unittest
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'adapters'))
from splitaligner_fast import (Library, Reference, Coverage, Mapping, Support,
                              InvalidArgument, NumericOverflow, validate_gene_id, _ValueBuffer)

class BindingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib = Library()

    def two(self):
        return self.lib.reference(2, 'two:A,B', [-1, 0, 0], [-1, 0, 1])

    def test_parent_close_and_full_queries(self):
        r = self.two()
        self.assertEqual(r.universe(), 'two:A,B')
        self.assertEqual(r.axes()[0]['sources'], (1, 2))
        c = r.coverage([0, 1])
        m = c.mapping([-1, 0, 0], [-1, 0, 1], gene_id=' g ')
        s = r.support()
        r.close()
        self.assertEqual(c.truth()['projections'][0]['members'], (0,))
        c.close()
        ledger = m.ledger()
        self.assertEqual(ledger['states'], 'D')
        self.assertEqual(ledger['gene_splits'][0]['sources'], (1, 2))
        self.assertEqual(ledger['projection_to_gene'], (0,))
        self.assertEqual(ledger['gene_to_projection'], (0,))
        self.assertEqual(ledger['reference_absent'], ())
        self.assertEqual(ledger['gene_unmatched'], ())
        s.update(m)
        self.assertEqual(s.counts()[0]['n_rec'], 1)
        m.close(); m.close()
        with self.assertRaises(RuntimeError): m.ledger()
        s.close()

    def test_values_skip_missing_root_and_reuse(self):
        with self.two() as r, r.coverage([0, 1]) as c, c.mapping([-1, 0, 0], [-1, 0, 1]) as m:
            first = m.apply([float('inf'), 0.25, 0.5])
            self.assertEqual(first['primitive'], ((0.75,), (True,)))
            self.assertEqual(m.apply([float('nan'), -0.5, 0.25])['gene_edge'], ((-0.25,), (True,)))
            missing = m.apply([0, 1e308, float('inf')])
            self.assertFalse(missing['projection'][1][0])
            self.assertTrue(math.isnan(missing['projection'][0][0]))
            self.assertFalse(m.apply([0, 1e308, 1e308], [True, False, True])['gene_edge'][1][0])
            with self.assertRaises(NumericOverflow): m.apply([0, 1e308, 1e308])
            self.assertEqual(m.apply([0, 0.25, 0.5])['primitive'], first['primitive'])
            self.assertEqual(m.apply(primitive=False, projection=False, gene_edge=False), {})
            self.assertEqual(set(m.apply([0, 0, 0], projection=False, gene_edge=False)), {'primitive'})
            sub = float.fromhex('0x0.0000000000001p-1022')
            self.assertEqual(m.apply([0, sub, 0])['gene_edge'][0][0], sub)
            with self.assertRaises(ValueError): m.apply([1, 2])
            self.assertEqual(m.ledger()['states'], 'D')

    def test_input_and_identity_errors(self):
        with self.assertRaises(ValueError): self.lib.reference(2, 'bad', [-1, 0], [-1, 0, 1])
        with self.assertRaises(InvalidArgument): self.lib.reference(2, 'bad', [-1, 0, 0], [-1, 0, 0])
        with self.two() as r:
            with self.assertRaises(InvalidArgument): r.coverage(mask_words=[4])
            with self.assertRaises(ValueError): r.coverage([0, 0])
            with r.coverage([0]) as c:
                self.assertTrue(all(c.truth()['structural']))
                with self.assertRaises(InvalidArgument): c.mapping([-1, 0, 0], [-1, 0, 1])
            with self.two() as other, other.coverage([0, 1]) as c, c.mapping([-1, 0, 0], [-1, 0, 1]) as m, r.support() as s:
                with self.assertRaises(InvalidArgument): s.update(m)
            with r.support(cohort_id='a') as a, r.support(cohort_id='b') as b:
                with self.assertRaises(InvalidArgument): a.merge(b)
        for gene in ('', ' ', 'a\tb', 'a\nb', 'a\x00b'):
            with self.assertRaises(ValueError): validate_gene_id(gene)
        self.assertEqual(validate_gene_id(' g '), ' g ')

    def test_capacity_does_not_write(self):
        with self.two() as r, r.coverage([0, 1]) as c, c.mapping([-1, 0, 0], [-1, 0, 1]) as m:
            count = C.c_size_t()
            sentinel = (C.c_size_t * 2)(999, 998)
            code = self.lib._c.saf_mapping_gene_sources(m._handle(), 0, sentinel, 1, C.byref(count))
            self.assertEqual((code, count.value, tuple(sentinel)), (2, 2, (999, 998)))
            values, valid = (C.c_double * 2)(123, 456), (C.c_uint8 * 2)(7, 8)
            lengths = (C.c_double * 3)(0, 0.25, 0.5)
            code = self.lib._c.saf_mapping_apply(m._handle(), lengths, None, 3, _ValueBuffer(values, valid, 0), _ValueBuffer(), _ValueBuffer())
            self.assertNotEqual(code, 0)
            self.assertEqual((tuple(values), tuple(valid)), ((123, 456), (7, 8)))
            code = self.lib._c.saf_mapping_apply(m._handle(), lengths, None, 3, _ValueBuffer(values, valid, 2), _ValueBuffer(), _ValueBuffer())
            self.assertEqual((code, tuple(values), tuple(valid)), (0, (0.75, 456), (1, 8)))

    def test_parallel_apply_and_mcmc_rebuild(self):
        p, t, proposal_t = [-1, 0, 1, 1, 0, 0], [-1, -1, 0, 1, 2, 3], [-1, -1, 0, 2, 1, 3]
        with self.lib.reference(4, 'ABCD', p, t) as r, r.coverage(range(4)) as c, c.mapping(p, t) as accepted, c.mapping(p, proposal_t) as proposal:
            before = accepted.ledger()
            self.assertNotEqual(before['states'], proposal.ledger()['states'])
            proposal.close()
            self.assertEqual(accepted.ledger(), before)
            with ThreadPoolExecutor(max_workers=2) as pool:
                outcomes = list(pool.map(lambda x: accepted.apply([0, x, x, x, x, x]), (0.25, 0.5)))
            self.assertEqual(tuple(2 * x for x in outcomes[0]['primitive'][0]), outcomes[1]['primitive'][0])
            perm = [4, 1, 5, 0, 3, 2]
            inv = {old: new for new, old in enumerate(perm)}
            pp = [-1 if p[old] == -1 else inv[p[old]] for old in perm]
            tt = [t[old] for old in perm]
            with c.mapping(pp, tt) as rebuilt:
                lengths = [0, 0.25, 0.5, 1, 2, 4]
                self.assertEqual(accepted.apply(lengths)['primitive'], rebuilt.apply([lengths[old] for old in perm])['primitive'])



class _FakeNative:
    """Safe ownership spy: no native library, dereference, or actual free occurs."""
    def __init__(self):
        self.calls, self.live, self.destroyed = [], set(), set()
        self.serial = 100

    def __getattr__(self, name):
        if not name.startswith('saf_'):
            raise AttributeError(name)
        def call(*args):
            self.calls.append(name)
            if name.endswith('_create'):
                self.serial += 1
                self.live.add(self.serial)
                C.cast(args[-1], C.POINTER(C.c_void_p))[0] = self.serial
                return 0
            if name.endswith('_destroy'):
                pointer = args[0].value
                if pointer not in self.live or pointer in self.destroyed:
                    raise AssertionError('fake backend detected repeated or invalid destruction')
                self.live.remove(pointer)
                self.destroyed.add(pointer)
                return None
            if name.endswith('_get_info'):
                value = args[1]._obj
                dimensions = {'n_taxa': 2, 'n_words': 1, 'n_axes': 1, 'n_nodes': 3,
                              'n_projections': 1, 'n_observed': 2, 'n_gene_splits': 1,
                              'n_reference_absent': 0, 'n_gene_unmatched': 0}
                for field, _ in value._fields_:
                    setattr(value, field, dimensions[field])
                return 0
            if name in ('saf_support_update', 'saf_support_merge'):
                return 0
            raise AssertionError('unexpected fake FFI call: ' + name)
        return call


class _PoisonObject:
    @property
    def _lib(self):
        raise AssertionError('wrong-kind object was inspected before its type check')

    def _handle(self):
        raise AssertionError('wrong-kind object reached handle access')


class BoundaryTests(unittest.TestCase):
    """Python-only safety regressions, executable without loading the real ABI."""
    def fake_library(self):
        library = Library.__new__(Library)
        library._c = _FakeNative()
        library.path = '<fake-backend-no-native-code>'
        self.libraries.append(library)
        return library

    def own(self, handle):
        self.owned.append(handle)
        return handle

    def setUp(self):
        self.owned = []
        self.libraries = []
        self.lib = self.fake_library()
        self.reference = self.own(self.lib.reference(2, 'A,B', [-1, 0, 0], [-1, 0, 1]))
        self.coverage = self.own(self.reference.coverage([0, 1]))
        self.mapping = self.own(self.coverage.mapping([-1, 0, 0], [-1, 0, 1]))
        self.support = self.own(self.reference.support())

    def tearDown(self):
        for handle in reversed(self.owned):
            handle.close()
        for library in self.libraries:
            self.assertFalse(library._c.live)

    def untouched(self, exception, operation, library=None):
        library = library or self.lib
        before = list(library._c.calls)
        with self.assertRaises(exception): operation()
        self.assertEqual(library._c.calls, before, 'rejected operation called the fake FFI')

    def test_owning_handles_reject_shallow_and_deep_copy(self):
        for handle in (self.reference, self.coverage, self.mapping, self.support):
            for copier in (copy.copy, copy.deepcopy):
                with self.subTest(kind=type(handle).__name__, copier=copier.__name__):
                    self.untouched(TypeError, lambda: copier(handle))
                    self.assertFalse(handle.closed)
        # Nested deepcopy must not manufacture a second owner either.
        self.untouched(TypeError, lambda: copy.deepcopy({'mapping': self.mapping}))
        self.mapping.close()
        self.untouched(TypeError, lambda: copy.copy(self.mapping))
        self.untouched(TypeError, lambda: copy.deepcopy(self.mapping))

    def test_wrong_kind_rejected_before_library_handle_or_ffi(self):
        values = [self.reference, self.coverage, self.mapping, self.support, None, object(), _PoisonObject()]
        constructors = [(Reference, lambda value: Reference(value, 2, 'A,B', [-1, 0, 0], [-1, 0, 1]), Library),
                        (Coverage, lambda value: Coverage(value, [3]), Reference),
                        (Mapping, lambda value: Mapping(value, [-1, 0, 0], [-1, 0, 1]), Coverage),
                        (Support, lambda value: Support(value, 'gene', 'cohort'), Reference)]
        for name, make, expected_type in constructors:
            for value in values:
                if isinstance(value, expected_type): continue
                with self.subTest(constructor=name.__name__, input_type=type(value).__name__):
                    self.untouched(TypeError, lambda: make(value))
        for value in values:
            if not isinstance(value, Mapping):
                self.untouched(TypeError, lambda: self.support.update(value))
            if not isinstance(value, Support):
                self.untouched(TypeError, lambda: self.support.merge(value))
        # Even closed receivers reject a wrong argument kind first.
        self.support.close()
        self.untouched(TypeError, lambda: self.support.update(self.reference))
        self.untouched(TypeError, lambda: self.support.merge(self.mapping))

    def test_closed_and_cross_library_handles_never_reach_ffi(self):
        other_lib = self.fake_library()
        other_ref = self.own(other_lib.reference(2, 'A,B', [-1, 0, 0], [-1, 0, 1]))
        other_cov = self.own(other_ref.coverage([0, 1]))
        other_mapping = self.own(other_cov.mapping([-1, 0, 0], [-1, 0, 1]))
        other_support = self.own(other_ref.support())
        foreign_before = list(other_lib._c.calls)
        self.untouched(ValueError, lambda: self.support.update(other_mapping))
        self.untouched(ValueError, lambda: self.support.merge(other_support))
        self.assertEqual(other_lib._c.calls, foreign_before)
        self.mapping.close()
        self.untouched(RuntimeError, lambda: self.support.update(self.mapping))
        other_support.close()
        self.untouched(RuntimeError, lambda: self.support.merge(other_support))
        self.coverage.close()
        self.untouched(RuntimeError, lambda: Mapping(self.coverage, [-1, 0, 0], [-1, 0, 1]))
        self.reference.close()
        self.untouched(RuntimeError, lambda: Coverage(self.reference, [3]))
        self.untouched(RuntimeError, lambda: Support(self.reference, 'gene', 'cohort'))
        self.support.close()
        self.untouched(RuntimeError, lambda: self.support.update(other_mapping))

    def test_valid_handles_still_create_update_merge_and_close_once(self):
        shard = self.own(self.reference.support())
        before = len(self.lib._c.calls)
        self.support.update(self.mapping)
        self.support.merge(shard)
        self.assertEqual(self.lib._c.calls[before:], ['saf_support_update', 'saf_support_merge'])
        pointer = self.mapping._ptr.value
        self.mapping.close()
        after = list(self.lib._c.calls)
        self.mapping.close()
        self.assertEqual(self.lib._c.calls, after)
        self.assertIn(pointer, self.lib._c.destroyed)
        self.assertEqual(self.reference.info['n_taxa'], 2)
        self.assertEqual(self.coverage.info['n_observed'], 2)

    def test_all_ascii_controls_and_del_rejected_ids_preserved(self):
        for code in list(range(32)) + [127]:
            with self.subTest(code=code):
                self.untouched(ValueError, lambda: self.coverage.mapping([-1, 0, 0], [-1, 0, 1], 'g' + chr(code) + 'x'))
                with self.assertRaises(ValueError): validate_gene_id('g' + chr(code) + 'x')
        for gene in (' g ', r'g\n', r' g\t ', '猫 gene', 'x\u00a0 y'):
            self.assertEqual(validate_gene_id(gene), gene)
            mapping = self.own(self.coverage.mapping([-1, 0, 0], [-1, 0, 1], gene))
            self.assertEqual(mapping.gene_id, gene)


if __name__ == '__main__':
    unittest.main(verbosity=2)
