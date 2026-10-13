'''
Tests for rapster/ns_formation.py and tools/make_ns_formation_table.py.

Run from the repository root with either
    python3 -m pytest tests/test_ns_formation.py
or
    python3 tests/test_ns_formation.py
'''

import os
import sys
import tempfile
import importlib.util

import numpy as np
from scipy.stats import ks_2samp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from rapster.ns_formation import NSFormationTable, sample_ns_formation_times

_spec = importlib.util.spec_from_file_location('mk', os.path.join(ROOT, 'tools', 'make_ns_formation_table.py'))
mk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mk)

M_SAMPLE, T_SAMPLE = mk.synthetic_sample(500_000, seed=3)
_TMPDIR = tempfile.mkdtemp()
TABLE_PATH = os.path.join(_TMPDIR, 'NS_formation_test.npz')
np.savez_compressed(TABLE_PATH, source=np.array('test'), **mk.build_table(M_SAMPLE, T_SAMPLE))
TABLE = NSFormationTable(TABLE_PATH)


def test_model0_is_legacy_and_draws_no_random_numbers():
    np.random.seed(123)
    t = sample_ns_formation_times(np.full(1000, 1.4), 0, 3.5, table=TABLE)
    after = np.random.rand()
    np.random.seed(123)
    assert np.all(t == 3.5)
    assert after == np.random.rand()


def test_marginal_reproduces_input_distribution():
    np.random.seed(1)
    t = sample_ns_formation_times(np.ones(50_000), 1, 3.5, table=TABLE)
    assert ks_2samp(np.log10(t), np.log10(T_SAMPLE)).pvalue > 1e-3
    assert t.min() >= TABLE.t_min and t.max() <= TABLE.t_max


def test_marginal_is_independent_of_mass():
    np.random.seed(2)
    t_light = sample_ns_formation_times(np.full(20_000, 1.3), 1, 3.5, table=TABLE)
    t_heavy = sample_ns_formation_times(np.full(20_000, 2.2), 1, 3.5, table=TABLE)
    assert ks_2samp(t_light, t_heavy).pvalue > 1e-3


def test_conditional_matches_input_within_each_bin():
    np.random.seed(4)
    for j in (0, TABLE.m_centers.size // 2, TABLE.m_centers.size - 1):
        lo, hi = TABLE.m_edges[j]
        in_bin = (M_SAMPLE >= lo) & (M_SAMPLE < hi)
        t = sample_ns_formation_times(np.full(20_000, TABLE.m_centers[j]), 2, 3.5, table=TABLE)
        assert ks_2samp(np.log10(t), np.log10(T_SAMPLE[in_bin])).pvalue > 1e-3, j


def test_conditional_heavier_NSs_form_earlier():
    np.random.seed(5)
    t_light = sample_ns_formation_times(np.full(20_000, 1.35), 2, 3.5, table=TABLE)
    t_heavy = sample_ns_formation_times(np.full(20_000, 2.2), 2, 3.5, table=TABLE)
    assert np.median(t_heavy) < np.median(t_light)


def test_conditional_out_of_range_masses_use_edge_bins():
    idx = TABLE.mass_bin_index([0.5, 1e3])
    assert idx[0] == 0 and idx[1] == TABLE.m_centers.size - 1


def test_reproducible_with_seed():
    m = np.random.default_rng(0).uniform(1.2, 2.3, 1000)
    np.random.seed(99); t1 = sample_ns_formation_times(m, 2, 3.5, table=TABLE)
    np.random.seed(99); t2 = sample_ns_formation_times(m, 2, 3.5, table=TABLE)
    assert np.array_equal(t1, t2)


def test_empty_input():
    assert sample_ns_formation_times(np.array([]), 2, 3.5, table=TABLE).size == 0


def test_bad_table_is_rejected():
    bad = os.path.join(_TMPDIR, 'bad.npz')
    np.savez(bad, q=np.linspace(0, 1, 5))
    try:
        NSFormationTable(bad)
    except ValueError:
        pass
    else:
        raise AssertionError('incomplete table was accepted')

    truncated = os.path.join(_TMPDIR, 'truncated.npz')
    with open(TABLE_PATH, 'rb') as f, open(truncated, 'wb') as g:
        g.write(f.read()[:1000])
    try:
        NSFormationTable(truncated)
    except ValueError:
        pass
    else:
        raise AssertionError('truncated table was accepted')


def test_unknown_model_is_rejected():
    try:
        sample_ns_formation_times([1.4], 3, 3.5, table=TABLE)
    except ValueError:
        pass
    else:
        raise AssertionError('model 3 was accepted')


if __name__ == '__main__':
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    for test in tests:
        test()
        print(f'PASS  {test.__name__}')
    print(f'\nall {len(tests)} tests passed')
