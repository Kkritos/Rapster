'''
Tests for the NS set-up in rapster/cluster_evolution.py:initialize_cluster
(-NS mass models combined with the -NSt formation-time models).

Run from the repository root with either
    python3 -m pytest tests/test_ns_initialization.py
or
    python3 tests/test_ns_initialization.py
'''

import io
import os
import sys
import contextlib

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from rapster.constants import tBH_form
from rapster.run_cluster import parse_args
from rapster.cluster_evolution import initialize_cluster
from rapster.stellar_evolution import init_stellar_mass_sampler
from rapster.functions import set_escape_factor

init_stellar_mass_sampler(0.0)
set_escape_factor(2.0)

# low kicks, so that many NSs are retained:
BASE_ARGS = ['-N', '2e5', '-S', '42', '-wK', '20']


def make_state(ns, nst, extra=()):
    sys.argv = ['rapster'] + BASE_ARGS + ['-NS', str(ns), '-NSt', str(nst)] + list(extra)
    return initialize_cluster(parse_args())


def n_ns_in_mBH(state):
    return int(np.sum(state['mBH'] < 3.0))


def test_tBH_form_is_10_Myr():
    assert tBH_form == 10.0


def test_NSt0_puts_retained_NSs_in_mBH_and_nothing_pending():
    for ns in (1, 2, 3):
        s = make_state(ns, 0)
        assert n_ns_in_mBH(s) > 0, ns
        assert s['mNS_pending'].size == 0 and s['tNS_pending'].size == 0 and s['vNS_pending'].size == 0
        assert s['mBH'].size == s['sBH'].size == s['gBH'].size == s['hBH'].size


def test_NSt_positive_defers_all_NSs():
    for ns in (1, 2, 3):
        for nst in (1, 2):
            s = make_state(ns, nst)
            assert n_ns_in_mBH(s) == 0, (ns, nst)
            m, t, v = s['mNS_pending'], s['tNS_pending'], s['vNS_pending']
            assert m.size > 0 and m.size == t.size == v.size
            assert np.all(np.diff(t) >= 0)                        # sorted by formation time
            assert np.all(t > 0) and np.all(v >= 0)


def test_BHs_unchanged_by_NSt():
    # BH masses are drawn before any NS random numbers, so they must not depend on -NSt:
    m0 = make_state(2, 0)['mBH']
    m1 = make_state(2, 1)['mBH']
    m2 = make_state(2, 2)['mBH']
    bh0 = m0[m0 >= 3.0]
    assert np.array_equal(bh0, m1) and np.array_equal(bh0, m2)


def test_pending_masses_follow_NS_model():
    s1 = make_state(1, 1)
    assert np.allclose(s1['mNS_pending'], 1.40)
    s3 = make_state(3, 2)
    assert s3['mNS_pending'].min() >= 1.0 and s3['mNS_pending'].std() > 0.1


def test_NSt2_heavier_NSs_form_earlier():
    s = make_state(3, 2, extra=['-N', '2e6'])
    m, t = s['mNS_pending'], s['tNS_pending']
    assert np.median(t[m > 2.0]) < np.median(t[m < 1.4])


def test_NSt_requires_NSs():
    sys.argv = ['rapster'] + BASE_ARGS + ['-NS', '0', '-NSt', '1']
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            parse_args()
    except SystemExit:
        pass
    else:
        raise AssertionError('-NS 0 with -NSt 1 was accepted')


if __name__ == '__main__':
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    for test in tests:
        test()
        print(f'PASS  {test.__name__}')
    print(f'\nall {len(tests)} tests passed')
