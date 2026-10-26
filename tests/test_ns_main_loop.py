'''
Tests for NS formation during the run (rapster/cluster_evolution.py:form_neutron_stars
and the BH-subsystem activation in compute_cluster_properties).

Run from the repository root with either
    python3 -m pytest tests/test_ns_main_loop.py
or
    python3 tests/test_ns_main_loop.py
'''

import os
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from rapster.constants import tBH_form, fBH_max
from rapster.run_cluster import parse_args
from rapster.cluster_evolution import (initialize_cluster, form_neutron_stars, compute_cluster_properties,
                                       compute_timescales, form_binaries, evolve_interactions,
                                       compute_external_params, evolve_tdes, record_evolution, update_cluster)
from rapster.compact_accretion import accrete_gas
from rapster.stellar_evolution import init_stellar_mass_sampler
from rapster.functions import set_escape_factor, v_esc

init_stellar_mass_sampler(0.0)
set_escape_factor(2.0)

BASE_ARGS = ['-N', '2e5', '-S', '42', '-wK', '20', '-NS', '2', '-P', '0']


def setup(nst, extra=()):
    sys.argv = ['rapster'] + BASE_ARGS + ['-NSt', str(nst)] + list(extra)
    config = parse_args()
    return initialize_cluster(config), config


def run(state, config, t_max):
    """Main loop as in run_cluster.main, checking the BH bookkeeping at every step."""
    while state['t'] < t_max and state['R_gal'] > 0 and state['Mcl'] > 0 and state['BH'].m.sum() < fBH_max * state['Mcl']:
        form_neutron_stars(state, config)
        if not compute_cluster_properties(state, config):
            break
        compute_timescales(state, config); form_binaries(state, config); evolve_interactions(state, config)
        compute_external_params(state, config); evolve_tdes(state, config); record_evolution(state)
        accrete_gas(state, config)
        if state['i_aux1'] == 1:
            n_single = state['N_BH'] - 2 * state['N_BBH'] - state['N_BHstar'] - 3 * state['N_Triples']
            assert n_single == state['BH'].m.size
        assert state['BH'].m.size == state['BH'].s.size == state['BH'].g.size == state['BH'].h.size
        if not update_cluster(state, config):
            break
    return state


def test_form_neutron_stars_takes_only_formed_NSs():
    s, c = setup(2)
    n0, t = s['tNS_pending'].size, np.median(s['tNS_pending'])
    s['t'] = t
    form_neutron_stars(s, c)
    k = s['N_NS_formed']
    assert abs(k - n0 / 2) <= 1                       # half have t_form <= median
    assert np.all(s['tNS_pending'] > t)                # the rest are still pending
    assert s['mNS_pending'].size == s['vNS_pending'].size == n0 - k
    assert s['N_NS_retained'] == np.sum(s['BH'].m < 3.0)


def test_retention_uses_current_escape_velocity():
    s, c = setup(1)
    s['t'] = 1e5                                       # everything has formed
    v = s['vNS_pending'].copy()
    expected = np.sum(v < v_esc(s['Mcl'] + s['M_gas'], s['rh']))
    form_neutron_stars(s, c)
    assert s['N_NS_retained'] == expected and s['tNS_pending'].size == 0


def test_NSs_added_to_N_BH_only_after_activation():
    s, c = setup(2)
    s['t'] = 5.0                                       # before tBH_form: i_aux1 == 0
    form_neutron_stars(s, c)
    assert s['N_BH'] == 0
    s['i_aux1'] = 1; s['N_BH'] = s['BH'].m.size; n = s['N_BH']; r = s['N_NS_retained']
    s['t'] = 100.0
    form_neutron_stars(s, c)
    assert s['N_BH'] == n + (s['N_NS_retained'] - r)


def test_NSt0_has_nothing_to_form():
    s, c = setup(0)
    m = s['BH'].m.copy()
    s['t'] = 1e5
    form_neutron_stars(s, c)
    assert np.array_equal(s['BH'].m, m) and s['N_NS_formed'] == 0


def test_full_runs_keep_bookkeeping_and_form_all_NSs():
    for nst in (1, 2):
        s, c = setup(nst)
        n0 = s['tNS_pending'].size
        run(s, c, 300.0)
        assert s['N_NS_formed'] + s['tNS_pending'].size == n0
        assert s['N_NS_retained'] > 0


def test_BH_free_cluster_waits_for_NSs():
    empty = os.path.join(tempfile.mkdtemp(), 'empty.npz')
    np.savez(empty, mBH_ini=np.array([]))
    s, c = setup(2, extra=['-BIi', '1', '-BIF', empty])
    assert s['BH'].m.size == 0
    run(s, c, 200.0)
    assert s['t'] >= 200.0                             # did not stop with 'BH SUBSYSTEM EVAPORATED'
    # the retained NSs activated the subsystem and interacted (whether some are still
    # in the cluster at 200 Myr depends on the realization):
    assert s['N_NS_retained'] > 0 and s['N_3bb'] + s['N_2cap'] > 0


if __name__ == '__main__':
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    for test in tests:
        test()
        print(f'PASS  {test.__name__}')
    print(f'\nall {len(tests)} tests passed')
