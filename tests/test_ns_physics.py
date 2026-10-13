'''
Physics checks of extended NS formation (-NSt 1/2) against analytic expectations.

1. Retention: in a full run, the number of retained NSs agrees with
       <N_ret> = sum_i P_Maxwell(v < v_esc(t_i)),
   where t_i is when NS i joins and v_esc(t) = v_esc(M_cl + M_gas, r_h) is
   taken from the run itself; the Maxwellian has 1D dispersion sqrt(3) wSN_kick.
2. Mass-time correlation: with -NSt 2 the sampled t_form of NSs in each COMPAS
   mass bin follows that bin's tabulated distribution, and with -NSt 1 it does not
   depend on mass.

Run from the repository root with either
    python3 -m pytest tests/test_ns_physics.py
or
    python3 tests/test_ns_physics.py
'''

import os
import sys

import numpy as np
from scipy.stats import maxwell, kstest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from rapster.constants import fBH_max
from rapster.run_cluster import parse_args
from rapster.cluster_evolution import (initialize_cluster, form_neutron_stars, compute_cluster_properties,
                                       compute_timescales, form_binaries, evolve_interactions,
                                       compute_external_params, evolve_tdes, record_evolution, update_cluster)
from rapster.compact_accretion import accrete_gas
from rapster.stellar_evolution import init_stellar_mass_sampler
from rapster.functions import set_escape_factor, v_esc
from rapster.ns_formation import NSFormationTable

init_stellar_mass_sampler(0.0)
set_escape_factor(2.0)


def run_with_expectation(args, t_max):
    """Run the main loop; return the state and the analytic mean and variance of N_ret."""
    sys.argv = ['rapster'] + list(args)
    config = parse_args()
    state = initialize_cluster(config)
    sigma = np.sqrt(3) * config['wSN_kick']
    mean, var, history = 0.0, 0.0, []
    while state['t'] < t_max and state['R_gal'] > 0 and state['Mcl'] > 0 and state['mBH'].sum() < fBH_max * state['Mcl']:
        # NSs about to join this step, and their retention probability:
        k = np.searchsorted(state['tNS_pending'], state['t'], side='right')
        if k > 0:
            p = maxwell.cdf(v_esc(state['Mcl'] + state['M_gas'], state['rh']), scale=sigma)
            mean += k * p
            var += k * p * (1 - p)
        form_neutron_stars(state, config)
        history.append((state['t'], state['N_NS_formed'], state['N_NS_retained'], mean, var))
        if not compute_cluster_properties(state, config):
            break
        compute_timescales(state, config); form_binaries(state, config); evolve_interactions(state, config)
        compute_external_params(state, config); evolve_tdes(state, config); record_evolution(state)
        accrete_gas(state, config)
        if not update_cluster(state, config):
            break
    return state, mean, var, np.array(history)


def test_retained_NSs_match_analytic_expectation():
    for nst in (1, 2):
        args = ['-N', '5e5', '-n', '2.65e5', '-S', '7', '-wK', '30', '-NS', '2', '-NSt', str(nst), '-P', '0']
        state, mean, var, _ = run_with_expectation(args, 2000.0)
        z = (state['N_NS_retained'] - mean) / np.sqrt(var)
        print(f'   NSt={nst}: retained {state["N_NS_retained"]}, expected {mean:.1f} +- {np.sqrt(var):.1f} (z = {z:+.2f})')
        assert abs(z) < 4


def test_late_NSs_are_retained_less():
    # the escape velocity drops as the cluster loses mass and expands, so the retention
    # probability of NSs joining late must be lower than that of NSs joining early:
    args = ['-N', '5e5', '-n', '2.65e5', '-S', '7', '-wK', '30', '-NS', '2', '-NSt', '1', '-P', '0']
    _, _, _, h = run_with_expectation(args, 2000.0)
    t, formed, retained = h[:, 0], h[:, 1], h[:, 2]
    i = np.searchsorted(t, 25.0)                       # ~ median COMPAS formation time
    f_early = retained[i] / formed[i]
    f_late = (retained[-1] - retained[i]) / (formed[-1] - formed[i])
    print(f'   retained fraction: {f_early:.3f} before 25 Myr, {f_late:.3f} after')
    assert f_late < f_early


def test_NSt2_times_follow_COMPAS_conditional_in_each_mass_bin():
    table = NSFormationTable()
    sys.argv = ['rapster', '-N', '2e7', '-S', '3', '-NS', '3', '-NSt', '2', '-P', '0']
    state = initialize_cluster(parse_args())
    m, t = state['mNS_pending'], state['tNS_pending']
    idx = table.mass_bin_index(m)
    n_tested = 0
    for j in np.unique(idx):
        sel = idx == j
        if sel.sum() < 500:
            continue
        cdf = lambda x, j=j: np.interp(np.log10(x), table.logt_conditional[j], table.q)
        assert kstest(t[sel], cdf).pvalue > 1e-4, j
        n_tested += 1
    assert n_tested >= 10


def test_NSt1_times_independent_of_mass():
    sys.argv = ['rapster', '-N', '2e7', '-S', '3', '-NS', '3', '-NSt', '1', '-P', '0']
    state = initialize_cluster(parse_args())
    m, t = state['mNS_pending'], state['tNS_pending']
    light, heavy = t[m < 1.4], t[m > 2.0]
    assert abs(np.median(light) / np.median(heavy) - 1) < 0.05


if __name__ == '__main__':
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    for test in tests:
        test()
        print(f'PASS  {test.__name__}')
    print(f'\nall {len(tests)} tests passed')
