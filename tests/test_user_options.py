"""Options of 2.15.0: -tge (gas expulsion time in Myr), -fCO (maximum mass fraction in
compact objects, counting all of them), -R 0 (cluster at the galactic center) and -GC
(continue at the galactic center after the inspiral)."""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from rapster.run_cluster import parse_args
from rapster.cluster_evolution import initialize_cluster, compact_object_mass
from rapster.constants import evolution_keys
from rapster.functions import set_escape_factor
from rapster.stellar_evolution import init_stellar_mass_sampler

init_stellar_mass_sampler(0.0)
set_escape_factor(2.0)


def _config(*args):
    argv = sys.argv
    sys.argv = ['rapster'] + list(args)
    try:
        return parse_args()
    finally:
        sys.argv = argv


def _run(d, *args):
    subprocess.run([sys.executable, '-m', 'rapster.run_cluster', '-P', '0', '-BOi', '0', '-RF', 'R'] + list(args),
                   cwd=d, check=True, stdout=subprocess.DEVNULL)
    res = os.path.join(d, 'R')
    log = open(os.path.join(res, 'log.txt')).read()
    ev = np.loadtxt(os.path.join(res, 'evolution.txt'), ndmin=2)
    return log, ev


K = {n: i for i, n in enumerate(evolution_keys)}


def test_gas_expulsion_time_in_myr():
    c = _config('-sfe', '0.5', '-tge', '3.5')
    assert c['t_ge'] == 3.5
    s = initialize_cluster(c)
    assert s['t_ge'] == 3.5 and s['M_gas0'] > 0
    assert _config()['t_ge'] == 10.0
    with pytest.raises(SystemExit):
        _config('-fge', '5')                  # removed in 2.15.0


def test_compact_object_mass_counts_all_compact_objects():
    s = initialize_cluster(_config('-N', '1e5'))
    single = float(s['BH'].m.sum())
    s['binaries'] = np.vstack([s['binaries'], np.r_[0, 3, 0.1, 0.5, 10.0, 20.0, np.zeros(9)]])
    s['triples'] = np.vstack([s['triples'], np.r_[0, 0, 0, 0, 5.0, 6.0, 7.0, np.zeros(17)]])
    s['pairs'] = np.vstack([s['pairs'], [0.1, 30.0, 0, 1, 0]])
    assert compact_object_mass(s) == pytest.approx(single + 30.0 + 18.0 + 30.0)


def test_fco_option_and_message():
    assert _config()['fCO_max'] == 0.2
    with tempfile.TemporaryDirectory() as d:          # small compact cluster: BH-dominated in Gyr
        log, ev = _run(d, '-N', '1e5', '-r', '0.316', '-S', '1')
        assert 'COMPACT OBJECTS REACHED A FRACTION fCO_max = 0.2 OF THE CLUSTER MASS' in log
        t_stop = ev[-1, K['t']]
    with tempfile.TemporaryDirectory() as d:          # a larger fraction allows a longer run
        log, ev = _run(d, '-N', '1e5', '-r', '0.316', '-S', '1', '-fCO', '0.5')
        assert ev[-1, K['t']] > t_stop


def test_cluster_at_galactic_center():
    with pytest.raises(SystemExit):
        _config('-R', '-1')
    with tempfile.TemporaryDirectory() as d:
        log, ev = _run(d, '-N', '1e6', '-R', '0', '-S', '2', '-tM', '3000')
        assert np.all(ev[:, K['R_gal']] == 0.0)
        assert 'GALAXY CENTER' not in log and 'MAXIMUM SIMULATION TIME REACHED' in log
        assert np.all(np.isfinite(ev[:, K['M_cl']]))


def test_continue_at_galactic_center():
    # a massive cluster close to the center spirals in by dynamical friction within ~1 Gyr
    args = ('-N', '1e7', '-R', '300', '-S', '3', '-tM', '3000')
    with tempfile.TemporaryDirectory() as d:
        log, ev = _run(d, *args)
        assert 'CLUSTER REACHED GALAXY CENTER' in log
        t_center = ev[-1, K['t']]
        assert t_center < 2900
    with tempfile.TemporaryDirectory() as d:
        log, ev2 = _run(d, *args, '-GC', '1')
        assert 'CLUSTER REACHED GALAXY CENTER' in log and 'CONTINUING AT THE GALACTIC CENTER' in log
        assert ev2[-1, K['t']] > t_center
        R = ev2[:, K['R_gal']]
        first = int(np.argmax(R == 0.0))                # first step at the center
        assert first > 0 and np.all(R[first:] == 0.0) and np.all(R[:first] > 0)
        assert len(R) - first > 10                      # and the cluster keeps evolving there
