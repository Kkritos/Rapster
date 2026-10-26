'''
Tests for the fixes made after version 2.13.1:

1. three-body binaries: masses follow the intended mass-weighted distribution;
2. -SD 2 (beta spin distribution) is reproducible for a given seed;
3. the hardening file records the termination condition on the current row;
4. every TDE row is written to tdes.txt;
5. outputBHs.pkl stores the TDE counts (hBH) and no placeholder entries;
6. the default central density follows the Plummer value for the chosen N and r_h.

Run from the repository root with
    python3 -m pytest tests/test_fixes.py
'''

import itertools
import os
import pickle
import subprocess
import sys
import tempfile

import numpy as np

from rapster import three_body_binary as tb


def _run(workdir, *args):
    subprocess.run([sys.executable, '-m', 'rapster.run_cluster', '-P', '0', *args],
                   cwd=workdir, check=True, capture_output=True)
    return os.path.join(workdir, 'Results')


def test_three_body_masses_follow_target_distribution():
    # exact target over ordered triples of distinct BHs vs. sampled frequencies;
    # includes an IMBH-dominated case:
    for masses in (np.array([1.3, 5.0, 8.0, 12.0, 30.0]),
                   np.array([2.0e4, 10.0, 15.0, 20.0, 30.0, 1.4])):
        triples = list(itertools.permutations(range(masses.size), 3))
        p = np.array([(masses[a]*masses[b])**4 * masses[c]**2.5
                      / (masses[a] + masses[b])**0.5 / (masses[a] + masses[b] + masses[c])**0.5
                      for a, b, c in triples])
        p /= p.sum()
        np.random.seed(7)
        n = 40000
        counts = dict.fromkeys(triples, 0)
        for _ in range(n):
            m1, m2, m3, a, b = tb.fast_sample_3bodyBinary(masses, return_indices=True)
            c = int(np.flatnonzero(masses == m3)[0])
            assert len({a, b, c}) == 3 and masses[a] == m1 and masses[b] == m2
            counts[(a, b, c)] += 1
        for (a, b, c), prob in zip(triples, p):
            expected = n * prob
            if expected > 50:
                assert abs(counts[(a, b, c)] - expected) < 5 * np.sqrt(expected), (masses, a, b, c)


def test_beta_spins_reproducible():
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
        args = ('-N', '2e5', '-SD', '2', '-s', '1', '-S', '11', '-tM', '200')
        r1, r2 = _run(d1, *args), _run(d2, *args)
        for name in ('evolution.txt', 'mergers.txt', 'tdes.txt'):
            with open(os.path.join(r1, name)) as f1, open(os.path.join(r2, name)) as f2:
                assert f1.read() == f2.read(), name
        with open(os.path.join(r1, 'outputBHs.pkl'), 'rb') as f:
            spins = pickle.load(f)['sBH']
        assert any(s.size and np.any((s > 0) & (s < 1)) for s in spins)


def test_hardening_condition_on_current_row():
    with tempfile.TemporaryDirectory() as d:
        res = _run(d, '-N', '2e6', '-r', '0.5', '-S', '1', '-BOi', '0')
        h = np.loadtxt(os.path.join(res, 'hardening.txt'), ndmin=2)
    cond = h[:, 10]
    # each BBH evolution that ends in a merger, disruption, ejection or capture
    # (condition 2-6) writes the condition on the row of the last interaction and
    # on the closing row that follows it:
    ends = [i for i in np.flatnonzero(cond >= 2)
            if i + 1 == len(cond) or not (h[i+1, 4] == h[i, 4] and cond[i+1] == cond[i])]
    assert ends
    assert all(i > 0 and h[i-1, 4] == h[i, 4] and cond[i-1] == cond[i] for i in ends)


def test_all_tde_rows_written():
    from rapster.cluster_evolution import write_output
    rows = [[0.0]*18] + [[7, 1.0, 2.0, typ] + [1.0]*14 for typ in (1, 22, 22)]
    state = {'mergers': [[0.0]*27], 'evolution': [[0.0]*71], 'hardening': [[0.0]*12], 'tdes': rows,
             'N_tdeBHWD': 0, 'N_tdeBHstar': 1, 'N_tdeBBHstar': 1,   # one whole-binary event -> 2 rows
             'N_me': 0, 'N_iter': 1, 'N_hardening': 0, 'Mcl0': 1.0, 'rh0': 1.0, 'Z': 0.001,
             'zCl_form': 3.0, 'R_gal0': 8000.0, 'Mcl': 1.0, 'rh': 1.0, 'R_gal': 8000.0}
    config = {'results_folder_name': 'Results', 'BOi': 0, 'Ti': 1, 'Mi': 0, 'Ei': 0, 'Hi': 0,
              'tdes_file': 'tdes'}
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        os.chdir(d)
        try:
            write_output(state, config)
            tdes = np.loadtxt(os.path.join(d, 'Results', 'tdes.txt'), ndmin=2)
        finally:
            os.chdir(cwd)
    assert len(tdes) == 3 and list(tdes[:, 3]) == [1, 22, 22]


def test_bh_snapshots_have_tde_counts_and_no_placeholders():
    with tempfile.TemporaryDirectory() as d:
        res = _run(d, '-N', '1e6', '-S', '3')
        with open(os.path.join(res, 'outputBHs.pkl'), 'rb') as f:
            out = pickle.load(f)
    assert len(out['hBH']) == len(out['mBH']) == len(out['t']) > 0
    for m, s, g, h in zip(out['mBH'], out['sBH'], out['gBH'], out['hBH']):
        assert m.shape == s.shape == g.shape == h.shape
        assert np.all(m > 0) and np.all(g >= 1) and np.all(h >= 0)
    assert max(h.max() for h in out['hBH'] if h.size) >= 1


def test_default_density_is_plummer_value():
    with tempfile.TemporaryDirectory() as d:
        res = _run(d, '-N', '2e5', '-r', '0.5', '-tM', '10', '-BOi', '0')
        with open(os.path.join(res, 'log.txt')) as f:
            log = f.read()
    assert 'Initial central stellar density: 848000 pc^-3' in log   # 0.53 * 2e5 / 0.5^3


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('passed', name)
