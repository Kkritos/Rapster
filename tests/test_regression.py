'''
Regression test: a default run with a fixed seed must reproduce stored results.

The run uses all default parameters (seed 1234567890) and takes a few seconds.
Counts must match exactly; floating-point quantities must match to a relative
tolerance of 1e-6, which absorbs round-off differences between NumPy versions
(results agree to ~10 significant digits across NumPy 1.26-2.5) but catches any
change in the physics.

If a change is MEANT to alter results (new physics, a bug fix), regenerate the
reference values with
    python3 tests/test_regression.py --update
paste the printed REFERENCE block below, and note the change in CHANGELOG.md.

Run from the repository root with either
    python3 -m pytest tests/test_regression.py
or
    python3 tests/test_regression.py
'''

import os
import subprocess
import sys
import tempfile

import numpy as np

REFERENCE = {
    'n_evolution_rows': 586,
    'n_mergers': 8,
    'n_tdes': 172,
    'N_BH_final': 691,
    'N_BHej_final': 300,
    't_final': 11526.43421710654,
    'M_cl_final': 209192.37607724348,
    'r_h_final': 9.831598127482364,
    'mBH_max_final': 40.661851175102434,
    'merger_remnant_mass_sum': 521.8577802355334,
}

COUNTS = ('n_evolution_rows', 'n_mergers', 'n_tdes', 'N_BH_final', 'N_BHej_final')
RTOL = 1e-6


def run_default(workdir):
    '''Run Rapster with default parameters in workdir and summarize the output.'''
    subprocess.run([sys.executable, '-m', 'rapster.run_cluster', '-P', '0'],
                   cwd=workdir, check=True, capture_output=True)
    res = os.path.join(workdir, 'Results')
    with open(os.path.join(res, 'evolution.txt')) as f:
        header = f.readline().lstrip('# ').split()
    evolution = np.loadtxt(os.path.join(res, 'evolution.txt'), ndmin=2)
    mergers = np.loadtxt(os.path.join(res, 'mergers.txt'), ndmin=2)
    tdes = np.loadtxt(os.path.join(res, 'tdes.txt'), ndmin=2)
    last = dict(zip(header, evolution[-1]))
    with open(os.path.join(res, 'mergers.txt')) as f:
        merger_cols = f.readline().lstrip('# ').split()
    return {
        'n_evolution_rows': len(evolution),
        'n_mergers': len(mergers),
        'n_tdes': len(tdes),
        'N_BH_final': int(last['N_BH']),
        'N_BHej_final': int(last['N_BHej']),
        't_final': last['t'],
        'M_cl_final': last['M_cl'],
        'r_h_final': last['r_h'],
        'mBH_max_final': last['mBH_max'],
        'merger_remnant_mass_sum': mergers[:, merger_cols.index('mRem')].sum(),
    }


def test_default_run_reproduces_reference():
    with tempfile.TemporaryDirectory() as workdir:
        result = run_default(workdir)
    mismatches = []
    for key, ref in REFERENCE.items():
        got = result[key]
        if key in COUNTS:
            ok = got == ref
        else:
            ok = np.isclose(got, ref, rtol=RTOL, atol=0)
        if not ok:
            mismatches.append(f'{key}: expected {ref}, got {got}')
    assert not mismatches, 'Results changed:\n  ' + '\n  '.join(mismatches)


if __name__ == '__main__':
    if '--update' in sys.argv:
        with tempfile.TemporaryDirectory() as workdir:
            result = run_default(workdir)
        print('REFERENCE = {')
        for key, val in result.items():
            val = int(val) if key in COUNTS else float(val)
            print(f"    '{key}': {val!r},")
        print('}')
    else:
        test_default_run_reproduces_reference()
        print('Regression test passed.')
