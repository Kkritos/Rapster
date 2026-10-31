"""-WT 0 (--with_tdes 0) switches off every TDE channel, including the disruptions of
stars by BBHs (tde_type 21: one component, 22: whole binary), which up to version
2.16.0 still happened (and changed the BH masses, spins and TDE counts)."""
import os
import subprocess
import sys
import tempfile

import numpy as np

from rapster.constants import evolution_keys
from rapster.tidal_disruptions import try_BBH_star_disruption

# a cluster with a 300 Msun seed BH in which BBH-star disruptions occur (58 of them with
# -WT 0 in version 2.16.0, one with -WT 1), in ~1.5 s:
ARGS = ['-N', '1e6', '-r', '0.8', '-mm', '0.08', '-mM', '200', '-Z', '0.001', '-z', '2', '-fb', '0.1',
        '-S', '4', '-s', '0.2', '-RP', '1', '-NS', '2', '-NSt', '2', '-fA', '0.3', '-mb', '2',
        '-MBH', '300', '-sBH', '0', '-fesc', '3.16227766', '-cxi', '2',
        '-Hi', '0', '-BOi', '0', '-P', '0', '-CKi', '0']


def _run(d, wt):
    folder = 'WT' + str(wt)
    subprocess.run([sys.executable, '-m', 'rapster.run_cluster'] + ARGS + ['-WT', str(wt), '-RF', folder],
                   cwd=d, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    res = os.path.join(d, folder)
    types = [float(line.split()[3]) for line in open(os.path.join(res, 'tdes.txt')) if not line.startswith('#')]
    ev = np.loadtxt(os.path.join(res, 'evolution.txt'), ndmin=2)
    mergers = np.loadtxt(os.path.join(res, 'mergers.txt'), ndmin=2)
    return types, ev, mergers


def test_WT0_has_no_TDEs_of_any_type():
    K = {n: i for i, n in enumerate(evolution_keys)}
    with tempfile.TemporaryDirectory() as d:
        types, ev, mergers = _run(d, 1)
        assert 21.0 in types                                  # (the BBH-star channel is active with -WT 1)
        types, ev, mergers = _run(d, 0)
        assert types == []                                    # no TDE rows of any type
        for key in ('N_tdeBBHstar', 'N_tdeBHstar', 'N_tdeBHWD'):
            assert not ev[:, K[key]].any(), key
        assert len(mergers) and not mergers[:, 25:27].any()   # no merging BH ever disrupted a star (h1, h2)


def test_try_BBH_star_disruption_with_tdes_0():
    """Passages inside the tidal radius disrupt the star with with_tdes=1 (one component,
    type 21, for rp < a; the whole binary, type 22, for a < rp < r_t), and do nothing with
    with_tdes=0 (inputs returned unchanged, no row, no random number used)."""
    common = dict(m1=30.0, m2=20.0, s1=0.1, s2=0.2, g1=1.0, g2=1.0, h1=0.0, h2=0.0, m_star=0.5,
                  R_star=1e-8, f_accreted=0.3, EoS='APR', seed=1, t=10.0, z=1.0, m_avg=0.6, vBH=10.0, v_star=12.0)
    for rp, a, tde_type in ((1e-10, 1e-4, 21.0), (2e-8, 1e-9, 22.0)):     # (r_t = 4.6e-8 pc)
        tdes = []
        np.random.seed(3)
        out = try_BBH_star_disruption(rp, a, tdes=tdes, with_tdes=1, **common)
        assert out[0] and tdes and all(row[3] == tde_type for row in tdes)
        tdes = []
        np.random.seed(3)
        out = try_BBH_star_disruption(rp, a, tdes=tdes, with_tdes=0, **common)
        u = np.random.rand()
        np.random.seed(3)
        assert out == (False, 30.0, 20.0, 0.1, 0.2, 0.0, 0.0, []) and tdes == [] and u == np.random.rand()
