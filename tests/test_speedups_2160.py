"""Optimizations of 2.16.0, which change the random numbers drawn but not the
distributions: the draws of large pools from their superblock sums, the TDEs of a step
(the pool updated at the end of the step), the GW captures of a step (captured BHs
excluded until the end of the step), and the background writer of the output files
(the same files)."""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from rapster.arraycache import RowLog
from rapster.bhpool import BHPool, FAM
import rapster.functions as fu
from rapster.tidal_disruptions import BH_TidalDisruptions, draw_TDE_BH, disrupt_star
from rapster.two_body_capture import two_body_capture, sample_2capture

fu.set_escape_factor(4.77)


def _changed_pool(n, seed):
    rng = np.random.default_rng(seed)
    P = BHPool(rng.uniform(5, 60, n), np.zeros(n), np.ones(n), np.zeros(n))
    for _ in range(2000):
        op = rng.integers(3)
        if op == 0:
            P.remove(int(rng.integers(P.n)))
        elif op == 1:
            P.add(float(rng.uniform(5, 90)), 0, 1, 0)
        else:
            k = int(rng.integers(P.n))
            P.set_mass(k, float(P.m[k] * rng.uniform(0.9, 1.3)))
    P.check()
    return P


@pytest.mark.parametrize('f', ['mass', 'p15', 'c13', 'w'])
def test_draws_follow_the_weights(f):
    """Slots drawn by draw and drawer, binned in 40 groups of equal probability: chi^2
    with 39 degrees of freedom."""
    P = _changed_pool(150000, 1)
    assert P.n > P.BLOCK * P.NB_SUPER                         # (superblocks)
    fi = FAM[f]
    w = P._W[:P.n, fi].copy()
    p = w / w.sum()
    edges = np.searchsorted(np.cumsum(p), np.linspace(0, 1, 41)[1:-1])
    expected = np.diff(np.concatenate([[0.0], [p[:e].sum() for e in edges], [1.0]]))
    np.random.seed(2)
    d = P.drawer(fi)
    for slots in (np.array([P.draw(fi) for _ in range(40000)]), np.array([d() for _ in range(40000)])):
        assert (w[slots] > 0).all()
        obs = np.bincount(np.searchsorted(edges, slots, side='right'), minlength=40)
        chi2 = ((obs - expected * len(slots))**2 / (expected * len(slots))).sum()
        assert chi2 < 80                       # (p ~ 1e-4 for 39 degrees of freedom)
    assert P.total(fi) == pytest.approx(w.sum(), rel=1e-12)


def test_superblock_sums_follow_large_pools():
    """Pools of more than NB_SUPER blocks: superblock sums and totals after additions,
    removals and mass changes, shrinking below NB_SUPER blocks (sums of 2.15.1) and
    growing above it again; pickling (also of a pool of version 2.15.x, without
    superblock sums)."""
    import pickle
    P = _changed_pool(140000, 3)
    assert P.n > P.BLOCK * P.NB_SUPER
    rng = np.random.default_rng(4)
    for step in range(3000):
        op = rng.integers(4)
        if op == 0:
            for _ in range(int(rng.integers(1, 30))):
                P.remove(int(rng.integers(P.n)))
        elif op == 1:
            P.add(float(rng.uniform(5, 90)), 0, 1, 0)
        elif op == 2:
            k = int(rng.integers(P.n)) if rng.random() < 0.8 else P.D
            P.set_mass(k, float(P.m[k] * rng.uniform(0.9, 1.3)))
        else:
            f = int(rng.integers(len(FAM)))
            P.draw(f); P.drawer(f)()
        if step % 300 == 0:
            P.check()
    P.check()
    assert P.n < P.BLOCK * P.NB_SUPER                         # (shrunk below NB_SUPER blocks)
    for _ in range(40000):
        P.add(float(rng.uniform(5, 90)), 0, 1, 0)
    P.check()
    Q = pickle.loads(pickle.dumps(P))
    Q.check()
    np.random.seed(3); x = [P.draw(FAM['c43']) for _ in range(200)]
    np.random.seed(3); y = [Q.draw(FAM['c43']) for _ in range(200)]
    assert x == y
    old = P.__getstate__()                                   # (a pool of 2.15.x)
    for key in ('_SST', '_sbver', '_sgen', '_sb_dirty', '_csb', '_wsb'):
        del old[key]
    R = BHPool.__new__(BHPool)
    R.__setstate__(old)
    R.check()
    np.random.seed(3); z = [R.draw(FAM['c43']) for _ in range(200)]
    assert z == x


M0 = np.array([3., 5., 8., 10., 12., 15., 20., 25., 30., 40., 50., 60., 80., 1.4, 1.8])


def _tdes_reference(P, k_tde, m_star, R_star, tdes):
    """The TDE loop of version 2.15.1: the pool updated after every TDE."""
    for _ in range(k_tde):
        k = draw_TDE_BH(P, m_star)
        m, s, g, h = P.get(k)
        m_new, s_new, h_new = disrupt_star(1, 10.0, 0.5, 1, 0.6, m_star, R_star, m, s, g, h, P.mean_mass(),
                                           10.0, 8.0, tdes, 0.5, 'APR')
        P.set_mass(k, m_new); P.set_spin(k, s_new); P.set_tde_count(k, h_new)


def _tde_outcomes(reference, reps):
    out = []
    for r in range(reps):
        np.random.seed(1000 + r)
        P = BHPool(M0.copy(), np.full(len(M0), 0.3), np.ones(len(M0)), np.zeros(len(M0)))
        tdes = RowLog(18, [0.0] * 18)
        for k_tde, m_star in ((60, 2.0), (40, 0.5)):
            if reference:
                _tdes_reference(P, k_tde, m_star, 1e-7, tdes)
            else:
                BH_TidalDisruptions(1, 10.0, 0.5, k_tde, 0, 1, 0.6, m_star, 1e-7, P, 10.0, 8.0, tdes,
                                    None, None, 0.5, 'APR')
        P.check()
        out.append(np.concatenate([P.m, P.s, P.h]))
    return np.array(out)


def test_tdes_of_a_step_have_the_same_distribution():
    """Final masses, spins and TDE counts of each BH after 100 TDEs (rich get richer):
    the same means as with the pool updated after every TDE."""
    a = _tde_outcomes(True, 1500)
    b = _tde_outcomes(False, 1500)
    z = (b.mean(0) - a.mean(0)) / np.sqrt(a.var(0) / len(a) + b.var(0) / len(b) + 1e-300)
    assert np.abs(z).max() < 4.5, z
    assert np.allclose(a[:, 2*len(M0):].sum(1), 100) and np.allclose(b[:, 2*len(M0):].sum(1), 100)


M2 = np.array([3., 5., 8., 10., 12., 15., 20., 25., 30., 40., 50., 60., 80., 150., 7., 9., 11., 13., 17., 19.])


def _capture_outcomes(reference, reps):
    out = []
    for r in range(reps):
        np.random.seed(5000 + r)
        P = BHPool(M2.copy(), np.zeros(len(M2)), np.arange(len(M2), dtype=float), np.zeros(len(M2)))
        taken = np.zeros(len(M2))
        if reference:          # (2.15.1: the captured BHs removed at once)
            for _ in range(6):
                m1, m2, k1, k2 = sample_2capture(P)
                taken[int(P.g[k1])] += 1; taken[int(P.g[k2])] += 1
                P.remove_two(k1, k2)
        else:
            mergers = RowLog(27, [0.0] * 27)
            res = two_body_capture(1, 10.0, 1.0, 0.5, 3.0, 6, 20.0, np.zeros((0, 15)), P, 30.0, 20.0,
                                   0, 40, 0, 0, 0, 0, mergers)
            for row in np.array(mergers)[1:]:
                taken[int(row[9])] += 1; taken[int(row[10])] += 1
            for row in res[7]:
                taken[int(row[8])] += 1; taken[int(row[9])] += 1
        out.append(taken)
    return np.array(out)


def test_captures_of_a_step_have_the_same_distribution():
    """Which BHs are captured in 6 GW captures among 20 BHs (the heaviest is often taken,
    which ends the exclusion): the same probabilities as with removal after each capture."""
    a = _capture_outcomes(True, 6000)
    b = _capture_outcomes(False, 6000)
    z = (b.mean(0) - a.mean(0)) / np.sqrt(a.var(0) / len(a) + b.var(0) / len(b) + 1e-300)
    assert np.abs(z).max() < 4.5, z
    assert (a.sum(1) == 12).all() and (b.sum(1) == 12).all()


HARNESS = '''
import os, sys
import rapster.cluster_evolution as ce
import rapster.arraycache as ac
ce.STREAM_CHUNK_ROWS = 50
ac.BACKGROUND_MIN_ROWS = 1
import rapster.run_cluster as rc
save = rc.save_checkpoint
count = [0]
def save_and_maybe_stop(state, config, version):
    save(state, config, version)
    count[0] += 1
    if count[0] == {kill}:
        os._exit(3)
rc.save_checkpoint = save_and_maybe_stop
sys.argv = ['rapster', '-P', '0'] + {argv!r}
rc.main()
'''


def _run(d, argv, kill=-1, env=None, check=True):
    r = subprocess.run([sys.executable, '-c', HARNESS.format(kill=kill, argv=list(argv))], cwd=d,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True,
                       env=dict(os.environ, **(env or {})))
    if check:
        assert r.returncode == 0, r.stderr
    return r


def test_background_writer_writes_the_same_files():
    """Rows written by the background writer (in chunks of 50 rows): the same files as
    written directly, also for a run stopped after a checkpoint, restarted and joined."""
    files = ('mergers.txt', 'evolution.txt', 'hardening.txt', 'tdes.txt', 'populations.txt')
    args = ['-N', '1e5', '-r', '0.6', '-S', '21', '-tM', '2500', '-CKdt', '300', '-POP', 'pops.txt']
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, 'pops.txt'), 'w') as f:
            f.write('400 4e4 0.002 0.7 0.5 3\n1300 2e4\n')
        _run(d, args + ['-RF', 'Direct'], env={'RAPSTER_BACKGROUND_WRITER': '0'})
        _run(d, args + ['-RF', 'Background'])
        r = _run(d, args + ['-RF', 'Stopped'], kill=4, check=False)
        assert r.returncode == 3
        _run(d, ['-RS', 'Stopped', '-RF', 'Cont'])
        subprocess.run([sys.executable, '-m', 'rapster.join_restart', '-o', 'Joined', 'Stopped', 'Cont'],
                       cwd=d, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for name in files:
            ref = open(os.path.join(d, 'Direct', name)).read()
            assert open(os.path.join(d, 'Background', name)).read() == ref, name
            assert open(os.path.join(d, 'Joined', name)).read() == ref, name
        with open(os.path.join(d, 'Direct', 'hardening.txt')) as f:
            assert sum(1 for _ in f) > 1000                      # (many chunks)
