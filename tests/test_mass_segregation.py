"""Mass segregation of the compact objects (-ms alpha, version 2.17.0): every single
compact object (BH or NS) takes part in core interactions in proportion to its core
occupation s(m) = m^alpha. The samplers are checked against their exact targets (the
targets of earlier versions times s(m) for every object drawn; chi-square tests with
fixed seeds, as in test_bhpool.py), the pool against a recomputation, and whole runs
(the option, restarts). alpha = 0 reproduces earlier versions bit for bit (checked on
whole runs by the comparison of 32 configurations with version 2.16.1, see CHANGELOG.md;
here: the weights and draws are those of the families without alpha)."""
import itertools
import json
import os
import pickle
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from rapster import bhpool
from rapster.bhpool import BHPool, FAM, set_occupation_power
from rapster import three_body_binary as tb
from rapster import binary_evolution as be
from rapster import two_body_capture as tc
from rapster import tidal_disruptions as td

from test_bhpool import _chi2_ok


@pytest.fixture(params=[0.6, 1.5])
def alpha(request):
    set_occupation_power(request.param)
    try:
        yield request.param
    finally:
        set_occupation_power(0.0)


def _pool(m):
    m = np.asarray(m, dtype=float)
    return BHPool(m, 0.1 * np.ones(m.size), np.ones(m.size), np.zeros(m.size))


def _populations():
    # BHs and NSs (1.3, 1.4 Msun), with and without an IMBH
    return [np.array([1.3, 5.0, 8.0, 12.0, 30.0, 9.0, 1.4]),
            np.array([2.0e4, 10.0, 15.0, 20.0, 30.0, 1.4, 25.0, 1.3]),
            np.array([60.0, 1.2, 1.3, 1.4, 1.5, 10, 11, 12, 13, 14, 15])]


def test_weights_and_pool_bookkeeping(alpha):
    rng = np.random.default_rng(1)
    m = np.concatenate([rng.uniform(1.1, 2.0, 600), rng.uniform(5, 40, 2400)])
    pool = _pool(m)
    pool.check()
    for name, base, k in (('w', 4.0, 1), ('wv', 6.5, 2), ('p15', 1.5, 1), ('c27', 2/7, 1),
                          ('c13', 1/3, 1), ('mass', 1.0, 0), ('m25', 2.5, 0), ('occ', 0.0, 1),
                          ('mocc', 1.0, 1)):
        assert bhpool._POWER_LIST[FAM[name]] == base + k * alpha
        w = pool._W[:pool.size, FAM[name]].copy()
        w_exact = m**(base + k * alpha)
        if name in ('w', 'wv'):
            w_exact[pool.D] = 0.0                       # (the heaviest BH is excluded)
        assert np.allclose(w, w_exact, rtol=1e-12, atol=0)
    for _ in range(400):                                # changes of the pool
        op = rng.integers(3)
        if op == 0:
            pool.remove(int(rng.integers(pool.n)))
        elif op == 1:
            pool.add(float(rng.uniform(1.1, 60)), 0.0, 1, 0)
        else:
            k = int(rng.integers(pool.n))
            pool.set_mass(k, float(pool.m[k] * rng.uniform(0.9, 1.3)))
    pool.check()


def test_occupation_draws(alpha):
    for masses in _populations():
        pool = _pool(masses)
        for name, p in (('occ', masses**alpha), ('mocc', masses**(1 + alpha))):
            np.random.seed(2)
            n = 60000
            counts = np.bincount([pool.draw(FAM[name]) for _ in range(n)], minlength=masses.size)
            ok, stat, dof = _chi2_ok(counts, p / p.sum(), n)
            assert ok, (name, masses, stat, dof)
        np.random.seed(2)
        counts = np.bincount([pool.draw_occupied() for _ in range(n)], minlength=masses.size)
        ok, stat, dof = _chi2_ok(counts, masses**alpha / np.sum(masses**alpha), n)
        assert ok, (masses, stat, dof)


def test_draws_without_the_heaviest_BH(alpha):
    """drawer_without_D (used where D would dominate the draws, e.g. an IMBH with
    s(m) = m^alpha) and draw_occupied_avoiding: the weights of the other BHs, also in
    pools with superblocks and after changes of the pool."""
    rng = np.random.default_rng(16)
    big = np.concatenate([rng.uniform(1.1, 40, 140000), [3e4]])
    for masses in _populations() + [big]:
        pool = _pool(masses)
        if masses.size > 1000:
            for _ in range(50):
                pool.set_mass(int(rng.integers(pool.n - 1)), float(rng.uniform(1.1, 40)))
            pool.remove(int(rng.integers(pool.n - 1)))
        m = pool.m.copy(); D = pool.D
        for name in ('c27', 'c127', 'occ'):
            w = m**bhpool._POWER_LIST[FAM[name]]; w[D] = 0.0
            dr = pool.drawer_without_D(FAM[name])
            assert np.isclose(dr.total, w.sum(), rtol=1e-12, atol=0)
            np.random.seed(17)
            n = 60000 if masses.size < 1000 else 200000
            counts = np.bincount([dr() for _ in range(n)], minlength=m.size)
            assert counts[D] == 0
            if masses.size > 1000:      # (binned by mass)
                edges = np.quantile(m[np.arange(m.size) != D], np.linspace(0, 1, 21)); edges[-1] += 1
                idx = np.clip(np.searchsorted(edges, m, side='right') - 1, 0, 19); idx[D] = 0
                counts = np.bincount(idx, weights=counts, minlength=20); w = np.bincount(idx, weights=w, minlength=20)
            ok, stat, dof = _chi2_ok(counts, w / w.sum(), n)
            assert ok, (name, masses.size, stat, dof)
        if masses.size < 1000:
            for excluded in ({D}, {D, 1 if D != 1 else 2}, (0,)):
                s = m**alpha; s[list(excluded)] = 0.0
                np.random.seed(18)
                counts = np.bincount([pool.draw_occupied_avoiding(excluded) for _ in range(40000)], minlength=m.size)
                assert not counts[list(excluded)].any()
                ok, stat, dof = _chi2_ok(counts, s / s.sum(), 40000)
                assert ok, (masses, excluded, stat, dof)


def test_draw_occupied_is_draw_uniform_without_segregation():
    pool = _pool(np.array([1.4, 5.0, 8.0, 30.0]))
    np.random.seed(3); a = [pool.draw_occupied() for _ in range(500)]
    np.random.seed(3); b = [pool.draw_uniform() for _ in range(500)]
    assert a == b


def _three_body_target(masses, alpha, remaining=None):
    idx = range(masses.size) if remaining is None else remaining
    triples = list(itertools.permutations(idx, 3))
    s = masses**alpha
    p = np.array([s[a] * s[b] * s[c] * (masses[a]*masses[b])**4 * masses[c]**2.5
                  / (masses[a] + masses[b])**0.5 / (masses[a] + masses[b] + masses[c])**0.5
                  for a, b, c in triples])
    return triples, p / p.sum()


def test_three_body_sampler_distribution(alpha):
    for masses in _populations():
        triples, p = _three_body_target(masses, alpha)
        index = {t: i for i, t in enumerate(triples)}
        n = 60000
        for sampler in ('pool', 'array'):
            pool = _pool(masses)
            np.random.seed(4)
            counts = np.zeros(len(triples))
            for _ in range(n):
                if sampler == 'pool':
                    m1, m2, m3, a, b = tb.sample_3bodyBinary(pool)
                else:
                    m1, m2, m3, a, b = tb.fast_sample_3bodyBinary(masses, return_indices=True)
                c = [k for k in range(masses.size) if masses[k] == m3 and k not in (a, b)][0]
                counts[index[(a, b, c)]] += 1
            ok, stat, dof = _chi2_ok(counts, p, n)
            assert ok, (sampler, masses, stat, dof)


def test_three_body_sampler_with_taken_BHs(alpha):
    masses = np.array([2.0e4, 10.0, 15.0, 20.0, 30.0, 1.4, 25.0, 12.0, 1.3])
    taken = {1, 4}
    pool = _pool(masses)
    remaining = [k for k in range(masses.size) if k not in taken]
    triples, p = _three_body_target(masses, alpha, remaining)
    index = {t: i for i, t in enumerate(triples)}
    np.random.seed(11)
    counts = np.zeros(len(triples))
    for _ in range(60000):
        sample = tb.sample_3bodyBinary(pool, exclude=taken, max_trials=30)
        if sample is None:
            continue
        m1, m2, m3, a, b = sample
        c = [k for k in remaining if masses[k] == m3 and k not in (a, b)][0]
        counts[index[(a, b, c)]] += 1
    ok, stat, dof = _chi2_ok(counts, p, int(counts.sum()))
    assert ok, (stat, dof)


def test_partner_sampler_distribution(alpha):
    for masses in _populations():
        pool = _pool(masses)
        for m12 in (3.0, 40.0, 3.0e3):
            p3 = (m12 + masses) / np.sqrt(m12**(-2/5) + masses**(-2/5)) * masses**(3/2) * masses**alpha
            np.random.seed(7)
            n = 40000
            counts = np.zeros(masses.size)
            for _ in range(n):
                m3, k3 = be.sample_BBH_partner(np.float64(m12), pool)
                counts[k3] += 1
            ok, stat, dof = _chi2_ok(counts, p3 / p3.sum(), n)
            assert ok, (masses, m12, stat, dof)


def test_tde_draw_distribution(alpha):
    for masses in _populations():
        pool = _pool(masses)
        for m_star in (0.6, 1e-2):
            p = (m_star + masses) * masses**(1/3) * masses**alpha
            np.random.seed(8)
            n = 40000
            counts = np.bincount([td.draw_TDE_BH(pool, m_star) for _ in range(n)], minlength=masses.size)
            ok, stat, dof = _chi2_ok(counts, p / p.sum(), n)
            assert ok, (masses, m_star, stat, dof)


def test_two_body_capture_sampler_distribution(alpha):
    # target: np.random.choice(size=2, replace=False) with p ~ m^(2/7) s(m), accepted
    # with probability ~ (m1 + m2)^(10/7); also with captured BHs excluded
    for masses, taken in ((np.array([1.3, 5.0, 8.0, 12.0, 30.0, 9.0, 1.4]), set()),
                          (np.array([2.0e4, 10.0, 15.0, 20.0, 30.0, 1.4, 25.0, 1.3]), set()),
                          (np.array([60.0, 1.2, 1.3, 1.4, 10, 11, 12, 13]), {3, 6})):
        pool = _pool(masses)
        remaining = [k for k in range(masses.size) if k not in taken]
        w = masses**(2/7 + alpha)
        p = w / w[remaining].sum()
        pairs = list(itertools.permutations(remaining, 2))
        target = np.array([p[i] * p[j] / (1 - p[i]) * (masses[i] + masses[j])**(10/7) for i, j in pairs])
        target /= target.sum()
        index = {t: k for k, t in enumerate(pairs)}
        np.random.seed(9)
        counts = np.zeros(len(pairs))
        for _ in range(60000):
            if taken:
                sample = tc.sample_2capture(pool, exclude=taken, w_excluded=float(w[list(taken)].sum()), max_trials=30)
                if sample is None:
                    continue
            else:
                sample = tc.sample_2capture(pool)
            m1, m2, i, j = sample
            counts[index[(i, j)]] += 1
        ok, stat, dof = _chi2_ok(counts, target, int(counts.sum()))
        assert ok, (masses, stat, dof)


def test_tdes_of_a_step(alpha):
    """The TDEs of a step (the pool updated at the end) against the pool updated after
    every TDE (both with the occupation weights): the same mean final masses."""
    m0 = np.array([1.3, 1.4, 1.5, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0])
    from rapster.arraycache import RowLog

    def run(reference, reps=600):
        out = []
        for r in range(reps):
            np.random.seed(100 + r)
            P = _pool(m0.copy())
            tdes = RowLog(18, [0.0] * 18)
            if reference:
                for _ in range(80):
                    k = td.draw_TDE_BH(P, 1.0)
                    m, s, g, h = P.get(k)
                    m_new, s_new, h_new = td.disrupt_star(1, 10.0, 0.5, 1, 0.6, 1.0, 1e-7, m, s, g, h, P.mean_mass(),
                                                          10.0, 8.0, tdes, 0.5, 'APR')
                    P.set_mass(k, m_new); P.set_spin(k, s_new); P.set_tde_count(k, h_new)
            else:
                td.BH_TidalDisruptions(1, 10.0, 0.5, 80, 0, 1, 0.6, 1.0, 1e-7, P, 10.0, 8.0, tdes, None, None, 0.5, 'APR')
            out.append(P.h.copy())
        return np.array(out)

    a = run(True); b = run(False)
    z = (b.mean(0) - a.mean(0)) / np.sqrt(a.var(0) / len(a) + b.var(0) / len(b) + 1e-300)
    assert np.abs(z).max() < 4.5, z


def test_pickled_pool_keeps_its_weights(alpha):
    rng = np.random.default_rng(12)
    pool = _pool(np.concatenate([rng.uniform(1.1, 2.0, 300), rng.uniform(5, 40, 900)]))
    for _ in range(100):
        pool.remove(int(rng.integers(pool.n))); pool.add(float(rng.uniform(1.1, 50)), 0.0, 1, 0)
    copy = pickle.loads(pickle.dumps(pool))
    copy.check()
    np.random.seed(13); x = [pool.draw(FAM['c13']) for _ in range(300)] + [pool.draw_occupied() for _ in range(300)]
    np.random.seed(13); y = [copy.draw(FAM['c13']) for _ in range(300)] + [copy.draw_occupied() for _ in range(300)]
    assert x == y


def test_pool_of_version_2161_without_the_new_families():
    """A pool pickled by version 2.15.1-2.16.1 (18 families, no m25, occ, mocc): the
    new families are computed from the masses."""
    rng = np.random.default_rng(14)
    pool = _pool(rng.uniform(3, 40, 2000))
    old = pool.__getstate__()
    keep = [FAM[f] for f in bhpool._FAMILY_NAMES_2151]
    old['_W'] = old['_W'][:, keep]
    old['_S'] = pool._S[:, keep]; old['_ST'] = pool._ST[keep, :]
    old['_families'] = bhpool._FAMILY_NAMES_2151
    for key in ('_SST', '_sbver', '_sgen', '_sb_dirty'):
        old.pop(key, None)
    q = BHPool.__new__(BHPool)
    q.__setstate__(old)
    q.check()
    assert np.array_equal(q._W, pool._W)
    np.random.seed(15); x = [pool.draw(FAM['c43']) for _ in range(200)]
    np.random.seed(15); y = [q.draw(FAM['c43']) for _ in range(200)]
    assert x == y


# ------------------------------------------------------------------- whole runs
HARNESS = '''
import os, sys
import rapster.cluster_evolution as ce
ce.STREAM_CHUNK_ROWS = 7
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


def _run(d, argv, kill=-1, check=True):
    r = subprocess.run([sys.executable, '-c', HARNESS.format(kill=kill, argv=list(argv))], cwd=d,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if check:
        assert r.returncode == 0, r.stderr
    return r


def test_runs_with_mass_segregation_and_restart():
    """A run with -ms 1.5: it states alpha in log.txt and the checkpoint; stopped after a
    checkpoint, restarted (alpha from the checkpoint) and joined, it gives the same files
    as without interruption; a restart cannot change alpha; alpha < 0 is refused."""
    args = ['-N', '3e5', '-r', '0.5', '-NS', '1', '-S', '5', '-tM', '3000', '-CKdt', '500', '-ms', '1.5']
    files = ('mergers.txt', 'evolution.txt', 'hardening.txt', 'tdes.txt')
    with tempfile.TemporaryDirectory() as d:
        _run(d, args + ['-RF', 'Ref'])
        assert 'alpha = 1.5' in open(os.path.join(d, 'Ref', 'log.txt')).read()
        with open(os.path.join(d, 'Ref', 'checkpoint.pkl'), 'rb') as f:
            assert pickle.load(f)['config']['mass_segregation_power'] == 1.5
        r = _run(d, args + ['-RF', 'Stopped'], kill=3, check=False)
        assert r.returncode == 3
        _run(d, ['-RS', 'Stopped', '-RF', 'Cont'])
        assert 'alpha = 1.5' in open(os.path.join(d, 'Cont', 'log.txt')).read()
        subprocess.run([sys.executable, '-m', 'rapster.join_restart', '-o', 'Joined', 'Stopped', 'Cont'],
                       cwd=d, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for name in files:
            assert open(os.path.join(d, 'Ref', name)).read() == open(os.path.join(d, 'Joined', name)).read(), name
        r = _run(d, ['-RS', 'Stopped', '-RF', 'Refused', '-ms', '0.5'], check=False)
        assert r.returncode != 0 and 'error' in r.stderr
        r = _run(d, ['-N', '1e4', '-ms', '-1', '-RF', 'Negative'], check=False)
        assert r.returncode != 0 and '-ms' in r.stderr
