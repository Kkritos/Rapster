'''
Tests of the pool of single BHs (rapster/bhpool.py) and of the samplers built on it.

The pool-based samplers draw from the same distributions as the array-based
samplers of versions <= 2.13.1, but not with the same random numbers, so they are
tested statistically: sampled frequencies are compared with the exact target
probabilities (chi-square test with a fixed seed; the thresholds correspond to
p-values ~1e-4, so a correct sampler fails only by an unlucky seed, which is fixed).
The pool's bookkeeping (weights, block sums, heaviest BH) is checked against a
recomputation after random sequences of additions, removals and mass changes.

Run from the repository root with
    python3 -m pytest tests/test_bhpool.py
'''

import itertools

import numpy as np
from scipy.stats import chi2

from rapster.bhpool import BHPool, FAM
from rapster import three_body_binary as tb
from rapster import binary_evolution as be
from rapster import two_body_capture as tc
from rapster import tidal_disruptions as td
from rapster import compact_accretion as ca


def _pool(m):
    m = np.asarray(m, dtype=float)
    return BHPool(m, 0.1 * np.ones(m.size), np.ones(m.size), np.zeros(m.size))


def _chi2_ok(counts, p, n, min_expected=20):
    """Chi-square test of counts against probabilities p (cells with small expectation pooled)."""
    counts = np.asarray(counts, dtype=float)
    expected = n * np.asarray(p, dtype=float)
    big = expected >= min_expected
    obs = list(counts[big]) + [counts[~big].sum()]
    exp = list(expected[big]) + [expected[~big].sum()]
    if exp[-1] < 5:
        obs, exp = obs[:-1], exp[:-1]
    obs = np.array(obs); exp = np.array(exp)
    stat = np.sum((obs - exp)**2 / exp)
    if len(obs) < 2:      # (one BH carries essentially all the weight)
        return abs(obs[0] - exp[0]) < 5 * np.sqrt(n - exp[0] + 1), stat, 0
    return chi2.sf(stat, len(obs) - 1) > 1e-4, stat, len(obs) - 1


def _populations():
    return [np.array([1.3, 5.0, 8.0, 12.0, 30.0, 9.0, 7.0]),
            np.array([2.0e4, 10.0, 15.0, 20.0, 30.0, 1.4, 25.0, 12.0]),   # IMBH
            np.array([60.0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22])]


def test_pool_bookkeeping_under_random_operations():
    rng = np.random.default_rng(1)
    m = rng.uniform(3, 40, 3000)
    pool = _pool(m)
    pool.check()
    for step in range(4000):
        op = rng.integers(0, 6)
        if op == 0 and pool.size > 3:
            pool.remove(int(rng.integers(0, pool.size)))
        elif op == 1:
            pool.add(rng.uniform(3, 80), 0.5, 2, 1)
        elif op == 2 and pool.size > 3:
            k = rng.choice(pool.size, 2, replace=False)
            pool.remove_two(int(k[0]), int(k[1]))
        elif op == 3:
            pool.set_mass(int(rng.integers(0, pool.size)), rng.uniform(3, 300))
        elif op == 4 and pool.size > 3:
            pool.remove(pool.D)                         # the heaviest BH leaves
        else:
            pool.replace(int(rng.integers(0, pool.size)), rng.uniform(3, 5000), 0.2, 3, 0)
        if step % 200 == 0:
            pool.check()
    pool.check()
    assert pool.m.size == pool.s.size == pool.g.size == pool.h.size == pool.size
    pool.set_all_masses(pool.m * 1.01)
    pool.check()


def test_pool_keeps_attributes_with_masses():
    pool = _pool([5.0, 6.0, 7.0, 8.0])
    pool.add(9.0, 0.9, 3, 2)
    pool.remove(0)                                      # the BH of slot 4 moves to slot 0
    assert pool.get(0) == (9.0, 0.9, 3.0, 2.0)
    assert sorted(pool.m) == [6.0, 7.0, 8.0, 9.0] and pool.D == 0


def test_weighted_draws():
    rng = np.random.default_rng(2)
    for m in _populations() + [rng.uniform(3, 40, 700)]:
        pool = _pool(m)
        for name in ('w', 'q', 'p15', 'c13', 'c127', 'mass'):
            f = FAM[name]
            w = pool._W[:pool.size, f]
            np.random.seed(3)
            n = 100000
            counts = np.bincount([pool.draw(f) for _ in range(n)], minlength=pool.size)
            assert counts[w == 0].sum() == 0
            ok, stat, dof = _chi2_ok(counts, w / w.sum(), n)
            assert ok, (name, m.size, stat, dof)


def test_three_body_sampler_distribution():
    for masses in _populations():
        pool = _pool(masses)
        triples = list(itertools.permutations(range(masses.size), 3))
        p = np.array([(masses[a]*masses[b])**4 * masses[c]**2.5
                      / (masses[a] + masses[b])**0.5 / (masses[a] + masses[b] + masses[c])**0.5
                      for a, b, c in triples])
        p /= p.sum()
        index = {t: i for i, t in enumerate(triples)}
        np.random.seed(4)
        n = 100000
        counts = np.zeros(len(triples))
        for _ in range(n):
            m1, m2, m3, a, b = tb.sample_3bodyBinary(pool)
            c = int(np.flatnonzero(masses == m3)[0])
            assert len({a, b, c}) == 3 and masses[a] == m1 and masses[b] == m2
            counts[index[(a, b, c)]] += 1
        ok, stat, dof = _chi2_ok(counts, p, n)
        assert ok, (masses, stat, dof)


def test_three_body_sampler_large_population_matches_exact_sampler():
    # marginal distributions of m1 + m2 and m3 against the exact array-based sampler:
    rng = np.random.default_rng(5)
    m = rng.uniform(3, 40, 5000); m[17] = 800.0
    pool = _pool(m)
    np.random.seed(6)
    a = np.array([tb.sample_3bodyBinary(pool)[:3] for _ in range(20000)])
    b = np.array([tb.fast_sample_3bodyBinary(m) for _ in range(20000)])
    from scipy.stats import ks_2samp
    for x, y in ((a[:, 0] + a[:, 1], b[:, 0] + b[:, 1]), (a[:, 2], b[:, 2])):
        assert ks_2samp(x, y).pvalue > 1e-4


def test_partner_sampler_distribution():
    for masses in _populations():
        pool = _pool(masses)
        for m12 in (3.0, 40.0, 3.0e3):
            p3 = (m12 + masses) / np.sqrt(m12**(-2/5) + masses**(-2/5)) * masses**(3/2)
            np.random.seed(7)
            n = 60000
            counts = np.zeros(masses.size)
            for _ in range(n):
                m3, k3 = be.sample_BBH_partner(np.float64(m12), pool)
                assert masses[k3] == m3
                counts[k3] += 1
            ok, stat, dof = _chi2_ok(counts, p3 / p3.sum(), n)
            assert ok, (masses, m12, stat, dof)


def test_tde_draw_distribution():
    for masses in _populations():
        pool = _pool(masses)
        for m_star in (0.6, 1e-2):
            p = (m_star + masses) * masses**(1/3)
            np.random.seed(8)
            n = 60000
            counts = np.bincount([td.draw_TDE_BH(pool, m_star) for _ in range(n)], minlength=masses.size)
            ok, stat, dof = _chi2_ok(counts, p / p.sum(), n)
            assert ok, (masses, m_star, stat, dof)


def test_two_body_capture_sampler_distribution():
    # target of the sampler of versions <= 2.13.1: np.random.choice(size=2, replace=False)
    # with p ~ m^(2/7), accepted with probability ~ (m1 + m2)^(10/7):
    for masses in _populations() + [np.array([3.0, 4.0]), np.array([1e5, 10.0, 20.0])]:
        pool = _pool(masses)
        p = masses**(2/7) / np.sum(masses**(2/7))
        pairs = list(itertools.permutations(range(masses.size), 2))
        target = np.array([p[i] * p[j] / (1 - p[i]) * (masses[i] + masses[j])**(10/7) for i, j in pairs])
        target /= target.sum()
        index = {t: k for k, t in enumerate(pairs)}
        np.random.seed(9)
        n = 60000
        counts = np.zeros(len(pairs))
        for _ in range(n):
            m1, m2, i, j = tc.sample_2capture(pool)
            assert i != j and masses[i] == m1 and masses[j] == m2
            counts[index[(i, j)]] += 1
        ok, stat, dof = _chi2_ok(counts, target, n)
        assert ok, (masses, stat, dof)


def test_bh_spin_after_accretion_matches_euler_loop():
    # bh_spin_after_accretion (RK4) against the Euler loop of evolve_spin_during_accretion,
    # for BHs growing by up to FAST_SPIN_MAX_GROWTH, prograde and retrograde, including
    # spins at 0 and at the Thorne limit:
    rng = np.random.default_rng(10)
    worst = 0.0
    for _ in range(1500):
        m = 10**rng.uniform(np.log10(3.0), np.log10(5000))
        s = rng.choice([0.0, rng.uniform(0, ca.THORNE_LIMIT), ca.THORNE_LIMIT, rng.uniform(0, 0.02)])
        dm = m * 10**rng.uniform(-6, np.log10(ca.FAST_SPIN_MAX_GROWTH))
        prograde = bool(rng.random() < 0.5)
        old = ca.evolve_spin_during_accretion(Mi=m, Mf=m+dm, NS=False, chi=s, dM=dm/100, prograde=prograde)['chi'][-1]
        new = ca.bh_spin_after_accretion(m, m + dm, s, dm/100, prograde)
        worst = max(worst, abs(new - old))
        assert 0.0 <= new <= ca.THORNE_LIMIT
    assert worst < 2e-6, worst


def test_spin_after_accretion_uses_euler_loop_otherwise():
    # NSs, and BHs growing by more than FAST_SPIN_MAX_GROWTH, use the original loop:
    for args in ((1.4, 1.5, True, 0.0), (5.0, 5.5, False, 0.3)):
        Mi, Mf, NS, chi = args
        ref = ca.evolve_spin_during_accretion(Mi=Mi, Mf=Mf, NS=NS, f=None, chi=chi, dM=(Mf-Mi)/100, eos='APR', prograde=True)['chi'][-1]
        assert ca.spin_after_accretion(Mi, Mf, NS, chi=chi, dM=(Mf-Mi)/100, eos='APR', prograde=True) == ref



def test_three_body_sampler_with_taken_BHs():
    # BHs taken by earlier binaries of the same step stay in the pool until the end of
    # the step and are excluded by rejection: the distribution over the remaining BHs
    # must be the exact target restricted to them (None = gave up, caller redraws)
    for masses, taken in ((np.array([1.3, 5.0, 8.0, 12.0, 30.0, 9.0, 7.0, 14.0, 6.0]), {2, 5}),
                          (np.array([2.0e4, 10.0, 15.0, 20.0, 30.0, 1.4, 25.0, 12.0, 18.0]), {1, 4, 8})):
        pool = _pool(masses)
        remaining = [k for k in range(masses.size) if k not in taken]
        triples = list(itertools.permutations(remaining, 3))
        p = np.array([(masses[a]*masses[b])**4 * masses[c]**2.5
                      / (masses[a] + masses[b])**0.5 / (masses[a] + masses[b] + masses[c])**0.5
                      for a, b, c in triples])
        p /= p.sum()
        index = {t: i for i, t in enumerate(triples)}
        np.random.seed(11)
        counts = np.zeros(len(triples))
        for _ in range(80000):
            sample = tb.sample_3bodyBinary(pool, exclude=taken, max_trials=30)
            if sample is None:
                continue
            m1, m2, m3, a, b = sample
            c = [k for k in remaining if masses[k] == m3][0]
            counts[index[(a, b, c)]] += 1
        ok, stat, dof = _chi2_ok(counts, p, int(counts.sum()))
        assert ok, (masses, stat, dof)


def test_zlk_smallest_root_matches_root_finder():
    # the dedicated solver gives the root that scipy.optimize.root (hybr, x0 = 1e-10)
    # converges to, and declines (None) when there is no positive root
    from scipy.optimize import root
    from rapster import functions as F
    rng = np.random.default_rng(12)
    n_fast = 0
    for _ in range(3000):
        m0, m1, m2 = 10**rng.uniform(0.5, 3, 3)
        a1 = 10**rng.uniform(-6, -2); a2 = a1 * 10**rng.uniform(0.5, 3)
        e1, e2 = rng.uniform(0, 0.99, 2); c1, c2 = rng.uniform(-1, 1, 2); w = rng.uniform(0, 2*np.pi)
        fun = F._DW_of_x(m0, m1, m2, a1, a2, e1, e2, c1, c2, w)
        x_fast = F._zlk_smallest_root(*fun.params)
        r = root(fun=fun, x0=1e-10, tol=1e-10)
        if x_fast is None:
            continue
        n_fast += 1
        assert r['success'] and abs(x_fast - r['x'][0]) <= 1e-10 * r['x'][0], (x_fast, r['x'][0])
    assert n_fast > 2000


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('passed', name)
