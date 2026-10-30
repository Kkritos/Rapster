"""Exact optimizations of 2.15.1: each faster routine gives the same values, and
consumes the same random numbers, as the code it replaces."""
import pickle

import numpy as np
import pytest

from rapster.bhpool import BHPool, FAM, N_FAMILIES, _FAMILY_NAMES, _FAMILY_NAMES_2150


def _pool(n, seed=1, imbh=False):
    rng = np.random.default_rng(seed)
    m = rng.uniform(3, 80, n)
    if imbh:
        m[n // 3] = 3e4
    return BHPool(m, rng.uniform(0, 1, n), np.ones(n), np.zeros(n)), rng


def _check_cdfs(P, fams=range(N_FAMILIES)):
    nb = (P.n + P.BLOCK - 1) // P.BLOCK
    for f in fams:
        assert np.array_equal(P._block_cdf(f), np.cumsum(P._ST[f, :nb]))


@pytest.mark.parametrize('n', [5, 300, 5000, 40000])
def test_batched_cumulative_sums_equal_cumsum(n):
    P, rng = _pool(n, imbh=n > 100)
    for step in range(400):
        op = rng.integers(4)
        if op == 0 and P.n > 3:
            P.remove(int(rng.integers(P.n)))
        elif op == 1:
            P.add(float(rng.uniform(3, 90)), 0.1, 1, 0)
        elif op == 2:
            k = int(rng.integers(P.n)) if rng.random() < 0.7 else P.D
            P.set_mass(k, float(P.m[k] * rng.uniform(0.9, 1.2)))
        else:
            f0, k = [(FAM['c13'], 3), (FAM['c27'], 2), (FAM['p15'], 4), (FAM['w'], 5), (FAM['mass'], 1)][int(rng.integers(5))]
            P.prepare(f0, k)
            _check_cdfs(P, range(f0, f0 + k))
        if step % 50 == 0:
            P.check()
    P.check()


def test_heavier_D_stays_D():
    P, rng = _pool(3000, imbh=True)
    for _ in range(200):
        P.set_mass(P.D, float(P.m[P.D]) + rng.uniform(0, 1))
        assert P.D == int(np.argmax(P.m))
    P.set_mass(P.D, 1.0)                  # D becomes lighter: another BH is the heaviest
    assert P.D == int(np.argmax(P.m))
    P.check()


@pytest.mark.parametrize('f', ['w', 'v', 'wC', 'q', 'r', 'p15', 'c13', 'mass'])
def test_drawer_matches_draw(f):
    P, rng = _pool(20000, imbh=True)
    fi = FAM[f]
    np.random.seed(5)
    a = [P.draw(fi) for _ in range(3000)]
    s_a = np.random.random()
    dr = P.drawer(fi)
    np.random.seed(5)
    b = [dr() for _ in range(3000)]
    s_b = np.random.random()
    assert a == b and s_a == s_b                  # same slots, same random numbers used
    P.add(10.0, 0.0, 1, 0)                        # a change of the pool: new drawer
    assert P.drawer(fi) is not dr


def test_pool_pickle_from_version_2150_order():
    """A pool pickled by 2.15.0 (other order of the weight families) is restored with
    the columns in this version's order."""
    P, rng = _pool(3000, imbh=True)
    state = P.__getstate__()
    perm = [_FAMILY_NAMES.index(f) for f in _FAMILY_NAMES_2150]
    old = dict(state)
    old['_W'] = state['_W'][:, perm]; old['_S'] = P._S[:, perm]; old['_ST'] = P._ST[perm, :]
    del old['_families']; old.pop('_CB'); old['_cb'] = {}
    Q = BHPool.__new__(BHPool)
    Q.__setstate__(old)
    assert np.array_equal(Q._W, P._W) and np.array_equal(Q._S, P._S) and np.array_equal(Q._ST, P._ST)
    Q.check()
    np.random.seed(3); x = [P.draw(FAM['c43']) for _ in range(100)]
    np.random.seed(3); y = [Q.draw(FAM['c43']) for _ in range(100)]
    assert x == y
    R = pickle.loads(pickle.dumps(P))
    R.check()


def test_remove_many_equals_removes():
    P, rng = _pool(5000, imbh=True)
    Q = pickle.loads(pickle.dumps(P))
    for _ in range(30):
        ks = sorted(set(rng.integers(0, P.n, 40).tolist()) | {P.D} if rng.random() < 0.3 else set(rng.integers(0, P.n, 40).tolist()), reverse=True)
        for k in ks:
            P.remove(k)
        Q.remove_many(ks)
        assert P.n == Q.n and P.D == Q.D
        for name in ('_m', '_s', '_g', '_h', '_W', '_S', '_ST'):
            P.totals(); Q.totals()
            assert np.array_equal(getattr(P, name), getattr(Q, name)), name
    Q.check()


def test_scalar_interpolation_equals_interp1d():
    """redshift_interp and lookback_interp of a scalar: the same 0-d array as the scipy
    interpolators, also at the grid points and the ends of the table."""
    import rapster.functions as fu
    rng = np.random.default_rng(8)
    for interp, ref, (lo, hi), table in ((fu.redshift_interp, fu._redshift_interp, (fu._T_MIN, fu._T_MAX), fu._redshift_np),
                                        (fu.lookback_interp, fu._lookback_interp, (fu._Z_MIN, fu._Z_MAX), fu._lookback_np)):
        xs = list(rng.uniform(lo, hi, 3000)) + [lo, hi] + (list(table[0][::97]) if table is not None else [])
        for x in xs:
            for v in (float(x), np.float64(x)):
                fu._lookback_memo.clear()
                a, b = interp(v), ref(v)
                assert isinstance(a, np.ndarray) and a.shape == b.shape == () and a.tobytes() == b.tobytes()
    with pytest.raises(ValueError):
        fu.redshift_interp(fu._T_MAX * 2)


def test_merger_remnant_equals_the_three_functions():
    """merger_remnant (one pass on Python floats) gives the values and types of
    remnant_mass, remnant_spin and remnant_kick, with the same random number."""
    import rapster.remnant as R
    rng = np.random.default_rng(9)
    for i in range(20000):
        m1 = np.float64(rng.uniform(3, 300)); m2 = np.float64(rng.uniform(3, 300)) if i % 7 else m1
        c1 = np.float64(rng.choice([0.0, 1.0, rng.random()])); c2 = np.float64(rng.choice([0.0, 0.998, rng.random()]))
        t1 = np.arccos(np.float64(rng.choice([1.0, -1.0, rng.uniform(-1, 1)]))); t2 = np.arccos(np.float64(rng.uniform(-1, 1)))
        dp = np.float64(rng.uniform(0, 2 * np.pi))
        if i % 3 == 0:
            m1 = float(m1); c1 = float(c1)
        np.random.seed(i)
        ref = (R.remnant_mass(m1, m2, c1, c2, t1, t2, dp), R.remnant_spin(m1, m2, c1, c2, t1, t2, dp),
               R.remnant_kick(m1, m2, c1, c2, t1, t2, dp))
        u_ref = np.random.rand()
        np.random.seed(i)
        got = R.merger_remnant(m1, m2, c1, c2, t1, t2, dp)
        assert np.random.rand() == u_ref
        for x, y in zip(ref, got):
            assert type(x) is type(y) and np.float64(x).tobytes() == np.float64(y).tobytes()


def test_delete_row_inplace_equals_np_delete():
    from rapster.functions import delete_row_inplace
    rng = np.random.default_rng(10)
    a = rng.random((700, 15)); ref = a.copy()
    for _ in range(650):
        k = int(rng.integers(0, a.shape[0]))
        a = delete_row_inplace(a, k); ref = np.delete(ref, k, axis=0)
        assert np.array_equal(a, ref)
    v = rng.random(50)
    assert np.array_equal(delete_row_inplace(v.copy(), 7), np.delete(v, 7))


def test_last_mass_of_the_euler_steps():
    """_last_mass_below (the last mass below Mf of the Euler steps, computed directly
    within a binade) equals the loop of additions, also for steps halfway between two
    multiples of the spacing (ties), near powers of 2 and for many steps."""
    import math
    from rapster.compact_accretion import _last_mass_below

    def loop(Mi, Mf, step):
        M = Mi; M_end = Mi
        while M < Mf:
            M_end = M; M += step
        return M_end

    rng = np.random.default_rng(11)
    cases = [(1.0, 1.0 + 300 * 2.5 * 2**-52, 2.5 * 2**-52), (1.0, 1.0 + 1e-13, 3.5 * 2**-52),
             (31.9, 32.3, 0.004), (63.99, 64.0, 1e-4), (10.0, 10.25, 0.0025)]
    for _ in range(30000):
        Mi = float(rng.choice([rng.uniform(1.1, 3), rng.uniform(3, 300), 2.0**int(rng.integers(1, 10)) * (1 - 1e-6)]))
        dm = float(rng.uniform(0.04, 50)) if rng.random() < 0.7 else float(rng.uniform(1e-6, 0.02)) * Mi
        cases.append((Mi, Mi + dm, dm / 100))
    for Mi, Mf, step in cases:
        assert _last_mass_below(Mi, Mf, step) == loop(Mi, Mf, step)
