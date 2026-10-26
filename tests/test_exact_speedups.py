'''
Checks that the performance helpers reproduce the original NumPy/SciPy calls bit for bit:
same returned values and same consumption of the global random stream (verified by
drawing one more number afterwards). If a future NumPy changes how
RandomState.choice works internally, these tests fail before any result changes.

Run from the repository root with
    python3 -m pytest tests/test_exact_speedups.py
'''

import io
import pickle
import warnings

import numpy as np
from scipy.optimize import root
from scipy.stats import poisson

from rapster import fastrandom as fr
from rapster import functions as F
from rapster.snapshots import SnapshotStream, dump_streaming


def _same(a, b):
    return np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True)


def test_choice_weighted_and_draw_index():
    rng = np.random.default_rng(1)
    for trial in range(300):
        n = int(rng.integers(1, 3000))
        a = rng.uniform(1, 100, n)
        if trial % 5 == 0:                       # repeated values
            a = np.round(a / 10) * 10 + 1
        w = a**rng.uniform(-1, 3)
        p = w / np.sum(w)
        for replace in (True, False):
            np.random.seed(trial)
            v_ref = np.random.choice(a, replace=replace, p=p)
            k_ref = int(np.atleast_1d(np.squeeze(np.where(a == v_ref)))[0])
            nxt_ref = np.random.rand()
            np.random.seed(trial)
            v, k = fr.choice_weighted(a, p)
            k = fr.first_index_of(a, k)
            nxt = np.random.rand()
            assert v == v_ref and k == k_ref and nxt == nxt_ref


def test_choice_pair_no_replace_including_repeated_draws():
    rng = np.random.default_rng(2)
    for trial in range(400):
        n = int(rng.integers(2, 400))
        w = rng.uniform(0, 1, n)**(1 + 30 * (trial % 2))   # peaked weights force repeated draws
        w[rng.integers(0, n)] += 1.0
        p = w / w.sum()
        a = rng.uniform(1, 50, n)
        cdf = fr.weights_cdf(p)
        np.random.seed(trial)
        ref = tuple(np.random.choice(a, size=2, replace=False, p=p)); nxt_ref = np.random.rand()
        np.random.seed(trial)
        got = fr.choice_pair_no_replace(a, p, cdf); nxt = np.random.rand()
        assert got == ref and nxt == nxt_ref


def test_permutation_equals_choice_without_replacement():
    for seed in range(50):
        n = 3 + seed * 37
        np.random.seed(seed); ref = np.random.choice(n, size=3, replace=False); nxt_ref = np.random.rand()
        np.random.seed(seed); got = np.random.permutation(n)[:3]; nxt = np.random.rand()
        assert _same(ref, got) and nxt == nxt_ref


def test_poisson_draw():
    for seed, mu in enumerate([0.0, 1e-9, 0.3, 2.5, 9.99, 10.0, 123.4, 5e4, 1e8]):
        np.random.seed(seed); ref = poisson.rvs(mu=mu); nxt_ref = np.random.rand()
        np.random.seed(seed); got = fr.poisson_draw(np.float64(mu)); nxt = np.random.rand()
        assert ref == got and nxt == nxt_ref


def test_array_helpers():
    rng = np.random.default_rng(3)
    a = rng.random(1000)
    for k in (0, 1, 500, 998, 999):
        assert _same(F.delete_at(a, k), np.delete(a, k))
        assert _same(F.delete_at(a, np.int64(k)), np.delete(a, k))
    for k1, k2 in ((0, 999), (999, 0), (10, 11), (500, 3)):
        assert _same(F.delete_at(a, [k1, k2]), np.delete(a, [k1, k2]))
    for v in (1.5, np.float64(2.5), 3, np.array([]), np.array([1.0, 2.0])):
        assert _same(F.append_to(a, v), np.append(a, v))
    b = rng.random((7, 15))
    for k in range(7):
        assert _same(F.delete_row(b, k), np.delete(b, k, axis=0))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for vals in ((1.0, 2.0), (2.0, 1.0), (3.0, 3.0), (np.nan, 1.0), (1.0, np.nan), (5.0, 1.0, 3.0), (1.0, np.inf, np.nan)):
            assert _same(F.min_of(*vals), np.min(list(vals)))
            assert _same(F.max_of(*vals), np.max(list(vals)))


def test_zlk_root_function():
    rng = np.random.default_rng(4)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(300):
            m0, m1, m2 = rng.uniform(5, 80, 3)
            a1 = 10**rng.uniform(-5, -2); a2 = a1 * 10**rng.uniform(0.5, 2)
            e1, e2 = np.sqrt(rng.random(2)); c1, c2 = rng.uniform(-1, 1, 2); w = rng.uniform(0, 2*np.pi)
            args = tuple(np.float64(x) for x in (m0, m1, m2, a1, a2, e1, e2, c1, c2)) + (w,)
            f = F._DW_of_x(*args)
            for x in (np.array([1e-10]), np.array([rng.random()])):
                assert _same(F.DW(x, *args), f(x))
            x_ref = root(fun=F.DW, x0=1e-10, args=args, tol=1e-10)['x'][0]
            x_new = root(fun=f, x0=1e-10, tol=1e-10)['x'][0]
            assert _same(x_ref, x_new)


def test_maxwell_table_and_lookback_memo():
    assert _same(np.vectorize(F.CDF_maxwell)(F.w_s), F.CDF_maxwell_s)
    for z in (3.0, np.float64(2.0), 3.0):
        got = F.lookback_interp(z)
        assert type(got) is np.ndarray and got.shape == () and got == F._lookback_interp(z)


def test_snapshot_stream_pickles_as_plain_list():
    rng = np.random.default_rng(5)
    arrays = [rng.random(int(rng.integers(0, 50))) for _ in range(30)]
    stream = SnapshotStream()
    for a in arrays:
        stream.append(a)
    assert len(stream) == 30 and _same(stream[-1], arrays[-1])
    buf = io.BytesIO()
    dump_streaming({'t': [0.1, 0.2], 'mBH': stream, 'hBH': []}, buf)
    data = pickle.loads(buf.getvalue())
    assert type(data['mBH']) is list and len(data['mBH']) == 30
    assert all(_same(a, b) and a.dtype == b.dtype for a, b in zip(arrays, data['mBH']))
    stream.close()



# ---------------------------------------------------------------- per-interaction helpers

def _ref_T_GW(m1, m2, a0, e0):
    Tc = 5 * F.c_light**5 * a0**4 / 256 / F.G_Newton**3 / m1 / m2 / (m1 + m2)
    factor_e = (1 + 0.27 * e0**10 + 0.33 * e0**20 + 0.2 * e0**1000) * (1 - e0**2)**(7/2)
    return Tc * factor_e


def test_T_GW_fast_path():
    rng = np.random.default_rng(20)
    for _ in range(20000):
        m1, m2 = 10**rng.uniform(0, 4, 2)
        a0 = 10**rng.uniform(-9, 0)
        e0 = rng.choice([0.0, rng.uniform(0, 1), 1 - 10**rng.uniform(-12, -1)])
        for args in ((np.float64(m1), np.float64(m2), np.float64(a0), np.float64(e0)), (m1, m2, a0, e0)):
            got, ref = F.T_GW(*args), _ref_T_GW(*args)
            assert got == ref and type(got) is np.float64
    x = np.array([1e-3, 2e-3])
    assert _same(F.T_GW(10.0, 20.0, x, 0.5), _ref_T_GW(10.0, 20.0, x, 0.5))   # arrays: NumPy path


def test_maxwell_sample_uses_same_interpolation():
    for seed in range(5):
        np.random.seed(seed)
        got = [F.get_maxwell_sample(np.float64(s)) for s in (0.3, 5.0, 40.0) * 200]
        nxt = np.random.rand()
        np.random.seed(seed)
        ref = [s * np.interp(np.random.rand(), F.CDF_maxwell_s, F.w_s) for s in (0.3, 5.0, 40.0) * 200]
        assert got == ref and np.random.rand() == nxt


def test_direct_minpack_root_matches_scipy_root():
    rng = np.random.default_rng(21)
    for _ in range(2000):
        m0, m1, m2 = 10**rng.uniform(0.5, 3, 3)
        a1 = 10**rng.uniform(-6, -2); a2 = a1 * 10**rng.uniform(0.5, 3)
        e1, e2 = rng.uniform(0, 0.99, 2); c1, c2 = rng.uniform(-1, 1, 2); w = rng.uniform(0, 2*np.pi)
        f = F._DW_of_x(m0, m1, m2, a1, a2, e1, e2, c1, c2, w)
        ref = root(fun=f, x0=1e-10, tol=1e-10)['x']
        got = F._root_hybr_x(f, 1e-10, 1e-10)
        assert np.array_equal(got.view(np.uint64), ref.view(np.uint64))


def test_ns_spin_helpers_match_original_helpers():
    from rapster import compact_accretion as ca
    rng = np.random.default_rng(22)
    eos = ca.get_eos('APR')
    for _ in range(20000):
        M = float(rng.uniform(1.0, eos.M_TOV)); chi = float(rng.uniform(0, 0.7))
        R = float(eos.radius(M))
        f, fK, chiK = ca._ns_spin_state(M, chi, R)
        assert f == ca.f_from_chi(M, chi, R=R, NS=True)
        assert fK == ca.f_kepler(M, NS=True, R=R) and chiK == ca.chi_kepler(M, NS=True, R=R)
        for pro in (True, False):
            for ff in (f, 0.0, float(rng.uniform(0, 2000))):
                assert ca._ns_E_J_isco(M, R, ff, pro) == ca._E_J_isco_NS(M, R, ff, prograde=pro)
        assert eos.radius(M) == np.interp(M, eos.M, eos.R)


def test_pool_total_equals_totals():
    from rapster.bhpool import BHPool, N_FAMILIES
    rng = np.random.default_rng(23)
    for n in (5, 300, 5000, 40000):
        m = rng.uniform(3, 60, n); m[0] = 3000.0
        pool = BHPool(m, 0 * m, 1 + 0 * m, 0 * m)
        for _ in range(50):
            pool.remove(int(rng.integers(0, pool.size)))
            pool.add(rng.uniform(3, 80), 0, 1, 0)
            tot = pool.totals()
            assert all(pool.total(f) == float(tot[f]) for f in range(N_FAMILIES))



def test_maxwell_draws_match_scipy():
    from scipy.stats import maxwell
    for seed in range(50):
        for scale in (np.sqrt(3) * 265.0, 34.6, 0.0):
            for size in (0, 1, 50):
                np.random.seed(seed); ref = maxwell.rvs(loc=0, scale=scale, size=size); nxt_ref = np.random.rand()
                np.random.seed(seed); got = fr.maxwell_draws(scale, size); nxt = np.random.rand()
                assert np.array_equal(got.view(np.uint64), ref.view(np.uint64)) and nxt == nxt_ref


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('passed', name)
