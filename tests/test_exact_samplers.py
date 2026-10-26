'''
Exactness of the array-based three-body sampler and its helpers.

fast_sample_3bodyBinary (used by sample_3bodyBinary when very few BHs carry
almost all the weight) uses cached powers of the BH masses (ElementwiseCache) and
locates random numbers in block sums instead of sequential cumulative sums
(CumulativeSearch). It must return exactly what the straightforward implementation
returns, with the same consumption of the global random stream. These tests
compare it against a reference implementation on BH populations large enough
(> CumulativeSearch.MIN_SIZE) to use the fast paths, including populations that
change between calls and an IMBH-dominated one, and probe random numbers placed
exactly on, and next to, the boundaries of the cumulative distributions.
(The pool-based samplers are tested statistically in tests/test_bhpool.py.)

Run from the repository root with
    python3 -m pytest tests/test_exact_samplers.py
'''

import numpy as np

from rapster import fastrandom as fr
from rapster import three_body_binary as tb
from rapster.arraycache import ElementwiseCache


# ---------------------------------------------------------------- reference code

def _ref_sums_of_others(w):
    total = np.sum(w)
    others = total - w
    j = int(np.argmax(w))
    if w[j] > 0.5 * total:
        others[j] = np.sum(w[:j]) + np.sum(w[j+1:])
    return others


def _ref_search(cdf, weights):
    j = int(np.searchsorted(cdf, np.random.rand() * cdf[-1], side='right'))
    if j >= cdf.size:
        j = int(np.flatnonzero(weights > 0)[-1])
    return j


def _ref_draw_excluding(weights, cdf, excluded, accept=None, efficiency=1.0):
    if efficiency > 0.25:
        while True:
            j = _ref_search(cdf, weights)
            if j not in excluded and (accept is None or np.random.rand() < accept(j)):
                return j
    w = weights.copy()
    w[list(excluded)] = 0.0
    if accept is not None:
        w *= np.clip(accept(np.arange(w.size)), 0.0, 1.0)
    return _ref_search(np.cumsum(w), w)


def _ref_3bb(mBH):
    m = np.asarray(mBH, dtype=float)
    D = int(np.argmax(m)); mD = m[D]; lm = np.log(m)

    def weights(power):
        w = np.exp(power * lm); w[D] = 0.0
        return w
    w = weights(4.0); v = weights(2.5); wC = weights(15/4); q = weights(43/12); r = weights(7/3)
    V = np.sum(v); wA = w * (V - v); Z_A = 2 * mD**4 / mD * np.sum(wA)
    wC_first = wC * _ref_sums_of_others(wC); Z_C = 4**(-0.25) * mD**2 * np.sum(wC_first)
    r_others = _ref_sums_of_others(r)
    wB_first = np.maximum(q * (r_others*_ref_sums_of_others(q) - _ref_sums_of_others(q*r)), 0.0)
    Z_B = np.sum(wB_first) / np.sqrt(6)
    cdf = {}

    def cdf_of(name, arr):
        if name not in cdf:
            cdf[name] = np.cumsum(arr)
        return cdf[name]
    Z = np.array([Z_A, Z_C, Z_B]); Z_cdf = np.cumsum(Z); Q = np.sum(q)
    while True:
        region = _ref_search(Z_cdf, Z)
        if region == 0:
            x = _ref_search(cdf_of('wA', wA), wA)
            c = _ref_draw_excluding(v, cdf_of('v', v), {x}, efficiency=(V - v[x]) / V)
            a, b = (D, x) if np.random.rand() < 0.5 else (x, D)
            bound = 1 / mD
        elif region == 1:
            a = _ref_search(cdf_of('wC_first', wC_first), wC_first)
            WC = cdf_of('wC', wC)[-1]
            b = _ref_draw_excluding(wC, cdf_of('wC', wC), {a}, efficiency=(WC - wC[a]) / WC)
            c = D
            bound = (4 * m[a] * m[b])**(-0.25) * mD**(-0.5)
        else:
            a = _ref_search(cdf_of('wB_first', wB_first), wB_first)
            Ra = r_others[a]
            b = _ref_draw_excluding(q, cdf_of('q', q), {a}, accept=lambda j: (Ra - r[j]) / Ra,
                                    efficiency=wB_first[a] / (q[a] * Q * Ra))
            c = _ref_draw_excluding(r, cdf_of('r', r), {a, b}, efficiency=(Ra - r[b]) / (Ra + r[a]))
            bound = (m[a]*m[b])**(-5/12) * m[c]**(-1/6) / np.sqrt(6)
        m1, m2, m3 = mBH[a], mBH[b], mBH[c]
        g = (m1 + m2)**(-0.5) * (m1 + m2 + m3)**(-0.5)
        if np.random.rand() * bound < g:
            return (m1, m2, m3, a, b)


# ---------------------------------------------------------------- helpers

def _population(rng, n, imbh=False):
    m = rng.uniform(3.0, 40.0, n) + 0.01 * rng.random(n)
    m[rng.integers(0, n, n // 20)] *= rng.uniform(1, 5, n // 20)   # some hierarchical remnants
    if imbh:
        m[rng.integers(0, n)] = 5000.0
    return m


def _evolve(rng, m):
    """Change a BH population as the simulation does: removals, appends, in-place changes."""
    op = rng.integers(0, 4)
    if op == 0:
        k = rng.choice(m.size, 2, replace=False)
        return np.delete(m, k)
    if op == 1:
        return np.append(np.delete(m, rng.integers(0, m.size)), rng.uniform(3, 60))
    if op == 2:
        return np.append(m, rng.uniform(3, 60, rng.integers(1, 3)))
    m = m.copy()
    m[rng.integers(0, m.size)] += 0.25
    return m


# ---------------------------------------------------------------- tests

def test_cumulative_search_matches_exact_cumsum_at_boundaries():
    rng = np.random.default_rng(10)
    for trial in range(40):
        n = int(rng.integers(fr.CumulativeSearch.MIN_SIZE, 60000))
        kind = trial % 4
        if kind == 0:
            w = rng.random(n)
        elif kind == 1:
            w = rng.random(n)**8 * (rng.random(n) < 0.7)        # many zeros
        elif kind == 2:
            w = rng.random(n); w[rng.integers(0, n)] = 1e6 * w.sum()   # one dominant weight
        else:
            w = 10**rng.uniform(-8, 8, n)                      # wide dynamic range
        search = fr.CumulativeSearch(w)
        exact_A = fr.weights_cdf(w / np.sum(w))               # np.random.choice (normalized p)
        c = np.cumsum(w)                                      # unnormalized (3bb draws)
        us = list(rng.random(300))
        for k in rng.integers(0, n, 100):
            v = exact_A[k]
            us += [v, np.nextafter(v, 0), np.nextafter(v, 2), v - 1e-12, v + 1e-12, v - 1e-9, v + 1e-9]
        for u in us:
            if not (0 <= u < 1):
                continue
            k = search.index(u)
            if k is None:
                continue
            kA = int(exact_A.searchsorted(u, side='right'))
            kB = int(np.searchsorted(c, u * c[-1], side='right'))
            # normalized probabilities as np.random.choice builds them from w/sum(w):
            assert k == kA, (trial, u, k, kA)
            if kB < n:
                assert k == kB, (trial, u, k, kB)


def test_weighted_index_matches_reference():
    rng = np.random.default_rng(11)
    for trial in range(30):
        n = int(rng.integers(100, 30000))
        w = rng.random(n)**int(rng.integers(1, 6))
        seed = int(rng.integers(0, 2**31))
        np.random.seed(seed)
        ref = [_ref_search(np.cumsum(w), w) for _ in range(40)]
        nxt_ref = np.random.rand()
        np.random.seed(seed)
        dist = fr.WeightedIndex(w)
        got = [dist.draw() for _ in range(40)]
        assert got == ref and np.random.rand() == nxt_ref


def test_three_body_sampler_matches_reference():
    rng = np.random.default_rng(12)
    for n, imbh in ((20000, False), (20000, True), (1500, True), (300, False)):
        m = _population(rng, n, imbh)
        for step in range(150):
            seed = int(rng.integers(0, 2**31))
            np.random.seed(seed); ref = _ref_3bb(m); nxt_ref = np.random.rand()
            np.random.seed(seed); got = tb.fast_sample_3bodyBinary(m, return_indices=True); nxt = np.random.rand()
            assert got == ref and nxt == nxt_ref, (imbh, step)
            m = _evolve(rng, m)


def test_elementwise_cache_with_zeroing():
    rng = np.random.default_rng(15)
    f = (lambda x: x**(-2/5), lambda x: np.exp(43/12 * np.log(x)))
    cache = ElementwiseCache(*f)
    m = _population(rng, 8000)
    for step in range(400):
        if step % 2:
            D = int(np.argmax(m))
            got = cache.get_zeroed(m, D)
            for g, fn in zip(got, f):
                ref = fn(m); ref[D] = 0.0
                assert np.array_equal(g.view(np.uint64), ref.view(np.uint64))
        else:
            got, _ = cache.get(m)
            for g, fn in zip(got, f):
                assert np.array_equal(g.view(np.uint64), fn(m).view(np.uint64))
        m = _evolve(rng, m)


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_'):
            fn()
            print('passed', name)
