'''
 Copyright (C) 2026  Konstantinos Kritos <kkritos1@jhu.edu>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''

# Bit-exact fast paths for the legacy global NumPy random stream (np.random.*).
#
# Every helper here returns exactly what the replaced np.random call returns and
# consumes exactly the same random numbers, so a run with a given seed produces
# identical output. The speed-up comes from building the cumulative distribution
# (CDF) of a weighted draw once and reusing it, instead of letting
# np.random.choice rebuild and re-validate it on every call, and, for large
# arrays, from locating the random number with block sums and a certified error
# bound instead of the sequential cumulative sum (CumulativeSearch). When an input
# is unusual (non-finite or negative weights, too few non-zero weights) the
# helpers call np.random.choice itself, so its error messages are unchanged.
#
# Reference: RandomState.choice in numpy/random/mtrand.pyx. For weights p it
# builds cdf = cumsum(p); cdf /= cdf[-1], draws uniforms u in [0, 1), and takes
# index = searchsorted(cdf, u, side='right'). Without replacement it zeroes the
# weights already drawn and repeats until it has enough distinct indices.

import numpy as np
from .constants import _LazyImport
_scipy_poisson = _LazyImport('scipy.stats', 'poisson')   # imported only if needed

# tolerance used by RandomState.choice when checking that p sums to 1:
_P_SUM_ATOL = np.sqrt(np.finfo(np.float64).eps)


def poisson_draw(mu):
    """Same value and same random-number consumption as scipy.stats.poisson.rvs(mu=mu).

    scipy draws through the global np.random state with RandomState.poisson, so
    for a finite mu >= 0 the direct NumPy call is identical and ~30x cheaper.
    Anything else is passed to scipy so its behaviour is unchanged.
    """
    mu_f = float(mu)
    if mu_f >= 0.0 and mu_f < np.inf:
        return np.random.poisson(mu_f)
    return _scipy_poisson.rvs(mu=mu)


def maxwell_draws(scale, size):
    """
    Same values and same random-number consumption as
    scipy.stats.maxwell.rvs(loc=0, scale=scale, size=size) (global np.random state).

    SciPy draws a Maxwell variate as chi(3) = sqrt(chi2(3)): RandomState.chisquare(3, size),
    then sqrt, then the scale (the intermediate *1 + 0 steps are exact). Done directly,
    this avoids importing scipy.stats (~0.6 s at start-up). As in SciPy, scale = 0
    returns zeros without drawing, and a negative scale is an error.
    """
    scale = float(scale)
    if not scale >= 0:
        raise ValueError('Domain error in arguments: the scale parameter must be positive.')
    if scale == 0:
        return 0 * np.ones(size, 'd')
    return np.sqrt(np.random.chisquare(3.0, size)) * scale + 0


def weights_cdf(p, nonnegative=False):
    """CDF of probabilities p, exactly as RandomState.choice builds it.

    Returns None when choice would reject p (NaN, negative entries, or a sum
    farther than its tolerance from 1); callers then fall back to np.random.choice.
    nonnegative=True: the caller guarantees that p has no negative entries (the
    check is skipped; NaNs still make the total NaN and are rejected).
    """
    cdf = np.cumsum(p)
    total = cdf[-1]
    if not np.isfinite(total) or abs(total - 1.0) > _P_SUM_ATOL or (not nonnegative and p.min() < 0):
        return None
    cdf /= total
    return cdf


def draw_index(cdf):
    """Index chosen by np.random.choice(a, p=p) (size=None, with or without
    replacement), given cdf = weights_cdf(p). Consumes one uniform."""
    return int(cdf.searchsorted(np.random.random_sample(), side='right'))


def first_index_of(a, k):
    """First index whose value equals a[k] (what np.where(a == a[k])[0][0] returns)."""
    v = a[k]
    if v != v:   # NaN: np.where finds nothing, so let it fail as before
        return int(np.flatnonzero(a == v)[0])
    if k > 0 and (a[:k] == v).any():
        return int(np.flatnonzero(a == v)[0])
    return k


def choice_weighted(a, p, nonnegative=False):
    """Exact, faster np.random.choice(a, p=p) for 1-D a and size=None.

    Valid for both replace=True and replace=False (one draw is the same either way).
    Returns (value, index). nonnegative: see weights_cdf.
    """
    cdf = weights_cdf(p, nonnegative)
    if cdf is None:
        value = np.random.choice(a, p=p)  # raises exactly as before
        return value, None
    k = draw_index(cdf)
    return a[k], k


def choice_pair_no_replace(a, p, cdf):
    """Exact np.random.choice(a, size=2, replace=False, p=p), with cdf = weights_cdf(p)
    built once by the caller. Returns the two sampled values in draw order.

    cdf may also be a CumulativeSearch of p (with p certified valid for choice),
    with which the exact cumulative sums are computed only if needed."""
    x = np.random.random_sample(2)
    if isinstance(cdf, CumulativeSearch):
        search = cdf
        search.n_draws += 1
        i1 = None
        if search._exact_normalized is None and search.n_draws <= CumulativeSearch.EXACT_AFTER:
            i0 = search.index(x[0])
            i1 = search.index(x[1]) if i0 is not None else None
        if i1 is None:
            i0, i1 = search.exact_normalized_cdf().searchsorted(x, side='right')
    else:
        i0, i1 = cdf.searchsorted(x, side='right')
    if i0 != i1:
        return a[i0], a[i1]

    # Both uniforms hit the same entry: numpy zeroes the weights already found and
    # draws the missing ones again from the renormalized distribution.
    found = [int(i0)]
    q = p.copy()
    while len(found) < 2:
        x = np.random.random_sample(2 - len(found))
        q[found] = 0
        c = np.cumsum(q)
        c /= c[-1]
        new = c.searchsorted(x, side='right')
        _, first = np.unique(new, return_index=True)
        first.sort()
        found.extend(int(v) for v in new.take(first))
    return a[found[0]], a[found[1]]


# unit roundoff of float64:
_UNIT_ROUNDOFF = 2.0**-53


class CumulativeSearch:
    """
    Locates a uniform number u in the cumulative distribution of nonnegative
    weights w, giving the index that searchsorted on the exact cumulative sums
    (np.cumsum) would give, without computing those sums in most cases.

    np.cumsum adds the entries one after the other (~3 ns per entry, the main cost
    of a draw from tens of thousands of BHs). Here the weights are summed in blocks
    with vectorized sums, which locates u within one block, and the partial sums of
    that block give approximations of the cumulative sums around u.

    Exactness: any way of adding m nonnegative numbers has a relative error below
    gamma_m = m u_r / (1 - m u_r), with u_r = 2^-53 (Higham, Accuracy and Stability
    of Numerical Algorithms, 2002, Sec. 4.2). This holds for the exact sequential
    sums as for the block sums, and the normalization by the total and the
    probabilities p = w/sum(w) used by np.random.choice add a few u_r. So, on the
    normalized scale [0, 1], the exact and the approximate cumulative values differ
    by less than the margin M = 16 (m + 4) u_r, with m >= n the number of additions
    (~1e-10 for n = 5e4). If u is farther than M from the approximate cumulative
    values on both sides of the located index, that index is provably the one the
    exact procedure returns. Otherwise (probability ~2M per draw) index() returns
    None and the caller computes the exact cumulative sums.

    The weights must be nonnegative (or NaN/inf, which make the search invalid).
    """

    BLOCK = 256        # block length
    MIN_SIZE = 4096    # below this the exact cumulative sums are cheap enough
    EXACT_AFTER = 16   # callers switch to the exact cumulative sums after this many draws

    def __init__(self, w):
        self.w = w
        n = w.shape[0]
        self.n = n
        B = self.BLOCK
        n_full = n // B
        main = n_full * B
        block_sums = w[:main].reshape(n_full, B).sum(axis=1) if n_full else np.empty(0)
        if main < n:
            block_sums = np.append(block_sums, w[main:].sum())
        self.cb = np.cumsum(block_sums)              # approximate cumulative sums at block ends
        self.total = self.cb[-1] if n else 0.0
        self.valid = bool(n > 0 and np.isfinite(self.total) and self.total > 0)
        n_additions = n + 2 * B + self.cb.shape[0] + 4
        self.margin = 16.0 * (n_additions + 4) * _UNIT_ROUNDOFF
        self._exact_unnormalized = None
        self._exact_normalized = None
        self._blocks = {}                            # block index -> approximate cumulative sums
        self.n_draws = 0                             # draws made by callers (see EXACT_AFTER)

    def index(self, u):
        """Index of u (in [0, 1), relative to the total) as searchsorted on the exact
        cumulative sums would give it, or None if u is too close to a boundary."""
        M = self.margin
        if not self.valid or not (M < u < 1.0 - M):
            return None
        T = self.total
        t = u * T
        b = int(self.cb.searchsorted(t, side='right'))
        if b >= self.cb.shape[0]:
            return None
        start = self.cb[b-1] if b > 0 else 0.0
        A = self._blocks.get(b)
        if A is None:
            B = self.BLOCK
            A = np.cumsum(self.w[b*B:(b+1)*B])
            A += start                               # approximate cumulative sums in this block
            self._blocks[b] = A
        i = int(A.searchsorted(t, side='right'))
        if i >= A.shape[0]:
            return None
        before = A[i-1] if i > 0 else start
        if before < (u - M) * T and A[i] > (u + M) * T:
            return b*self.BLOCK + i
        return None

    def exact_cumsum(self):
        """np.cumsum(w) (computed once, when needed)."""
        if self._exact_unnormalized is None:
            self._exact_unnormalized = np.cumsum(self.w)
        return self._exact_unnormalized

    def exact_normalized_cdf(self):
        """weights_cdf(w): the normalized cumulative distribution np.random.choice uses."""
        if self._exact_normalized is None:
            self._exact_normalized = weights_cdf(self.w)
        return self._exact_normalized

    def certified_choice_valid(self):
        """True if w (probabilities summing to ~1) is certainly accepted by
        np.random.choice, i.e. weights_cdf(w) would not return None."""
        return self.valid and abs(self.total - 1.0) < 1e-9


class WeightedIndex:
    """
    Exact, faster replacement of
        cdf = np.cumsum(w)
        j = int(np.searchsorted(cdf, np.random.rand() * cdf[-1], side='right'))
        if j >= cdf.size: j = int(np.flatnonzero(w > 0)[-1])
    for repeated draws from nonnegative weights w (same index, same random number).
    """

    def __init__(self, w):
        self.w = w
        self.search = CumulativeSearch(w) if w.shape[0] >= CumulativeSearch.MIN_SIZE else None
        self.cdf = None

    def exact_cumsum(self):
        if self.cdf is None:
            self.cdf = self.search.exact_cumsum() if self.search is not None else np.cumsum(self.w)
        return self.cdf

    def total(self):
        """cdf[-1] of the exact cumulative sums (computed only if needed)."""
        return self.exact_cumsum()[-1]

    def draw(self):
        u = np.random.rand()
        if self.search is not None and self.cdf is None:
            # after many draws from the same weights, computing the exact cumulative
            # sums once is cheaper than locating each draw with the block sums:
            self.search.n_draws += 1
            if self.search.n_draws <= CumulativeSearch.EXACT_AFTER:
                j = self.search.index(u)
                if j is not None:
                    return j
        cdf = self.exact_cumsum()
        j = int(np.searchsorted(cdf, u * cdf[-1], side='right'))
        if j >= cdf.size:  # u*total rounded up to the total: take the last positive weight
            j = int(np.flatnonzero(self.w > 0)[-1])
        return j

class NormalizedChoice:
    """
    Exact, faster np.random.choice(a, p=w/np.sum(w)) (size=None) for nonnegative
    weights w, for one or more draws from the same weights: the uniform number is
    located with a CumulativeSearch of w (no sequential cumulative sum in most
    cases), with the exact cumulative distribution of p = w/np.sum(w) used for
    small arrays, for numbers too close to a boundary, and after many draws.
    """

    def __init__(self, w):
        self.w = w
        self.p = None
        self.cdf = None
        self.search = None
        if w.shape[0] >= CumulativeSearch.MIN_SIZE:
            search = CumulativeSearch(w)
            # w >= 0 (or NaN, which makes the search invalid). With a finite, positive
            # total far below the overflow limit, np.sum(w) is finite and positive too,
            # so p = w/np.sum(w) is certainly accepted by np.random.choice:
            if search.valid and search.total < 1e300:
                self.search = search

    def probabilities(self):
        if self.p is None:
            self.p = self.w/np.sum(self.w)
        return self.p

    def exact_cdf(self):
        if self.cdf is None:
            self.cdf = weights_cdf(self.probabilities())
        return self.cdf

    def draw(self):
        """Index drawn by np.random.choice(a, p=w/np.sum(w)), or None if numpy would
        reject the probabilities (no random number consumed; the caller then lets
        np.random.choice raise as before)."""
        if self.search is not None:
            u = np.random.random_sample()
            self.search.n_draws += 1
            k = None
            if self.cdf is None and self.search.n_draws <= CumulativeSearch.EXACT_AFTER:
                k = self.search.index(u)
            if k is None:
                k = int(self.exact_cdf().searchsorted(u, side='right'))
            return k
        cdf = self.exact_cdf()
        if cdf is None:
            return None
        return draw_index(cdf)

# End of file.
