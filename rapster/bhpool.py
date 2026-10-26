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

# Pool of single (unbound) compact objects (BHs and retained NSs) with fast
# weighted draws.
#
# The single BHs used to be four NumPy arrays (masses, spins, generations, TDE
# counts) that were copied whenever a BH was added or removed, and every weighted
# draw (three-body binaries, BBH-single partners, TDEs, two-body captures,
# exchanges) recomputed weights such as m^4 for all BHs and their cumulative sum.
# In compact clusters, with ~1e4-1e5 single BHs and up to ~1e6 such events, these
# passes over all BHs dominated the run time.
#
# BHPool keeps the BHs in arrays with spare capacity ("slots" 0..n-1). Removing a
# BH moves the last BH into its slot; adding a BH appends it (both O(1), no copy).
# For each weight family (m^4, m^(5/2), ...) the weight of every BH is stored,
# together with the sum over each block of BLOCK consecutive slots. A change of
# one BH updates its weights and the sums of its block. A weighted draw locates
# the random number among the block sums and then inside one block, which costs
# O(n/BLOCK + BLOCK) (a few microseconds) instead of O(n).
#
# Draws follow the same probabilities as before; the samplers built on them (in
# three_body_binary.py, binary_evolution.py, two_body_capture.py,
# tidal_disruptions.py and exchanges.py) use rejection sampling where the old
# weights do not factorize. Because the order in which BHs are stored and the
# random numbers consumed differ, a given seed gives a different (equally
# probable) realization than versions <= 2.13.1; ensemble statistics agree
# (tests/test_bhpool.py, and the validation described in CHANGELOG.md).
#
# The families used by the three-body sampler exclude the most massive single BH
# (called D in three_body_binary.py): its weights are kept at zero there.

import math

import numpy as np

# weight families: (name, power of m, excludes the heaviest BH D)
#   three-body binaries (three_body_binary.py):
#     w = m^4, v = m^(5/2), wC = m^(15/4), q = m^(43/12), r = m^(7/3), and the
#     products wv, wCwC, qq, qr (probabilities that proposed BHs coincide);
#   BBH-single partners (binary_evolution.py): p15, p25, p17, p27 = m^1.5, m^2.5, m^1.7, m^2.7
#     (p25 also gives sum(m^(5/2)) for the BH relaxation factor);
#   BH-star TDEs (tidal_disruptions.py): c13 = m^(1/3), c43 = m^(4/3);
#   two-body captures (two_body_capture.py): c27 = m^(2/7), c127 = m^(12/7);
#   exchanges and mean mass: mass = m.
FAMILIES = (
    ('w',    4.0,          True),
    ('v',    2.5,          True),
    ('wC',   15/4,         True),
    ('q',    43/12,        True),
    ('r',    7/3,          True),
    ('wv',   6.5,          True),
    ('wCwC', 7.5,          True),
    ('qq',   43/6,         True),
    ('qr',   43/12 + 7/3,  True),
    ('p15',  1.5,          False),
    ('p25',  2.5,          False),
    ('p17',  1.7,          False),
    ('p27',  2.7,          False),
    ('c13',  1/3,          False),
    ('c43',  4/3,          False),
    ('c27',  2/7,          False),
    ('c127', 12/7,         False),
    ('mass', 1.0,          False),
)
FAM = {name: i for i, (name, _, _) in enumerate(FAMILIES)}
N_FAMILIES = len(FAMILIES)
_POWER_LIST = [p for _, p, _ in FAMILIES]
_POWERS = np.array(_POWER_LIST)[:, None]
_N_EXCL = sum(e for _, _, e in FAMILIES)       # the families excluding D come first
assert all(e for _, _, e in FAMILIES[:_N_EXCL]) and not any(e for _, _, e in FAMILIES[_N_EXCL:])


def family_weights(m):
    """Weights of all families for the masses m (1-D array): shape (N_FAMILIES, len(m))."""
    m = np.asarray(m, dtype=float)
    W = np.exp(_POWERS * np.log(m))
    W[FAM['mass']] = m
    return W


def _weights_of(m):
    """Weights of all families for one mass (list)."""
    m = float(m)
    return [m**p for p in _POWER_LIST]

_ZEROS_EXCL = [0.0] * _N_EXCL


class BHPool:
    """
    Single compact objects: masses m, spins s, generations g, TDE counts h
    (views of the first n slots, not to be modified directly), with weighted draws
    over the weight families of FAMILIES (index FAM[name]).

    Slots change when BHs are removed (the last BH moves into the freed slot), so a
    slot number is valid only until the next removal.
    """

    BLOCK = 128

    def __init__(self, m=(), s=(), g=(), h=()):
        m = np.asarray(m, dtype=float)
        n = m.shape[0]
        B = self.BLOCK
        cap = (max(1024, 2 * n) // B + 1) * B
        self._m = np.zeros(cap); self._s = np.zeros(cap); self._g = np.zeros(cap); self._h = np.zeros(cap)
        self._m[:n] = m; self._s[:n] = s; self._g[:n] = g; self._h[:n] = h
        self.n = n
        self._W = np.zeros((cap, N_FAMILIES))          # weights, one row per slot
        self._S = np.zeros((cap // B, N_FAMILIES))     # sums over blocks of BLOCK slots
        self._ST = np.zeros((N_FAMILIES, cap // B))    # the same, one contiguous row per family
        self._ones = np.ones(B)
        self.D = -1
        self._tot = None   # totals of all families, until the next change
        self._cb = {}      # cumulative block sums per family, until the next change
        self._bver = [0] * (cap // B)   # version of each block, increased whenever it changes
        self._wcum = {}    # (family, block) -> (block version, weights, cumulative weights)
        self._rebuild()

    # ------------------------------------------------------------------ views
    @property
    def m(self):
        return self._m[:self.n]

    @property
    def s(self):
        return self._s[:self.n]

    @property
    def g(self):
        return self._g[:self.n]

    @property
    def h(self):
        return self._h[:self.n]

    @property
    def size(self):
        return self.n

    def __len__(self):
        return self.n

    def get(self, k):
        """(m, s, g, h) of slot k, as Python floats."""
        return float(self._m[k]), float(self._s[k]), float(self._g[k]), float(self._h[k])

    def max_mass(self):
        return float(self._m[self.D]) if self.n else 0.0

    def mean_mass(self):
        return self.total(FAM['mass']) / self.n

    # --------------------------------------------------------------- internals
    def _changed(self):
        self._tot = None
        if self._cb:
            self._cb = {}

    def _rebuild(self):
        """Recompute all weights and block sums from the masses (O(n))."""
        n = self.n
        self._W[:] = 0.0
        if n:
            self._W[:n] = family_weights(self._m[:n]).T
            self.D = int(np.argmax(self._m[:n]))
            self._W[self.D, :_N_EXCL] = 0.0
        else:
            self.D = -1
        self._S[:] = self._W.reshape(-1, self.BLOCK, N_FAMILIES).sum(axis=1)
        self._ST[:] = self._S.T
        self._wcum = {}
        self._changed()

    def _update_block(self, k):
        B = self.BLOCK
        b = k // B
        np.dot(self._ones, self._W[b*B:(b+1)*B], out=self._S[b])   # (BLAS: faster than .sum(axis=0))
        self._ST[:, b] = self._S[b]
        self._bver[b] += 1

    def _set_weights(self, k):
        """Weights of slot k from its mass (zero in the D-excluding families if k is D)."""
        w = _weights_of(self._m[k])
        if k == self.D:
            w[:_N_EXCL] = _ZEROS_EXCL
        self._W[k] = w

    def _find_D(self):
        """Make the heaviest BH D (after the old D was removed or became lighter)."""
        old = self.D
        new = int(np.argmax(self._m[:self.n])) if self.n else -1
        if old == new:
            return
        self.D = new
        if 0 <= old < self.n:
            self._set_weights(old)
            self._update_block(old)
        if new >= 0:
            self._W[new, :_N_EXCL] = 0.0
            self._update_block(new)

    def _make_D(self, k):
        """Slot k (whose mass now exceeds the current D's) becomes D."""
        old = self.D
        self.D = k
        self._W[k, :_N_EXCL] = 0.0
        self._update_block(k)
        if 0 <= old < self.n and old != k:
            self._set_weights(old)
            self._update_block(old)

    def _grow(self):
        B = self.BLOCK
        cap = self._m.shape[0]
        new_cap = (2 * cap // B + 1) * B
        for name in ('_m', '_s', '_g', '_h'):
            b = np.zeros(new_cap); b[:cap] = getattr(self, name)
            setattr(self, name, b)
        W = np.zeros((new_cap, N_FAMILIES)); W[:cap] = self._W; self._W = W
        S = np.zeros((new_cap // B, N_FAMILIES)); S[:cap // B] = self._S; self._S = S
        ST = np.zeros((N_FAMILIES, new_cap // B)); ST[:, :cap // B] = self._ST; self._ST = ST
        self._bver += [0] * (new_cap // B - cap // B)

    # ----------------------------------------------------------------- updates
    def add(self, m, s, g, h):
        """Add a compact object; returns its slot."""
        if self.n == self._m.shape[0]:
            self._grow()
        k = self.n
        self._m[k] = m; self._s[k] = s; self._g[k] = g; self._h[k] = h
        self.n += 1
        self._changed()
        if self.D < 0 or self._m[k] > self._m[self.D]:
            self._W[k] = _weights_of(m)
            self._make_D(k)
        else:
            self._W[k] = _weights_of(m)
            self._update_block(k)
        return k

    def add_many(self, m, s, g, h):
        for x in zip(np.atleast_1d(m), np.atleast_1d(s), np.atleast_1d(g), np.atleast_1d(h)):
            self.add(*x)

    def remove(self, k):
        """Remove slot k; the last BH moves into slot k."""
        n = self.n
        if not 0 <= k < n:
            raise IndexError(f'slot {k} out of range for {n} single BHs')
        self._changed()
        last = n - 1
        removed_D = (k == self.D)
        if k != last:
            self._m[k] = self._m[last]; self._s[k] = self._s[last]
            self._g[k] = self._g[last]; self._h[k] = self._h[last]
            self._W[k] = self._W[last]
            if self.D == last:
                self.D = k
        self._m[last] = 0.0; self._s[last] = 0.0; self._g[last] = 0.0; self._h[last] = 0.0
        self._W[last] = 0.0
        self.n = last
        if removed_D:
            self.D = -1
            if k != last:
                self._set_weights(k)          # the moved BH, with all its weights
        self._update_block(k)
        if last // self.BLOCK != k // self.BLOCK:
            self._update_block(last)
        if removed_D:
            self._find_D()

    def remove_two(self, k1, k2):
        """Remove two distinct slots."""
        if k1 == k2:
            raise ValueError('slots must differ')
        self.remove(max(k1, k2))
        self.remove(min(k1, k2))

    def replace(self, k, m, s, g, h):
        """Put another compact object in slot k (e.g. after an exchange)."""
        self._s[k] = s; self._g[k] = g; self._h[k] = h
        self.set_mass(k, m)

    def set_mass(self, k, m):
        """Change the mass of slot k (e.g. after accretion)."""
        self._changed()
        self._m[k] = m
        if k == self.D:
            self._set_weights(k)
            self._update_block(k)
            self._find_D()                    # D may have become lighter than another BH
        elif m > self._m[self.D]:
            self._W[k] = _weights_of(m)
            self._make_D(k)
        else:
            self._W[k] = _weights_of(m)
            self._update_block(k)

    def set_spin(self, k, s):
        self._s[k] = s

    def set_tde_count(self, k, h):
        self._h[k] = h

    def set_all_masses(self, m):
        """Replace all masses (gas accretion), keeping spins, generations and TDE counts."""
        self._m[:self.n] = m
        self._rebuild()

    # ------------------------------------------------------------------- draws
    def totals(self):
        """Sums of the weights of all families over all single BHs (array, do not modify)."""
        tot = self._tot
        if tot is None:
            nb = (self.n + self.BLOCK - 1) // self.BLOCK
            tot = self._tot = self._S[:nb].sum(axis=0)
        return tot

    def _block_cdf(self, f):
        """Cumulative block sums of family f (sequential sums; cached until the next change)."""
        cb = self._cb.get(f)
        if cb is None:
            nb = (self.n + self.BLOCK - 1) // self.BLOCK
            cb = self._cb[f] = np.cumsum(self._ST[f, :nb])
        return cb

    def total(self, f):
        """Sum of the weights of family f over all single BHs (= totals()[f]: both are
        the sequential sum of the block sums)."""
        cb = self._block_cdf(f)
        return float(cb[-1]) if cb.shape[0] else 0.0

    def weight(self, f, k):
        return float(self._W[k, f])

    def draw(self, f):
        """Slot drawn with probability proportional to its weight in family f."""
        cb = self._cb.get(f)
        if cb is None:
            cb = self._block_cdf(f)
        B = self.BLOCK
        nb = cb.shape[0]
        total = cb[-1]
        wcum = self._wcum
        while True:
            t = np.random.rand() * total
            b = int(cb.searchsorted(t, side='right'))
            if b >= nb:
                continue                       # (rounding at the end: draw again)
            if b:
                t -= cb[b-1]
            # cumulative weights within block b (reused while the block is unchanged):
            ver = self._bver[b]
            entry = wcum.get((f, b))
            if entry is None or entry[0] != ver:
                w = self._W[b*B:(b+1)*B, f]
                entry = wcum[(f, b)] = (ver, w, np.cumsum(w))
            w = entry[1]
            i = int(entry[2].searchsorted(t, side='right'))
            if i < B and w[i] > 0.0:
                return b*B + i
            # (t fell beyond the block's own cumulative sum, by rounding: draw again)

    def draw_uniform(self):
        """Slot drawn uniformly."""
        return int(np.random.rand() * self.n) % self.n

    def check(self):
        """Consistency check (used by the tests): weights, block sums, D."""
        n = self.n
        W = np.zeros_like(self._W)
        if n:
            W[:n] = family_weights(self._m[:n]).T
            assert self._m[self.D] == self._m[:n].max()
            W[self.D, :_N_EXCL] = 0.0
        assert np.allclose(self._W, W, rtol=1e-12, atol=0)
        S = self._W.reshape(-1, self.BLOCK, N_FAMILIES).sum(axis=1)
        assert np.allclose(self._S, S, rtol=1e-12, atol=0)
        assert np.array_equal(self._ST, self._S.T)
        assert not np.any(self._m[n:]) and not np.any(self._W[n:])
        # total(f) (from the cumulative block sums) equals totals()[f] bit for bit:
        tot = self.totals()
        assert all(self.total(f) == float(tot[f]) for f in range(N_FAMILIES))

# End of file.
