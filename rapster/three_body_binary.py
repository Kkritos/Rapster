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

from .constants import *
from .functions import *
from .arraycache import ElementwiseCache
from .fastrandom import WeightedIndex
from .bhpool import FAM

def _sums_of_others(w):
    """
    For each entry j, the sum of all other entries of w (w >= 0). Computed as
    total - w_j; at most one entry can exceed half the total, and for it (where
    total - w_j would lose precision) the other entries are summed directly.
    """
    total = np.sum(w)
    others = total - w
    j = int(np.argmax(w))
    if w[j] > 0.5 * total:
        others[j] = np.sum(w[:j]) + np.sum(w[j+1:])
    return others

# Draws of an index j with probability weights[j] / sum(weights) use
# WeightedIndex(weights).draw() (fastrandom.py): the uniform number u is placed in
# the cumulative sums, j = searchsorted(cumsum(weights), u * total, side='right'),
# with the sequential cumulative sum computed only when needed (exact).

def _draw_excluding(dist, excluded, accept=None, efficiency=1.0):
    """
    Index j, not in `excluded`, drawn with probability ~ w[j] * accept(j), where
    dist = WeightedIndex(w) and accept(j) in [0, 1] (default 1; accept must also
    take an index array). Uses rejection from the full distribution when its
    expected acceptance `efficiency` exceeds 0.25, and otherwise builds the
    conditional distribution directly.
    """
    if efficiency > 0.25:
        while True:
            j = dist.draw()
            if j not in excluded and (accept is None or np.random.rand() < accept(j)):
                return j
    w = dist.w.copy()
    w[list(excluded)] = 0.0
    if accept is not None:
        w *= np.clip(accept(np.arange(w.size)), 0.0, 1.0)
    return WeightedIndex(w).draw()

def _efficiency_others(dist, a):
    """
    (T - w[a]) / T > 0.25, with T = np.cumsum(w)[-1] the sequential total of
    dist = WeightedIndex(w). Decided from the block-sum total when it certainly
    gives the same answer (the two totals differ by a relative ~n*2^-53, far
    below the margin used), otherwise from the exact total. Returns the
    efficiency to pass to _draw_excluding (1.0 or 0.0 when decided approximately).
    """
    search = dist.search
    if search is not None and search.valid:
        T = search.total
        eff = (T - dist.w[a]) / T
        if abs(eff - 0.25) > search.margin:
            return 1.0 if eff > 0.25 else 0.0
    T = dist.total()
    return (T - dist.w[a]) / T

# exp(power*log(m)) for the weight powers of fast_sample_3bodyBinary, one row per
# power 4, 5/2, 15/4, 43/12, 7/3 (ElementwiseCache: exact, reuses unchanged BHs;
# each row equals np.exp(power * np.log(m)) computed on its own):
_WEIGHT_POWERS = np.array([4.0, 2.5, 15/4, 43/12, 7/3])[:, None]
_weight_powers = ElementwiseCache(lambda x: np.exp(_WEIGHT_POWERS * np.log(x)), min_size=6000)


def fast_sample_3bodyBinary(mBH, n_samples=1, return_indices=False):
    """
    Sample the three single BHs of a three-body binary formation: m1 and m2 form
    the binary, m3 is the third body that carries away the binding energy.

    The masses are drawn, as three distinct BHs, from the joint distribution
        p(m1, m2, m3) ~ (m1 m2)^4 m3^(5/2) g,   g = (m1 + m2)^(-1/2) (m1 + m2 + m3)^(-1/2)
    (the p_3bodyBinary density of earlier versions), by rejection sampling with
    mass-dependent bounds on g. With D the most massive BH, the triples split into

      A) D in the binary:      g <= 1/m_D
      C) D is the third body:  g <= (4 m1 m2)^(-1/4) m_D^(-1/2)   [m1+m2 >= 2 sqrt(m1 m2)]
      B) D not involved:       g <= (m1 m2)^(-5/12) m3^(-1/6) / sqrt(6)
                               [AM-GM: m1+m2 >= 2 sqrt(m1 m2), m1+m2+m3 >= 3 (m1 m2 m3)^(1/3)]

    In each region the bound times (m1 m2)^4 m3^(5/2) factorizes, so candidates are
    drawn exactly from it (region chosen with its total weight, then the BHs one
    after another), and accepted with probability g / bound <= 1. The bounds are
    tight for equal masses (B) and when D dominates (A, C), so the acceptance stays
    ~0.6-1 both for ordinary BH populations and with an IMBH or a massive seed BH;
    a single bound over all triples would accept only a fraction ~(m/m_D) when
    one BH of mass m_D dominates.

    @in mBH: array of single BH masses [Msun] (at least 3, all > 0)
    @in n_samples: number of samples (m1, m2, m3)
    @in return_indices: if True, each sample also carries the indices in mBH of m1 and m2

    @out: (m1, m2, m3), or (m1, m2, m3, index of m1, index of m2); a list if n_samples > 1
    """

    m = np.asarray(mBH, dtype=float)
    D = int(np.argmax(m))
    mD = m[D]

    # weights m**power = exp(power*log(m)) for the five powers below, with the most
    # massive BH D removed (set to zero). The powers are evaluated only for BHs that
    # are new or changed since the last call; the array is read-only and valid
    # until the next call (see ElementwiseCache.get_zeroed).
    (W,) = _weight_powers.get_zeroed(m, D)
    w, v, wC, q, r = W
    # w:  binary members, regions A: m^4
    # v:  third body, region A: m^(5/2)
    # wC: binary members, region C: m^4 (4 m)^(-1/4)
    # q:  binary members, region B: m^4 m^(-5/12)
    # r:  third body, region B: m^(5/2) m^(-1/6)
    # (Rows of 2-D arrays are processed with one NumPy call for all rows; each row's
    # sum, argmax and elementwise results are identical to those of the row alone.)
    totals = W.sum(axis=1)                        # np.sum of each row
    V = totals[1]
    Q = totals[3]

    # sums over all BHs except one (as _sums_of_others) for wC, q and r:
    others = totals[2:5, None] - W[2:5]
    j_max = W[2:5].argmax(axis=1)
    for k in range(3):
        j = int(j_max[k])
        row = W[2+k]
        if row[j] > 0.5 * totals[2+k]:
            others[k, j] = np.sum(row[:j]) + np.sum(row[j+1:])
    wC_others, q_others, r_others = others        # r_others[a] = sum_{c != a, D} r_c

    # first-member weights of the three regions, one row each:
    F = np.empty((3, m.shape[0]))
    wA, wC_first, wB_first = F

    # region A: (D, x, c) or (x, D, c); weight w_D w_x v_c / m_D; x ~ w_x (V - v_x), c ~ v_c:
    # wA = w * (V - v)
    np.subtract(V, v, out=wA)
    np.multiply(w, wA, out=wA)

    # region C: (a, b, D); weight 4^(-1/4) m_D^(5/2) m_D^(-1/2) wC_a wC_b:
    # wC_first = wC * sum_{b != a} wC_b
    np.multiply(wC, wC_others, out=wC_first)

    # region B: (a, b, c), none of them D; weight q_a q_b r_c / sqrt(6):
    # wB_first = max(q * (r_others*sum_{b != a} q_b - sum_{b != a} q_b r_b), 0)
    qr_others = _sums_of_others(q*r)
    np.multiply(r_others, q_others, out=wB_first)
    np.subtract(wB_first, qr_others, out=wB_first)
    np.multiply(q, wB_first, out=wB_first)
    np.maximum(wB_first, 0.0, out=wB_first)

    sums = F.sum(axis=1)
    Z_A = 2 * mD**4 / mD * sums[0]
    Z_C = 4**(-0.25) * mD**2 * sums[1]
    Z_B = sums[2] / np.sqrt(6)

    # weighted-draw helpers, built when first needed:
    dists = {}
    def dist_of(name, arr):
        if name not in dists:
            dists[name] = WeightedIndex(arr)
        return dists[name]

    Z = np.array([Z_A, Z_C, Z_B])
    Z_dist = WeightedIndex(Z)

    results = []
    while len(results) < n_samples:

        region = Z_dist.draw()

        if region == 0:   # A: D is a binary member
            x = dist_of('wA', wA).draw()
            c = _draw_excluding(dist_of('v', v), {x}, efficiency=(V - v[x]) / V)
            a, b = (D, x) if np.random.rand() < 0.5 else (x, D)
            bound = 1 / mD

        elif region == 1: # C: D is the third body
            a = dist_of('wC_first', wC_first).draw()
            # efficiency (WC - wC[a]) / WC, with WC the total of wC:
            b = _draw_excluding(dist_of('wC', wC), {a}, efficiency=_efficiency_others(dist_of('wC', wC), a))
            c = D
            bound = (4 * m[a] * m[b])**(-0.25) * mD**(-0.5)

        else:             # B: D not involved
            a = dist_of('wB_first', wB_first).draw()
            Ra = r_others[a]   # sum of r over all BHs but a (and D)
            # b ~ q_b (Ra - r_b), b != a:
            b = _draw_excluding(dist_of('q', q), {a}, accept=lambda j: (Ra - r[j]) / Ra,
                                efficiency=wB_first[a] / (q[a] * Q * Ra))
            # c ~ r_c, c != a, b:
            c = _draw_excluding(dist_of('r', r), {a, b}, efficiency=(Ra - r[b]) / (Ra + r[a]))
            bound = (m[a]*m[b])**(-5/12) * m[c]**(-1/6) / np.sqrt(6)

        m1, m2, m3 = mBH[a], mBH[b], mBH[c]
        g = (m1 + m2)**(-0.5) * (m1 + m2 + m3)**(-0.5)
        if np.random.rand() * bound < g:
            results.append((m1, m2, m3, a, b) if return_indices else (m1, m2, m3))

    return results[0] if n_samples == 1 else results

# indices of the pool weight families used below (bhpool.FAMILIES):
_I_W, _I_V, _I_WC, _I_Q, _I_R, _I_WV, _I_WCWC, _I_QQ, _I_QR = (FAM[f] for f in ('w', 'v', 'wC', 'q', 'r', 'wv', 'wCwC', 'qq', 'qr'))

def sample_3bodyBinary(BH, exclude=None, max_trials=None):
    """
    Sample the three single BHs of a three-body binary formation from the pool BH
    (bhpool.BHPool), with the same distribution as fast_sample_3bodyBinary:
    three distinct BHs with
        p(m1, m2, m3) ~ (m1 m2)^4 m3^(5/2) g,   g = (m1 + m2)^(-1/2) (m1 + m2 + m3)^(-1/2).

    The regions A, C, B and their bounds on g are those of fast_sample_3bodyBinary
    (D = the most massive single BH). Here the BHs of a region are proposed
    independently of each other from the pool's weight families (with D excluded),
    i.e. including proposals in which the same BH appears twice; such proposals are
    rejected, together with the usual acceptance g / bound. This keeps each trial
    to a few draws of O(n/BLOCK + BLOCK) operations, with no pass over all BHs.

    The fraction of proposals with a repeated BH is estimated beforehand; it is
    small unless very few BHs carry almost all the weight, and if it exceeds 1/2
    the exact sampler fast_sample_3bodyBinary is used on the pool's mass array.

    exclude: optional set of slots (not containing the heaviest BH D) that are still in
    the pool but no longer available (taken by earlier binaries of the same step);
    proposals that contain one of them are rejected, which gives the same
    distribution over the remaining BHs as removing them first.
    max_trials: with exclude, give up (return None) after this many proposals, so that
    the caller can remove the excluded BHs from the pool and draw again.

    @in BH: pool of single BHs (at least 3)

    @out: (m1, m2, m3, slot of m1, slot of m2), or None (see max_trials)
    """
    D = BH.D
    m = BH.m
    mD = float(m[D])

    T = BH.totals()
    W_, V_, WC_, Q_, R_ = T[_I_W], T[_I_V], T[_I_WC], T[_I_Q], T[_I_R]

    # region totals of the proposal weights (sums over all non-D BHs, repeats included):
    Z_A = 2 * mD**3 * W_ * V_                 # (D, x, c) and (x, D, c): m_D^4 w_x v_c / m_D
    Z_C = 4**(-0.25) * mD**2 * WC_**2         # (a, b, D): 4^(-1/4) m_D^(5/2) m_D^(-1/2) wC_a wC_b
    Z_B = Q_**2 * R_ / np.sqrt(6)             # (a, b, c): q_a q_b r_c / sqrt(6)
    Z = Z_A + Z_C + Z_B

    # probability that a proposal repeats a BH (for region B an upper bound):
    p_rep = (Z_A * T[_I_WV] / (W_ * V_)
             + Z_C * T[_I_WCWC] / WC_**2
             + Z_B * (T[_I_QQ] / Q_**2 + 2 * T[_I_QR] / (Q_ * R_))) / Z
    if not p_rep < 0.5:
        if exclude:
            return None          # (the caller removes the excluded BHs and draws again)
        return fast_sample_3bodyBinary(m, n_samples=1, return_indices=True)

    fw, fv, fwC, fq, fr = _I_W, _I_V, _I_WC, _I_Q, _I_R
    trials = 0
    while True:
        if exclude:
            trials += 1
            if trials > max_trials:
                return None
        u = np.random.rand() * Z
        if u < Z_A:     # A: D is a binary member
            x = BH.draw(fw)
            c = BH.draw(fv)
            if x == c:
                continue
            a, b = (D, x) if np.random.rand() < 0.5 else (x, D)
            bound = 1 / mD
        elif u < Z_A + Z_C:   # C: D is the third body
            a = BH.draw(fwC)
            b = BH.draw(fwC)
            if a == b:
                continue
            c = D
            bound = (4 * m[a] * m[b])**(-0.25) * mD**(-0.5)
        else:           # B: D not involved
            a = BH.draw(fq)
            b = BH.draw(fq)
            c = BH.draw(fr)
            if a == b or c == a or c == b:
                continue
            bound = (m[a] * m[b])**(-5/12) * m[c]**(-1/6) / np.sqrt(6)

        if exclude and (a in exclude or b in exclude or c in exclude):
            continue

        m1, m2, m3 = m[a], m[b], m[c]
        g = (m1 + m2)**(-0.5) * (m1 + m2 + m3)**(-0.5)
        if np.random.rand() * bound < g:
            return m1, m2, m3, a, b

def three_body_binary(t, z, k_3bb, mBH_avg, binaries, BH, vBH, N_3bb, N_BBH, random_pairing=False):
    """
    @in t: simulation time
    @in z: simulation redshift
    @in k_3bb: number of 3bbs in current step
    @in mBH_avg: average BH mass
    @in binaries: array of BBHs
    @in BH: pool of single BHs (bhpool.BHPool: masses, spins, generations, TDE counts)
    @in vBH: 3D BH velocity dispersion
    @in N_3bb: number of 3bbs
    @in N_BBH: number of BBHs
    @in random_pairing: if True, use uniform random pairing instead of the mass-weighted
                       distribution of sample_3bodyBinary

    @out: all inputs
    """
    
    if k_3bb>0:

        new_rows = []

        # BHs that form binaries are removed from the pool at the end of the step (or
        # when the heaviest BH is taken, or when many proposals hit taken BHs): until
        # then they are excluded from the draws, which leaves the pool's sums unchanged
        # and the draws cheap (see sample_3bodyBinary):
        taken = set()

        def flush():
            for k in sorted(taken, reverse=True):   # (descending: no taken BH is moved)
                BH.remove(k)
            taken.clear()
        
        for i in range(k_3bb):

            if BH.size - len(taken) < 3:  # not enough single BHs left to form a 3bb, exit loop
                break

            N_3bb+=1 # update number of three-body-binaries

            # sample the BHs that form the 3bb (slots k1 != k2 in the pool):
            if random_pairing:
                k1 = BH.draw_uniform()
                while k1 in taken:
                    k1 = BH.draw_uniform()
                k2 = BH.draw_uniform()
                while k2 == k1 or k2 in taken:
                    k2 = BH.draw_uniform()
                m1, m2 = BH.m[k1], BH.m[k2]
                m3 = mBH_avg
            else:
                # mass-weighted pairing:
                if BH.D in taken:
                    flush()
                sample = sample_3bodyBinary(BH, exclude=taken, max_trials=30) if taken else None
                if sample is None:
                    flush()
                    sample = sample_3bodyBinary(BH)
                m1, m2, m3, k1, k2 = sample

            # initial hardness of newly formed 3bb:
            eta = sample_hardness()
            
            # semimajor axis:
            sma = G_Newton * m1 * m2 / eta / m3 / vBH**2
            
            # eccentricity (thermal):
            eccen = np.sqrt(np.random.rand())
            
            # new binary (appended to binaries after the loop):
            s1, s2, g1, g2, h1, h2 = BH._s[k1], BH._s[k2], BH._g[k1], BH._g[k2], BH._h[k1], BH._h[k2]
            new_rows.append([np.random.randint(0, 999999999), 3, sma, eccen, m1, m2, s1, s2, g1, g2, t, z, 0, h1, h2])
            
            # update number of in-cluster BBHs:
            N_BBH+=1
            
            # the single BHs leave the pool (at the next flush):
            taken.add(k1)
            taken.add(k2)

        flush()

        if new_rows:
            binaries = np.append(binaries, new_rows, axis=0)
            
    return t, z, k_3bb, mBH_avg, binaries, BH, vBH, N_3bb, N_BBH

# End of file.
