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

# Multiple stellar populations.
#
# Rapster describes the stars of a cluster by a few global quantities: their total mass
# M_cl, the mean stellar mass <m>, the half-mass radius r_h and the central density. With
# a single population born at t = 0 these follow from the age t alone: <m> =
# <m>_0 (t/t_sev)^(-nu_sev), stellar-evolution mass loss dM = -nu_sev M dt/t, stars up to
# the turn-off mass (1e4 Myr/t)^(1/2.5) Msun, white dwarfs forming from stars at the
# turn-off, and BHs and NSs forming at t ~ tBH_form = 10 Myr.
#
# With several populations k, born at times t_k, each keeps its own stellar mass M_k and
# age tau_k = t - t_k, and the cluster quantities are sums over them:
#   M_cl = sum_k M_k,   N = sum_k N_k with N_k = M_k / <m>_k(tau_k),   <m> = M_cl / N,
#   dM_k(stellar evolution) = -nu_sev M_k dt / tau_k   (tau_k > t_sev),
# and stars for encounters and tidal disruptions are drawn from population k with
# probability proportional to its number of stars in the cluster center (see below), from
# the IMF up to that population's turn-off mass. All populations share the IMF (-mm, -mM,
# -mb) and the binary fraction (-fb); each has its own birth time, number of stars,
# metallicity and star formation efficiency. A population born after t = 0 forms its BHs
# and NSs like the first one (BHs at tau = tBH_form, NSs as set by -NSt), with natal kicks
# compared with the escape velocity of the cluster at that time; its residual gas, if any,
# joins the gas of the cluster and is expelled on the same timescale -tge.
#
# Spatial segregation (optional). Rapster has no radial structure, so by default all
# populations are distributed like the cluster as a whole ("well mixed"): each loses the
# same fraction of its mass by escape and contributes to the central density in proportion
# to its number of stars. A population can instead be born more (or less) concentrated,
# with half-mass radius r_h,k = x_k r_h (x_k = 1: well mixed). Then, for Plummer profiles
# of half-mass radius x_k r_h:
#   - its central density is larger by x_k^(-3): stars in the cluster center (encounters,
#     tidal disruptions) are drawn from population k with weight N_k x_k^(-3), and the
#     central stellar density n_star is multiplied by sum_k (N_k / N) x_k^(-3);
#   - it loses mass by escape in proportion to its mass outside 2 r_h (where escapers come
#     from): weight M_k f_out(x_k), with f_out the Plummer mass fraction outside 2 r_h;
#   - it mixes with the rest of the cluster: x_k relaxes to 1 exponentially on the
#     timescale n_mix,k t_rh (several half-mass relaxation times in N-body models of
#     multiple-population globular clusters, e.g. Vesperini et al. 2013, MNRAS 429, 1913).
# The cluster r_h keeps its own evolution (it is not recomputed from the populations), and
# the BHs follow the cluster as before. These are approximations within Rapster's one-zone
# model; x_k and n_mix,k are choices of the user (README, "Multiple stellar populations").
#
# With a single population the code takes the original expressions exactly (results are
# unchanged, bit for bit); the population bookkeeping starts when a second population is
# born.

import math

import numpy as np

from .constants import *

# ratio of the Plummer scale radius to the half-mass radius: r_h = 1.305 a
_PLUMMER_RH_OVER_A = 1.0 / math.sqrt(2.0**(2.0/3.0) - 1.0)

# escapers are taken to come from outside this many cluster half-mass radii:
ESCAPE_RADIUS_RH = 2.0

# default mixing timescale of a segregated population, in units of the half-mass
# relaxation time (used only for populations born with x != 1):
DEFAULT_N_MIX = 5.0

# columns of populations.txt (one row per population and time step, written when there
# are several populations): time [Myr], index k of the population (0: initial one), its
# birth time [Myr], number of stars, stellar mass [Msun], mean stellar mass [Msun], and
# half-mass radius relative to that of the cluster:
POPULATION_KEYS = ['t', 'k', 't_birth', 'N_k', 'M_k', 'm_avg_k', 'x_k']

# populations born within this time [Myr] of the current time are born now (time steps
# are shortened to land on birth times, up to rounding):
_BIRTH_TOLERANCE = 1e-9


def read_populations_file(path, Z_default):
    """Stellar populations born after t = 0, from a text file (-POP): one population per
    line, with columns

        t_birth  N  [Z  [SFE  [x  [n_mix]]]]

    t_birth: birth time [Myr] (> 0, from the formation of the cluster), N: number of stars,
    Z: metallicity (default: that of the cluster, -Z), SFE: star formation efficiency
    (default 1: no residual gas), x: half-mass radius relative to that of the cluster at
    birth (default 1: well mixed), n_mix: mixing time of a segregated population in
    half-mass relaxation times (default DEFAULT_N_MIX). Text after '#' is a comment.
    Returns a list of dicts, sorted by birth time; raises ValueError on invalid input."""
    names = ('t_birth', 'N', 'Z', 'SFE', 'x', 'n_mix')
    pops = []
    with open(path) as f:
        for i, line in enumerate(f, 1):
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            fields = line.replace(',', ' ').split()
            if not 2 <= len(fields) <= 6:
                raise ValueError(f'{path}, line {i}: expected 2 to 6 columns (t_birth N [Z] [SFE] [x] [n_mix]), got {len(fields)}')
            try:
                values = [float(v) for v in fields]
            except ValueError:
                raise ValueError(f'{path}, line {i}: not a number')
            if not all(math.isfinite(v) for v in values):
                raise ValueError(f'{path}, line {i}: values must be finite')
            p = dict(zip(names, [None, None, float(Z_default), 1.0, 1.0, DEFAULT_N_MIX]))
            p.update(zip(names, values))
            if p['t_birth'] <= 0:
                raise ValueError(f'{path}, line {i}: t_birth must be positive (the first population is born at t = 0)')
            if p['N'] < 1:
                raise ValueError(f'{path}, line {i}: N must be at least 1')
            if p['Z'] < 0:
                raise ValueError(f'{path}, line {i}: Z must be non-negative')
            if not 0 < p['SFE'] <= 1:
                raise ValueError(f'{path}, line {i}: SFE must be in (0, 1]')
            if p['x'] <= 0:
                raise ValueError(f'{path}, line {i}: x must be positive')
            if p['n_mix'] <= 0:
                raise ValueError(f'{path}, line {i}: n_mix must be positive')
            pops.append(p)
    pops.sort(key=lambda p: p['t_birth'])
    return pops


def plummer_mass_fraction_outside(x):
    """Fraction of the mass of a Plummer sphere of half-mass radius x r_h that lies
    outside ESCAPE_RADIUS_RH r_h."""
    y = ESCAPE_RADIUS_RH * _PLUMMER_RH_OVER_A / x
    return 1.0 - y**3 / (y**2 + 1.0)**1.5


def mean_stellar_mass(m_avg0, tau):
    """Mean stellar mass of a population of age tau [Myr] (same law as the cluster's
    <m> in update_cluster): <m>_0 (tau/t_sev)^(-nu_sev) after t_sev."""
    if tau > t_sev:
        return m_avg0 * (tau / t_sev)**(-nu_sev)
    return m_avg0


class StellarPopulations:
    """The stellar populations of a cluster: those already born (lists, index k, k = 0 is
    the initial population) and those scheduled to be born later (sorted by birth time)."""

    def __init__(self, M0, Z, SFE, scheduled=()):
        # born populations:
        self.t_birth = [0.0]
        self.N_birth = [None]        # number of stars at birth (None: initial population, see -N)
        self.M = [float(M0)]         # current stellar mass [Msun]; for k = 0 kept up to date
                                     # only once a second population is born
        self.Z = [float(Z)]
        self.SFE = [float(SFE)]
        self.x = [1.0]               # r_h,k / r_h
        self.n_mix = [DEFAULT_N_MIX]
        # populations still to be born: list of dicts (t_birth, N, Z, SFE, x, n_mix)
        self.scheduled = sorted((dict(p) for p in scheduled), key=lambda p: p['t_birth'])

    # ------------------------------------------------------------------ queries
    @property
    def size(self):
        return len(self.t_birth)

    @property
    def multi(self):
        """True once more than one population has been born."""
        return len(self.t_birth) > 1

    @property
    def segregated(self):
        return any(xk != 1.0 for xk in self.x)

    def next_birth_time(self):
        return self.scheduled[0]['t_birth'] if self.scheduled else math.inf

    def due(self, t):
        """Scheduled populations to be born at time t (removed from the schedule)."""
        out = []
        while self.scheduled and self.scheduled[0]['t_birth'] <= t + _BIRTH_TOLERANCE:
            out.append(self.scheduled.pop(0))
        return out

    def ages(self, t):
        return [t - tk for tk in self.t_birth]

    def mean_masses(self, t, m_avg0):
        return [mean_stellar_mass(m_avg0, t - tk) for tk in self.t_birth]

    def numbers(self, t, m_avg0):
        """Numbers of stars N_k = M_k / <m>_k."""
        return [Mk / mk for Mk, mk in zip(self.M, self.mean_masses(t, m_avg0))]

    def central_weights(self, t, m_avg0):
        """Weights of the populations among the stars of the cluster center: N_k x_k^-3."""
        return [Nk / xk**3 for Nk, xk in zip(self.numbers(t, m_avg0), self.x)]

    def density_factor(self, t, m_avg0):
        """Factor sum_k (N_k/N) x_k^-3 of the central stellar density (1 if well mixed)."""
        Ns = self.numbers(t, m_avg0)
        N = sum(Ns)
        if N <= 0:
            return 1.0
        return sum(Nk / xk**3 for Nk, xk in zip(Ns, self.x)) / N

    def escape_weights(self):
        """Shares of the escaping mass: M_k f_out(x_k), normalized (= M_k / M if well mixed)."""
        if self.segregated:
            w = [Mk * plummer_mass_fraction_outside(xk) for Mk, xk in zip(self.M, self.x)]
        else:
            w = list(self.M)
        total = sum(w)
        return [wk / total for wk in w] if total > 0 else [0.0 for _ in w]

    def rows(self, t, m_avg0, Mcl):
        """Rows of populations.txt at time t (columns POPULATION_KEYS). Before a second
        population is born, the initial population is the whole cluster (stellar mass Mcl)."""
        M = list(self.M) if self.multi else [float(Mcl)]
        out = []
        for k, (tk, Mk, xk) in enumerate(zip(self.t_birth, M, self.x)):
            mk = mean_stellar_mass(m_avg0, t - tk)
            out.append([t, k, tk, Mk / mk, Mk, mk, xk])
        return out

    # ------------------------------------------------------------------ updates
    def add(self, t_birth, N, M, Z, SFE, x, n_mix):
        self.t_birth.append(float(t_birth))
        self.N_birth.append(float(N))
        self.M.append(float(M))
        self.Z.append(float(Z))
        self.SFE.append(float(SFE))
        self.x.append(float(x))
        self.n_mix.append(float(n_mix))

    def evolve(self, t, dt, m_avg0, dM_rlx, t_rlx):
        """Advance the populations from t to t + dt: stellar-evolution mass loss of each
        population at its own age, the escaping mass dM_rlx (< 0) shared according to
        escape_weights, and the mixing of segregated populations. Returns the total
        stellar-evolution mass change sum_k dM_k(sev) (<= 0)."""
        shares = self.escape_weights()
        dM_sev_total = 0.0
        for k in range(len(self.M)):
            tau = t - self.t_birth[k]
            dM_sev = - nu_sev * self.M[k] * dt / tau if tau > t_sev else 0.0
            dM_sev_total += dM_sev
            self.M[k] = max(self.M[k] + dM_sev + shares[k] * dM_rlx, 0.0)
            if self.x[k] != 1.0 and t_rlx > 0:
                self.x[k] = 1.0 + (self.x[k] - 1.0) * math.exp(-dt / (self.n_mix[k] * t_rlx))
        return dM_sev_total

    def white_dwarf_formation_rate(self, t, m_avg0, Kroupa_norm, IMF):
        """Rate [1/Myr] at which white dwarfs form: stars of each population reaching the
        end of their life (turn-off mass), as in evolve_tdes for one population."""
        solar_life = 1.0e4
        m_WN = 8.0
        t_WN = solar_life / m_WN**2.5
        rate = 0.0
        for Nk, tk in zip(self.numbers(t, m_avg0), self.t_birth):
            tau = t - tk
            if tau > t_WN:
                m_to = (solar_life / tau)**(1 / 2.5)
                rate += Nk / tau / 2.5 * Kroupa_norm * IMF(m_to) * m_to
        return rate

# End of file.
