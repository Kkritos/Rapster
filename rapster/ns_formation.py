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

# Neutron-star formation times.
#
# NS formation-time models (-NSt):
#   0 : all NSs form at tBH_form (legacy behaviour; draws no random numbers,
#       so runs are bit-identical to earlier versions)
#   1 : t_form drawn from the marginal distribution p(t), independent of mass
#   2 : t_form drawn from p(t | m_NS) in the NS-mass bin that contains m_NS
#       (NS masses themselves are set by -NS 1/2/3, as before)
#
# The distributions are tabulated in Data/NS_formation.npz, built from a
# binary population-synthesis sample by tools/make_ns_formation_table.py.
# Times are measured from ZAMS (= cluster formation) in Myr. Sampling is by
# inverse transform in log10(t), using quantile tables.

import os
import numpy as np
from zipfile import BadZipFile

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, '..', 'Data')
DEFAULT_NS_FORMATION_TABLE = os.path.join(DATA_DIR, 'NS_formation.npz')

NS_TIME_MODELS = (0, 1, 2)

_REQUIRED_KEYS = ('q', 'logt_marginal', 'm_centers', 'logt_conditional')


class NSFormationTable:
    """Quantile tables of NS formation times, marginal and conditional on mass.

    Arrays (all float):
        q                 (n_q,)        cumulative probabilities, 0 ... 1
        logt_marginal     (n_q,)        log10(t/Myr) quantiles of p(t)
        m_centers         (n_m,)        centres of populated NS-mass bins [Msun]
        m_edges           (n_m, 2)      lower/upper edge of each bin [Msun]
        logt_conditional  (n_m, n_q)    log10(t/Myr) quantiles of p(t | m-bin)
        counts            (n_m,)        number of NSs in each bin
    """

    def __init__(self, path=DEFAULT_NS_FORMATION_TABLE):
        try:
            archive = np.load(path, allow_pickle=False)
        except (BadZipFile, EOFError, OSError, ValueError) as error:
            raise ValueError(f'could not read NS formation-time table {path!r}') from error

        try:
            missing = [k for k in _REQUIRED_KEYS if k not in archive.files]
            if missing:
                raise ValueError(f'NS formation-time table {path!r} is missing {missing}')
            self.q = np.asarray(archive['q'], dtype=float)
            self.logt_marginal = np.asarray(archive['logt_marginal'], dtype=float)
            self.m_centers = np.asarray(archive['m_centers'], dtype=float)
            self.logt_conditional = np.asarray(archive['logt_conditional'], dtype=float)
            self.m_edges = np.asarray(archive['m_edges'], dtype=float) if 'm_edges' in archive.files else None
            self.counts = np.asarray(archive['counts'], dtype=float) if 'counts' in archive.files else None
            self.source = str(archive['source']) if 'source' in archive.files else 'unknown'
        finally:
            archive.close()

        self.path = path
        self._validate()

    def _validate(self):
        n_q = self.q.size
        if self.q.ndim != 1 or n_q < 2 or self.q[0] != 0 or self.q[-1] != 1 or np.any(np.diff(self.q) <= 0):
            raise ValueError('q must be strictly increasing from 0 to 1')
        if self.logt_marginal.shape != (n_q,):
            raise ValueError('logt_marginal must have the same length as q')
        if self.m_centers.ndim != 1 or self.m_centers.size < 1 or np.any(np.diff(self.m_centers) <= 0):
            raise ValueError('m_centers must be a non-empty, strictly increasing array')
        if self.logt_conditional.shape != (self.m_centers.size, n_q):
            raise ValueError('logt_conditional must have shape (len(m_centers), len(q))')
        for name, arr in (('logt_marginal', self.logt_marginal), ('logt_conditional', self.logt_conditional)):
            if not np.all(np.isfinite(arr)):
                raise ValueError(f'{name} must be finite')
            if np.any(np.diff(arr, axis=-1) < 0):
                raise ValueError(f'{name} quantiles must be non-decreasing')

    @property
    def t_min(self):
        return 10**self.logt_marginal[0]

    @property
    def t_max(self):
        return 10**self.logt_marginal[-1]

    def sample_marginal(self, n):
        """Draw n formation times [Myr] from p(t)."""
        u = np.random.rand(int(n))
        return 10**np.interp(u, self.q, self.logt_marginal)

    def mass_bin_index(self, m):
        """Index of the populated mass bin whose centre is nearest to each m.

        Masses outside the tabulated range map to the first/last bin.
        """
        m = np.atleast_1d(np.asarray(m, dtype=float))
        j = np.searchsorted(self.m_centers, m)
        j = np.clip(j, 1, self.m_centers.size - 1) if self.m_centers.size > 1 else np.zeros_like(j)
        if self.m_centers.size > 1:
            left_closer = (m - self.m_centers[j - 1]) <= (self.m_centers[j] - m)
            j = np.where(left_closer, j - 1, j)
        return j

    def sample_conditional(self, m):
        """Draw one formation time [Myr] from p(t | m) for each mass in m."""
        m = np.atleast_1d(np.asarray(m, dtype=float))
        out = np.empty(m.size)
        if m.size == 0:
            return out
        idx = self.mass_bin_index(m)
        u = np.random.rand(m.size)
        # loop over bins (n_m ~ tens), vectorised within each bin:
        for j in np.unique(idx):
            sel = idx == j
            out[sel] = 10**np.interp(u[sel], self.q, self.logt_conditional[j])
        return out


def sample_ns_formation_times(mNS, model, t_default, table=None):
    """Formation times [Myr] for NSs of masses mNS under formation-time model `model`.

    Args:
        mNS (array): NS masses [Msun].
        model (int): 0 (all at t_default), 1 (marginal p(t)), 2 (p(t | m)).
        t_default (float): formation time for model 0 [Myr] (tBH_form).
        table (NSFormationTable, optional): pre-loaded table; loaded from the
            default path when needed and not given.

    Returns:
        array: formation times [Myr], same length as mNS.
    """
    mNS = np.atleast_1d(np.asarray(mNS, dtype=float))
    if model not in NS_TIME_MODELS:
        raise ValueError(f'unknown NS formation-time model {model!r}; choose from {NS_TIME_MODELS}')

    if model == 0:
        return np.full(mNS.size, float(t_default))

    if table is None:
        table = NSFormationTable()

    if model == 1:
        return table.sample_marginal(mNS.size)
    return table.sample_conditional(mNS)

# End of file.
