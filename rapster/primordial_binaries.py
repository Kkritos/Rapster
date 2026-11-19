'''
Primordial BH binaries (channel 11): a fraction pbf of the single BHs, at the BH
formation time, is paired up with semi-major axes drawn from bbh_df.csv
(column SemiMajorAxis_AU).
'''

import os
import numpy as np
import pandas as pd

_dir = os.path.dirname(os.path.abspath(__file__))
_sma_cache = {}


def _sma_samples_AU():
    """Semi-major axes [AU] from bbh_df.csv (next to this file, or in the working directory)."""
    if 'sma' not in _sma_cache:
        for path in (os.path.join(_dir, 'bbh_df.csv'), 'bbh_df.csv'):
            if os.path.exists(path):
                _sma_cache['sma'] = pd.read_csv(path)['SemiMajorAxis_AU'].values
                break
        else:
            raise FileNotFoundError('bbh_df.csv (primordial binaries, -pbf > 0) not found in the package folder or the working directory')
    return _sma_cache['sma']


def extract_primordial_binaries(BH, binaries, pbf, R_sun, G_Newton, mBH_avg_ini, vBH_ini, tBH_form, z_i_0):
    """
    Remove pairs of single BHs from the pool BH and append them to binaries as
    primordial binaries. Returns (binaries, number of binaries created).
    """
    if BH.size < 2:
        return binaries, 0

    num_prim_bi = min(int(pbf * BH.size) // 2, BH.size // 2)
    if num_prim_bi <= 0:
        return binaries, 0

    ind_prim = np.random.choice(BH.size, size=num_prim_bi * 2, replace=False)
    sma_samples_rsun = _sma_samples_AU() * 215.032 * R_sun
    m, s, g, h = BH.m.copy(), BH.s.copy(), BH.g.copy(), BH.h.copy()

    used = []
    new_rows = []
    for i in range(0, len(ind_prim), 2):
        i1, i2 = int(ind_prim[i]), int(ind_prim[i + 1])
        m1, m2 = m[i1], m[i2]

        aBHb_min = 2 * R_sun
        aBHb_max = G_Newton * (m1 + m2 + mBH_avg_ini) * m1 * m2 / \
                   ((m1 + m2) * mBH_avg_ini * 2 * vBH_ini**2)
        if aBHb_max <= aBHb_min:
            continue

        a_prim = np.clip(np.random.choice(sma_samples_rsun), aBHb_min, aBHb_max)
        # row layout: id, channel, a, e, m1, m2, s1, s2, g1, g2, t, z, 0, h1, h2
        new_rows.append([np.random.randint(0, 999999999), 11, a_prim, 0.0,
                         m1, m2, s[i1], s[i2], g[i1], g[i2], tBH_form, z_i_0, 0, h[i1], h[i2]])
        used += [i1, i2]

    if used:
        binaries = np.append(binaries, np.array(new_rows), axis=0)
        BH.remove_many(sorted(used, reverse=True))
    return binaries, len(used) // 2
