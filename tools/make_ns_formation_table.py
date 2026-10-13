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

# Build Data/NS_formation.npz (the NS formation-time table read by
# rapster/ns_formation.py) from a population-synthesis sample of NSs.
#
# Usage (real data):
#   python3 tools/make_ns_formation_table.py ns_formation_fiducial.npz \
#       --mass-key remnant_mass --time-key sn_time -o Data/NS_formation.npz
#
# Usage (synthetic test sample, no input file):
#   python3 tools/make_ns_formation_table.py --synthetic 2000000 -o Data/NS_formation.npz
#
# Run without --mass-key/--time-key to list the arrays in the input file.

import argparse
import numpy as np


def quantile_grid(n_uniform=2001, n_tail=60, q_tail_min=1e-6):
    """Cumulative-probability grid: uniform in q, refined in both tails."""
    tail = np.logspace(np.log10(q_tail_min), -3, n_tail)
    q = np.unique(np.concatenate(([0.0, 1.0], np.linspace(0, 1, n_uniform), tail, 1 - tail)))
    return q


def synthetic_sample(n, seed=1):
    """Toy NS sample with a COMPAS-like two-branch mass-time structure (testing only)."""
    rng = np.random.default_rng(seed)
    # ZAMS masses from a Salpeter-like slope in [8, 25] Msun:
    a = -2.3
    u = rng.random(n)
    M = (u * (25**(a + 1) - 8**(a + 1)) + 8**(a + 1))**(1 / (a + 1))
    # lifetime ~ 10 Myr (M/20)^-2.5; 30% on a delayed branch (x2.8):
    t = 10 * (M / 20)**(-2.5) * np.where(rng.random(n) < 0.3, 2.8, 1.0)
    t *= 10**rng.normal(0, 0.03, n)
    # heavier progenitors -> heavier NSs, floor at 1.26 Msun:
    m = 1.26 + 1.2 * ((M - 8) / 17)**1.5 + rng.normal(0, 0.02, n)
    m = np.clip(m, 1.26, 2.5)
    # 0.5% late (~Gyr) channel at the ECSN mass:
    late = rng.random(n) < 0.005
    m[late] = 1.26
    t[late] = 10**rng.uniform(np.log10(2e3), np.log10(1.4e4), late.sum())
    return m, t


def build_table(m, t, dm=0.05, min_count=1000, q=None):
    """Quantile tables of log10(t): marginal and conditional on fixed-width mass bins.

    Bins with fewer than min_count NSs are dropped; at sampling time their
    masses map to the nearest populated bin.
    """
    if q is None:
        q = quantile_grid()

    logt = np.log10(t)
    logt_marginal = np.quantile(logt, q)

    lo = np.floor(m.min() / dm) * dm
    hi = np.ceil(m.max() / dm) * dm
    edges = np.arange(lo, hi + dm / 2, dm)
    if edges.size < 2:
        edges = np.array([lo, lo + dm])
    ibin = np.clip(np.digitize(m, edges) - 1, 0, edges.size - 2)

    centers, bin_edges, counts, cond = [], [], [], []
    for j in range(edges.size - 1):
        sel = ibin == j
        nj = int(sel.sum())
        if nj < min_count:
            continue
        centers.append(np.mean(m[sel]))          # mass-weighted centre of the bin's NSs
        bin_edges.append((edges[j], edges[j + 1]))
        counts.append(nj)
        cond.append(np.quantile(logt[sel], q))

    if not centers:
        raise ValueError('no mass bin has at least min_count NSs; lower --min-count or raise --dm')

    return {
        'q': q,
        'logt_marginal': logt_marginal,
        'm_centers': np.array(centers),
        'm_edges': np.array(bin_edges),
        'counts': np.array(counts, dtype=float),
        'logt_conditional': np.array(cond),
    }


def main():
    p = argparse.ArgumentParser(description='Build the Rapster NS formation-time table')
    p.add_argument('input', nargs='?', help='input .npz with NS masses and formation times')
    p.add_argument('--mass-key', help='name of the NS-mass array [Msun]')
    p.add_argument('--time-key', help='name of the formation-time array (from ZAMS)')
    p.add_argument('--time-unit', type=float, default=1.0, help='multiply times by this to get Myr (default 1: input in Myr)')
    p.add_argument('--synthetic', type=int, default=0, help='ignore input; build from N synthetic NSs (testing only)')
    p.add_argument('--dm', type=float, default=0.05, help='NS-mass bin width [Msun] (default 0.05)')
    p.add_argument('--min-count', type=int, default=1000, help='minimum NSs per mass bin (default 1000)')
    p.add_argument('-o', '--output', default='Data/NS_formation.npz', help='output table')
    args = p.parse_args()

    if args.synthetic > 0:
        m, t = synthetic_sample(args.synthetic)
        source = f'synthetic test sample (N={args.synthetic}); NOT physical'
    else:
        if args.input is None:
            p.error('give an input .npz or --synthetic N')
        d = np.load(args.input, allow_pickle=False)
        if args.mass_key is None or args.time_key is None:
            print('Arrays in', args.input)
            for k in d.files:
                a = d[k]
                print(f'  {k:24s} {str(a.shape):16s} {a.dtype}')
            p.error('choose --mass-key and --time-key from the list above')
        m = np.asarray(d[args.mass_key], dtype=float)
        t = np.asarray(d[args.time_key], dtype=float) * args.time_unit
        if m.shape != t.shape or m.ndim != 1:
            p.error('mass and time arrays must be 1D and of equal length')
        source = f'{args.input} [{args.mass_key}, {args.time_key}]'

    good = np.isfinite(m) & np.isfinite(t) & (m > 0) & (t > 0)
    if not good.all():
        print(f'dropping {np.sum(~good)} NSs with non-finite or non-positive mass/time')
    m, t = m[good], t[good]

    table = build_table(m, t, dm=args.dm, min_count=args.min_count)
    n_cond = int(table['counts'].sum())

    np.savez_compressed(args.output, source=np.array(source), **table)

    print(f'source: {source}')
    print(f'N_NS = {m.size}, mass range [{m.min():.3f}, {m.max():.3f}] Msun, '
          f't range [{t.min():.3g}, {t.max():.3g}] Myr')
    print(f'marginal t quantiles [Myr]: 5%={np.quantile(t, .05):.3g}  50%={np.median(t):.3g}  95%={np.quantile(t, .95):.3g}')
    print(f'{table["m_centers"].size} populated mass bins (dm={args.dm}), '
          f'containing {n_cond} NSs ({100 * n_cond / m.size:.2f}%)')
    print(f'written: {args.output}')


if __name__ == '__main__':
    main()

# End of file.
