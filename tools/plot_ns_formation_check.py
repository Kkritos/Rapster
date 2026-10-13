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

# Diagnostic figures for extended NS formation (-NSt):
#   ns_check_a_mass_time.png   sampled (m_NS, t_form) with -NS 3 -NSt 2 vs the COMPAS table
#   ns_check_b_cdf.png         cumulative distribution of NS formation times per model
#   ns_check_c_retention.png   retained NSs vs time in a run vs the analytic expectation
#
# Usage (from the repository root):
#   python3 tools/plot_ns_formation_check.py [output_folder]

import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tests'))

from rapster.constants import tBH_form
from rapster.run_cluster import parse_args
from rapster.cluster_evolution import initialize_cluster
from rapster.ns_formation import NSFormationTable
from test_ns_physics import run_with_expectation      # also initializes the samplers

# Wong (2011) colorblind-safe palette:
BLUE, VERMILION, GREEN, BLACK = '#0072B2', '#D55E00', '#009E73', '#000000'
FIGSIZE = (3.5, 3.06)


def style(ax):
    ax.tick_params(which='both', direction='in', top=True, right=True)


def pending(args):
    sys.argv = ['rapster'] + list(args) + ['-P', '0']
    s = initialize_cluster(parse_args())
    return s['mNS_pending'], s['tNS_pending']


def main(out):
    os.makedirs(out, exist_ok=True)
    table = NSFormationTable()
    q = table.q

    # (a) sampled (m, t) with uniform NS masses and -NSt 2, vs COMPAS conditional quantiles
    m, t = pending(['-N', '2e7', '-S', '3', '-NS', '3', '-NSt', '2'])
    fig, ax = plt.subplots(figsize=FIGSIZE, constrained_layout=True)
    h = ax.hist2d(m, np.log10(t), bins=[60, 60], range=[[1.0, 2.25], [0.6, 2.2]], cmap='Blues', cmin=1)
    cb = fig.colorbar(h[3], ax=ax, pad=0.02)
    cb.set_label('sampled NSs per bin', fontsize=8)
    cb.ax.tick_params(direction='in', labelsize=7)
    med = [np.interp(0.50, q, row) for row in table.logt_conditional]
    lo = [np.interp(0.05, q, row) for row in table.logt_conditional]
    hi = [np.interp(0.95, q, row) for row in table.logt_conditional]
    ax.errorbar(table.m_centers, med, yerr=[np.subtract(med, lo), np.subtract(hi, med)],
                fmt='o', ms=4, color=VERMILION, elinewidth=1, capsize=1.5, label='COMPAS table: median, 5-95%')
    # sampled median in each COMPAS mass bin (only masses inside the bin edges):
    ms_, md_ = [], []
    for j, (a, b) in enumerate(table.m_edges):
        sel = (m >= a) & (m < b)
        if sel.sum() > 100:
            ms_.append(table.m_centers[j]); md_.append(np.median(np.log10(t[sel])))
    ax.plot(ms_, md_, 'x', ms=5, mew=1.2, color=BLACK, label='sampled median')
    ax.set_xlim(1.0, 2.25); ax.set_ylim(0.6, 2.2)
    ax.set_xlabel(r'$m_{\rm NS}\ (M_\odot)$')
    ax.set_ylabel(r'$\log_{10}(t_{\rm form}/{\rm Myr})$')
    ax.set_title(r'-NS 3 -NSt 2 (uniform $m_{\rm NS}$)', fontsize=9)
    ax.legend(fontsize=7, loc='upper right', frameon=False)
    style(ax)
    fig.savefig(os.path.join(out, 'ns_check_a_mass_time.png'), dpi=200)
    plt.close(fig)

    # (b) cumulative distributions of NS formation times
    fig, ax = plt.subplots(figsize=FIGSIZE, constrained_layout=True)
    ax.plot(10**table.logt_marginal, q, color=BLACK, lw=2, label='COMPAS table (all NSs)')
    for args, color, ls, lab in [(['-NS', '2', '-NSt', '1'], BLUE, '--', '-NS 2 -NSt 1'),
                                 (['-NS', '2', '-NSt', '2'], VERMILION, '-', '-NS 2 -NSt 2'),
                                 (['-NS', '3', '-NSt', '2'], GREEN, '-', '-NS 3 -NSt 2')]:
        _, tt = pending(['-N', '5e6', '-S', '5'] + args)
        ax.plot(np.sort(tt), np.arange(1, tt.size + 1) / tt.size, color=color, ls=ls, lw=1.5, label=lab)
    ax.axvline(tBH_form, color='0.5', lw=1, ls=':', label=r'-NSt 0 (all at $t_{\rm BH,form}$)')
    ax.set_xscale('log'); ax.set_xlim(5, 2e2); ax.set_ylim(0, 1)
    ax.set_xlabel(r'$t_{\rm form}\ ({\rm Myr})$')
    ax.set_ylabel('cumulative fraction of NSs')
    ax.legend(fontsize=7, loc='lower right', frameon=False)
    style(ax)
    fig.savefig(os.path.join(out, 'ns_check_b_cdf.png'), dpi=200)
    plt.close(fig)

    # (c) retained NSs vs time, compared with sum_i P_Maxwell(v < v_esc(t_i))
    fig, ax = plt.subplots(figsize=FIGSIZE, constrained_layout=True)
    base = ['-N', '5e5', '-n', '2.65e5', '-S', '7', '-wK', '30', '-NS', '2', '-P', '0']
    for nst, color in [(1, BLUE), (2, VERMILION)]:
        _, mean, var, hh = run_with_expectation(base + ['-NSt', str(nst)], 2000.0)
        tt, ret, ex, sd = hh[:, 0], hh[:, 2], hh[:, 3], np.sqrt(hh[:, 4])
        ax.fill_between(tt, ex - 2 * sd, ex + 2 * sd, color=color, alpha=0.2, lw=0, step='post')
        ax.step(tt, ex, where='post', color=color, lw=1, ls='--')
        ax.step(tt, ret, where='post', color=color, lw=2, label=f'-NSt {nst}')
    sys.argv = ['rapster'] + base + ['-NSt', '0']
    n0 = int(np.sum(initialize_cluster(parse_args())['mBH'] < 3.0))
    ax.axhline(n0, color='0.4', lw=1, ls=':', label=f'-NSt 0: {n0} retained at {tBH_form:g} Myr')
    ax.plot([], [], color='0.3', lw=2, label='solid: run')
    ax.plot([], [], color='0.3', lw=1, ls='--', label=r'dashed, band: expected $\pm2\sigma$')
    ax.set_xscale('log'); ax.set_xlim(5, 2e3); ax.set_ylim(0, 1.6 * n0)
    ax.set_xlabel(r'$t\ ({\rm Myr})$')
    ax.set_ylabel('retained NSs (cumulative)')
    ax.set_title(r'$N=5\times10^5$, $w_{\rm SN}=30$ km/s', fontsize=9)
    ax.legend(fontsize=6.5, loc='upper left', frameon=False)
    style(ax)
    fig.savefig(os.path.join(out, 'ns_check_c_retention.png'), dpi=200)
    plt.close(fig)

    print('figures written to', out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'Results_ns_check')

# End of file.
