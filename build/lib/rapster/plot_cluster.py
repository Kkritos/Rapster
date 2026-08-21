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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os

from .constants import evolution_keys, merger_keys, hardening_keys, tdes_keys

# Unit conversions:
PC_TO_AU = 3.086e+16 / 1.496e+11
MYR_TO_SEC = 1e6 * 365 * 24 * 60 * 60


def load_results(results_dir):
    """Load all output files from a Rapster simulation.

    Args:
        results_dir (str): Path to the Results directory.

    Returns:
        dict: Dictionary with keys 'evolution', 'mergers', 'hardening', 'tdes', 'outBHs'.
            Each value is a pandas DataFrame or dict (for outBHs).
            Missing files are set to None.
    """
    data = {}

    evo_path = os.path.join(results_dir, 'evolution.txt')
    if os.path.exists(evo_path):
        data['evolution'] = pd.DataFrame(np.loadtxt(evo_path), columns=evolution_keys)
    else:
        data['evolution'] = None

    mer_path = os.path.join(results_dir, 'mergers.txt')
    if os.path.exists(mer_path) and os.path.getsize(mer_path) > 0:
        data['mergers'] = pd.DataFrame(np.loadtxt(mer_path), columns=merger_keys)
    else:
        data['mergers'] = None

    hard_path = os.path.join(results_dir, 'hardening.txt')
    if os.path.exists(hard_path) and os.path.getsize(hard_path) > 0:
        data['hardening'] = pd.DataFrame(np.loadtxt(hard_path), columns=hardening_keys)
    else:
        data['hardening'] = None

    tde_path = os.path.join(results_dir, 'tdes.txt')
    if os.path.exists(tde_path) and os.path.getsize(tde_path) > 0:
        data['tdes'] = pd.DataFrame(np.loadtxt(tde_path), columns=tdes_keys)
    else:
        data['tdes'] = None

    bh_path = os.path.join(results_dir, 'outputBHs.pkl')
    if os.path.exists(bh_path):
        with open(bh_path, 'rb') as f:
            data['outBHs'] = pickle.load(f)
    else:
        data['outBHs'] = None

    return data


def plot_cluster_evolution(data, save_dir):
    """Plot number of stars, BHs, BBHs, mergers, and TDEs vs time."""
    evolution = data['evolution']
    if evolution is None:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(evolution['t'], evolution['M_cl'] / evolution['m_avg'], lw=3, label='stars')
    ax.plot(evolution['t'], evolution['N_BH'], lw=2, label='BHs')
    ax.plot(evolution['t'], evolution['N_BBH'], lw=1, label='BBHs')
    ax.plot(evolution['t'], evolution['N_me'], lw=2, ls='--', label='mergers')
    ax.plot(evolution['t'], evolution['N_tdeBHWD'] + evolution['N_tdeBHstar'], lw=2, ls='--', label='TDEs')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 1e4)
    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel('Number')
    ax.set_title('Cluster evolution')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'cluster_evolution.png'), dpi=150)
    plt.close(fig)


def plot_radii_evolution(data, save_dir):
    """Plot half-mass radius, BH half-mass radius, and BH core radius vs time."""
    evolution = data['evolution']
    if evolution is None:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(evolution['t'], evolution['r_h'], lw=3, label='half-mass radius')
    ax.plot(evolution['t'], evolution['r_hBH'], lw=2, label='BH half-mass radius')
    ax.plot(evolution['t'], evolution['r_cBH'], lw=1, label='BH core radius')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 1e4)
    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel('Radius (pc)')
    ax.set_title('Radii evolution')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'radii_evolution.png'), dpi=150)
    plt.close(fig)


def plot_cluster_mass(data, save_dir):
    """Plot cluster mass evolution with time."""
    evolution = data['evolution']
    if evolution is None:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(evolution['t'], evolution['M_cl'], lw=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 1e4)
    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel(r'Cluster mass ($M_\odot$)')
    ax.set_title('Cluster mass evolution')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'cluster_mass.png'), dpi=150)
    plt.close(fig)


def plot_bh_mass_function(data, save_dir):
    """Plot BH mass function at several snapshots."""
    outBHs = data['outBHs']
    if outBHs is None:
        return

    n_snapshots = len(outBHs['t'])
    if n_snapshots < 5:
        return

    # Pick ~5 evenly spaced snapshots after the first few:
    indices = np.linspace(max(1, n_snapshots // 10), n_snapshots - 1, 5, dtype=int)
    bins = np.linspace(0, 80)

    fig, ax = plt.subplots(figsize=(8, 5))
    for idx in indices:
        masses = outBHs['mBH'][idx]
        t_val = outBHs['t'][idx]
        ax.hist(masses, histtype='step', bins=bins, lw=2, label=f't={t_val:.3g} Myr')
    ax.legend()
    ax.set_yscale('log')
    ax.set_xlabel(r'Black-hole mass ($M_\odot$)')
    ax.set_ylabel('Number')
    ax.set_title('BH mass function evolution')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'bh_mass_function.png'), dpi=150)
    plt.close(fig)


def plot_merger_masses(data, save_dir):
    """Plot m1 vs m2 for BBH mergers."""
    mergers = data['mergers']
    if mergers is None or len(mergers) == 0:
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(mergers['m1'], mergers['m2'], s=50, marker='o')
    m_range = np.linspace(0, mergers['m1'].max() * 1.1)
    ax.plot(m_range, m_range, color='k', ls='--', label=r'$m_1=m_2$')
    ax.legend()
    ax.set_xlabel(r'Primary mass, $m_1$ ($M_\odot$)')
    ax.set_ylabel(r'Secondary mass, $m_2$ ($M_\odot$)')
    ax.set_title('BBH merger masses')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'merger_masses.png'), dpi=150)
    plt.close(fig)


def plot_merger_channels(data, save_dir):
    """Plot pie chart of merger formation channels."""
    mergers = data['mergers']
    if mergers is None or len(mergers) == 0:
        return

    channels = mergers['channel']
    cap2 = (channels == 2)
    cap3 = (channels == 6)
    zlk = (channels == 4)
    eject = (channels < 0)
    incluster = ~(cap2 | cap3 | eject | zlk)

    sizes = [
        cap2.sum(),
        cap3.sum(),
        incluster.sum(),
        eject.sum(),
        zlk.sum(),
    ]
    labels = ['2-capture', '3-capture', 'in-cluster', 'ejected', 'Lidov-Kozai']

    # Filter out zero-count channels:
    nonzero = [(s, l) for s, l in zip(sizes, labels) if s > 0]
    if not nonzero:
        return
    sizes, labels = zip(*nonzero)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title('Merger channels')
    ax.pie(sizes, labels=labels, autopct='%1.1f%%')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'merger_channels.png'), dpi=150)
    plt.close(fig)


def plot_tdes(data, save_dir):
    """Plot TDE BH mass vs penetration factor and available time column."""
    tdes = data['tdes']
    if tdes is None or len(tdes) == 0:
        return

    if 't_fb' in tdes:
        y_values = tdes['t_fb'] * MYR_TO_SEC / 60
        y_label = 'Fallback time (min)'
    else:
        y_values = tdes['t']
        y_label = 'Cluster time (Myr)'

    fig, ax = plt.subplots(figsize=(8, 5))
    sc = ax.scatter(tdes['beta'], y_values, c=tdes['m_BH'], cmap='gnuplot')
    plt.colorbar(sc, ax=ax, label=r'Black-hole mass ($M_\odot$)')
    if np.all(y_values > 0):
        ax.set_yscale('log')
    ax.set_xlabel(r'Penetration factor ($r_t/r_p$)')
    ax.set_ylabel(y_label)
    ax.set_title('Tidal disruption events')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'tdes.png'), dpi=150)
    plt.close(fig)


def plot_hardening(data, save_dir):
    """Plot BBH semimajor axis vs time."""
    hardening = data['hardening']
    if hardening is None or len(hardening) == 0:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(hardening['t'] + hardening['t_local'] + hardening['dt_local'],
               hardening['a'] * PC_TO_AU, s=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel('Semimajor axis (AU)')
    ax.set_title('BBH hardening')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'hardening.png'), dpi=150)
    plt.close(fig)


def plot_merger_spins(data, save_dir):
    """Plot histograms of merger spin parameters (chi1, chi2, chiEff, chiRem)."""
    mergers = data['mergers']
    if mergers is None or len(mergers) == 0:
        return

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    axes[0, 0].hist(mergers['chi1'], bins=30, histtype='step', lw=2)
    axes[0, 0].set_xlabel(r'$\chi_1$')
    axes[0, 0].set_ylabel('Number')
    axes[0, 0].set_title('Primary spin')

    axes[0, 1].hist(mergers['chi2'], bins=30, histtype='step', lw=2)
    axes[0, 1].set_xlabel(r'$\chi_2$')
    axes[0, 1].set_ylabel('Number')
    axes[0, 1].set_title('Secondary spin')

    axes[1, 0].hist(mergers['chiEff'], bins=30, histtype='step', lw=2)
    axes[1, 0].set_xlabel(r'$\chi_{\rm eff}$')
    axes[1, 0].set_ylabel('Number')
    axes[1, 0].set_title('Effective spin')

    axes[1, 1].hist(mergers['chiRem'], bins=30, histtype='step', lw=2)
    axes[1, 1].set_xlabel(r'$\chi_{\rm rem}$')
    axes[1, 1].set_ylabel('Number')
    axes[1, 1].set_title('Remnant spin')

    fig.suptitle('Merger spin distributions')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'merger_spins.png'), dpi=150)
    plt.close(fig)


def plot_eccentricity_by_channel(data, save_dir):
    """Plot eccentricity histograms for different merger formation channels."""
    mergers = data['mergers']
    if mergers is None or len(mergers) == 0:
        return

    channel_labels = {
        1: 'Exchange (ch 1)',
        2: '2-body capture (ch 2)',
        3: '3-body binary (ch 3)',
        4: 'Lidov-Kozai (ch 4)',
        5: 'Triple-induced (ch 5)',
        6: '3-body capture (ch 6)',
    }

    # Also handle ejected mergers (negative channels):
    channels_present = sorted(mergers['channel'].unique())

    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.linspace(0, 1, 30)

    for ch in channels_present:
        mask = mergers['channel'] == ch
        if mask.sum() == 0:
            continue
        ch_int = int(ch)
        if ch_int < 0:
            label = f'Ejected (ch {ch_int})'
        else:
            label = channel_labels.get(ch_int, f'Channel {ch_int}')
        ax.hist(mergers['e'][mask], bins=bins, histtype='step', lw=2, label=label)

    ax.legend()
    ax.set_xlabel('Eccentricity')
    ax.set_ylabel('Number')
    ax.set_title('Eccentricity distribution by merger channel')
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, 'eccentricity_by_channel.png'), dpi=150)
    plt.close(fig)


def generate_all_plots(results_dir):
    """Load results and generate all diagnostic plots.

    Saves PNG files to a 'plots/' subdirectory inside results_dir.

    Args:
        results_dir (str): Path to the Results directory.
    """
    plots_dir = os.path.join(results_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    data = load_results(results_dir)

    plot_cluster_evolution(data, plots_dir)
    plot_radii_evolution(data, plots_dir)
    plot_cluster_mass(data, plots_dir)
    plot_bh_mass_function(data, plots_dir)
    plot_merger_masses(data, plots_dir)
    plot_merger_channels(data, plots_dir)
    plot_tdes(data, plots_dir)
    plot_hardening(data, plots_dir)
    plot_merger_spins(data, plots_dir)
    plot_eccentricity_by_channel(data, plots_dir)

    print(f'Plots saved to {plots_dir}/')


# End of file.
