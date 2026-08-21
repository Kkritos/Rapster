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
import pickle
import os

from .constants import merger_keys, evolution_keys, tdes_keys


CHANNEL_NAMES = {
    1: 'Exchange',
    2: '2-body capture',
    3: '3-body binary',
    4: 'Lidov-Kozai',
    5: 'Triple-induced',
    6: '3-body capture',
}


def analyze_cluster(results_dir):
    """Print a summary of simulation results.

    Reports merger statistics (total, in-cluster, ejected, per channel),
    retained BHs after mergers, maximum dynamically-formed BH mass,
    generational BH counts, and TDE summary.

    Args:
        results_dir (str): Path to the Results directory.
    """
    print('=' * 70)
    print('CLUSTER ANALYSIS SUMMARY')
    print('=' * 70)

    # --- Mergers ---
    mergers_path = os.path.join(results_dir, 'mergers.txt')
    if os.path.exists(mergers_path) and os.path.getsize(mergers_path) > 0:
        mergers = pd.DataFrame(np.loadtxt(mergers_path), columns=merger_keys)
        _print_merger_summary(mergers)
    else:
        print('\nNo mergers recorded.')

    # --- BH populations ---
    bh_path = os.path.join(results_dir, 'outputBHs.pkl')
    if os.path.exists(bh_path):
        with open(bh_path, 'rb') as f:
            outBHs = pickle.load(f)
        _print_bh_summary(outBHs)
    else:
        print('\nNo BH output file found.')

    # --- TDEs ---
    tdes_path = os.path.join(results_dir, 'tdes.txt')
    if os.path.exists(tdes_path) and os.path.getsize(tdes_path) > 0:
        tdes = pd.DataFrame(np.loadtxt(tdes_path), columns=tdes_keys)
        _print_tde_summary(tdes)
    else:
        print('\nNo TDEs recorded.')

    # --- Evolution ---
    evo_path = os.path.join(results_dir, 'evolution.txt')
    if os.path.exists(evo_path):
        evolution = pd.DataFrame(np.loadtxt(evo_path), columns=evolution_keys)
        _print_evolution_summary(evolution)

    print('\n' + '=' * 70)


def _print_merger_summary(mergers):
    """Print merger statistics."""
    N_total = len(mergers)
    in_cluster = mergers[mergers['channel'] > 0]
    ejected = mergers[mergers['channel'] < 0]
    retained = mergers[mergers['vGW'] < mergers['v_esc']]

    print('\n--- MERGERS ---')
    print(f'  Total mergers:        {N_total}')
    print(f'  In-cluster mergers:   {len(in_cluster)}')
    print(f'  Ejected mergers:      {len(ejected)}')
    print(f'  Retained after merger (vGW < v_esc): {len(retained)}')

    # Per-channel breakdown:
    print('\n  Mergers by channel:')
    channels = sorted(mergers['channel'].unique())
    for ch in channels:
        count = (mergers['channel'] == ch).sum()
        ch_int = int(ch)
        if ch_int < 0:
            name = f'Ejected (ch {ch_int})'
        else:
            name = CHANNEL_NAMES.get(ch_int, f'Channel {ch_int}')
        print(f'    {name:25s}: {count}')

    # Maximum dynamically-formed BH mass:
    if len(mergers) > 0:
        max_rem_mass = mergers['mRem'].max()
        max_row = mergers.loc[mergers['mRem'].idxmax()]
        print(f'\n  Max remnant mass (dynamical): {max_rem_mass:.2f} Msun')
        print(f'    from m1={max_row["m1"]:.2f} + m2={max_row["m2"]:.2f} Msun '
              f'(channel {int(max_row["channel"])}, gen {int(max_row["gRem"])})')

    # Mass and spin ranges:
    if len(mergers) > 0:
        print(f'\n  Primary mass range:   [{mergers["m1"].min():.2f}, {mergers["m1"].max():.2f}] Msun')
        print(f'  Secondary mass range: [{mergers["m2"].min():.2f}, {mergers["m2"].max():.2f}] Msun')
        print(f'  chiEff range:         [{mergers["chiEff"].min():.4f}, {mergers["chiEff"].max():.4f}]')
        print(f'  Remnant spin range:   [{mergers["chiRem"].min():.4f}, {mergers["chiRem"].max():.4f}]')


def _print_bh_summary(outBHs):
    """Print BH population summary from the final snapshot."""
    if len(outBHs['t']) == 0:
        print('\nNo BH snapshots available.')
        return

    # Use the final snapshot:
    final_masses = outBHs['mBH'][-1]
    final_spins = outBHs['sBH'][-1]
    final_gens = outBHs['gBH'][-1]
    t_final = outBHs['t'][-1]

    print(f'\n--- BH POPULATION (t = {t_final:.1f} Myr) ---')
    print(f'  Total BHs remaining:  {len(final_masses)}')
    if len(final_masses) > 0:
        print(f'  Mass range:           [{final_masses.min():.2f}, {final_masses.max():.2f}] Msun')
        print(f'  Mean mass:            {final_masses.mean():.2f} Msun')

    # Generation breakdown:
    print('\n  BHs by generation:')
    unique_gens = sorted(np.unique(final_gens))
    for g in unique_gens:
        count = (final_gens == g).sum()
        gen_masses = final_masses[final_gens == g]
        print(f'    Generation {int(g):d}: {count:5d} BHs  '
              f'(mass range [{gen_masses.min():.2f}, {gen_masses.max():.2f}] Msun)')


def _print_tde_summary(tdes):
    """Print TDE summary statistics."""
    N_total = len(tdes)
    bh_wd = tdes[tdes['type'] == 11]
    bh_star = tdes[tdes['type'] == 1]

    print(f'\n--- TIDAL DISRUPTION EVENTS ---')
    print(f'  Total TDEs:           {N_total}')
    print(f'  BH-WD TDEs:          {len(bh_wd)}')
    print(f'  BH-star TDEs:        {len(bh_star)}')
    if N_total > 0:
        print(f'  BH mass range:       [{tdes["m_BH"].min():.2f}, {tdes["m_BH"].max():.2f}] Msun')
        print(f'  Penetration range:   [{tdes["beta"].min():.2f}, {tdes["beta"].max():.2f}]')


def _print_evolution_summary(evolution):
    """Print final cluster state."""
    final = evolution.iloc[-1]
    print(f'\n--- FINAL CLUSTER STATE ---')
    print(f'  Time:                 {final["t"]:.1f} Myr (z = {final["z"]:.3f})')
    print(f'  Cluster mass:         {final["M_cl"]:.0f} Msun')
    print(f'  Half-mass radius:     {final["r_h"]:.3f} pc')
    print(f'  Galactocentric R:     {final["R_gal"]:.0f} pc')
    print(f'  N_BH:                 {int(final["N_BH"])}')
    print(f'  N_BBH:                {int(final["N_BBH"])}')
    print(f'  N_Triples:            {int(final["N_triples"])}')


# End of file.
