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

# Checkpoints and restarts.
#
# A run writes checkpoint.pkl to its results folder every -CKdt Myr of simulated time
# and when it ends by reaching t_max (or another condition of the main loop: -fCO, the
# cluster's mass). A checkpoint holds the whole state of the simulation at the start of
# a time step (cluster, compact objects, binaries, populations, counters, the state of
# the random number generator, the configuration), except the rows of the output files:
# those written so far are flushed to the files of the run, and the checkpoint keeps
# only the rows that may still change, with the number of rows written before them.
#
# A restart (-RS <folder or checkpoint file>) continues from a checkpoint in a new
# results folder (-RF), whose files hold the rows from the checkpoint on. Its restart.json
# says how many rows of each file come before (in the folder restarted from), and
# join_restart.py joins the files of a run and its restarts. A restart without changes
# continues the run exactly: the joined files are the same, byte for byte, as those of a
# run without interruption. A restart may change the maximum time, the seed of the random
# numbers, the half-mass radius, and add BHs or a new stellar population (README,
# "Checkpoints and restarts").

import json
import os
import pickle

import numpy as np

CHECKPOINT_FILE = 'checkpoint.pkl'
RESTART_INFO_FILE = 'restart.json'
CHECKPOINT_FORMAT = 1

# state entries kept outside the checkpoint's state: the records (rows of output files)
# and the per-step BH snapshots (-BOi 1), which start again in the restarted run:
RECORD_KEYS = ('mergers', 'evolution', 'hardening', 'tdes', 'populations')
SNAPSHOT_KEYS = ('simulation_times', 'black_hole_masses', 'black_hole_spins',
                 'black_hole_generations', 'black_hole_tdes')

# options that a restart may set (the others are those of the run restarted from):
RESTART_OPTIONS = {
    'restart_from', 'results_folder_name', 'maximum_time', 'seed', 'half_mass_radius',
    'print_information', 'generate_plots', 'analyze_results', 'checkpoint_indicator',
    'checkpoint_interval', 'max_compact_object_fraction', 'continue_at_galactic_center',
    'populations_file', 'add_black_holes_file', 'add_number', 'add_metallicity',
    'add_star_formation_efficiency', 'add_radius_ratio', 'add_mixing_time',
}


def checkpoint_path(path):
    """The checkpoint file of a results folder (or the file itself)."""
    return os.path.join(path, CHECKPOINT_FILE) if os.path.isdir(path) else path


def save_checkpoint(state, config, version):
    """Write checkpoint.pkl in the results folder of the run (atomically: a run stopped
    while writing keeps its previous checkpoint). Writes the rows of the output files
    that can no longer change to their files first (StreamedRowLog.sync)."""
    results_dir = os.path.join(os.getcwd(), config['results_folder_name'])
    records = {}
    for key in RECORD_KEYS:
        record = state.get(key)
        records[key] = None if record is None else record.sync()
    data = {
        'format': CHECKPOINT_FORMAT,
        'version': version,
        'config': config,
        'state': {k: v for k, v in state.items() if k not in RECORD_KEYS and k not in SNAPSHOT_KEYS},
        'records': records,
        'n_snapshots': len(state['simulation_times']) + state.get('n_snapshots_before', 0),
        'rng': np.random.get_state(),
    }
    path = os.path.join(results_dir, CHECKPOINT_FILE)
    tmp = path + '.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)
    print('CHECKPOINT WRITTEN AT t = {:.6g} Myr'.format(state['t']))


def load_checkpoint(path):
    path = checkpoint_path(path)
    with open(path, 'rb') as f:
        data = pickle.load(f)
    if not isinstance(data, dict) or data.get('format') != CHECKPOINT_FORMAT:
        raise ValueError(f'{path} is not a Rapster checkpoint (format {CHECKPOINT_FORMAT})')
    return data


def rows_before(records):
    """Number of data rows of each output file that are in the folder restarted from."""
    return {key: (None if rec is None else max(int(rec[0]) - 1, 0)) for key, rec in records.items()}


def write_restart_info(results_dir, info):
    with open(os.path.join(results_dir, RESTART_INFO_FILE), 'w') as f:
        json.dump(info, f, indent=2)
        f.write('\n')


def read_restart_info(results_dir):
    with open(os.path.join(results_dir, RESTART_INFO_FILE)) as f:
        return json.load(f)

# End of file.
