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

# Join the output files of a run and of its restarts into one results folder:
#
#     python -m rapster.join_restart -o Joined Results Results_restart [Results_restart2 ...]
#
# The folders are given in order: each one after the first is a restart (-RS) of the
# previous one. Its restart.json says how many rows of each file come from the folders
# before it: the joined file has those rows of the previous folder (its finished file, or
# the partial file .<name>.txt.part of a run that was stopped), then the rows of the
# restart. Merger lines end with the cluster properties at the end of the last run. With
# -BOi 1, the BH snapshots of outputBHs.pkl are joined the same way (all folders must have
# the file). log.txt files are concatenated.

import argparse
import os
import pickle
import shutil
import sys

from .restart import read_restart_info
from .snapshots import dump_streaming

MERGER_COLUMNS = 27


class _SnapshotList:
    """A list of arrays pickled as SnapshotStream pickles them (list() plus its items),
    so that a joined outputBHs.pkl has the same bytes as that of an uninterrupted run."""

    def __init__(self, items):
        self.items = items

    def __reduce__(self):
        return (list, (), None, iter(self.items))


def _data_file(folder, name):
    """Path of an output file of a folder: name.txt, or the partial file of a stopped run."""
    for path in (os.path.join(folder, name + '.txt'), os.path.join(folder, '.' + name + '.txt.part')):
        if os.path.isfile(path):
            return path
    return None


def _join_text(prev_path, new_path, n_before, out_path, merger_suffix=None):
    with open(out_path, 'w') as out:
        header = None
        if prev_path is not None and n_before:
            with open(prev_path) as f:
                header = f.readline()
                out.write(header)
                n = 0
                for line in f:
                    if n == n_before:
                        break
                    if merger_suffix is not None:
                        line = ' '.join(line.split()[:MERGER_COLUMNS]) + merger_suffix + '\n'
                    out.write(line)
                    n += 1
            if n < n_before:
                raise ValueError(f'{prev_path} has {n} rows, {n_before} expected before the restart')
        if new_path is not None:
            with open(new_path) as f:
                first = f.readline()
                if header is None:
                    out.write(first)
                shutil.copyfileobj(f, out)


def join(folders, output):
    if len(folders) < 2:
        raise ValueError('give a run and at least one restart of it')
    if os.path.exists(output) and os.listdir(output):
        raise ValueError(f'the output folder {output} exists and is not empty')
    infos = [read_restart_info(f) for f in folders[1:]]
    for prev, info in zip(folders, infos):
        source = os.path.basename(os.path.dirname(info['restarted_from']))
        if source != os.path.basename(os.path.normpath(prev)):
            print(f'join_restart: warning: {prev} is not the folder restarted from ({info["restarted_from"]}); '
                  'check the order of the folders', file=sys.stderr)
    os.makedirs(output, exist_ok=True)
    names = infos[-1]['file_names']
    suffix = infos[-1].get('merger_suffix')

    prev = folders[0]
    tmp_dirs = []
    for i, (folder, info) in enumerate(zip(folders[1:], infos)):
        last = i == len(infos) - 1
        dest = output if last else os.path.join(output, '.join_tmp{}'.format(i))
        os.makedirs(dest, exist_ok=True)
        if not last:
            tmp_dirs.append(dest)
        for key in ('mergers', 'evolution', 'hardening', 'tdes', 'populations'):
            name = names[key]
            prev_path = _data_file(prev, name)
            new_path = _data_file(folder, name)
            n_before = info['rows_before'].get(key) or 0
            if prev_path is None and new_path is None:
                continue
            if prev_path is None and n_before:
                raise ValueError(f'{name}.txt missing in {prev}')
            _join_text(prev_path, new_path, n_before, os.path.join(dest, name + '.txt'),
                       merger_suffix=suffix if key == 'mergers' else None)
        # BH snapshots (-BOi 1):
        bh_prev = os.path.join(prev, names['outputBHs'] + '.pkl')
        bh_new = os.path.join(folder, names['outputBHs'] + '.pkl')
        if os.path.isfile(bh_prev) and os.path.isfile(bh_new):
            with open(bh_prev, 'rb') as f:
                a = pickle.load(f)
            with open(bh_new, 'rb') as f:
                b = pickle.load(f)
            n = info['snapshots_before']
            if len(a['t']) < n:
                raise ValueError(f'{bh_prev} has {len(a["t"])} snapshots, {n} expected before the restart')
            joined = {k: (list(a[k][:n]) + list(b[k]) if k == 't' else _SnapshotList(list(a[k][:n]) + list(b[k])))
                      for k in b}
            with open(os.path.join(dest, names['outputBHs'] + '.pkl'), 'wb') as f:
                dump_streaming(joined, f)
        elif os.path.isfile(bh_prev) or os.path.isfile(bh_new):
            missing = prev if not os.path.isfile(bh_prev) else folder
            print(f'join_restart: warning: {names["outputBHs"]}.pkl not joined: missing in {missing} '
                  '(a run stopped abruptly writes none)', file=sys.stderr)
        # logs:
        with open(os.path.join(dest, 'log.txt'), 'w') as out:
            for path in (os.path.join(prev, 'log.txt'), os.path.join(folder, 'log.txt')):
                if os.path.isfile(path):
                    with open(path) as f:
                        shutil.copyfileobj(f, out)
        prev = dest

    for d in tmp_dirs:
        shutil.rmtree(d)


def main():
    parser = argparse.ArgumentParser(description='Join the output files of a Rapster run and of its restarts')
    parser.add_argument('folders', nargs='+', help='results folders: the run, then its restarts in order')
    parser.add_argument('-o', '--output', required=True, help='new folder for the joined files')
    args = parser.parse_args()
    try:
        join(args.folders, args.output)
    except (ValueError, OSError, KeyError) as err:
        parser.error(str(err))
    print('Joined files written to', args.output)


if __name__ == '__main__':
    main()

# End of file.
