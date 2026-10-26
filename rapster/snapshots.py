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

# Disk-backed storage for the per-step BH snapshots written to outputBHs.pkl (-BOi 1).
#
# The snapshots hold every compact object at every step. For heavy clusters they
# reach tens of GB (e.g. ~20 GB for N = 2.5e7, r_h = 0.4 pc: ~8e4 BHs x 1.1e4 steps
# x 3 arrays), which used to be kept in RAM and then pickled. Here each snapshot is
# appended to a temporary file as soon as it is taken, and outputBHs.pkl is written
# by streaming the snapshots back one at a time. The file loads exactly as before:
# a dict whose 'mBH', 'sBH', 'gBH' values are plain lists of the same float arrays.

import os
import pickle
import tempfile

import numpy as np


class SnapshotStream:
    """Append-only, disk-backed list of NumPy arrays.

    Supports append, len, iteration and indexing like a list. When pickled it is
    written as an ordinary list (list() plus its items), so pickle.load returns a
    plain list of arrays and readers of the output file need no changes.
    """

    def __init__(self, directory=None):
        self._directory = directory if directory and os.path.isdir(directory) else None
        self._file = None       # temporary file, created at the first append
        self._offsets = []      # byte offset of every stored array

    def append(self, array):
        if self._file is None:
            # deleted automatically when closed or when the process exits:
            try:
                self._file = tempfile.TemporaryFile(dir=self._directory, prefix='rapster_snapshots_')
            except OSError:
                # results folder not writable: use the system temporary directory
                self._file = tempfile.TemporaryFile(prefix='rapster_snapshots_')
        self._file.seek(0, os.SEEK_END)
        self._offsets.append(self._file.tell())
        np.save(self._file, np.asarray(array), allow_pickle=False)

    def __len__(self):
        return len(self._offsets)

    def _read(self, offset):
        self._file.seek(offset)
        return np.load(self._file, allow_pickle=False)

    def __iter__(self):
        for offset in list(self._offsets):
            yield self._read(offset)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self._read(o) for o in self._offsets[index]]
        return self._read(self._offsets[index])

    def __reduce__(self):
        # pickled as list() followed by its items, read back lazily from disk:
        return (list, (), None, iter(self))

    def close(self):
        if self._file is not None:
            self._file.close()
            self._file = None
            self._offsets = []


def dump_streaming(obj, file):
    """pickle.dump(obj, file) without keeping every pickled array alive.

    The pickler's memo normally holds a reference to every object written, which
    would pull all streamed snapshots back into memory. Fast mode disables the memo;
    it is safe here because the output is a tree (dict -> lists -> arrays) with no
    shared or self-referencing objects.
    """
    pickler = pickle.Pickler(file)
    pickler.fast = True
    pickler.dump(obj)

# End of file.
