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

# Reuse of elementwise functions of the BH-mass array between calls.
#
# Several samplers need powers, logarithms and exponentials of every single-BH
# mass (e.g. mBH**(-2/5), exp(4 log mBH)). Between two calls typically only a few
# BHs have been added, removed or changed, out of ~1e4-1e5. ElementwiseCache keeps
# the previous mass array and the function values, matches the new array against it
# in long runs of identical entries, and evaluates the functions only on entries
# that are new or changed. Values are reused only where the masses are verified to
# be equal, and the functions act on each entry independently, so the result is
# identical, bit for bit, to evaluating them on the whole new array
# (tests/test_exact_samplers.py checks this for the functions used).
#
# RowLog stores the hardening and merger records compactly (array of doubles).

from array import array

import os
import numpy as np


def _align(new, old, max_segments=8):
    """
    Pieces that rebuild `new` from runs of identical entries of `old`:
    a list of (reused, start, stop), meaning old[start:stop] when reused is True,
    and new[start:stop] (to be evaluated) when it is False. Pieces follow the order
    of `new`. Runs of equal entries are found by direct comparison, so any piece
    marked reused is exactly equal to the corresponding part of `new`; the search
    only assumes that entries are mostly removed, appended at the end, or changed
    in place, which keeps the number of pieces small.
    """
    n1 = new.shape[0]
    n0 = old.shape[0]
    pieces = []
    i = j = 0    # current positions in new and old
    for _ in range(max_segments):
        L = min(n1 - i, n0 - j)
        if L <= 0:
            break
        diff = new[i:i+L] != old[j:j+L]
        d = int(diff.argmax())
        if not diff[d]:                    # equal over the whole overlap
            pieces.append((True, j, j + L))
            i += L
            j += L
            break
        if d:
            pieces.append((True, j, j + d))
            i += d
            j += d
        # new[i] differs from old[j]:
        v = new[i]
        if j + 1 < n0 and old[j+1] == v:     # old[j] was removed
            j += 1
        elif j + 2 < n0 and old[j+2] == v:   # old[j], old[j+1] were removed
            j += 2
        else:                                # new[i] is a new or changed value
            pieces.append((False, i, i + 1))
            i += 1
            j += 1
    if i < n1:
        pieces.append((False, i, n1))

    # merge neighbouring pieces of the same kind that are contiguous (in old for
    # reused pieces, in new for evaluated ones):
    merged = []
    for p in pieces:
        if merged and merged[-1][0] == p[0] and merged[-1][2] == p[1]:
            merged[-1] = (p[0], merged[-1][1], p[2])
        else:
            merged.append(p)
    return merged


class ElementwiseCache:
    """
    Values f(x) of elementwise functions f for a 1-D float array x, reused from
    the previous call wherever the entries of x are unchanged.

    `funcs` must act on every entry independently (NumPy ufunc expressions such as
    x**p or np.exp(c*np.log(x))), so that f(x[a:b]) == f(x)[..., a:b] exactly;
    a function may return several values per entry, stacked along the first axes
    (the entries along the last axis). The returned arrays are read-only and must
    not be modified by the caller. For arrays shorter than `min_size` the functions
    are evaluated directly (cheaper than matching entries).
    """

    def __init__(self, *funcs, max_segments=8, min_size=0):
        self.funcs = funcs
        self.max_segments = max_segments
        self.min_size = min_size
        self.x = None
        self.x_object = None
        self.values = None
        self._zeroed = None    # (index, saved values) while get_zeroed's arrays are in use

    def _restore(self):
        """Undo the zeroing done by get_zeroed (before the stored values are reused)."""
        if self._zeroed is not None:
            index, saved = self._zeroed
            for v, value in zip(self.values, saved):
                v.flags.writeable = True
                v[..., index] = value
                v.flags.writeable = False
            self._zeroed = None

    def get_zeroed(self, x, index):
        """
        As get(x)[0], but with entry `index` of every returned array set to 0.0,
        without copying the arrays: the stored arrays are modified and restored
        at the next call. The returned arrays are valid only until then.
        """
        values, _ = self.get(x)
        if self.values is not values:
            # not stored: these arrays are private (copies if they are read-only)
            values = tuple(v if v.flags.writeable else np.array(v) for v in values)
            for v in values:
                v[..., index] = 0.0
            return values
        saved = []
        for v in values:
            saved.append(np.array(v[..., index]))
            v.flags.writeable = True
            v[..., index] = 0.0
            v.flags.writeable = False
        self._zeroed = (index, saved)
        return values

    def get(self, x, same_object_unmodified=False):
        """
        Return (values, changed): the tuple (f(x) for f in funcs), identical to
        evaluating each f on x, and whether x differs from the previous call
        (always True for arrays shorter than min_size, which are not stored).

        same_object_unmodified: the caller guarantees that, if x is the very array
        object passed in the previous call, it has not been modified in place since;
        the comparison of its entries is then skipped.
        """
        self._restore()
        cacheable = (isinstance(x, np.ndarray) and x.ndim == 1 and x.dtype == np.float64
                     and x.flags.c_contiguous)
        if not cacheable or x.shape[0] < self.min_size:
            return tuple(np.asarray(f(x)) for f in self.funcs), True
        if same_object_unmodified and self.x is not None and x is self.x_object:
            return self.values, False
        old = self.x
        if old is not None:
            pieces = _align(x, old, self.max_segments)
            if x.shape == old.shape and pieces == [(True, 0, old.shape[0])]:
                self.x_object = x
                return self.values, False      # unchanged
            values = tuple(self._assemble(pieces, x, ov, f) for ov, f in zip(self.values, self.funcs))
        else:
            values = tuple(np.asarray(f(x)) for f in self.funcs)

        for v in values:
            v.flags.writeable = False
        self.x = x.copy()
        self.x_object = x
        self.values = values
        return values, True

    @staticmethod
    def _assemble(pieces, x, old_values, f):
        parts = [old_values[..., a:b] if reused else f(x[a:b]) for reused, a, b in pieces]
        if len(parts) == 1:
            out = np.array(parts[0], copy=True)
        elif not parts:
            out = f(x)
        else:
            out = np.concatenate(parts, axis=-1)
        return out

class RowLog:
    """
    Append-only table of float rows, stored compactly (8 bytes per value instead of
    a Python list of Python floats, ~40 bytes per value). Used for the hardening,
    merger, TDE and evolution records, which reach millions of rows in compact clusters.
    np.array(log, dtype=float) gives the same array as np.array(list_of_rows, dtype=float).
    """

    def __init__(self, ncols, first_row=None):
        self.ncols = ncols
        self.data = array('d')
        if first_row is not None:
            self.append(first_row)

    def append(self, row):
        if len(row) != self.ncols:
            raise ValueError(f'row of length {len(row)}, expected {self.ncols}')
        self.data.extend(row)

    def set_last(self, col, value):
        """Set column `col` of the last row."""
        self.data[len(self.data) - self.ncols + col] = value

    def __len__(self):
        return len(self.data) // self.ncols

    def __array__(self, dtype=None, copy=None):
        out = np.array(np.frombuffer(self.data, dtype=np.float64)).reshape(-1, self.ncols)
        return out if dtype is None else out.astype(dtype, copy=False)


class StreamedRowLog(RowLog):
    """
    RowLog whose rows are written to their output file during the run, in chunks,
    instead of all at the end, so that memory stays bounded however many rows a run
    produces (millions of hardening and TDE rows in compact clusters).

    The rows go to a temporary file `.<name>.part` next to the output file, which is
    renamed to the output file by finish() when the run ends normally; a run that is
    interrupted leaves only the .part file. The text is the same as writing all rows
    at the end: the same writer (and `prepare` step) is applied chunk by chunk.

    - Row 0 is the placeholder row and is never written.
    - The last `keep` rows always stay in memory: they may still change (set_last)
      or be cut at the end of the run (finish(n_rows) writes only the first n_rows
      data rows).
    - If `enabled` is False the rows are discarded instead of written.
    - len() counts all rows appended so far, including those already written.
    """

    def __init__(self, ncols, first_row, path, header, writer, keep=0, chunk=131072,
                 prepare=None, enabled=True):
        super().__init__(ncols, first_row)
        self.path = path
        self.part_path = os.path.join(os.path.dirname(path), '.' + os.path.basename(path) + '.part')
        self.header = header
        self.writer = writer
        self.keep = keep
        self.chunk = chunk
        self.prepare = prepare
        self.enabled = enabled
        self.removed = 0      # rows (counting the placeholder row) no longer in memory
        self._f = None

    def __len__(self):
        return self.removed + len(self.data) // self.ncols

    def __array__(self, dtype=None, copy=None):
        if self.removed:
            raise RuntimeError('rows of this log were already written to its output file')
        return super().__array__(dtype, copy)

    def _open(self):
        if self._f is None:
            os.makedirs(os.path.dirname(self.path) or '.', exist_ok=True)
            self._f = open(self.part_path, 'w')
            self._f.write(self.header)

    def _emit(self, a, b):
        """Write in-memory rows a..b-1 (the placeholder row excluded) and drop
        in-memory rows 0..b-1."""
        nc = self.ncols
        if self.removed == 0:
            a = max(a, 1)                      # in-memory row 0 is the placeholder row
        if self.enabled and b > a:
            rows = np.frombuffer(self.data, dtype=np.float64, count=b*nc).reshape(b, nc)[a:b].copy()
            if self.prepare is not None:
                rows = self.prepare(rows)
            self._open()
            self.writer(self._f, rows)
        del self.data[:b*nc]
        self.removed += b

    def flush(self):
        """Write the rows that can no longer change, once there are `chunk` of them."""
        n = len(self.data) // self.ncols - self.keep
        if n >= self.chunk:
            self._emit(0, n)

    def finish(self, n_rows=None):
        """Write the data rows not yet written, up to data row n_rows (all rows if
        None; the placeholder row is not counted), close the file and rename it to the
        output file. In-memory row i is row removed + i of the whole log (row 0: the
        placeholder), so data row n_rows is in-memory row n_rows - removed."""
        m = len(self.data) // self.ncols
        if n_rows is None:
            n_rows = self.removed + m - 1
        if n_rows < self.removed - 1:
            raise RuntimeError('rows to be cut at the end were already written')
        end = min(max(n_rows + 1 - self.removed, 0), m)
        self._emit(0, end)
        del self.data[:]
        self.removed += m - end
        if self.enabled:
            self._open()
            self._f.close()
            self._f = None
            os.replace(self.part_path, self.path)

# End of file.
