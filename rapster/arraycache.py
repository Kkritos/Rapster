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

import atexit
import os
import pickle
import queue
import struct
import subprocess
import sys
import threading
from itertools import islice
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


def write_rows(f, rows, suffix=''):
    """
    Write a 2-D float array, one row per line, values separated by spaces and
    `suffix` appended to every line. repr of a Python float is the same text as
    str of a NumPy float64 (shortest round-trip representation), so the file is
    identical to writing str(rows[i][j]) value by value, which it replaces
    (~10x faster for the large hardening and evolution files).
    """
    chunk = 10000   # rows converted to Python floats at a time (bounded memory)
    for start in range(0, len(rows), chunk):
        f.writelines(' '.join(map(repr, row)) + suffix + '\n' for row in rows[start:start+chunk].tolist())


# ---------------------------------------------------------------- background writer
#
# Writing the rows of a StreamedRowLog to its file (write_rows: the text of each value,
# ~0.5 us per value) is done by a separate process, while the simulation continues:
# the rows are sent to it through a pipe (by a thread of this process, so that a send
# never waits for the writer), and it appends them to the file, in order. The text is
# the same. StreamedRowLog.sync (checkpoints) and finish wait until the writer has
# written (and flushed) every row sent before. Used for chunks of at least
# BACKGROUND_MIN_ROWS rows (the large hardening, TDE and evolution files of long runs),
# unless the environment variable RAPSTER_BACKGROUND_WRITER is 0. If the writer cannot
# be started, the rows are written by this process as before.

BACKGROUND_MIN_ROWS = 20000


def _read_exact(f, n):
    data = f.read(n)
    if len(data) != n:
        raise EOFError
    return data


def _send_message(f, header, payload=b''):
    h = pickle.dumps(header, protocol=pickle.HIGHEST_PROTOCOL)
    f.write(struct.pack('<QQ', len(h), len(payload)))
    f.write(h)
    if payload:
        f.write(payload)


def _read_message(f):
    nh, npl = struct.unpack('<QQ', _read_exact(f, 16))
    header = pickle.loads(_read_exact(f, nh))
    return header, (_read_exact(f, npl) if npl else b'')


def _writer_main():
    """The writer process: append the rows received on stdin to their files, until stdin
    is closed; answer each sync request on stdout (None, or the error met so far)."""
    inp = sys.stdin.buffer
    out = sys.stdout.buffer
    files = {}
    error = None
    while True:
        try:
            header, payload = _read_message(inp)
        except EOFError:
            break
        op = header[0]
        if op == 'w':
            _, path, ncols = header
            if error is None:
                try:
                    f = files.get(path)
                    if f is None:
                        f = files[path] = open(path, 'a')
                    write_rows(f, np.frombuffer(payload, dtype=np.float64).reshape(-1, ncols))
                except Exception as err:          # (reported at the next sync)
                    error = repr(err)
        elif op == 's':
            _, path, close = header
            if error is None:
                try:
                    f = files.get(path)
                    if f is not None:
                        f.flush()
                        if close:
                            f.close()
                            del files[path]
                except Exception as err:
                    error = repr(err)
            _send_message(out, (error,))
            out.flush()
    for f in files.values():
        try:
            f.close()
        except Exception:
            pass


# (the writer runs this module by its file, without importing the rapster package)
_WRITER_BOOT = ("import sys, importlib.util as u; s = u.spec_from_file_location('_rapster_writer', sys.argv[1]); "
                "m = u.module_from_spec(s); s.loader.exec_module(m); m._writer_main()")


class _BackgroundWriter:
    PIECE_ROWS = 10000

    def __init__(self):
        self.proc = subprocess.Popen([sys.executable, '-c', _WRITER_BOOT, os.path.abspath(__file__)],
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        # (at most two chunks wait to be sent: if the writer falls behind, the simulation
        # waits for it rather than keeping more rows in memory)
        self.queue = queue.Queue(maxsize=2)
        self.error = None
        self.thread = threading.Thread(target=self._send_loop, daemon=True, name='rapster-writer-send')
        self.thread.start()
        atexit.register(self.close)

    def alive(self):
        return self.proc.poll() is None

    def _send_loop(self):
        while True:
            item = self.queue.get()
            try:
                if item is None:
                    return
                if self.error is None:
                    _send_message(self.proc.stdin, *item)
                    self.proc.stdin.flush()
            except Exception as err:
                self.error = err
            finally:
                self.queue.task_done()

    def write(self, path, rows):
        # (in pieces of PIECE_ROWS rows: at most a few pieces are in memory on their way
        # to the writer, instead of whole chunks)
        rows = np.ascontiguousarray(rows, dtype=np.float64)
        for start in range(0, len(rows), self.PIECE_ROWS):
            self.queue.put((('w', path, rows.shape[1]), rows[start:start + self.PIECE_ROWS].tobytes()))

    def wait(self, path, close=False):
        """Wait until every row sent for path is written and flushed (and the file closed)."""
        self.queue.join()                     # (the send thread is idle from here)
        if self.error is None:
            try:
                _send_message(self.proc.stdin, ('s', path, close))
                self.proc.stdin.flush()
                (err,), _ = _read_message(self.proc.stdout)
                if err is not None:
                    self.error = RuntimeError(err)
            except Exception as err:
                self.error = err
        if self.error is not None:
            raise RuntimeError(f'rows of {path} could not be written by the background writer: '
                               f'{self.error!r}')

    def close(self):
        """Stop the writer once it has written everything (at exit)."""
        try:
            self.queue.put(None, timeout=60)
            self.thread.join(timeout=60)
            self.proc.stdin.close()
            self.proc.wait(timeout=60)
        except Exception:
            pass


_background = None
_background_failed = False


def _background_writer():
    """The background writer of this process (started at the first use), or None."""
    global _background, _background_failed
    if _background is not None and _background.alive():
        return _background
    if _background_failed or os.environ.get('RAPSTER_BACKGROUND_WRITER', '1') == '0':
        return None
    if _background is not None:            # (the writer stopped: rows may be lost)
        raise RuntimeError('the background writer of the output files stopped')
    try:
        _background = _BackgroundWriter()
    except Exception as err:
        _background_failed = True
        print(f'(background writer not available, rows written directly: {err!r})', file=sys.stderr)
        return None
    return _background


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
        self.offset = 0       # data rows of the log written by an earlier run (restarts)
        self._f = None
        self._bg = None       # background writer of the file, once a chunk was sent to it

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
            bg = self._bg
            if bg is None and self.writer is write_rows and len(rows) >= BACKGROUND_MIN_ROWS:
                bg = _background_writer()
                if bg is not None:
                    self._f.flush()               # (the header and rows written so far)
                    self._bg = bg
            if bg is not None:
                bg.write(self.part_path, rows)
            else:
                self.writer(self._f, rows)
        del self.data[:b*nc]
        self.removed += b

    def _wait_background(self):
        """Wait until the background writer has written this file's rows, and closed it."""
        if self._bg is not None:
            self._f.flush()
            self._bg.wait(self.part_path, close=True)
            self._bg = None

    def flush(self):
        """Write the rows that can no longer change, once there are `chunk` of them."""
        n = len(self.data) // self.ncols - self.keep
        if n >= self.chunk:
            self._emit(0, n)

    def sync(self):
        """Write every row but the last `keep` to the file now (whatever the chunk size),
        and flush the file; used at checkpoints. Returns (removed, rows in
        memory): the state from which restore() continues the log in another run.
        (No fsync: the rows are in the operating system's buffers, which a stopped process
        does not lose; an fsync per checkpoint would slow down short runs noticeably.)"""
        n = len(self.data) // self.ncols - self.keep
        if n > 0:
            self._emit(0, n)
        if self.enabled:
            self._open()
            self._f.flush()
            if self._bg is not None:
                self._bg.wait(self.part_path)
        rows = np.frombuffer(self.data, dtype=np.float64).reshape(-1, self.ncols).copy()
        return self.removed, rows

    def restore(self, removed, rows):
        """Continue a log saved by sync() in another run (a restart): its rows before
        the saved in-memory rows are in the earlier run's file, not in this one. Row
        numbers (len, finish) stay those of the whole log."""
        if self.removed or self._f is not None:
            raise RuntimeError('restore() must be called on a new log')
        rows = np.asarray(rows, dtype=np.float64).reshape(-1, self.ncols)
        del self.data[:]
        self.data.frombytes(rows.tobytes())
        self.removed = int(removed)
        self.offset = max(self.removed - 1, 0)

    def finish(self, n_rows=None, suffix=''):
        """Write the data rows not yet written, up to data row n_rows (all rows if
        None; the placeholder row is not counted), close the file and rename it to the
        output file. In-memory row i is row removed + i of the whole log (row 0: the
        placeholder), so data row n_rows is in-memory row n_rows - removed.

        With a suffix (text known only at the end, e.g. the final cluster properties
        on every merger line), all rows are written first and the file is then copied
        line by line, keeping the first n_rows data lines and appending the suffix to
        each: the same text as writing rows + suffix at the end."""
        if suffix:
            self._finish_with_suffix(n_rows, suffix)
            return
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
            self._wait_background()
            self._f.close()
            self._f = None
            os.replace(self.part_path, self.path)

    def _finish_with_suffix(self, n_rows, suffix):
        m = len(self.data) // self.ncols
        self._emit(0, m)
        if not self.enabled:
            return
        self._open()
        self._wait_background()
        self._f.close()
        self._f = None
        tmp_path = self.part_path + '2'
        with open(self.part_path) as src, open(tmp_path, 'w') as dst:
            dst.write(src.readline())                       # header line
            lines = src if n_rows is None else islice(src, max(n_rows - self.offset, 0))
            dst.writelines(line[:-1] + suffix + '\n' for line in lines)
        os.remove(self.part_path)
        os.replace(tmp_path, self.path)

# End of file.
