"""Output records written during the run (StreamedRowLog) give the same files as
writing all rows at the end."""
import filecmp
import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from rapster.arraycache import StreamedRowLog


def _writer(f, rows):
    f.writelines(' '.join(map(repr, r)) + '\n' for r in rows.tolist())


def _expected(rows, n_rows):
    data = rows[1:] if n_rows is None else rows[1:1 + n_rows]
    return '# h\n' + ''.join(' '.join(map(repr, r)) + '\n' for r in data)


@pytest.mark.parametrize('chunk', [1, 2, 3, 5, 100])
@pytest.mark.parametrize('keep, cut', [(0, 0), (1, 0), (2, 1), (2, 2)])
def test_streamed_rows_match_rows_written_at_the_end(chunk, keep, cut):
    rng = np.random.default_rng(chunk + 10 * keep + cut)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'out', 'x.txt')
        log = StreamedRowLog(3, [0.0] * 3, path, '# h\n', _writer, keep=keep, chunk=chunk)
        rows = [[0.0] * 3]
        for i in range(37):
            row = rng.normal(size=3).tolist()
            log.append(row); rows.append(row)
            if keep and rng.random() < 0.5:            # the last row may still change
                log.set_last(2, float(i)); rows[-1][2] = float(i)
            if rng.random() < 0.6:
                log.flush()
        assert len(log) == len(rows)
        n_rows = None if cut == 0 else len(rows) - 1 - cut
        log.finish(n_rows)
        with open(path) as f:
            assert f.read() == _expected(rows, n_rows)
        assert os.listdir(os.path.dirname(path)) == ['x.txt']   # no .part file left


@pytest.mark.parametrize('chunk', [1, 4, 100])
@pytest.mark.parametrize('n_rows', [None, 0, 5, 37, 50])
def test_suffix_added_at_the_end(chunk, n_rows):
    """Rows streamed during the run get a suffix known only at the end (merger rows:
    the final cluster properties), and only the first n_rows are kept."""
    rng = np.random.default_rng(chunk)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'm.txt')
        log = StreamedRowLog(3, [0.0] * 3, path, '# h\n', _writer, chunk=chunk)
        rows = [[0.0] * 3]
        for i in range(37):
            row = rng.normal(size=3).tolist()
            log.append(row); rows.append(row)
            log.flush()
        log.finish(n_rows, suffix=' 1.5 2.5')
        data = rows[1:] if n_rows is None else rows[1:1 + n_rows]
        with open(path) as f:
            assert f.read() == '# h\n' + ''.join(' '.join(map(repr, r)) + ' 1.5 2.5\n' for r in data)
        assert os.listdir(d) == ['m.txt']


def test_disabled_and_empty_records():
    with tempfile.TemporaryDirectory() as d:
        off = StreamedRowLog(2, [0.0, 0.0], os.path.join(d, 'off.txt'), '# h\n', _writer,
                             chunk=1, enabled=False)
        for i in range(10):
            off.append([i, i]); off.flush()
        off.finish()
        empty = StreamedRowLog(2, [0.0, 0.0], os.path.join(d, 'empty.txt'), '# h\n', _writer)
        empty.finish()
        assert sorted(os.listdir(d)) == ['empty.txt']
        with open(os.path.join(d, 'empty.txt')) as f:
            assert f.read() == '# h\n'


def test_cutting_rows_already_written_raises():
    with tempfile.TemporaryDirectory() as d:
        log = StreamedRowLog(1, [0.0], os.path.join(d, 'x.txt'), '', _writer, keep=0, chunk=1)
        for i in range(5):
            log.append([i]); log.flush()
        with pytest.raises(RuntimeError):
            log.finish(2)


def _run(directory, chunk, args):
    code = ('import sys, rapster.cluster_evolution as ce; ce.STREAM_CHUNK_ROWS = %d; '
            'sys.argv = ["rapster"] + %r; from rapster.run_cluster import main; main()' % (chunk, args))
    subprocess.run([sys.executable, '-c', code], cwd=directory, check=True,
                   stdout=subprocess.DEVNULL, env=dict(os.environ))


@pytest.mark.parametrize('args', [
    ['-N', '2e6', '-r', '0.4', '-S', '5', '-P', '0', '-BOi', '0'],
    ['-N', '1e5', '-r', '0.3', '-S', '6', '-P', '0', '-BOi', '0'],           # ends early
    ['-N', '1e6', '-S', '7', '-P', '0', '-BOi', '0', '-Hi', '0', '-tM', '3000'],
])
def test_run_files_independent_of_chunk_size(args):
    with tempfile.TemporaryDirectory() as d:
        os.mkdir(os.path.join(d, 'a')); os.mkdir(os.path.join(d, 'b'))
        _run(os.path.join(d, 'a'), 3, args)              # flushes at almost every step
        _run(os.path.join(d, 'b'), 10**9, args)          # everything written at the end
        ra, rb = os.path.join(d, 'a', 'Results'), os.path.join(d, 'b', 'Results')
        names = sorted(os.listdir(ra))
        assert names == sorted(os.listdir(rb))
        assert not any(n.endswith('.part') for n in names)
        for n in names:
            if n.endswith('.txt') and n != 'log.txt':
                assert filecmp.cmp(os.path.join(ra, n), os.path.join(rb, n), shallow=False), n
