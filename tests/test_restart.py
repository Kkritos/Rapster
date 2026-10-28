"""Checkpoints and restarts (restart.py, -RS, join_restart.py): a run stopped after a
checkpoint and restarted, or continued from its final checkpoint to a larger t_max, gives
after joining the same files, byte for byte, as a run without interruption; changes at a
restart; options a restart refuses."""
import json
import os
import pickle
import subprocess
import sys
import tempfile

import numpy as np
import pytest

FILES = ('mergers.txt', 'evolution.txt', 'hardening.txt', 'tdes.txt', 'populations.txt')

# runs rapster in-process with a small chunk of streamed rows (so that rows are written
# during the run), stopping the process abruptly (os._exit) after checkpoint KILL:
HARNESS = '''
import os, sys
import rapster.cluster_evolution as ce
ce.STREAM_CHUNK_ROWS = 7
import rapster.run_cluster as rc
save = rc.save_checkpoint
count = [0]
def save_and_maybe_stop(state, config, version):
    save(state, config, version)
    count[0] += 1
    if count[0] == {kill}:
        os._exit(3)
rc.save_checkpoint = save_and_maybe_stop
sys.argv = ['rapster'] + {argv!r}
rc.main()
'''


def _run(d, argv, kill=-1, check=True):
    r = subprocess.run([sys.executable, '-c', HARNESS.format(kill=kill, argv=['-P', '0'] + list(argv))],
                       cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if check:
        assert r.returncode == 0, r.stderr
    return r


def _join(d, out, *folders):
    subprocess.run([sys.executable, '-m', 'rapster.join_restart', '-o', out] + list(folders),
                   cwd=d, check=True, stdout=subprocess.DEVNULL)


def _same(d, a, b, bhs=True):
    for name in FILES:
        pa, pb = os.path.join(d, a, name), os.path.join(d, b, name)
        assert os.path.isfile(pa) == os.path.isfile(pb), name
        if os.path.isfile(pa):
            assert open(pa).read() == open(pb).read(), name
    if bhs:
        assert open(os.path.join(d, a, 'outputBHs.pkl'), 'rb').read() == open(os.path.join(d, b, 'outputBHs.pkl'), 'rb').read()
        with open(os.path.join(d, a, 'outputBHs.pkl'), 'rb') as f:
            x = pickle.load(f)
        with open(os.path.join(d, b, 'outputBHs.pkl'), 'rb') as f:
            y = pickle.load(f)
        assert x.keys() == y.keys()
        for k in x:
            assert len(x[k]) == len(y[k]) and all(np.array_equal(u, v) for u, v in zip(x[k], y[k])), k


ARGS = ['-N', '1e5', '-r', '0.6', '-S', '21', '-tM', '2500', '-CKdt', '300', '-POP', 'pops.txt']


@pytest.fixture(scope='module')
def reference():
    d = tempfile.mkdtemp()
    with open(os.path.join(d, 'pops.txt'), 'w') as f:
        f.write('400 4e4 0.002 0.7 0.5 3\n1300 2e4\n')
    _run(d, ARGS + ['-RF', 'Ref'])
    return d


def test_stopped_run_restarted_and_joined_is_identical(reference):
    d = reference
    r = _run(d, ARGS + ['-RF', 'Stopped'], kill=4, check=False)
    assert r.returncode == 3
    assert not os.path.isfile(os.path.join(d, 'Stopped', 'evolution.txt'))      # only .part files
    _run(d, ['-RS', 'Stopped', '-RF', 'Cont'])
    info = json.load(open(os.path.join(d, 'Cont', 'restart.json')))
    assert info['t_restart'] >= 1200 and info['changes'] == {} and info['rows_before']['evolution'] > 0
    _join(d, 'Joined', 'Stopped', 'Cont')
    _same(d, 'Ref', 'Joined', bhs=False)                   # (the stopped run has no outputBHs.pkl)


def test_chain_of_restarts_and_continuation_after_the_end(reference):
    d = reference
    # a run to 900 Myr, continued to 1700 Myr, then to 2500 Myr (each from its final checkpoint):
    _run(d, ARGS[:-4] + ['-tM', '900', '-CKdt', '300', '-POP', 'pops.txt', '-RF', 'P1'])
    _run(d, ['-RS', 'P1', '-RF', 'P2', '-tM', '1700'])
    _run(d, ['-RS', os.path.join('P2', 'checkpoint.pkl'), '-RF', 'P3', '-tM', '2500'])
    _join(d, 'Chain', 'P1', 'P2', 'P3')
    _same(d, 'Ref', 'Chain')


def test_restart_with_changes(reference):
    d = reference
    np.savez(os.path.join(d, 'bhs.npz'), mBH=np.array([25.0, 60.0]), sBH=np.array([0.2, 0.7]))
    with open(os.path.join(d, 'late.txt'), 'w') as f:
        f.write('2000 1e4 0.0003\n')
    _run(d, ARGS + ['-tM', '1000', '-RF', 'Base'])
    _run(d, ['-RS', 'Base', '-RF', 'Changed', '-tM', '2200', '-S', '5', '-r', '2.5', '-addBH', 'bhs.npz',
             '-addN', '3e4', '-addZ', '0.0001', '-addSFE', '0.6', '-addx', '0.4', '-POP', 'late.txt'])
    info = json.load(open(os.path.join(d, 'Changed', 'restart.json')))
    t0 = info['t_restart']
    assert info['changes']['seed'] == 5 and info['changes']['r_h'][1] == 2.5 and info['changes']['BHs_added'] == 2
    births = [p['t_birth'] for p in info['changes']['populations_added']]
    assert births == [t0, 2000.0]
    log = open(os.path.join(d, 'Changed', 'log.txt')).read()
    assert 'BORN AT t = {:.6g} Myr: N = 30000'.format(t0) in log and 'BORN AT t = 2000 Myr' in log
    ev = np.loadtxt(os.path.join(d, 'Changed', 'evolution.txt'), ndmin=2)
    from rapster.constants import evolution_keys
    K = {n: i for i, n in enumerate(evolution_keys)}
    # (the file starts with the rows of the last steps before the checkpoint that the run
    # restarted from had not written yet)
    assert (ev[:, K['t']] < t0).sum() <= 2
    first = ev[ev[:, K['t']] == t0][0]
    assert first[K['r_h']] == 2.5 and first[K['seed']] == 5
    assert first[K['M_gas']] == pytest.approx(0.4 / 0.6 * 3e4 * 0.5859340259680028, rel=1e-9)
    _join(d, 'JoinedChanged', 'Base', 'Changed')
    with open(os.path.join(d, 'JoinedChanged', 'evolution.txt')) as f:
        assert sum(1 for line in f if not line.startswith('#')) == info['rows_before']['evolution'] + len(ev)


@pytest.mark.parametrize('argv', [
    ['-N', '1e4'],                       # a cluster parameter
    ['-RF', 'Base'],                     # the folder restarted from
    ['-addZ', '0.001'],                  # a population without -addN
    ['-POP', 'early.txt'],               # a population born before the restart
])
def test_restart_refuses(reference, argv):
    d = reference
    if not os.path.isdir(os.path.join(d, 'Base')):
        _run(d, ARGS + ['-tM', '1000', '-RF', 'Base'])
    with open(os.path.join(d, 'early.txt'), 'w') as f:
        f.write('100 1e4\n')
    r = _run(d, ['-RS', 'Base', '-RF', 'Refused'] + argv if '-RF' not in argv else ['-RS', 'Base'] + argv, check=False)
    assert r.returncode != 0 and 'error' in r.stderr


def test_restart_options_need_a_restart():
    with tempfile.TemporaryDirectory() as d:
        r = _run(d, ['-addN', '1e4'], check=False)
        assert r.returncode != 0 and 'restart' in r.stderr


def test_no_checkpoint_with_CKi_0():
    with tempfile.TemporaryDirectory() as d:
        _run(d, ['-N', '2e4', '-tM', '1200', '-CKdt', '200', '-CKi', '0', '-RF', 'R'])
        assert not os.path.exists(os.path.join(d, 'R', 'checkpoint.pkl'))


def test_bh_pool_pickle_round_trip_continues_identically():
    from rapster.bhpool import BHPool, FAM
    rng = np.random.RandomState(4)
    pool = BHPool(rng.uniform(5, 50, 3000), np.zeros(3000), np.ones(3000), np.zeros(3000))
    for _ in range(500):                                     # changes, draws (fill the caches)
        pool.remove(rng.randint(pool.n))
        pool.set_mass(rng.randint(pool.n), rng.uniform(5, 80))
        np.random.seed(rng.randint(1000)); pool.draw(FAM['mass'])
    copy = pickle.loads(pickle.dumps(pool))
    assert len(pickle.dumps(pool)) < 0.7 * pool._W.nbytes        # trailing empty slots not stored
    for _ in range(800):
        k, m = rng.randint(pool.n), rng.uniform(5, 300)
        seed = rng.randint(10**6)
        out = []
        for p in (pool, copy):
            p.set_mass(k, m)
            if k % 3 == 0:
                p.add(m / 2, 0.1, 2, 0)
            np.random.seed(seed)
            out.append((p.draw(FAM['mass']), p.draw(FAM['p25']), p.total(FAM['mass']), p.max_mass(), p.n))
        assert out[0] == out[1]
    assert np.array_equal(pool._W, copy._W) and np.array_equal(pool._S, copy._S)


def test_restart_adding_the_first_population():
    with tempfile.TemporaryDirectory() as d:
        _run(d, ['-N', '5e4', '-S', '2', '-tM', '1200', '-RF', 'A'])
        with open(os.path.join(d, 'A', 'checkpoint.pkl'), 'rb') as f:
            t0 = float(pickle.load(f)['state']['t'])
        r = _run(d, ['-RS', 'A', '-RF', 'B', '-addN', '1e4'], check=False)      # t_max not later
        assert r.returncode != 0 and 'later -tM' in r.stderr
        for name, t_max in (('C', t0 + 0.05), ('D', t0 + 400)):                # (ends within a step, or later)
            _run(d, ['-RS', 'A', '-RF', name, '-addN', '1e4', '-tM', repr(t_max)])
            ev = [l.split() for l in open(os.path.join(d, name, 'evolution.txt')) if not l.startswith('#')]
            po = [l.split() for l in open(os.path.join(d, name, 'populations.txt')) if not l.startswith('#')]
            times = sorted({r[0] for r in po}, key=float)
            assert times == [r[1] for r in ev]
            M = {}
            for row in po:
                M[row[0]] = M.get(row[0], 0.0) + float(row[4])
            assert all(M[row[1]] == pytest.approx(float(row[5]), rel=1e-12) for row in ev)
