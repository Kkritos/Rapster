"""Central massive BH (-eta, version 2.18.0; massive_bh.py): the promotion criterion, the
loss-cone rates in their limits, the escape velocities, the 3bbs of the other BHs, the
massive BH taken out of the pool, and whole runs (a run without promotion is the same with
-eta 0 and -eta 3; a -MBH seed run; restarts; the stopping conditions). -eta 0 reproduces
version 2.17.0 bit for bit (checked on 32 configurations, see CHANGELOG.md)."""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from rapster import functions as fu
from rapster import massive_bh as mb
from rapster.bhpool import BHPool
from rapster.constants import G_Newton, R_sun
from rapster.three_body_binary import three_body_binary


def _pool(m):
    m = np.asarray(m, dtype=float)
    return BHPool(m, 0.1 * np.ones(m.size), np.ones(m.size), np.zeros(m.size))


def test_promotion_ratio_and_threshold():
    """r_infl / r_wand = (G M / sigma^2) / sqrt(9 m sigma^2 / (4 pi G rho_c M)); the
    critical mass scales as eta^(2/3) (ratio ~ M^(3/2))."""
    m, sigma2, rho = 22.7, 8.4, 1.2e7
    for M in (50.0, 130.0, 1e4):
        r_infl = G_Newton * M / sigma2
        r_wand = np.sqrt(3 * (m / M) * sigma2 / (4 * np.pi * G_Newton * rho / 3))
        assert np.isclose(mb.promotion_ratio(M, m, sigma2, rho), r_infl / r_wand, rtol=1e-12)
    M1 = (9 * m * sigma2**3 / (4 * np.pi * G_Newton**3 * rho))**(1/3)        # ratio 1
    assert np.isclose(mb.promotion_ratio(M1, m, sigma2, rho), 1.0, rtol=1e-12)
    assert np.isclose(mb.promotion_ratio(M1 * 3**(2/3), m, sigma2, rho), 3.0, rtol=1e-12)


def test_loss_cone_rate_limits():
    """1 / (1 / rate_full + t_r L / N): the full loss cone when its rate is the smaller,
    the diffusive rate N / (t_r L) otherwise."""
    assert np.isclose(mb.loss_cone_rate(1e-6, 1e3, 1.0, 10.0), 1e-6, rtol=1e-6)
    assert np.isclose(mb.loss_cone_rate(1e6, 1e3, 1.0, 10.0), 100.0, rtol=1e-3)
    assert np.isclose(mb.loss_cone_rate(100.0, 1e3, 1.0, 10.0), 50.0, rtol=1e-12)
    assert mb.loss_cone_rate(0.0, 1e3, 1.0, 10.0) == 0.0
    assert mb.loss_cone_rate(3.0, 0.5, 1.0, 10.0) == 3.0      # (< 1 bound object: full)


def test_imri_rate_full_loss_cone_is_the_capture_rate():
    """In the full loss cone (diffusion not limiting), the GW capture rate of BHs by the
    massive BH is the two-body capture rate of functions.Rate_cap for the pair: for M = m,
    Rate_cap(m, n, v / sqrt 2) / n (Rate_cap uses the relative velocity sqrt(2) v)."""
    m, n, v = 20.0, 5e5, 6.0
    full = mb.imri_rate(m, m, n, v, 1e30)       # (the diffusive limit is ~N: far larger)
    assert np.isclose(full, fu.Rate_cap(m, n, v / np.sqrt(2)) / n, rtol=1e-6)
    # a heavier massive BH captures faster (cross section ~ M^(12/7) at M >> m):
    assert mb.imri_rate(1e3, m, n, v, 1e30) > 100 * full


def test_tde_rate_full_loss_cone_is_rapsters_micro_tde_rate():
    """In the full loss cone, the TDE rate of the massive BH is Rapster's micro-TDE rate of
    one BH of mass M at rest (cluster_evolution.evolve_tdes)."""
    M, m_avg, n_star, v_star = 300.0, 0.6, 1e5, 30.0
    full = 2 * np.sqrt(2 * (3 * np.pi - 8)) * G_Newton * (M + m_avg) * n_star * R_sun * (M / m_avg)**(1/3) / v_star
    assert np.isclose(mb.tde_rate(M, m_avg, R_sun, n_star, v_star, m_avg, n_star), full, rtol=1e-3)
    # (a heavy BH in a hot, sparse system empties its loss cone: the diffusive rate)
    M, n_star, v_star = 1e6, 1e3, 300.0
    full = 2 * np.sqrt(2 * (3 * np.pi - 8)) * G_Newton * (M + m_avg) * n_star * R_sun * (M / m_avg)**(1/3) / v_star
    assert mb.tde_rate(M, m_avg, R_sun, n_star, v_star, m_avg, n_star) < 0.5 * full


def test_escape_velocities():
    """No massive BH: v_esc_cl = F_ESC sqrt(v_star^2 + vBH^2) (bit for bit); with one:
    sqrt(F^2 (v_star^2 + vBH^2) + 2 G M / max(r_cBH, r_infl)), except for the retention of
    a merger remnant of the massive BH itself."""
    F = fu.F_ESC
    try:
        fu.set_escape_factor(2.0)
        fu.set_massive_bh(0.0, 0.0)
        v0 = fu.v_esc_cl(30.0, 6.0)
        assert v0 == fu.F_ESC * np.sqrt(30.0**2 + 6.0**2)
        assert fu.v_esc_remnant(v0, 30.0, 6.0, 20.0, 30.0) == v0
        fu.set_massive_bh(1e3, 500.0)
        v1 = fu.v_esc_cl(30.0, 6.0)
        assert np.isclose(v1, np.sqrt(v0**2 + 500.0), rtol=1e-12)
        assert fu.v_esc_remnant(v1, 30.0, 6.0, 20.0, 30.0) == v1           # other BHs
        assert fu.v_esc_remnant(v1, 30.0, 6.0, 1e3, 30.0) == v0            # the massive BH
    finally:
        fu.set_massive_bh(0.0, 0.0)
        fu.F_ESC = F


def test_take_out_and_put_back():
    m = np.array([10.0, 2e3, 15.0, 30.0])
    state = {'BH': _pool(m), 'M_MBH': 2e3}
    cfg = {'eta_MBH': 3.0, 'MBH_physics': 1}
    assert mb.slot_in_pool(state) == state['BH'].D
    assert mb.take_out(state, {'eta_MBH': 3.0, 'MBH_physics': 0}) is None       # (Phase 1 only)
    rec = mb.take_out(state, cfg)
    assert rec[0] == 2e3 and state['BH'].size == 3 and state['BH'].max_mass() == 30.0
    mb.put_back(state, rec)
    assert state['BH'].size == 4 and state['BH'].max_mass() == 2e3
    state['M_MBH'] = 0.0
    assert mb.take_out(state, cfg) is None and state['BH'].size == 4


def test_three_body_binaries_of_the_other_BHs_suppressed():
    """f_others = 0: only 3bbs with the massive BH form; f_others = 1: as before (no random
    number used)."""
    m = np.concatenate([[1e3], np.linspace(10, 40, 60)])
    for f, expect_all_with_MBH in ((0.0, True), (1.0, False)):
        BH = _pool(m)
        np.random.seed(3)
        out = three_body_binary(0.0, 0.0, 20, 25.0, np.zeros((1, 15)), BH, 5.0, 0, 0, f_others=f, M_MBH=1e3)
        binaries = out[4][1:]
        assert len(binaries) >= 1
        with_MBH = (binaries[:, 4] == 1e3) | (binaries[:, 5] == 1e3)
        assert with_MBH.sum() <= 1                      # (it is taken by its first binary)
        if expect_all_with_MBH:
            assert with_MBH.all()
    BH1 = _pool(m); BH2 = _pool(m)
    np.random.seed(4); a = three_body_binary(0.0, 0.0, 5, 25.0, np.zeros((1, 15)), BH1, 5.0, 0, 0)[4]
    u1 = np.random.rand()
    np.random.seed(4); b = three_body_binary(0.0, 0.0, 5, 25.0, np.zeros((1, 15)), BH2, 5.0, 0, 0, f_others=1.0, M_MBH=1e3)[4]
    assert np.array_equal(a, b) and np.random.rand() == u1


def test_promotion_bookkeeping():
    """promote: a seed stays promoted before the BH subsystem is active; once active, a
    promoted BH stays promoted while the heaviest BH is at least as heavy; -eta 0: never."""
    from rapster.cluster_evolution import _bh_mean_field
    m = np.concatenate([[5e3], np.full(500, 20.0)])
    state = {'M_MBH': 5e3, 'triples': np.zeros((1, 24)), 'Mcl': 6e5, 'rh': 1.0}
    cfg = {'eta_MBH': 3.0, 'c_xi': 1.0}
    args = (m.size, m.sum(), (m**2.5).sum(), m.max(), 1e5, 0.6, np.log(np.e + 0.02 * 1e6), _bh_mean_field)
    assert mb.promote(state, cfg, 0, *args) == (5e3, None)
    M, mf = mb.promote(state, cfg, 1, *args)
    assert M == 5e3 and np.isclose(mf['mBH_avg'], 20.0) and np.isclose(mf['psi_BH'], 1.0)
    assert mb.promote(state, {'eta_MBH': 0.0, 'c_xi': 1.0}, 1, *args) == (0.0, None)
    # (a BH only a few times heavier than the others is not promoted, even if r_infl >=
    # eta r_wand, MIN_MASS_RATIO)
    m3 = np.full(30, 25.0); m3[0] = 77.0
    st = {'M_MBH': 0.0, 'triples': np.zeros((1, 24)), 'Mcl': 2e4, 'rh': 0.3}
    args3 = (m3.size, m3.sum(), (m3**2.5).sum(), m3.max(), 1e6, 0.6, np.log(np.e + 0.02 * 3e4), _bh_mean_field)
    assert mb.promote(st, {'eta_MBH': 1e-6, 'c_xi': 1.0}, 1, *args3) == (0.0, None)
    # (ejected: the heaviest BH is lighter, and an ordinary BH is not promoted)
    m2 = np.full(500, 20.0); m2[0] = 25.0
    args2 = (m2.size, m2.sum(), (m2**2.5).sum(), m2.max(), 1e5, 0.6, np.log(np.e + 0.02 * 1e6), _bh_mean_field)
    assert mb.promote(state, cfg, 1, *args2) == (0.0, None)


# ------------------------------------------------------------------- whole runs
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
sys.argv = ['rapster', '-P', '0'] + {argv!r}
rc.main()
'''


def _run(d, argv, kill=-1, check=True):
    r = subprocess.run([sys.executable, '-c', HARNESS.format(kill=kill, argv=list(argv))], cwd=d,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if check:
        assert r.returncode == 0, r.stderr
    return r


def _read(d, folder, name):
    return open(os.path.join(d, folder, name)).read()


def test_run_without_promotion_is_unchanged():
    """A run in which no BH is promoted: the same mergers, TDEs and hardening, and the same
    first 71 evolution columns, with -eta 0 (71 columns) and -eta 3 (9 more)."""
    args = ['-N', '3e5', '-r', '1', '-S', '2', '-tM', '2000', '-CKi', '0']
    with tempfile.TemporaryDirectory() as d:
        _run(d, args + ['-eta', '0', '-RF', 'A'])
        _run(d, args + ['-RF', 'B'])
        for name in ('mergers.txt', 'tdes.txt', 'hardening.txt'):
            assert _read(d, 'A', name) == _read(d, 'B', name), name
        a = np.loadtxt(os.path.join(d, 'A', 'evolution.txt')); b = np.loadtxt(os.path.join(d, 'B', 'evolution.txt'))
        assert a.shape[1] == 71 and b.shape[1] == 71 + len(mb.EVOLUTION_KEYS_MBH)
        assert np.array_equal(a, b[:, :71]) and not b[:, 71].any()
        assert 'CENTRAL MASSIVE BH' not in _read(d, 'B', 'log.txt')


def test_seed_run_and_restart():
    """-MBH 1e4 -MBHph 1 in a 6e5 Msun cluster: the seed is promoted from the start (psi_BH
    of the other BHs stays ~1, it captures BHs, merger channel 7, and/or disrupts stars,
    tde_type 31); stopped after a checkpoint, restarted and joined, the run gives the same
    files; a restart cannot change -eta."""
    args = ['-N', '1e6', '-r', '1', '-S', '7', '-tM', '400', '-CKdt', '100', '-MBH', '1e4', '-MBHph', '1']
    files = ('mergers.txt', 'evolution.txt', 'hardening.txt', 'tdes.txt')
    with tempfile.TemporaryDirectory() as d:
        _run(d, args + ['-RF', 'Ref'])
        with open(os.path.join(d, 'Ref', 'evolution.txt')) as f:
            keys = f.readline().lstrip('# ').split()
        E = np.loadtxt(os.path.join(d, 'Ref', 'evolution.txt'))
        K = {k: i for i, k in enumerate(keys)}
        active = E[:, K['N_BH']] > 0
        assert (E[active, K['M_MBH']] >= 1e4).all()
        assert E[active, K['psi_BH']].max() < 3
        M = np.loadtxt(os.path.join(d, 'Ref', 'mergers.txt'), ndmin=2)
        T = np.loadtxt(os.path.join(d, 'Ref', 'tdes.txt'), ndmin=2)
        assert (M.size and (M[:, 2] == mb.CHANNEL_IMRI).any()) or (T.size and (T[:, 3] == mb.TDE_TYPE_STAR).any())
        r = _run(d, args + ['-RF', 'Stopped'], kill=2, check=False)
        assert r.returncode == 3
        _run(d, ['-RS', 'Stopped', '-RF', 'Cont'])
        subprocess.run([sys.executable, '-m', 'rapster.join_restart', '-o', 'Joined', 'Stopped', 'Cont'],
                       cwd=d, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for name in files:
            assert _read(d, 'Ref', name) == _read(d, 'Joined', name), name
        r = _run(d, ['-RS', 'Stopped', '-RF', 'Refused', '-eta', '1'], check=False)
        assert r.returncode != 0 and 'error' in r.stderr


def test_stopping_conditions():
    """A seed heavier than the cluster stops the run (outside the model); a promoted seed
    does not count towards -fCO (1e5 Msun in a 6e5 Msun cluster no longer stops at t = 0);
    -eta < 0 is refused."""
    with tempfile.TemporaryDirectory() as d:
        _run(d, ['-N', '1e5', '-r', '1', '-MBH', '1e5', '-tM', '200', '-CKi', '0', '-RF', 'Heavy'])
        assert 'HEAVIER THAN THE CLUSTER' in _read(d, 'Heavy', 'log.txt')
        _run(d, ['-N', '1e6', '-r', '1', '-MBH', '1e5', '-tM', '30', '-CKi', '0', '-RF', 'Seed'])
        assert 'COMPACT OBJECTS REACHED' not in _read(d, 'Seed', 'log.txt')
        assert len(np.loadtxt(os.path.join(d, 'Seed', 'evolution.txt'), ndmin=2)) > 2
        _run(d, ['-N', '1e6', '-r', '1', '-MBH', '1e5', '-tM', '30', '-CKi', '0', '-eta', '0', '-RF', 'Legacy'])
        assert 'COMPACT OBJECTS REACHED' in _read(d, 'Legacy', 'log.txt')
        r = _run(d, ['-N', '1e4', '-eta', '-1', '-RF', 'Negative'], check=False)
        assert r.returncode != 0 and '-eta' in r.stderr


def test_phase1_without_phase2():
    """-MBHph 0 (the default): the seed is promoted (psi_BH ~1, the column M_MBH), but it
    has no loss-cone captures or TDEs, no escape-velocity term and f_3bb = 1; -MBHph must
    be 0 or 1."""
    with tempfile.TemporaryDirectory() as d:
        _run(d, ['-N', '1e6', '-r', '1', '-S', '7', '-tM', '300', '-CKi', '0', '-MBH', '1e4', '-RF', 'P1'])
        with open(os.path.join(d, 'P1', 'evolution.txt')) as f:
            keys = f.readline().lstrip('# ').split()
        E = np.loadtxt(os.path.join(d, 'P1', 'evolution.txt'))
        K = {k: i for i, k in enumerate(keys)}
        active = E[:, K['N_BH']] > 0
        assert (E[active, K['M_MBH']] >= 1e4).all() and E[active, K['psi_BH']].max() < 3
        assert not E[:, K['N_IMRI']].any() and not E[:, K['N_tdeMBH']].any() and (E[:, K['f_3bb']] == 1).all()
        M = np.loadtxt(os.path.join(d, 'P1', 'mergers.txt'), ndmin=2)
        T = np.loadtxt(os.path.join(d, 'P1', 'tdes.txt'), ndmin=2)
        assert not (M.size and (M[:, 2] == mb.CHANNEL_IMRI).any()) and not (T.size and (T[:, 3] >= 31).any())
        r = _run(d, ['-N', '1e4', '-MBHph', '2', '-RF', 'Bad'], check=False)
        assert r.returncode != 0 and '-MBHph' in r.stderr
