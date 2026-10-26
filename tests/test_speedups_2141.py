"""Exact (bit-identical) speedups of version 2.14.1."""
import math
import numpy as np

import rapster.compact_accretion as ca
from rapster.bhpool import BHPool, FAM, N_FAMILIES
from rapster.arraycache import RowLog


def test_kept_cumulative_block_sums_are_bit_identical():
    """After every pool change the cumulative block sums kept between draws equal
    np.cumsum of all block sums bit for bit (BHPool.check)."""
    rng = np.random.default_rng(4)
    np.random.seed(4)
    pool = BHPool(rng.uniform(3, 40, 700), np.zeros(700), np.ones(700), np.zeros(700))
    for step in range(6000):
        op = rng.integers(0, 6)
        if op == 0 and pool.size > 3:
            pool.remove(int(rng.integers(0, pool.size)))
        elif op == 1:
            for _ in range(2):
                pool.add(rng.uniform(3, 80), 0.5, 2, 0)
        elif op == 2 and pool.size > 3:
            pool.remove(pool.size - 1)                  # last slot (may empty the last block)
        else:
            pool.set_mass(int(rng.integers(0, pool.size)), rng.uniform(3, 300))
        f = int(rng.integers(N_FAMILIES))
        T = pool.total(f)
        nb = (pool.size + pool.BLOCK - 1) // pool.BLOCK
        assert T == float(np.cumsum(pool._ST[f, :nb])[-1])
        if T > 0:
            pool.draw(f)
        pool.total(FAM['c13']); pool.mean_mass()
        if step % 50 == 0:
            pool.check()
    pool.check()


def _ns_euler_reference(Mi, Mf, chi, dM, eos, prograde):
    """NS branch of evolve_spin_during_accretion as in version 2.14.0 (helpers
    _ns_spin_state and _ns_E_J_isco called at every step); final chi."""
    t = ca.get_eos(eos)
    M, M_TOV = float(Mi), t.M_TOV
    Mg = M * ca.Msun_to_km
    chi = float(min(max(chi, 0.0), ca.THORNE_LIMIT))
    J = chi * Mg**2
    pro = prograde
    last = None
    while M < Mf:
        assert M < t.M_max_rot
        Mg = M * ca.Msun_to_km
        chi = float(min(max(J / Mg**2, 0.0), ca.THORNE_LIMIT))
        R = float(t.radius(min(M, M_TOV)))
        f, fK, chiK = ca._ns_spin_state(M, chi, R)
        last = chi
        if f < fK:
            e, j = ca._ns_E_J_isco(M, R, f, prograde=pro)
            J_new = J + (j / e) * dM * ca.Msun_to_km
            if J_new < 0.0 and not pro:
                J_new, pro = 0.0, True
            J = J_new
            M += dM
        else:
            M += dM
            Mg = M * ca.Msun_to_km
            J = chiK * Mg**2
    return last


def test_ns_spin_evolution_matches_previous_loop():
    rng = np.random.default_rng(7)
    for eos in sorted(ca.EOS_TABLES):
        M_rot = ca.get_eos(eos).M_max_rot
        for _ in range(150):
            dm = 10**rng.uniform(-3, -0.5)
            Mi = rng.uniform(1.0, M_rot - 1.2 * dm)
            chi = float(rng.choice([0.0, rng.uniform(0, 0.7)]))
            pro = bool(rng.integers(2))
            new = ca.evolve_spin_during_accretion(Mi, Mi + dm, True, chi=chi, dM=dm / 100,
                                                  eos=eos, prograde=pro)['chi'][-1]
            assert new == _ns_euler_reference(Mi, Mi + dm, chi, dm / 100, eos, pro)


def test_rowlog_array_equals_list_of_rows():
    rows = [[0.0] * 18] + [list(np.random.default_rng(i).normal(size=18)) for i in range(50)]
    log = RowLog(18, rows[0])
    for r in rows[1:]:
        log.append(r)
    assert np.array_equal(np.array(log, dtype=float), np.array(rows, dtype=float))
