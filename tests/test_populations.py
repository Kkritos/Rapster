"""Multiple stellar populations (populations.py, -POP): input file, bookkeeping of the
populations, segregation limits, field-star mixture, and whole runs (population rows in
step with the evolution rows, masses adding up to the cluster mass, compact objects of
later populations, unchanged results while only one population is born)."""
import math
import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest

from rapster import populations as P
from rapster.constants import t_sev, nu_sev, evolution_keys
from rapster.stellar_evolution import init_stellar_mass_sampler, get_star, set_field_star_mixture

init_stellar_mass_sampler(0.0)
K = {n: i for i, n in enumerate(evolution_keys)}
M_AVG0 = 0.5859340259680028


# ----------------------------------------------------------------------- input file
def _write(d, text):
    path = os.path.join(d, 'pops.txt')
    with open(path, 'w') as f:
        f.write(text)
    return path


def test_read_populations_file_defaults_comments_and_order():
    with tempfile.TemporaryDirectory() as d:
        pops = P.read_populations_file(_write(d, '# t N Z SFE x n_mix\n\n900 1e4 # late\n300, 2e5, 0.002, 0.5, 0.4, 3\n'), 0.001)
    assert [p['t_birth'] for p in pops] == [300.0, 900.0]
    assert pops[0] == dict(t_birth=300.0, N=2e5, Z=0.002, SFE=0.5, x=0.4, n_mix=3.0)
    assert pops[1] == dict(t_birth=900.0, N=1e4, Z=0.001, SFE=1.0, x=1.0, n_mix=P.DEFAULT_N_MIX)


@pytest.mark.parametrize('line', ['0 1e4', '-5 1e4', '10', '10 0.5', '10 1e4 -0.1', '10 1e4 0.001 0',
                                  '10 1e4 0.001 1.5', '10 1e4 0.001 1 0', '10 1e4 0.001 1 1 0',
                                  '10 1e4 0.001 1 1 5 7', '10 abc', '10 nan', '10 inf'])
def test_read_populations_file_rejects_invalid_lines(line):
    with tempfile.TemporaryDirectory() as d:
        with pytest.raises(ValueError):
            P.read_populations_file(_write(d, line + '\n'), 0.001)


# ----------------------------------------------------------------------- bookkeeping
def test_plummer_mass_fraction_outside():
    # x = 1: Plummer sphere with the cluster's r_h; mass within r of a Plummer sphere of scale
    # radius a: r^3 / (r^2 + a^2)^(3/2), and r_h = 1.305 a
    a = 1 / math.sqrt(2**(2/3) - 1)
    r = 2 * a
    assert P.plummer_mass_fraction_outside(1.0) == pytest.approx(1 - r**3 / (r**2 + 1)**1.5)
    assert P.plummer_mass_fraction_outside(0.5) == pytest.approx(1 - (2 * r)**3 / ((2 * r)**2 + 1)**1.5)
    assert P.plummer_mass_fraction_outside(1e-3) < 1e-5            # very concentrated: no escapers
    assert P.plummer_mass_fraction_outside(1e3) > 0.999           # very extended: all outside
    # half of the mass lies outside the half-mass radius:
    assert (r / 2)**3 / ((r / 2)**2 + 1)**1.5 == pytest.approx(0.5)


def _two_populations(x=1.0, n_mix=5.0):
    pops = P.StellarPopulations(4e5, 0.001, 1.0, [dict(t_birth=500.0, N=1e5, Z=0.002, SFE=1.0, x=x, n_mix=n_mix)])
    assert not pops.multi and pops.next_birth_time() == 500.0
    assert pops.due(499.0) == []
    (p,) = pops.due(500.0 - 1e-12)
    pops.add(500.0, p['N'], p['N'] * M_AVG0, p['Z'], p['SFE'], p['x'], p['n_mix'])
    assert pops.multi and pops.scheduled == [] and pops.next_birth_time() == math.inf
    return pops


def test_numbers_and_mean_masses():
    pops = _two_populations()
    t = 700.0
    m = pops.mean_masses(t, M_AVG0)
    assert m[0] == pytest.approx(M_AVG0 * (t / t_sev)**(-nu_sev))
    assert m[1] == pytest.approx(M_AVG0 * (200.0 / t_sev)**(-nu_sev))
    N = pops.numbers(t, M_AVG0)
    assert N[1] == pytest.approx(1e5 * M_AVG0 / m[1])
    assert P.mean_stellar_mass(M_AVG0, 1.0) == M_AVG0          # before t_sev


def test_well_mixed_populations_have_no_segregation_effects():
    pops = _two_populations(x=1.0)
    assert not pops.segregated
    assert pops.density_factor(800.0, M_AVG0) == 1.0
    w = pops.escape_weights()
    assert w == pytest.approx([pops.M[0] / sum(pops.M), pops.M[1] / sum(pops.M)])
    assert pops.central_weights(800.0, M_AVG0) == pytest.approx(pops.numbers(800.0, M_AVG0))


def test_segregated_population_weights_and_density_factor():
    x = 0.5
    pops = _two_populations(x=x)
    assert pops.segregated
    t = 800.0
    N = pops.numbers(t, M_AVG0)
    assert pops.central_weights(t, M_AVG0) == pytest.approx([N[0], N[1] / x**3])
    assert pops.density_factor(t, M_AVG0) == pytest.approx((N[0] + N[1] / x**3) / sum(N))
    w = pops.escape_weights()
    f0, f1 = P.plummer_mass_fraction_outside(1.0), P.plummer_mass_fraction_outside(x)
    expected = np.array([pops.M[0] * f0, pops.M[1] * f1])
    assert w == pytest.approx(list(expected / expected.sum()))
    assert w[1] < pops.M[1] / sum(pops.M)          # concentrated: loses less than its share


def test_evolve_conserves_mass_and_mixes():
    pops = _two_populations(x=0.4, n_mix=3.0)
    t, dt, t_rlx = 900.0, 10.0, 300.0
    M_before = list(pops.M)
    shares = pops.escape_weights()
    dM_rlx = -0.01 * sum(M_before)
    dM_sev = pops.evolve(t, dt, M_AVG0, dM_rlx, t_rlx)
    sev = [-nu_sev * M_before[0] * dt / t, -nu_sev * M_before[1] * dt / (t - 500.0)]
    assert dM_sev == pytest.approx(sum(sev))
    for k in range(2):
        assert pops.M[k] == pytest.approx(M_before[k] + sev[k] + shares[k] * dM_rlx)
    assert sum(pops.M) == pytest.approx(sum(M_before) + dM_sev + dM_rlx)
    assert pops.x[0] == 1.0
    assert pops.x[1] == pytest.approx(1 - 0.6 * math.exp(-dt / (3.0 * t_rlx)))


def test_evolve_young_population_has_no_stellar_evolution_mass_loss_before_t_sev():
    pops = _two_populations()
    dM = pops.evolve(501.0, 0.5, M_AVG0, 0.0, 100.0)
    assert dM == pytest.approx(-nu_sev * 4e5 * 0.5 / 501.0)
    assert pops.M[1] == 1e5 * M_AVG0


def test_white_dwarf_formation_rate_is_the_sum_over_populations():
    from rapster.stellar_evolution import IMF_kroupa
    from scipy import integrate
    Kn = 1 / integrate.quad(IMF_kroupa, 0.08, 150)[0]
    def single(N, age):
        m = (1e4 / age)**(1 / 2.5)
        return N / age / 2.5 * Kn * IMF_kroupa(m) * m if age > 1e4 / 8**2.5 else 0.0
    pops = _two_populations()
    t = 2000.0
    N = pops.numbers(t, M_AVG0)
    assert pops.white_dwarf_formation_rate(t, M_AVG0, Kn, IMF_kroupa) == pytest.approx(single(N[0], t) + single(N[1], t - 500.0))
    # a newborn population adds no WDs:
    assert pops.white_dwarf_formation_rate(510.0, M_AVG0, Kn, IMF_kroupa) == pytest.approx(single(pops.numbers(510.0, M_AVG0)[0], 510.0))


def test_rows_before_and_after_the_second_birth():
    pops = P.StellarPopulations(4e5, 0.001, 1.0, [dict(t_birth=500.0, N=1e5, Z=0.002, SFE=1.0, x=0.5, n_mix=5.0)])
    (row,) = pops.rows(100.0, M_AVG0, 3.9e5)            # before: the cluster's mass
    assert row[0] == 100.0 and row[1] == 0 and row[4] == 3.9e5
    pops = _two_populations(x=0.5)
    rows = pops.rows(600.0, M_AVG0, -1.0)
    assert [r[1] for r in rows] == [0, 1] and rows[1][6] == 0.5
    assert len(P.POPULATION_KEYS) == len(rows[0])


# ----------------------------------------------------------------------- field stars
def test_field_star_mixture():
    tBH = 10.0
    try:
        # all weight on a 1000 Myr old population: no star above its turn-off mass
        set_field_star_mixture([1000.0, 5.0], [1.0, 0.0])
        np.random.seed(1)
        m = np.array([get_star(5.0, tBH, 0.08, 150)[0] for _ in range(3000)])
        assert m.max() <= (1e4 / 1000.0)**(1 / 2.5)
        # all weight on a young population: massive stars appear
        set_field_star_mixture([1000.0, 5.0], [0.0, 1.0])
        np.random.seed(1)
        m = np.array([get_star(1000.0, tBH, 0.08, 150)[0] for _ in range(20000)])
        assert m.max() > 3.0
        # equal weights: about half of the stars from each population (the young one alone
        # gives stars above 2.51 Msun at a rate p; the mixture at a rate p / 2)
        np.random.seed(2)
        p = np.mean([get_star(5.0, tBH, 0.08, 150)[0] > 2.6 for _ in range(40000)])
        set_field_star_mixture([1000.0, 5.0], [1.0, 1.0])
        np.random.seed(3)
        q = np.mean([get_star(1000.0, tBH, 0.08, 150)[0] > 2.6 for _ in range(40000)])
        assert q == pytest.approx(p / 2, rel=0.15)
    finally:
        set_field_star_mixture(None, None)


# ----------------------------------------------------------------------- whole runs
def _run(d, name, *args, chunk=None):
    argv = ['-P', '0', '-BOi', '0', '-RF', name] + list(args)
    if chunk is None:
        subprocess.run([sys.executable, '-m', 'rapster.run_cluster'] + argv, cwd=d, check=True, stdout=subprocess.DEVNULL)
    else:
        code = ('import sys, rapster.cluster_evolution as ce; ce.STREAM_CHUNK_ROWS = %d; '
                'sys.argv = ["rapster"] + %r; from rapster.run_cluster import main; main()' % (chunk, argv))
        subprocess.run([sys.executable, '-c', code], cwd=d, check=True, stdout=subprocess.DEVNULL)
    return os.path.join(d, name)


def _lines(path):
    with open(path) as f:
        return [line for line in f if not line.startswith('#')]


def _check_population_rows(res):
    ev = [l.split() for l in _lines(os.path.join(res, 'evolution.txt'))]
    po = [l.split() for l in _lines(os.path.join(res, 'populations.txt'))]
    times = []
    for r in po:
        if not times or times[-1] != r[0]:
            times.append(r[0])
    assert times == [r[1] for r in ev]                       # one block per evolution row
    M = {}
    for r in po:
        M[r[0]] = M.get(r[0], 0.0) + float(r[4])
    for r in ev:
        assert M[r[1]] == pytest.approx(float(r[K['M_cl']]), rel=1e-12)   # sum_k M_k = M_cl
    return ev, po


def test_run_with_populations():
    with tempfile.TemporaryDirectory() as d:
        _write(d, '300 2e4 0.002 0.5 0.5 2\n700 1e4\n')
        res = _run(d, 'R', '-N', '5e4', '-r', '0.8', '-tM', '1500', '-S', '3', '-POP', 'pops.txt')
        log = open(os.path.join(res, 'log.txt')).read()
        assert 'STELLAR POPULATION 1 BORN AT t = 300 Myr' in log and 'STELLAR POPULATION 2 BORN AT t = 700 Myr' in log
        ev, po = _check_population_rows(res)
        t = np.array([float(r[1]) for r in ev])
        assert 300.0 in t and 700.0 in t                       # steps land on the birth times
        k_rows = [r for r in po if r[1] == '1.0']
        assert float(k_rows[0][0]) == 300.0 and float(k_rows[0][6]) == 0.5
        x = [float(r[6]) for r in k_rows]
        assert all(b >= a for a, b in zip(x, x[1:])) and x[-1] < 1   # mixing towards x = 1
        gas = np.array([float(r[K['M_gas']]) for r in ev])
        i = int(np.searchsorted(t, 300.0))
        assert gas[i - 1] == 0 and gas[i] == pytest.approx(2e4 * M_AVG0, rel=1e-9)   # SFE 0.5


def test_streamed_population_rows_match_unstreamed():
    with tempfile.TemporaryDirectory() as d:
        _write(d, '200 3e4 0.001 1 0.6\n')
        args = ('-N', '5e4', '-r', '0.7', '-tM', '1200', '-S', '5', '-POP', 'pops.txt')
        a = _run(d, 'A', *args)
        b = _run(d, 'B', *args, chunk=7)
        for name in ('populations.txt', 'evolution.txt', 'mergers.txt', 'hardening.txt', 'tdes.txt'):
            assert open(os.path.join(a, name)).read() == open(os.path.join(b, name)).read(), name
        assert not [f for f in os.listdir(b) if f.endswith('.part')]


def test_population_born_after_t_max_changes_nothing():
    with tempfile.TemporaryDirectory() as d:
        _write(d, '1e6 1e5 0.001 0.5 0.3\n')
        args = ('-N', '1e5', '-r', '0.6', '-tM', '3000', '-S', '8')
        a = _run(d, 'A', *args)
        b = _run(d, 'B', *args, '-POP', 'pops.txt')
        for name in ('evolution.txt', 'mergers.txt', 'hardening.txt', 'tdes.txt'):
            assert open(os.path.join(a, name)).read() == open(os.path.join(b, name)).read(), name
        _check_population_rows(b)


def test_compact_objects_of_a_later_population():
    # metal-rich first population in a small cluster: no BH is retained; the BH subsystem
    # waits for the BHs of the metal-poor population born at 300 Myr (they form at 310 Myr)
    with tempfile.TemporaryDirectory() as d:
        _write(d, '300 1e5 0.0001\n')
        res = _run(d, 'R', '-N', '5e3', '-r', '1', '-Z', '0.02', '-tM', '1000', '-S', '1', '-POP', 'pops.txt')
        ev = np.loadtxt(os.path.join(res, 'evolution.txt'), ndmin=2)
        t, NBH = ev[:, K['t']], ev[:, K['N_BH']]
        assert (NBH[t < 310] == 0).all() and NBH[-1] > 0
        assert 'BH SUBSYSTEM EVAPORATED' not in open(os.path.join(res, 'log.txt')).read()


def test_split_population_matches_single_population():
    # N = 2e5 stars born at t = 0, or 1e5 at t = 0 and 1e5 (well mixed) at t = 0.1 Myr: the
    # same stellar mass and mean stellar mass to within a percent
    with tempfile.TemporaryDirectory() as d:
        _write(d, '0.1 1e5\n')
        a = np.loadtxt(os.path.join(_run(d, 'A', '-N', '2e5', '-tM', '2000', '-S', '11'), 'evolution.txt'))
        b = np.loadtxt(os.path.join(_run(d, 'B', '-N', '1e5', '-tM', '2000', '-S', '11', '-POP', 'pops.txt'), 'evolution.txt'))
        for t0 in (1, 5, 100, 500, 1000, 1900):
            ra = a[np.searchsorted(a[:, K['t']], t0)]
            rb = b[np.searchsorted(b[:, K['t']], t0)]
            assert rb[K['M_cl']] == pytest.approx(ra[K['M_cl']], rel=0.02)
            assert rb[K['m_avg']] == pytest.approx(ra[K['m_avg']], rel=0.02)
            assert rb[K['n_star']] == pytest.approx(ra[K['n_star']], rel=0.25)


def test_populations_option_validation():
    from rapster.run_cluster import parse_args
    with tempfile.TemporaryDirectory() as d:
        bad = _write(d, '0 1e4\n')
        argv = sys.argv
        try:
            sys.argv = ['rapster', '-POP', bad]
            with pytest.raises(SystemExit):
                parse_args()
            sys.argv = ['rapster', '-POP', os.path.join(d, 'missing.txt')]
            with pytest.raises(SystemExit):
                parse_args()
            good = _write(d, '50 1e4 # Z from -Z\n')
            sys.argv = ['rapster', '-POP', good, '-Z', '0.003']
            c = parse_args()
            assert c['populations'] == [dict(t_birth=50.0, N=1e4, Z=0.003, SFE=1.0, x=1.0, n_mix=P.DEFAULT_N_MIX)]
        finally:
            sys.argv = argv
