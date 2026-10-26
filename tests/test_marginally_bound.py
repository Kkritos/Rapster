"""Marginally bound radius of TDE orbits (x_mb, R_mb, fill_R_mb in functions.py)."""
import numpy as np
import pytest
from scipy.optimize import brentq

from rapster.functions import x_mb, R_mb, fill_R_mb, R_MB_PENDING
from rapster.constants import G_Newton, c_light


def _x_mb_direct(a, iota):
    """Independent reference: the marginally bound radius as the root of
        x^4 - 4x^3 - b^2 (1 - 3 cos^2 psi) x^2 + b^4 cos^2 psi + 4 b sin(psi) sqrt(x^5 - b^2 x^3 cos^2 psi) = 0
    (the equation of the former lookup table, Data/.marginally_bound.py, in which
    psi = pi/2 is the equatorial plane and b < 0 a retrograde orbit), i.e.
    psi = pi/2 - iota, b = a for prograde orbits (iota <= pi/2) and
    psi = iota - pi/2, b = -a for retrograde ones; bracketed between the
    equatorial radii and solved with brentq."""
    if iota <= np.pi / 2:
        psi, b = np.pi / 2 - iota, a
    else:
        psi, b = iota - np.pi / 2, -a
    c2, sp = np.cos(psi) ** 2, np.sin(psi)

    def f(x):
        return (x**4 - 4 * x**3 - b**2 * (1 - 3 * c2) * x**2 + b**4 * c2
                + 4 * b * sp * np.sqrt(x**5 - b**2 * x**3 * c2))

    lo = 2 - a + 2 * np.sqrt(1 - a) - 1e-9
    hi = 2 + a + 2 * np.sqrt(1 + a) + 1e-9
    return brentq(f, lo, hi, xtol=1e-14, rtol=1e-14)


@pytest.mark.parametrize('a', [0.1, 0.3, 0.7, 0.9, 0.998])
@pytest.mark.parametrize('iota', [0.05, 0.4, 0.8, 1.2, 1.5, 1.65, 2.0, 2.6, 3.09])
def test_against_direct_solution(a, iota):
    assert x_mb(a, iota) == pytest.approx(_x_mb_direct(a, iota), rel=1e-8)


def test_limits():
    a = np.array([0.1, 0.5, 0.9, 0.998])
    # equatorial orbits: r_mb = 2 -+ a + 2 sqrt(1 -+ a) (Bardeen, Press & Teukolsky 1972)
    assert np.allclose(x_mb(a, 0.0), 2 - a + 2 * np.sqrt(1 - a), rtol=1e-12)
    assert np.allclose(x_mb(a, np.pi), 2 + a + 2 * np.sqrt(1 + a), rtol=1e-12)
    # Schwarzschild: 4 GM/c^2 for every inclination
    assert np.all(x_mb(np.zeros(5), np.linspace(0, np.pi, 5)) == 4.0)
    assert x_mb(1e-9, 1.0) == pytest.approx(4.0, abs=1e-6)
    # prograde orbits are tighter than retrograde ones, increasing monotonically with iota
    iota = np.linspace(0, np.pi, 200)
    assert np.all(np.diff(x_mb(np.full(200, 0.7), iota)) > 0)


def test_scalar_and_array_agree():
    rng = np.random.default_rng(2)
    a = rng.uniform(0, 0.998, 500)
    iota = np.arccos(rng.uniform(-1, 1, 500))
    xs = np.array([x_mb(a[i], iota[i]) for i in range(500)])
    assert np.array_equal(xs, x_mb(a, iota))
    assert np.all(np.isfinite(xs))
    assert isinstance(x_mb(0.5, 1.0), float)
    assert R_mb(10.0, 0.5, 1.0) == pytest.approx(x_mb(0.5, 1.0) * G_Newton * 10.0 / c_light**2)


def test_fill_R_mb():
    rng = np.random.default_rng(3)
    t = np.zeros((6, 18))
    t[:, 6] = rng.uniform(5, 50, 6)
    t[:, 7] = rng.uniform(0, 0.99, 6)
    t[:, 12] = [0.0, 0.5, 1.5, 1.6, 2.5, np.pi]       # prograde and retrograde orbits
    t[:, 13] = R_MB_PENDING
    t[5, 13] = 0.0                                    # row of a BBH-star TDE (no r_mb)
    out = fill_R_mb(t.copy())
    assert out[5, 13] == 0.0
    assert np.all(np.isfinite(out[:5, 13])) and np.all(out[:5, 13] > 0)
    assert np.allclose(out[:5, 13], R_mb(t[:5, 6], t[:5, 7], t[:5, 12]), rtol=0, atol=0)
