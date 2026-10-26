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

from .constants import *
import math

# the compiled interpolation that np.interp calls for real-valued data (same results):
try:
    from numpy._core.multiarray import interp as _compiled_interp      # NumPy >= 2
except ImportError:                                                    # NumPy 1.x
    from numpy.core.multiarray import interp as _compiled_interp

# Get the directory where constants.py lives (the 'rapster' folder):
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up one level to the root, then into the Data folder:
DATA_PATH_Planck18 = os.path.join(BASE_DIR, '..', 'Data', 'Planck18_lookup_table.npz')

# load Planck18 lookup tables:
Planck18_lookup_table = np.load(DATA_PATH_Planck18)


# interpolate redshift-lookback and lookback-redshift relations:
# Build the interpolators once over the Planck18 lookup table; the public
# wrappers below add a domain check so out-of-range inputs raise a clear error
# instead of silently extrapolating.
_lookback_interp = interpolate.interp1d(Planck18_lookup_table['z'], Planck18_lookup_table['lookback'])
_redshift_interp = interpolate.interp1d(Planck18_lookup_table['t'], Planck18_lookup_table['redshift'])

# tabulated domains (set by the Planck18 lookup table):
_Z_MIN, _Z_MAX = float(Planck18_lookup_table['z'].min()), float(Planck18_lookup_table['z'].max())
_T_MIN, _T_MAX = float(Planck18_lookup_table['t'].min()), float(Planck18_lookup_table['t'].max())


# memo of scalar lookback times: the main loop asks for lookback_interp(zCl_form),
# always with the same redshift, tens of thousands of times per run:
_lookback_memo = {}

def lookback_interp(z):
    """
    Lookback time [Myr] at redshift z, via interpolation over the Planck18
    lookup table (fast replacement for lookback()).

    @in z: redshift (scalar or array)

    @out: lookback time [Myr]
    """
    scalar = type(z) is float or type(z) is np.float64
    if scalar:
        cached = _lookback_memo.get((type(z), z))
        if cached is not None:
            return cached.copy()   # same 0-d array the interpolator returns
    if np.any(np.asarray(z) < _Z_MIN) or np.any(np.asarray(z) > _Z_MAX):
        raise ValueError(
            f"lookback_interp: redshift z={z} is outside the Planck18 lookup "
            f"table range z in [{_Z_MIN}, {_Z_MAX}]."
        )
    out = _lookback_interp(z)
    if scalar and z == z:          # skip NaN, which is never equal to itself
        _lookback_memo[(type(z), z)] = out.copy()
    return out


def redshift_interp(t):
    """
    Redshift at lookback time t [Myr], via interpolation over the Planck18
    lookup table (fast replacement for redshift()).

    @in t: lookback time [Myr] (scalar or array)

    @out: redshift
    """
    if np.any(np.asarray(t) < _T_MIN) or np.any(np.asarray(t) > _T_MAX):
        raise ValueError(
            f"redshift_interp: lookback time t={t} Myr is outside the Planck18 "
            f"lookup table range t in [{_T_MIN}, {_T_MAX}] Myr."
        )
    return _redshift_interp(t)

# Marginally bound (parabolic, E = 1) spherical orbits around a Kerr BH:
#
# With G = c = M = 1, E = 1 and x = L - a (L: axial angular momentum, Q: Carter
# constant), the radial potential is R(r) = (r^2 - a x)^2 - Delta (r^2 + x^2 + Q),
# Delta = r^2 - 2r + a^2. The conditions R = R' = 0 of a marginally bound spherical
# orbit give, for each r,
#     u = r^2 - a x = Delta sqrt(r) / (sqrt(r) - 1),
#     x = (r^2 - u) / a,   L = x + a,   Q = u^2 / Delta - r^2 - x^2,
# and the orbital inclination cos(iota) = L / sqrt(L^2 + Q) (iota = 0: prograde
# equatorial, iota = pi: retrograde equatorial). cos(iota) decreases monotonically
# from 1 at the prograde equatorial radius 2 - a + 2 sqrt(1 - a) to -1 at the
# retrograde one 2 + a + 2 sqrt(1 + a), so r_mb(a, iota) is the unique root
# between the two (bisection, then Newton steps inside the bracket).

def _cos_inclination_mb(x, a):
    """cos(iota) of the marginally bound spherical orbit of radius x (units GM/c^2)
    around a Kerr BH of spin 0 < a < 1."""
    sq = np.sqrt(x)
    Delta = x*x - 2.0*x + a*a
    u = Delta * sq / (sq - 1.0)
    xL = (x*x - u) / a
    L = xL + a
    Q = np.maximum(u*u / Delta - x*x - xL*xL, 0.0)   # (>= 0; rounding at the equatorial ends)
    return L / np.sqrt(L*L + Q)

def _x_mb_spinning(a, c):
    """x_mb for arrays of spins 0 < a < 1 and cos(iota) = c: bisection to a bracket
    of ~1e-5 GM/c^2, then Newton steps kept inside the bracket (converged to
    rounding, ~1e-15 relative)."""
    lo = 2.0 - a + 2.0*np.sqrt(1.0 - a)            # prograde equatorial (cos(iota) = 1)
    hi = 2.0 + a + 2.0*np.sqrt(1.0 + a)            # retrograde equatorial (cos(iota) = -1)
    for _ in range(19):
        mid = 0.5*(lo + hi)
        above = _cos_inclination_mb(mid, a) > c
        lo = np.where(above, mid, lo)
        hi = np.where(above, hi, mid)
    x = 0.5*(lo + hi)
    for _ in range(4):
        h = 1e-7 * x
        f = _cos_inclination_mb(x, a) - c
        df = (_cos_inclination_mb(x + h, a) - _cos_inclination_mb(x - h, a)) / (2.0*h)
        with np.errstate(divide='ignore', invalid='ignore'):
            step = np.where(df < 0.0, f / df, 0.0)
        x = np.clip(x - step, lo, hi)
    return x

def x_mb(a_BH, iota):
    """
    Radius of the marginally bound orbit in units of the gravitational radius GM/c^2,
    for spin magnitude a_BH (0 <= a_BH < 1) and orbital inclination iota (radians,
    0 <= iota <= pi, measured from the BH spin axis: cos(iota) > 0 is prograde).
    Exact Kerr result (see the notes above); arrays are evaluated elementwise.
    """
    a = np.asarray(a_BH, dtype=float)
    c = np.cos(np.asarray(iota, dtype=float))
    a, c = np.broadcast_arrays(a, c)
    out = np.full(a.shape, 4.0)                    # Schwarzschild: 4 GM/c^2 for any iota
    spin = a > 1e-12
    if np.any(spin):
        aa = a[spin]; cc = c[spin]
        res = np.empty(aa.shape)
        for k in range(0, aa.shape[0], 65536):     # (chunks: faster for millions of rows)
            res[k:k+65536] = _x_mb_spinning(aa[k:k+65536], cc[k:k+65536])
        out[spin] = res
    return out if out.ndim else float(out)

def R_mb(M_BH, a_BH, iota):
    """
    Marginally bound orbit radius (pc) of a BH of mass M_BH (Msun), spin a_BH and
    orbital inclination iota (see x_mb).
    """
    rg = G_Newton*np.asarray(M_BH, dtype=float)/c_light**2   # gravitational radius
    return x_mb(a_BH, iota) * rg

# placeholder for r_mb in TDE rows, filled in by fill_R_mb when the output is written:
R_MB_PENDING = -1.0

def fill_R_mb(tdes):
    """Replace the R_MB_PENDING entries of column r_mb (13) of a TDE array by
    R_mb(m_BH, s_BH, iota) (columns 6, 7, 12), evaluated for all rows at once."""
    rows = tdes[:, 13] == R_MB_PENDING
    if np.any(rows):
        tdes[rows, 13] = R_mb(tdes[rows, 6], tdes[rows, 7], tdes[rows, 12])
    return tdes

def E_cosmo(z):
    """
    @in z: redshift
    
    @out: auxiliary cosmological function
    """
    
    return np.sqrt(Omega_R * (1 + z)**4 + Omega_M * (1 + z)**3 + Omega_K * (1 + z)**2 + Omega_V)

def lookback(z):
    """
    @in z: redshift
    
    @out: lookback time [Myr]
    """

    tL = 0.0
    zz = 0.0
    dzz = 1e-3
    while zz < z:
        tL += dzz / (1 + zz) / E_cosmo(zz) * t_Hubble
        zz += dzz
    
    return tL

def redshift(t):
    """
    @in t: lookback time [Myr]
    
    @out: redshift
    """
    
    z = 0.0
    
    dz = 1e-3

    while lookback(z) < t:
        z += dz
    
    return z

def min_of(*values):
    """
    Same value as np.min([...]) for a few real scalars, without building an array
    (called several times per BBH interaction). NaN is handed to np.min, so it
    propagates exactly as before.
    """
    m = values[0]
    for v in values[1:]:
        if v < m:
            m = v
        elif not v >= m:   # a NaN is involved
            return np.min(list(values))
    return m

def max_of(*values):
    """
    Same value as np.max([...]) for a few real scalars (see min_of).
    """
    m = values[0]
    for v in values[1:]:
        if v > m:
            m = v
        elif not v <= m:   # a NaN is involved
            return np.max(list(values))
    return m

def delete_at(arr, index):
    """
    Same result as np.delete(arr, index) for a 1-D array and an int index or a
    list of two distinct ints, but built from slices, without np.delete's
    overhead (the BH arrays are edited thousands of times per step).
    """
    if arr.ndim == 1:
        n = arr.shape[0]
        if isinstance(index, (int, np.integer)):
            k = int(index)
            if 0 <= k < n:
                return np.concatenate((arr[:k], arr[k+1:]))
        elif isinstance(index, list) and len(index) == 2:
            a, b = int(index[0]), int(index[1])
            if 0 <= a < n and 0 <= b < n and a != b:
                a, b = min(a, b), max(a, b)
                return np.concatenate((arr[:a], arr[a+1:b], arr[b+1:]))
    return np.delete(arr, index)

def delete_row(arr, index):
    """
    Same result as np.delete(arr, index, axis=0) for an int index (row arrays
    such as binaries, triples and pairs), built from slices.
    """
    if isinstance(index, (int, np.integer)) and arr.ndim >= 1:
        k = int(index)
        if 0 <= k < arr.shape[0]:
            return np.concatenate((arr[:k], arr[k+1:]), axis=0)
    return np.delete(arr, index, axis=0)

def append_to(arr, values):
    """
    Same result as np.append(arr, values) for a 1-D array.
    """
    return np.concatenate((arr, np.ravel(values)))

def append2_to(arr, value1, value2):
    """
    Same result as append_to(append_to(arr, value1), value2), with one copy of arr.
    """
    return np.concatenate((arr, np.ravel(value1), np.ravel(value2)))

def delete_append(arr, index, value):
    """
    Same result as append_to(delete_at(arr, index), value), with one copy of arr.
    """
    if arr.ndim == 1 and isinstance(index, (int, np.integer)) and 0 <= index < arr.shape[0]:
        k = int(index)
        return np.concatenate((arr[:k], arr[k+1:], np.ravel(value)))
    return append_to(delete_at(arr, index), value)

def T_GW(m1, m2, a0, e0):
    """
    I. Mandel (2021) fit to Peters timescale.
    
    @in m1: primary mass of BBH [Msun]
    @in m2: secondary mass of BBH [Msun]
    @in a0: initial semimajor axis of BBH [pc]
    @in e0: initial eccentricity of BBH
    
    @out: GW coalescence timescale [Myr]
    """
    
    # Called for every BBH interaction: for ordinary values the expressions below are
    # evaluated with Python floats (the same IEEE operations, and the same pow, as with
    # NumPy scalars, so the value is identical) and returned as a NumPy float; other
    # inputs (arrays, e0 outside [0, 1), non-positive masses) take the NumPy path.
    if (type(m1) in _FLOAT_TYPES and type(m2) in _FLOAT_TYPES and type(a0) in _FLOAT_TYPES
            and type(e0) in _FLOAT_TYPES and m1 > 0.0 and m2 > 0.0 and 0.0 <= e0 < 1.0
            and 0.0 <= a0 < 1e60):
        m1 = float(m1); m2 = float(m2); a0 = float(a0); e0 = float(e0)
        Tc = 5 * c_light**5 * a0**4 / 256 / G_Newton**3 / m1 / m2 / (m1 + m2)
        factor_e = (1 + 0.27 * e0**10 + 0.33 * e0**20 + 0.2 * e0**1000) * (1 - e0**2)**(7/2)
        return np.float64(Tc * factor_e)

    # Coalescence time for circular orbits:
    Tc = 5 * c_light**5 * a0**4 / 256 / G_Newton**3 / m1 / m2 / (m1 + m2)
    
    # Eccentricity factor:
    factor_e = (1 + 0.27 * e0**10 + 0.33 * e0**20 + 0.2 * e0**1000) * (1 - e0**2)**(7/2)
    
    return Tc * factor_e

_FLOAT_TYPES = (float, np.float64)

def t_relax(Mcl, rh, m_avg, psi, logL):
    """
    @in Mcl: cluster's mass [Msun]
    @in rh: half-mass radius [pc]
    @in m_avg: average mass [Msun]
    @in psi: multimass relaxation factor
    @in logL: Coulomb logarithm
    
    @out: half-mass relaxation timescale [Myr]
    """
    
    return 0.138 * np.sqrt(Mcl * rh**3 / G_Newton) / m_avg / psi / logL

# Escape-velocity prefactor f, set by the user at runtime (see set_escape_factor):
#   v_esc = f * sqrt(<v_star^2> + <v_BH^2>),  with <v^2> = 0.4 G M / r.
# f = 2 (default): mass-weighted rms escape speed, v_esc^2 = 1.6 G M_cl / r_h.
# f = sqrt(10):    v_esc^2 = 4 G M_cl / r_h (Mai et al. 2025, Eq. 1).

F_ESC = None   # set at runtime by set_escape_factor(); run_cluster.py passes the user's -fesc value

def set_escape_factor(f):
    """
    @in f: escape-velocity prefactor (> 0)
    """
    global F_ESC
    if f <= 0:
        raise ValueError('escape-velocity prefactor must be positive')
    F_ESC = float(f)

def v_esc_cl(v_star, vBH):
    """
    @in v_star: 3D rms velocity of stars [km/s]
    @in vBH: 3D rms velocity of BHs [km/s]
    
    @out: escape velocity used for retention/ejection of BHs and BBHs [km/s]
    """
    
    return F_ESC * np.sqrt(v_star**2 + vBH**2)

def v_esc(Mcl, rh):
    """
    @in Mcl: cluster mass [Msun]
    @in rh: half-mass radius [pc]
    
    @out: escape velocity used for natal-kick retention [km/s]
    """
    
    return F_ESC * np.sqrt(0.4 * G_Newton * Mcl / rh)

def sample_angles():
    """
    @out theta1: [rad]
    @out theta2: [rad]
    @out dPhi: [rad]
    """
    
    cos_theta1 = np.random.uniform(-1, 1)
    cos_theta2 = np.random.uniform(-1, 1)
    
    theta1 = np.arccos(cos_theta1)
    theta2 = np.arccos(cos_theta2)
    
    dPhi = np.random.uniform(0, 2*np.pi)
    
    return theta1, theta2, dPhi

def Rate_3bb(m, n, v):
    """
    @in m: mass scale [Msun]
    @in n: number density [pc^-3]
    @in v: velocity dispersion [km/s]
    
    @out: three-body binary formation rate [1/Myr]
    """
    
    return 8 * np.pi / np.sqrt(3) * n**3 * (G_Newton * m)**5 / v**9 * eta_min**(-11/2) * (1 + 3 * eta_min) * (1 + 6 * eta_min) * P_3bb

def Rate_cap(m, n, v):
    """
    @in m: mass [Msun]
    @in n: number density [pc^-3]
    @in v: velocity dispersion [km/s]
    
    @out: 2-body capture rate [1/Myr]
    """

    # relative velocity:
    v_rel = np.sqrt(2) * v
    
    # capture cross section:
    Sigma_cap = 2 * np.pi * (85 * np.pi / 6 / np.sqrt(2))**(2/7) * G_Newton**2 * (2 * m)**(10/7) * m**(4/7) / c_light**(10/7) / v_rel**(18/7)
    
    return n**2 * Sigma_cap * v_rel

def Rate_int(m, n, v, rp):
    """
    @in m: total mass [Msun]
    @in n: number density [pc^-3]
    @in v: velocity dispersion [km/s]
    @in rp: pericenter of interaction [pc]
    
    @out: interaction rate [1/Myr]
    """
    
    # relative velocity:
    v_rel = np.sqrt(2) * v
    
    # interaction cross section:
    Sigma_int = 2 * np.pi * G_Newton * m * rp / v_rel**2
    
    return n * Sigma_int * v_rel

def Rate_exc(m1, m2, m3, n3, v_rel, a):
    """
    @in m1: to-be-exchanged mass of binary
    @in m2: primary mass of binary
    @in m3: incoming mass
    @in n3: number density of incoming mass
    @in v_rel: binary-single relative velocity before interaction
    @in a: initial binary semimajor axis
    
    @out: 1-2 -> 3-2 exchange rate [1/Myr]
    """

    m12 = m1 + m2
    m123 = m12 + m3
    m13 = m1 + m3
    m23 = m2 + m3
    mu1 = m1 / m12
    mu2 = m3 / m123
    
    # exchange cross section [Heggie et al. (1996)]:
    Sigma_exc = 6.76e-5 * a * (10 / v_rel)**2 * m123 * (m23/m123)**(1/6) * (m3/m13)**(7/2) * (m123/m12)**(1/3) * (m13/m123) \
    * np.exp(3.70 + 7.49*mu1 - 1.89*mu2 - 15.49*mu1**2 - 2.93*mu1*mu2 - 2.92*mu2**2 + 3.07*mu1**3 + 13.15*mu1**2*mu2 - 5.23*mu1*mu2**2 + 3.12*mu2**3)
    
    return n3 * Sigma_exc * v_rel
    
def sample_hardness():
    """
    @out: 3bb hardness parameter
    """
    
    return eta_min * (1 - np.random.rand())**(-2/7)

def CDF_maxwell(w):
    """
    Returns cumulative density function of Maxwellian distribution.

    @in w: velocity value, normalized to the 1D velocity dispersion
    """

    return erf(w/np.sqrt(2)) - np.sqrt(2/np.pi)*w*np.exp(-w**2/2)

w_s = np.linspace(0, 10000, 10**6) # normalized velocities (assuming 1D velo. disp. =1)

def _maxwell_cdf_table(w):
    """CDF_maxwell evaluated on the grid w, identical to np.vectorize(CDF_maxwell)(w).

    Above w ~ 9 the CDF is exactly 1.0 in double precision: erf(w/sqrt(2)) rounds
    to 1 and w*exp(-w^2/2) is below half an ulp of 1 (it underflows to 0 for
    w > 38.6). So the function is evaluated point by point only up to w = 50, and
    the remaining ~995,000 grid points are set to 1.0. This replaces ~4 s of
    start-up time with ~10 ms. The shortcut is used only after checking that the
    function really returns 1.0 at the end of the evaluated range and at sample
    points across the filled range; otherwise every point is evaluated as before.
    """
    n_eval = int(np.searchsorted(w, 50.0, side='right'))
    head = np.vectorize(CDF_maxwell)(w[:n_eval])
    tail_samples = w[n_eval::997]
    tail_is_one = (n_eval < w.size
                   and np.all(head[-500:] == 1.0)
                   and np.all(np.vectorize(CDF_maxwell)(tail_samples) == 1.0)
                   and CDF_maxwell(w[-1]) == 1.0)
    if not tail_is_one:
        return np.vectorize(CDF_maxwell)(w)
    table = np.ones(w.size, dtype=head.dtype)
    table[:n_eval] = head
    return table

CDF_maxwell_s = _maxwell_cdf_table(w_s)

def get_maxwell_sample(sigma):
    """
    Returns a sample from the Maxwellian with the inverse sampling method.

    @in sigma: 1D velocity dispersion parameter
    """
    
    # (the compiled function np.interp calls, without its Python-level checks; same value)
    return sigma * _compiled_interp(np.random.rand(), CDF_maxwell_s, w_s, None, None)

def DW(x, m0, m1, m2, a1, a2, e1, e2, cosi1, cosi2, w):
    """
    Change in the double averaged Hamiltonian of a hierarchical triple as inner eccentricity is varied.
    Notation from M.C.Miller & D.P.Hamilton ApJ 576, 894 (2002).
    
    @in x: 1-e1^2 (see definition of e1 below)
    @in m0: primary mass of inner binary
    @in m1: secondary mass of inner binary
    @in m2: tertiary mass
    @in a1: inner binary semimajor axis
    @in a2: outer binary semimajor axis
    @in e1: inner binary eccentricity
    @in e2: outer binary eccentricity
    @in cosi1: cosine of inner binary's inclination relative to invariant plane
    @in cosi2:    "   of outer    "          "          "     "     "       "
    @in w: argument of periapsis for inner binary
    """
 
    E = 1 - e1**2

    # inclinations:
    i1, i2 = np.arccos(cosi1), np.arccos(cosi2)
 
    sinw = np.sin(w)
    
    M1 = m0 + m1 # total inner binary mass
    M2 = M1 + m2 # total outer binary mass
    
    mu1 = m0*m1/M1 # reduced mass of inner binary
    mu2 = M1*m2/M2 # reduced mass of outer binary

    # integrals of motion:
    B = mu2/mu1*np.sqrt(M1*a1/M2/a2*(1-e2**2))
    A = np.sqrt(E)*np.cos(i1) + B*np.cos(i2)

    # general relativistic correction term at zero inner eccentricity:
    tPN = 8*M1/m2*(a2/a1*np.sqrt(1-e2**2))**3*M1/232/a1/(3e5)**2
    
    cosi = (A**2-B**2-E)/2/B/np.sqrt(E)
    W = - 2*E + E*cosi**2 + 5*(1-E)*sinw**2*(cosi**2-1) + tPN/np.sqrt(E)
    
    cosI = (A**2-B**2-x)/2/B/np.sqrt(x)
    
    return -2*E + E*cosI**2 + 5*(1-E)*(cosI**2-1) + tPN/np.sqrt(x) - W

def _DW_of_x(m0, m1, m2, a1, a2, e1, e2, cosi1, cosi2, w):
    """
    DW as a function of x alone, for the root finder in find_eMAX.

    Every quantity that does not depend on x is computed once here, with the same
    expressions as in DW, instead of at each of the ~50 evaluations the root finder
    makes; the returned function gives values identical to DW(x, ...).
    """

    E = 1 - e1**2
    i1, i2 = np.arccos(cosi1), np.arccos(cosi2)
    sinw = np.sin(w)
    M1 = m0 + m1
    M2 = M1 + m2
    mu1 = m0*m1/M1
    mu2 = M1*m2/M2
    B = mu2/mu1*np.sqrt(M1*a1/M2/a2*(1-e2**2))
    A = np.sqrt(E)*np.cos(i1) + B*np.cos(i2)
    tPN = 8*M1/m2*(a2/a1*np.sqrt(1-e2**2))**3*M1/232/a1/(3e5)**2
    cosi = (A**2-B**2-E)/2/B/np.sqrt(E)
    W = - 2*E + E*cosi**2 + 5*(1-E)*sinw**2*(cosi**2-1) + tPN/np.sqrt(E)

    # x-independent pieces of the expressions in DW (same operation order):
    A2_B2 = A**2-B**2       # (A**2-B**2-x) == (A**2-B**2) - x
    minus_2E = -2*E
    five_1mE = 5*(1-E)      # 5*(1-E)*(cosI**2-1) == (5*(1-E)) * (cosI**2-1)

    def DW_x_numpy(x):
        cosI = (A2_B2-x)/2/B/np.sqrt(x)
        return minus_2E + E*cosI**2 + five_1mE*(cosI**2-1) + tPN/np.sqrt(x) - W

    # The root finder calls this ~50 times per triple with a 1-element array x.
    # The same IEEE operations on Python floats (cosI**2 of an array is cosI*cosI,
    # and sqrt is correctly rounded in both) give identical values without NumPy's
    # per-operation overhead. Division by zero, which NumPy turns into inf/nan
    # but Python raises on, and anything unusual are handed to the NumPy version.
    c_A2_B2 = float(A2_B2); c_B = float(B); c_E = float(E); c_minus_2E = float(minus_2E)
    c_five_1mE = float(five_1mE); c_tPN = float(tPN); c_W = float(W)
    scalars_ok = all(type(v) in (float, np.float64) for v in (A2_B2, B, E, minus_2E, five_1mE, tPN, W))

    def DW_x(x):
        if scalars_ok and type(x) is np.ndarray and x.shape == (1,) and x.dtype == np.float64:
            xv = float(x[0])
            try:
                sx = math.sqrt(xv) if xv >= 0.0 else math.nan
                cosI = (c_A2_B2-xv)/2/c_B/sx
                cosI2 = cosI*cosI
                val = c_minus_2E + c_E*cosI2 + c_five_1mE*(cosI2-1) + c_tPN/sx - c_W
                return np.array([val])
            except (ZeroDivisionError, OverflowError, ValueError):
                pass
        return DW_x_numpy(x)

    # x-independent quantities, as Python floats, for _zlk_smallest_root (None if unusual):
    DW_x.params = (c_E, c_B, c_A2_B2, c_tPN, c_W) if scalars_ok else None
    return DW_x

# root(fun, x0, tol=tol) (default method 'hybr') hands the problem to MINPACK's hybrd
# through several Python layers (argument checks, a call of fun to find the output
# shape, a wrapper counting the evaluations). _root_hybr_x makes the same MINPACK call
# with the same parameters directly (identical result, less overhead per triple);
# if SciPy's internal module is not available, it calls root.
try:
    from scipy.optimize import _minpack as _scipy_minpack
    _hybrd = _scipy_minpack._hybrd
except (ImportError, AttributeError):
    _hybrd = None
_EPS_FLOAT64 = float(np.finfo(np.float64).eps)

def _root_hybr_x(fun, x0, tol):
    """Same as root(fun=fun, x0=x0, tol=tol)['x'] for a function returning float64 values."""
    if _hybrd is None:
        return root(fun=fun, x0=x0, tol=tol)['x']
    x0 = np.asarray(x0).flatten()
    # _root_hybr: xtol=tol, maxfev=200*(n+1), ml=mu=-10 (no band), epsfcn=eps of the
    # output dtype, factor=100, diag=None:
    return _hybrd(fun, x0, (), 1, tol, 200 * (len(x0) + 1), -10, -10, _EPS_FLOAT64, 100, None)[0]

def _zlk_smallest_root(E, B, C, tPN, W):
    """
    Smallest positive root x of DW(x) (see _DW_of_x), or None.

    With y = sqrt(x) and cosI = (C - x) / (2 B y), C = A^2 - B^2,
        DW = K1 cosI^2 + c0 + tPN / y,   K1 = 5 - 4E,   c0 = -2E - 5(1-E) - W,
    and 4 B^2 y^2 DW = P(y) = K1 y^4 + b2 y^2 + b1 y + b0, with b2 = -2 C K1 + 4 B^2 c0,
    b1 = 4 B^2 tPN >= 0, b0 = K1 C^2 >= 0. By Descartes' rule P has two positive roots
    or none; if the local minimum y* of P (found by Newton's method on P') is below
    zero, the smaller root lies in (0, y*), where it is found by safeguarded Newton
    iterations and then polished with Newton steps on DW(x) itself.

    This is the root that the MINPACK root finder (scipy.optimize.root, method 'hybr')
    started at x0 = 1e-10 converges to: for ~40,000 triples from compact-cluster runs
    the two agree to <= 1e-13 (relative), ~4x faster. When there is no positive root
    (P(y*) >= 0), or anything is unusual, None is returned and find_eMAX uses the root
    finder as before (which then stops where it can make no further progress).
    """
    K1 = 5.0 - 4.0*E
    c0 = -2.0*E - 5.0*(1.0 - E) - W
    B2 = B*B
    b2 = -2.0*C*K1 + 4.0*B2*c0
    b1 = 4.0*B2*tPN
    b0 = K1*C*C
    if not (K1 > 0.0 and b2 < 0.0 and b0 > 0.0 and b1 >= 0.0 and B > 0.0):
        return None
    try:
        # local minimum of P (P' = 4 K1 y^3 + 2 b2 y + b1): Newton from the minimum of
        # K1 y^4 + b2 y^2, where P' = b1 >= 0 and P'' > 0, moving left:
        y = math.sqrt(-b2/(2.0*K1))
        for _ in range(60):
            d1 = (4.0*K1*y*y + 2.0*b2)*y + b1
            d2 = 12.0*K1*y*y + 2.0*b2
            if not d2 > 0.0:
                return None
            step = d1/d2
            y_new = y - step
            if not y_new > 0.0:
                return None
            done = abs(step) <= 1e-15*y
            y = y_new
            if done:
                break
        else:
            return None
        if not ((K1*y*y + b2)*y + b1)*y + b0 < 0.0:
            return None          # no positive root
        # the root in (0, y*): P(lo) > 0 > P(hi)
        lo, hi = 0.0, y
        v = 0.5*y
        for _ in range(200):
            pv = ((K1*v*v + b2)*v + b1)*v + b0
            if pv > 0.0:
                lo = v
            else:
                hi = v
            dp = (4.0*K1*v*v + 2.0*b2)*v + b1
            v_new = v - pv/dp if dp != 0.0 else 0.5*(lo + hi)
            if not lo < v_new < hi:
                v_new = 0.5*(lo + hi)
            done = abs(v_new - v) <= 1e-15*v
            v = v_new
            if done:
                break
        x = v*v
        # polish on DW(x) (Newton steps in x):
        for _ in range(2):
            sx = math.sqrt(x)
            cI = (C - x)/2/B/sx
            f = K1*cI*cI + c0 + tPN/sx
            df = 2.0*K1*cI*(-(C + x)/(4.0*B*x*sx)) - tPN/(2.0*x*sx)
            if not df != 0.0:
                break
            x_new = x - f/df
            if not x_new > 0.0:
                break
            x = x_new
        return x if 0.0 < x < math.inf else None
    except (ZeroDivisionError, OverflowError, ValueError):
        return None

def find_eMAX(m0, m1, m2, a1, a2, e1, e2, cosi1, cosi2, w, x0=1e-10, tol=1e-10):
    """
    Maximum eccentricity during a von Zeipel-Lidov-Kozai cycle.
    Same inputs as function DW in this file.
    """
 
    fun = _DW_of_x(m0, m1, m2, a1, a2, e1, e2, cosi1, cosi2, w)
    # the smallest positive root, which the root finder started at x0 = 1e-10 finds
    # (see _zlk_smallest_root); the root finder itself when there is none:
    x_min = None
    if x0 == 1e-10 and tol == 1e-10 and fun.params is not None:
        x_min = _zlk_smallest_root(*fun.params)
    if x_min is None:
        x_min = _root_hybr_x(fun, x0, tol)[0]

    # physically motivated floor: pericenter must exceed 4*GM/c^2 (a few Schwarzschild radii),
    # giving a minimum value of x_min = epsilon_min = 8*G*(m0+m1)/(c^2*a1):
    x_min_floor = 8 * G_Newton * (m0 + m1) / c_light**2 / a1
    x_min = np.clip(x_min, x_min_floor, 1.0)

    eMAX = np.sqrt(1 - x_min)
    
    return eMAX

def X_isco(x, s):
    """
    Returns the isco in units of gravitational radius (G*m/c^2).
    Inputs:
    x: dimensionless spin parameter in [0, 1]
    s: 1 for prograde and -1 for retrograde orbits.
    """
    
    x = np.clip(x, 0.0, THORNE_SPIN_LIMIT)
    Z1 = 1 + (1 - x**2)**(1/3)*((1+x)**(1/3) + (1-x)**(1/3))
    Z2 = np.sqrt(3*x**2 + Z1**2)
    XX = 3 + Z2 - s*np.sqrt((3-Z1)*(3+Z1+2*Z2))
    return XX

def dxdM(M, x, s):
    """
    Bardeen formula; Bardeen 1970 ("Kerr Metric Black Holes")
    M: mass
    x: dimensionless spin parameter in [0, 1]
    """
    
    return 2/3/np.sqrt(3)*s/M*(1 + 2*np.sqrt(3*X_isco(x, s) - 2))/np.sqrt(1 - 2/3/X_isco(x, s)) - 2*x/M

def evolve_spin_RungeKutta(Mi, Mf, xi, s, dM):
    """
    Returns the final spin of an accretion episode. Uses Runge-Kutta's 4th-order method.
    Mi: initial compact object mass
    Mf: final compact object mass so that (Mf-Mi)=accreted mass
    xi: initial compact object dimensionless spin parameter in [0, 1]
    s: +1 for prograde and -1 for retrograde accretion
    dM: numerical integration mass-step
    """
    
    M = Mi
    x = xi
    while M<Mf:
        k1 = dM*dxdM(M, x, s)
        k2 = dM*dxdM(M+dM/2, x+k1/2, s)
        k3 = dM*dxdM(M+dM/2, x+k2/2, s)
        k4 = dM*dxdM(M+dM, x+k3, s)
        dx = (k1 + 2*k2 + 2*k3 + k4)/6
        if x + dx > THORNE_SPIN_LIMIT:
            x = THORNE_SPIN_LIMIT
            break
        elif x + dx < 0 and s==-1:
            s = +1
        else:
            x = x + dx
            M = M + dM
    xf = x;
    return xf

# End of file.
