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

import math

# the compiled interpolation that np.interp calls for real-valued data (same results):
try:
    from numpy._core.multiarray import interp as _compiled_interp      # NumPy >= 2
except ImportError:                                                    # NumPy 1.x
    from numpy.core.multiarray import interp as _compiled_interp

from .constants import *

# ---------------------------------------------------------------------------
# Constants (SI) used internally for physics computations
# ---------------------------------------------------------------------------
Msun = 1.989e30   # kg
c    = 2.998e8    # m/s
km   = 1e3        # m  (conversion: 1 km = 1e3 m)
G    = 6.674e-11  # m^3 kg^-1 s^-2

# Geometric mass conversion: G*Msun/c^2 in km
Msun_to_km = G * Msun / c**2 / km    # about 1.477 km / M_sun  (= G*Msun/c^2 in km)

# Thorne (1974) spin limit: radiation captured by the BH prevents chi=1
THORNE_LIMIT = THORNE_SPIN_LIMIT

# ---------------------------------------------------------------------------
# EOS tables  (mass in M_sun, radius in km)
# ---------------------------------------------------------------------------
# Get the directory where this module lives (the 'rapster' folder):
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up one level to the root, then into the Data folder:
DATA_DIR = os.path.join(BASE_DIR, '..', 'Data')

# supported NS equations of state: name -> two-column table (radius [km], mass
# [Msun]) in Data/, stable branch only, mass strictly increasing:
EOS_TABLES = {
    'APR':  'APR.txt',
    'AU':   'AU.txt',
    'APR4': 'APR4.txt',
    'MPA1': 'MPA1.txt',
    'DD2':  'DD2.txt',
}

# ---------------------------------------------------------------------------
# Supramassive limit  -  Breu & Rezzolla (2016)
#   M_max,rot ~ 1.203 * M_TOV  (EOS-insensitive to ~2%)
# ---------------------------------------------------------------------------
BREU_REZZOLLA_FACTOR = 1.203

# lowest NS mass accepted, shared by every EoS: a physical floor rather than
# wherever a given table happens to start:
M_NS_FLOOR = 1.0

class EquationOfState:
    """Mass-radius relation and mass limits for one nuclear equation of state."""
    __slots__ = ('name', 'M', 'R', 'M_min', 'M_TOV', 'M_max_rot', 'R_TOV')

    def __init__(self, name, filename):
        R, M = np.loadtxt(os.path.join(DATA_DIR, filename), unpack=True)
        if np.any(np.diff(M) <= 0):
            raise ValueError(f"EoS table '{filename}' is not strictly increasing in mass.")
        self.name      = name
        self.M         = M
        self.R         = R
        self.M_min     = M_NS_FLOOR
        self.M_TOV     = float(M.max())
        self.M_max_rot = BREU_REZZOLLA_FACTOR * self.M_TOV
        self.R_TOV     = float(np.interp(self.M_TOV, M, R))

    def radius(self, M):
        """Circumferential radius [km] at gravitational mass M [Msun]."""
        # (the compiled function np.interp calls, without its Python-level checks; same value)
        return _compiled_interp(M, self.M, self.R, None, None)

EOS_REGISTRY = {name: EquationOfState(name, f) for name, f in EOS_TABLES.items()}

def get_eos(name):
    """Return the EquationOfState registered under `name`."""
    try:
        return EOS_REGISTRY[name]
    except KeyError:
        raise ValueError(
            f"Unknown EoS '{name}'. Available: {', '.join(sorted(EOS_TABLES))}."
        ) from None

# ---------------------------------------------------------------------------
# Geometric mass
# ---------------------------------------------------------------------------

def geom_mass(M):
    """
    Geometric mass  M^ = G M / c^2  in km.

    Parameters
    ----------
    M : float - gravitational mass [M_sun]

    Returns
    -------
    Mg : float [km]
    """
    return M * Msun_to_km

# ---------------------------------------------------------------------------
# Moment of inertia  -  Lattimer & Schutz (2005)
#
#   I / (M R^2) = 0.237 (1 + 0.674 C + 4.48 C^4),   C = M/R
#
#   Corrected coefficients: a1 = 0.674, a4 = 4.48
# ---------------------------------------------------------------------------

def moment_of_inertia(M, R):
    """
    Moment of inertia of a neutron star.
    Lattimer & Schutz (2005).

    Parameters
    ----------
    M : float - mass [M_sun]
    R : float - radius [km]

    Returns
    -------
    I : float [M_sun km^2]
    """
    C = M / R
    return 0.237 * M * R**2 * (1.0 + 0.674*C + 4.48*C**4)


def dIdM_partial(M, R):
    """
    Partial dI/dM at fixed R.
    Expanding I = 0.237 [M R^2 + 0.674 M^2 R + 4.48 M^5/R^2]:
        dI/dM = 0.237 [R^2 + 2*0.674 M R + 5*4.48 M^4/R^2]

    Returns
    -------
    dI/dM : float [km^2]
    """
    return 0.237 * (
        R**2
        + 2.0 * 0.674 * M * R
        + 5.0 * 4.48  * M**4 / R**2
    )


def dIdR_partial(M, R):
    """
    Partial dI/dR at fixed M.
    Expanding I = 0.237 [M R^2 + 0.674 M^2 R + 4.48 M^5/R^2]:
        dI/dR = 0.237 M [2 R + 0.674 M - 4*4.48 M^4/R^3]

    Returns
    -------
    dI/dR : float [M_sun km]
    """
    return 0.237 * M * (
        2.0 * R
        + 0.674 * M
        - 4.0 * 4.48 * M**4 / R**3
    )

# ---------------------------------------------------------------------------
# chi <-> frequency
# ---------------------------------------------------------------------------

def chi_from_frequency(M, R, f):
    """
    Dimensionless spin parameter chi = Jc / (GM^2) from spin frequency.

    Parameters
    ----------
    M : float - mass [M_sun]
    R : float - radius [km]
    f : float - spin frequency [Hz]

    Returns
    -------
    chi : float
    """
    I_SI = moment_of_inertia(M, R) * Msun * km**2   # kg m^2
    M_SI = M * Msun                                  # kg
    return c * I_SI * 2.0 * np.pi * f / (G * M_SI**2)


def f_from_chi(M, chi, R=None, NS=True):
    """
    Spin frequency [Hz] from chi.

    NS=True  : f = chi G M^2 / (2pi c I)                         R required
    NS=False : f_H = c^3 chi / (4pi G M (1 + sqrt(1 - chi^2)))

    Parameters
    ----------
    M   : float - mass [M_sun]
    chi : float - dimensionless spin parameter
    R   : float - radius [km], required if NS=True
    NS  : bool  - True for neutron star, False for black hole

    Returns
    -------
    f : float [Hz]
    """
    M_SI = M * Msun
    if NS:
        if R is None:
            raise ValueError("R must be provided for NS (NS=True).")
        I_SI = moment_of_inertia(M, R) * Msun * km**2
        return chi * G * M_SI**2 / (2.0 * np.pi * c * I_SI)
    else:
        chi = min(max(chi, 0.0), THORNE_LIMIT)
        return c**3 * chi / (4.0 * np.pi * G * M_SI
                             * (1.0 + np.sqrt(1.0 - chi**2)))

# ---------------------------------------------------------------------------
# ISCO radius
# ---------------------------------------------------------------------------

def R_ISCO(M, NS=True, R=None, f=None, chi=None, prograde=True):
    """
    ISCO radius [km].

    NS=True  : Luk & Lin (2018) Eq. (6) universal NS relation.  f required.
               No Schwarzschild floor applied - the fit is self-consistent.
    NS=False : Bardeen, Press & Teukolsky (1972) exact Kerr.     chi required.
    R        : accepted but unused; kept for uniform call signature.

    Parameters
    ----------
    M        : float - mass [M_sun]
    NS       : bool  - True for NS, False for BH
    R        : float - radius [km], unused
    f        : float - spin frequency [Hz], required if NS=True
    chi      : float - dimensionless spin, required if NS=False
    prograde : bool  - prograde orbit (default True)

    Returns
    -------
    r_isco : float [km]
    """
    if NS:
        if f is None:
            raise ValueError("f must be provided for NS (NS=True).")
        a1 =  8.809
        a2 = -9.166e-4
        a3 =  8.787e-8
        a4 = -6.019e-12
        if f <= 0.0:
            # f -> 0 limit of the polynomial: y/f -> a1*M
            return a1 * M
        x  = M * f
        y  = a1*x + a2*x**2 + a3*x**3 + a4*x**4
        return y / f
    else:
        if chi is None:
            raise ValueError("chi must be provided for BH (NS=False).")
        chi  = min(max(chi, 0.0), THORNE_LIMIT)
        Mg   = geom_mass(M)
        Z1   = 1.0 + (1.0 - chi**2)**(1.0/3.0) * (
                   (1.0 + chi)**(1.0/3.0) + (1.0 - chi)**(1.0/3.0)
               )
        Z2   = np.sqrt(3.0 * chi**2 + Z1**2)
        sign = +1.0 if prograde else -1.0
        return Mg * (3.0 + Z2 - sign * np.sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0*Z2)))

# ---------------------------------------------------------------------------
# Kepler / Thorne-limit frequency and chi
# ---------------------------------------------------------------------------

def f_kepler(M, NS=True, R=None):
    """
    Mass-shedding (NS) or Thorne-limit (BH) frequency [Hz].

    NS=True  : Lattimer & Prakash (2007)  f_K ~ 1045 sqrt(M) (10/R)^1.5
               R required.
    NS=False : horizon frequency at chi = THORNE_LIMIT.

    Parameters
    ----------
    M  : float - mass [M_sun]
    NS : bool  - True for NS, False for BH
    R  : float - radius [km], required if NS=True

    Returns
    -------
    f : float [Hz]
    """
    if NS:
        if R is None:
            raise ValueError("R must be provided for NS (NS=True).")
        return 1045.0 * np.sqrt(M) * (10.0 / R)**1.5
    else:
        return f_from_chi(M, THORNE_LIMIT, NS=False)


def chi_kepler(M, NS=True, R=None):
    """
    Spin parameter at mass-shedding (NS) or THORNE_LIMIT (BH).

    Parameters
    ----------
    M  : float - mass [M_sun]
    NS : bool  - True for NS, False for BH
    R  : float - radius [km], required if NS=True

    Returns
    -------
    chi : float
    """
    if NS:
        if R is None:
            raise ValueError("R must be provided for NS (NS=True).")
        return chi_from_frequency(M, R, f_kepler(M, NS=True, R=R))
    else:
        return THORNE_LIMIT

# ---------------------------------------------------------------------------
# Specific energy and angular momentum at the ISCO
#
# NS:  Bejger, Zdunik & Haensel (2010)  [A&A 520, A16]
#      Frame-dragging correction to Schwarzschild orbital velocity.
#      All intermediate lengths in km; only J_SI in SI to compute Nphi.
#
# BH:  Bardeen, Press & Teukolsky (1972)  [ApJ 178, 347]  Eqs. (2.12)-(2.13)
#
#      Numerically stable simplified form (avoids cancellation near chi->1):
#
#        rt   = r_isco / Mg           (dimensionless)
#        sign = +1 prograde, -1 retrograde
#
#        E~ = (1 - 2/rt + sign*chi/rt^1.5) / sqrt(1 - 3/rt + 2*sign*chi/rt^1.5)
#
#        J~ [km] = sign * Mg^0.5 * (r^2 - 2a*sqrt(Mg*r) + a^2)
#                  / (r^0.75 * sqrt(r^1.5 - 3*Mg*r^0.5 + 2*sign*a*Mg^0.5))
#
#      Sign convention:
#        sign = +1  prograde  -> E~ < sqrt(8/9),  J~ > 0  (spin-up)
#        sign = -1  retrograde -> E~ > sqrt(8/9), J~ < 0  (spin-down)
#
#      J~/Mg = sign*(rt^2 - 2*sign*chi*sqrt(rt) + chi^2)
#              / (rt^0.75 * sqrt(rt^1.5 - 3*sqrt(rt) + 2*sign*chi))
# ---------------------------------------------------------------------------

def _E_J_isco_NS(M, R, f, prograde=True):
    """
    Specific energy (dimensionless) and angular momentum [km] at NS ISCO.
    Bejger, Zdunik & Haensel (2010), extended to retrograde orbits.

    For prograde orbits frame-dragging aids the orbital motion:
        v = (Omega_Schw - N^phi) * r / N
    For retrograde orbits it opposes it:
        v = (Omega_Schw + N^phi) * r / N
    and J~ < 0 (angular momentum anti-aligned with spin).

    The Luk & Lin (2018) ISCO fit is for prograde orbits; for retrograde
    we use the same ISCO radius as a reasonable approximation (the
    retrograde ISCO is larger, but Luk & Lin do not provide a retrograde
    fit for NSs, so this is the best available approximation).
    """
    sign = +1.0 if prograde else -1.0

    Mg   = geom_mass(M)                                # km
    r    = R_ISCO(M, NS=True, f=f)                     # km  (prograde fit)
    r_S  = 2.0 * Mg                                    # km

    # guard: ISCO must be outside Schwarzschild radius
    if r <= r_S:
        r = r_S * 1.001

    I_SI = moment_of_inertia(M, R) * Msun * km**2      # kg m^2
    J_SI = I_SI * 2.0 * np.pi * f                      # kg m^2 s^-1

    r_m      = r * km                                   # m
    Nphi     = 2.0 * G * J_SI / (r_m**3 * c**3)        # m^-1 (geometric)
    Omega_S  = np.sqrt(Mg * km) / r_m**1.5             # m^-1 (Schwarzschild)
    N        = np.sqrt(1.0 - r_S / r)                  # dimensionless lapse

    # prograde: frame-dragging aids orbit (subtract N^phi)
    # retrograde: frame-dragging opposes orbit (add N^phi)
    v     = min(max((Omega_S - sign * Nphi) * r_m / N, 0.0), 0.9999)
    gamma = 1.0 / np.sqrt(1.0 - v**2)

    E = (N + Nphi * v * r_m) * gamma                   # dimensionless (always > 0)
    J = sign * v * r_m * gamma / km                    # km  (positive prograde,
                                                        #      negative retrograde)
    return E, J


# The NS branch of evolve_spin_during_accretion evaluates the helpers above once per
# mass step (~100 times per TDE onto a NS). _ns_spin_state and _ns_E_J_isco repeat
# their expressions with Python floats and math.sqrt, operation by operation (same
# IEEE results as NumPy scalars and np.sqrt), without the function-call and
# NumPy-scalar overhead; tests/test_exact_speedups.py checks that they agree exactly.

def _ns_spin_state(M, chi, R):
    """(f, fK, chiK) of a NS = (f_from_chi(M, chi, R=R, NS=True), f_kepler(M, NS=True, R=R),
    chi_kepler(M, NS=True, R=R)), for Python floats M, chi, R > 0."""
    C = M / R
    I_SI = 0.237 * M * R**2 * (1.0 + 0.674*C + 4.48*C**4) * Msun * km**2   # moment_of_inertia * Msun * km**2
    M_SI = M * Msun
    f = chi * G * M_SI**2 / (2.0 * math.pi * c * I_SI)
    fK = 1045.0 * math.sqrt(M) * (10.0 / R)**1.5
    chiK = c * I_SI * 2.0 * math.pi * fK / (G * M_SI**2)
    return f, fK, chiK

def _ns_E_J_isco(M, R, f, prograde=True):
    """_E_J_isco_NS(M, R, f, prograde) for Python floats M, R > 0 and f."""
    sign = +1.0 if prograde else -1.0
    Mg = M * Msun_to_km
    # R_ISCO(M, NS=True, f=f):
    if f <= 0.0:
        r = 8.809 * M
    else:
        x = M * f
        r = (8.809*x + -9.166e-4*x**2 + 8.787e-8*x**3 + -6.019e-12*x**4) / f
    r_S = 2.0 * Mg
    if r <= r_S:
        r = r_S * 1.001
    C = M / R
    I_SI = 0.237 * M * R**2 * (1.0 + 0.674*C + 4.48*C**4) * Msun * km**2
    J_SI = I_SI * 2.0 * math.pi * f
    r_m = r * km
    Nphi = 2.0 * G * J_SI / (r_m**3 * c**3)
    Omega_S = math.sqrt(Mg * km) / r_m**1.5
    N = math.sqrt(1.0 - r_S / r)
    v = min(max((Omega_S - sign * Nphi) * r_m / N, 0.0), 0.9999)
    gamma = 1.0 / math.sqrt(1.0 - v**2)
    E = (N + Nphi * v * r_m) * gamma
    J = sign * v * r_m * gamma / km
    return E, J

def E_isco(M, NS=True, R=None, f=None, chi=None, prograde=True):
    """
    Specific energy at the ISCO  (dimensionless, E~ = E/mc^2).

    NS=True  : Bejger, Zdunik & Haensel (2010).  R and f required.
    NS=False : Bardeen, Press & Teukolsky (1972). chi required.
               Uses numerically stable form: E~ = (1-2/rt+sign*chi/rt^1.5)
                                                  / sqrt(1-3/rt+2*sign*chi/rt^1.5)
               where rt = r_isco/Mg.

    Parameters
    ----------
    M        : float - mass [M_sun]
    NS       : bool  - True for NS, False for BH
    R        : float - radius [km], required if NS=True
    f        : float - spin frequency [Hz], required if NS=True
    chi      : float - dimensionless spin, required if NS=False
    prograde : bool  - prograde orbit (default True)

    Returns
    -------
    E : float - dimensionless specific energy
    """
    if NS:
        if R is None or f is None:
            raise ValueError("R and f must be provided for NS (NS=True).")
        E, _ = _E_J_isco_NS(M, R, f, prograde=prograde)
        return E
    else:
        if chi is None:
            raise ValueError("chi must be provided for BH (NS=False).")
        chi  = min(max(chi, 0.0), 1.0 - 1e-10)
        Mg   = geom_mass(M)                             # km
        r    = R_ISCO(M, NS=False, chi=chi,
                      prograde=prograde)                 # km
        rt   = r / Mg                                   # dimensionless r/Mg
        sign = +1.0 if prograde else -1.0
        # BPT (1972) simplified — stable near chi -> 1
        return ((1.0 - 2.0/rt + sign*chi/rt**1.5)
                / np.sqrt(1.0 - 3.0/rt + 2.0*sign*chi/rt**1.5))


def J_isco(M, NS=True, R=None, f=None, chi=None, prograde=True):
    """
    Specific angular momentum at the ISCO [km]  (J~ = L/(mc) geometric).

    Positive for prograde, negative for retrograde.

    NS=True  : Bejger, Zdunik & Haensel (2010).  R and f required.
    NS=False : Bardeen, Press & Teukolsky (1972). chi required.

    Parameters
    ----------
    M        : float - mass [M_sun]
    NS       : bool  - True for NS, False for BH
    R        : float - radius [km], required if NS=True
    f        : float - spin frequency [Hz], required if NS=True
    chi      : float - dimensionless spin, required if NS=False
    prograde : bool  - prograde orbit (default True)

    Returns
    -------
    J : float - specific angular momentum [km]
    """
    if NS:
        if R is None or f is None:
            raise ValueError("R and f must be provided for NS (NS=True).")
        _, J = _E_J_isco_NS(M, R, f, prograde=prograde)
        return J
    else:
        if chi is None:
            raise ValueError("chi must be provided for BH (NS=False).")
        chi  = min(max(chi, 0.0), 1.0 - 1e-10)
        Mg   = geom_mass(M)                             # km
        r    = R_ISCO(M, NS=False, chi=chi,
                      prograde=prograde)                 # km
        rt   = r / Mg                                   # dimensionless r/Mg
        sign = +1.0 if prograde else -1.0
        # BPT (1972) eq. (2.13) fully in dimensionless units (rt = r/Mg, a/Mg = chi).
        # Using rt throughout avoids the dimensional inconsistency in the factored
        # form that gives wrong results for retrograde orbits at chi > 0.
        #   J~/Mg = sign*(rt^2 - 2*sign*chi*sqrt(rt) + chi^2)
        #           / (rt^0.75 * sqrt(rt^1.5 - 3*sqrt(rt) + 2*sign*chi))
        num   = sign * (rt**2 - 2.0*sign*chi*rt**0.5 + chi**2)
        denom = rt**0.75 * np.sqrt(rt**1.5 - 3.0*rt**0.5 + 2.0*sign*chi)
        return Mg * num / denom

# ---------------------------------------------------------------------------
# Spin evolver
#
# State variable: J [km^2]  (geometric angular momentum = chi * Mg^2)
# J is continuous across NS -> BH transition - no state variable switch.
# chi = J/Mg^2 is always >= 0 (spin magnitude). Orbit direction is tracked
# separately via prograde_now (mutable flag inside evolve()). When retrograde
# accretion drives chi to 0, prograde_now flips to True, matching Bardeen (1970).
#
# Evolution equations
# -------------------
#   NS phase  (M < M_max_rot):
#     dJ/dM_acc = (J_isco / E_isco) * Msun_to_km   [km^2 / M_sun]
#     M_grav   += dM_acc   (gravitational mass ~ rest mass for NS)
#
#   BH phase  (M >= M_max_rot):
#     dJ/dM_grav = (J_isco / E_isco) * Msun_to_km   [km^2 / M_sun]
#     M_grav    += dM   (same as NS - M steps in gravitational mass)
#
#     Bardeen (1970) gives dchi/dM_grav directly; the E_isco factor
#     relating dM_grav to dM_acc is already absorbed into dJ/dM_grav.
#     Applying M += E*dM on top of j/e would double-count 1/E.
#
#   At Kepler/Thorne limit:
#     J is clamped to chi_kepler * Mg^2 at the next mass step.
#
# Units of J update:
#   J_isco [km] * dM [M_sun] * Msun_to_km [km/M_sun] / E_isco [dimensionless]
#   = km^2 
# ---------------------------------------------------------------------------

def evolve_spin_during_accretion(Mi, Mf, NS, f=None, chi=0.0, dM=1e-3, eos='APR', prograde=True):
    """
    Evolve compact object spin via disk accretion.

    The object starts as a NS (NS=True) or BH (NS=False).  If NS=True and the
    mass exceeds M_max_rot the evolution switches to the Bardeen (1970) BH regime.

    Returns a dict with keys NS, M, J, R, f, chi, fK, chiK; only the FINAL
    per-step state is stored (as length-1 arrays), which is all the callers use.

    Implementation notes (performance, physics unchanged):
      - math.sqrt (Python-float, no numpy scalar dispatch) for the sqrt calls,
      - the BH-branch isco/f_from_chi math is inlined into the loop, and
      - the BH ISCO radius is computed once per step (shared by E and J).
    The rarer NS branch uses the helper functions (E_isco/J_isco/f_from_chi/...).
    """
    # --- input validation ---
    if f is not None and chi != 0.0:
        raise ValueError("Provide either f or chi as initial condition, not both.")
    if Mi == Mf:
        raise ValueError("Mi and Mf must differ.")
    if NS and eos is None:
        raise ValueError(
            f"NS=True requires an EoS, one of: {', '.join(sorted(EOS_TABLES))}."
        )

    # --- EOS selection ---
    if NS:
        eos_table = get_eos(eos)
        Radius    = eos_table.radius
        M_TOV     = eos_table.M_TOV
        M_max_rot = eos_table.M_max_rot
        if Mi > M_max_rot:
            raise ValueError(
                f"Mi={Mi:.3f} M_sun exceeds M_max_rot={M_max_rot:.3f} M_sun "
                f"for EOS={eos}. Use NS=False to treat as a BH."
            )
        R_TOV = Radius(M_TOV)
    else:
        Radius    = None
        M_TOV     = None
        M_max_rot = -math.inf
        R_TOV     = None

    # dM sign follows direction (np.sign identical to this for nonzero) ---
    # (as a Python float: with a NumPy scalar dM, M += dM would turn M and the whole
    # loop into slower NumPy-scalar arithmetic, with identical IEEE results)
    step = abs(float(dM))
    dM = step if (Mf - Mi) > 0 else (-step if (Mf - Mi) < 0 else 0.0)

    # --- initialise ---
    M  = float(Mi)
    Mg = M * Msun_to_km

    if NS:
        R = Radius(min(M, M_TOV))
        if f is not None:
            chi = chi_from_frequency(M, R, f)
        chi = float(min(max(chi, 0.0), THORNE_LIMIT))
    else:
        chi = float(min(max(chi, 0.0), THORNE_LIMIT))
        R   = Mg * (1.0 + math.sqrt(1.0 - chi**2))

    J = chi * Mg**2

    decreasing = dM < 0

    # constant factors of the BH-branch expressions below (same operations, hoisted):
    c3 = c**3                                    # in f = c**3 * chi / (...)
    c3_thorne = c**3 * THORNE_LIMIT              # numerator of fK
    four_pi_G = 4.0 * math.pi * G                # leading factor of both denominators
    thorne_term = 1.0 + math.sqrt(1.0 - THORNE_LIMIT**2)

    prograde_now = prograde
    last = None   # only the final per-step state is needed by callers (evo[...][-1])

    if NS:
        # The masses of the steps, M_0 = Mi, M_{i+1} = M_i + dM, do not depend on the
        # spin, so the NS radii R(min(M_i, M_TOV)) are interpolated for all steps in
        # one call (np.interp gives the same value for each point as one by one):
        Ms = []
        Mx = M
        while (Mx > Mf) if decreasing else (Mx < Mf):
            Ms.append(Mx if Mx < M_TOV else M_TOV)    # = min(Mx, M_TOV)
            Mx += dM
        Rs = Radius(np.array(Ms)).tolist() if Ms else []
        istep = 0

    while (M > Mf) if decreasing else (M < Mf):

        is_NS = NS and (M < M_max_rot)

        Mg  = M * Msun_to_km
        chi = float(min(max(J / Mg**2, 0.0), THORNE_LIMIT))

        if is_NS:
            # NS path: f_from_chi, f_kepler and chi_kepler, inlined from _ns_spin_state
            # (same operations; C and I_SI are reused by the ISCO quantities below)
            R    = Rs[istep]
            C    = M / R
            I_SI = 0.237 * M * R**2 * (1.0 + 0.674*C + 4.48*C**4) * Msun * km**2
            M_SI = M * Msun
            f    = chi * G * M_SI**2 / (2.0 * math.pi * c * I_SI)
            fK   = 1045.0 * math.sqrt(M) * (10.0 / R)**1.5
            chiK = c * I_SI * 2.0 * math.pi * fK / (G * M_SI**2)
        else:
            # BH path, inlined with math.sqrt (mirrors f_from_chi / R_ISCO / E_isco / J_isco)
            chi2 = chi**2
            horizon_term = 1.0 + math.sqrt(1.0 - chi2)
            R    = Mg * horizon_term
            M_SI = M * Msun
            f    = c3 * chi / (four_pi_G * M_SI * horizon_term)
            fK   = c3_thorne / (four_pi_G * M_SI * thorne_term)
            chiK = THORNE_LIMIT

        last = (is_NS, M, J, R, f, chi, fK, chiK)   # current pre-advance state

        if f < fK:
            if is_NS:
                # E_isco and J_isco for NS=True (_ns_E_J_isco inlined, same operations;
                # Mg, C and I_SI as computed above)
                sgn = +1.0 if prograde_now else -1.0
                if f <= 0.0:
                    r = 8.809 * M
                else:
                    x = M * f
                    r = (8.809*x + -9.166e-4*x**2 + 8.787e-8*x**3 + -6.019e-12*x**4) / f
                r_S = 2.0 * Mg
                if r <= r_S:
                    r = r_S * 1.001
                J_SI = I_SI * 2.0 * math.pi * f
                r_m = r * km
                Nphi = 2.0 * G * J_SI / (r_m**3 * c**3)
                Omega_S = math.sqrt(Mg * km) / r_m**1.5
                N_ = math.sqrt(1.0 - r_S / r)
                v = min(max((Omega_S - sgn * Nphi) * r_m / N_, 0.0), 0.9999)
                gamma = 1.0 / math.sqrt(1.0 - v**2)
                e = (N_ + Nphi * v * r_m) * gamma
                j = sgn * v * r_m * gamma / km
            else:
                # BH: R_ISCO computed once, shared by E~ and J~ (mirrors E_isco/J_isco, NS=False)
                # (chi**2, rt**1.5 and rt**0.5 are each computed once and reused)
                sign = +1.0 if prograde_now else -1.0
                Z1   = 1.0 + (1.0 - chi2)**(1.0/3.0) * (
                           (1.0 + chi)**(1.0/3.0) + (1.0 - chi)**(1.0/3.0))
                Z2   = math.sqrt(3.0 * chi2 + Z1**2)
                r    = Mg * (3.0 + Z2 - sign * math.sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0*Z2)))
                rt   = r / Mg
                rt15 = rt**1.5
                rt05 = rt**0.5
                e    = ((1.0 - 2.0/rt + sign*chi/rt15)
                        / math.sqrt(1.0 - 3.0/rt + 2.0*sign*chi/rt15))
                num   = sign * (rt**2 - 2.0*sign*chi*rt05 + chi2)
                denom = rt**0.75 * math.sqrt(rt15 - 3.0*rt05 + 2.0*sign*chi)
                j    = Mg * num / denom

            J_new = J + (j / e) * dM * Msun_to_km
            if J_new < 0.0 and not prograde_now:
                J_new        = 0.0
                prograde_now = True
            J  = J_new
            M += dM
        else:
            M  += dM
            Mg  = M * Msun_to_km
            J   = chiK * Mg**2
        if NS:
            istep += 1

    # only the final per-step state is returned (as length-1 arrays so callers'
    # evo[key][-1] indexing still works); intermediate history is not stored.
    is_NS, M, J, R, f, chi, fK, chiK = last
    return {
        'NS'  : np.array([is_NS]),
        'M'   : np.array([M]),
        'J'   : np.array([J]),
        'R'   : np.array([R]),
        'f'   : np.array([f]),
        'chi' : np.array([chi]),
        'fK'  : np.array([fK]),
        'chiK': np.array([chiK]),
    }

# ---------------------------------------------------------------------------
# Fast evaluation of the final spin of a BH (not NS) for small relative growth
# ---------------------------------------------------------------------------
#
# For a BH (NS=False), the loop of evolve_spin_during_accretion integrates, with
# Euler steps of dM in M,
#     dJ/dM = (j~/e~) Msun_to_km,   j~ = Mg J~/Mg  (Bardeen 1970),
# which with J = chi Mg^2 is the autonomous equation
#     dchi/dlnM = (J~/Mg)(chi) / E~(chi) - 2 chi
# (spin-up for prograde, spin-down for retrograde accretion; chi capped at
# THORNE_LIMIT; a retrograde disk that spins the BH down to chi = 0 continues as
# a prograde disk). The value returned there, evo['chi'][-1], is the spin at the
# start of the last Euler step, i.e. at M_{N-1}, the last mass of the sequence
# M_0 = Mi, M_{i+1} = M_i + dM below Mf.
#
# bh_spin_after_accretion integrates the same equation from Mi to M_{N-1} with
# fourth-order Runge-Kutta steps in ln M (at most 0.01 each) instead of the ~100
# Euler steps, ~10x faster. The two differ by the error of the Euler steps, which
# relative to the change of chi is ~ 1.5 (Mf - Mi)/(100 Mi): for growth up to
# FAST_SPIN_MAX_GROWTH = 2% (used by spin_after_accretion) below 3e-4 of the
# change of chi, i.e. < 2e-5 in chi. (A retrograde disk reaching chi = 0 switches
# to prograde at the next mass M_i of the sequence, as in the Euler loop; the
# Euler steps of a BH at the Thorne limit are repeated exactly, since whether it
# leaves the limit under retrograde accretion depends on the rounding of J/Mg^2.)

FAST_SPIN_MAX_GROWTH = 0.02

def _dchi_dlnM(chi, sign):
    """dchi/dlnM for a thin disk at the BH ISCO (sign = +1 prograde, -1 retrograde);
    same expressions as the BH branch of evolve_spin_during_accretion."""
    chi2 = chi * chi
    Z1 = 1.0 + (1.0 - chi2)**(1.0/3.0) * ((1.0 + chi)**(1.0/3.0) + (1.0 - chi)**(1.0/3.0))
    Z2 = math.sqrt(3.0 * chi2 + Z1**2)
    rt = 3.0 + Z2 - sign * math.sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0*Z2))
    rt15 = rt**1.5
    rt05 = rt**0.5
    e = (1.0 - 2.0/rt + sign*chi/rt15) / math.sqrt(1.0 - 3.0/rt + 2.0*sign*chi/rt15)
    j = sign * (rt**2 - 2.0*sign*chi*rt05 + chi2) / (rt**0.75 * math.sqrt(rt15 - 3.0*rt05 + 2.0*sign*chi))
    return j / e - 2.0 * chi

def _rk4_spin(chi, x, x_end, sign):
    """
    Integrate dchi/dlnM from x = ln M to x_end with RK4 steps <= 0.01.
    Returns (chi, x, crossed): for retrograde accretion, crossed is True if chi
    would become negative within the next step, which then starts at x.
    """
    L = x_end - x
    if L <= 0.0:
        return chi, x, False
    nsub = max(1, int(math.ceil(L / 0.01)))
    h = L / nsub
    for _ in range(nsub):
        if chi >= THORNE_LIMIT and sign > 0:
            return THORNE_LIMIT, x_end, False   # prograde: stays at the limit
        k1 = _dchi_dlnM(chi, sign)
        c = chi + 0.5*h*k1
        if c < 0.0:
            return chi, x, True
        k2 = _dchi_dlnM(min(c, THORNE_LIMIT), sign)
        c = chi + 0.5*h*k2
        if c < 0.0:
            return chi, x, True
        k3 = _dchi_dlnM(min(c, THORNE_LIMIT), sign)
        c = chi + h*k3
        if c < 0.0:
            return chi, x, True
        k4 = _dchi_dlnM(min(c, THORNE_LIMIT), sign)
        c = chi + h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
        if c < 0.0:
            return chi, x, True
        chi = min(c, THORNE_LIMIT)
        x += h
    return chi, x_end, False

def bh_spin_after_accretion(Mi, Mf, chi=0.0, dM=1e-3, prograde=True):
    """
    Final spin of a BH accreting from Mi to Mf > Mi, as
    evolve_spin_during_accretion(Mi, Mf, NS=False, chi=chi, dM=dM, prograde=prograde)['chi'][-1],
    computed with RK4 steps instead of Euler steps (see the notes above).
    """
    step = abs(float(dM))
    Mi = float(Mi)
    if not Mf > Mi:
        raise ValueError('bh_spin_after_accretion requires Mf > Mi')

    # last mass M_{N-1} < Mf of the sequence M_{i+1} = M_i + dM (as in the Euler loop):
    M = Mi
    M_end = Mi
    while M < Mf:
        M_end = M
        M += step

    # A BH at the Thorne limit stays there in the Euler loop as long as J/Mg^2 does
    # not round to below THORNE_LIMIT (then f = fK and J is reset to the limit);
    # these steps are repeated exactly, also for retrograde accretion:
    chi = float(min(max(chi, 0.0), THORNE_LIMIT))
    M = Mi
    J = chi * (M * Msun_to_km)**2
    while True:
        chi = float(min(max(J / (M * Msun_to_km)**2, 0.0), THORNE_LIMIT))
        if chi < THORNE_LIMIT:
            break
        if not M + step < Mf:
            return chi                         # (last step of the Euler loop)
        M += step
        J = THORNE_LIMIT * (M * Msun_to_km)**2

    x_end = math.log(M_end)
    if prograde:
        return _rk4_spin(chi, math.log(M), x_end, 1.0)[0]

    chi, x, crossed = _rk4_spin(chi, math.log(M), x_end, -1.0)
    if not crossed:
        return chi
    # chi reaches 0 within the next RK step: approach the crossing with RK4 steps of
    # 0.9 times the linear estimate of the distance (each reduces chi ~10-fold), then
    # extrapolate linearly:
    while True:
        k = _dchi_dlnM(chi, -1.0)
        d = chi / (-k)
        if chi < 1e-8:
            x_c = x + d
            break
        x_next = min(x + 0.9 * d, x_end)
        chi, x, crossed = _rk4_spin(chi, x, x_next, -1.0)
        if x >= x_end and not crossed:
            return chi
    # the Euler loop sets J = 0 and turns prograde at the first mass of the
    # sequence above the crossing:
    M = Mi
    while math.log(M) <= x_c:
        M += step
    if M > M_end:
        return 0.0
    return _rk4_spin(0.0, math.log(M), x_end, 1.0)[0]

def spin_after_accretion(Mi, Mf, NS, chi=0.0, dM=1e-3, eos='APR', prograde=True):
    """
    Final spin evo['chi'][-1] of evolve_spin_during_accretion(Mi, Mf, NS, chi=chi,
    dM=dM, eos=eos, prograde=prograde), with bh_spin_after_accretion for BHs that
    grow by at most FAST_SPIN_MAX_GROWTH (relative).
    """
    if not NS and Mf > Mi and Mf - Mi <= FAST_SPIN_MAX_GROWTH * Mi:
        return bh_spin_after_accretion(Mi, Mf, chi=chi, dM=dM, prograde=prograde)
    return evolve_spin_during_accretion(Mi=Mi, Mf=Mf, NS=NS, f=None, chi=chi, dM=dM, eos=eos, prograde=prograde)['chi'][-1]

def radiative_efficiency(m, chi, is_NS, eos):
    """eta = 1 - E_isco using the object's spin magnitude (prograde disk). Spin not evolved."""
    if is_NS:
        R = np.interp(m, eos.M, eos.R)                    # NS radius [km]
        f = f_from_chi(m, chi, R=R, NS=True)              # spin freq [Hz] from chi
        return 1.0 - E_isco(m, NS=True, R=R, f=f)
    return 1.0 - E_isco(m, NS=False, chi=chi)


def gas_accretion_rate(m, rho_gas, v_eff2, eta, f_Edd, lam=1.0):
    """Three-way-min gas accretion rate [Msun/Myr]: Bondi, dynamical, Eddington.
    Rapster units (Msun, pc, Myr, km/s; G_Newton=1/232). Vectorized over m/eta."""
    v_eff3    = v_eff2**1.5
    bondi     = 4*np.pi*lam*G_Newton**2 * m**2 * rho_gas / v_eff3
    dynamical = v_eff3 / G_Newton
    eddington = f_Edd * m / (eta * t_Edd)
    return np.minimum(np.minimum(bondi, dynamical), eddington)


def accrete_gas(state, config):
    """Grow compact-object masses (singles/binaries/pairs/triples; BHs+NSs) by gas
    accretion using each object's spin for the efficiency, and drain the reservoir.
    Mass-only (spin not evolved)."""
    M_gas = state['M_gas']; dt = state['dt']
    if M_gas <= 0 or dt <= 0:
        return
    rh = state['rh']; vBH = state['vBH']
    c_s = config['c_s']; f_Edd = config['f_Edd']
    eos = get_eos(config['EoS']); M_TOV = eos.M_TOV

    rho_gas = 3*M_gas / (8*np.pi*rh**3)                   # mean within r_h [Msun/pc^3]
    v_eff2  = c_s**2 + vBH**2

    binaries = state['binaries']; pairs = state['pairs']; triples = state['triples']
    BH = state['BH']
    mass_cols = [BH.m, binaries[:, 4], binaries[:, 5], pairs[:, 1],
                 triples[:, 4], triples[:, 5], triples[:, 6]]
    spin_cols = [BH.s, binaries[:, 6], binaries[:, 7], pairs[:, 2],
                 triples[:, 7], triples[:, 8], triples[:, 9]]

    dms = []
    for m, sp in zip(mass_cols, spin_cols):
        dm = np.zeros_like(m, dtype=float)
        valid = m > 0
        if valid.any():
            mv = m[valid]; sv = np.abs(sp[valid]); is_NS = mv < M_TOV
            eta = np.array([radiative_efficiency(mv[i], sv[i], is_NS[i], eos)
                            for i in range(mv.size)])
            dm[valid] = gas_accretion_rate(mv, rho_gas, v_eff2, eta, f_Edd) * dt
        dms.append(dm)

    total = float(sum(d.sum() for d in dms))
    if total > M_gas and total > 0:                       # conserve: can't accrete more than exists
        dms = [d*(M_gas/total) for d in dms]; total = M_gas

    BH.set_all_masses(mass_cols[0] + dms[0])
    binaries[:, 4] += dms[1]; binaries[:, 5] += dms[2]
    pairs[:, 1]    += dms[3]
    triples[:, 4]  += dms[4]; triples[:, 5] += dms[5]; triples[:, 6] += dms[6]
    state['binaries'] = binaries; state['pairs'] = pairs; state['triples'] = triples
    state['M_gas'] = M_gas - total

# End of file.
