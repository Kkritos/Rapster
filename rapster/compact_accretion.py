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
# Get the directory where constants.py lives (the 'rapster' folder):
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up one level to the root, then into the Data folder:
DATA_PATH_APR = os.path.join(BASE_DIR, '..', 'Data', 'APR.txt')
DATA_PATH_AU = os.path.join(BASE_DIR, '..', 'Data', 'AU.txt')

R_APR, M_APR = np.loadtxt(DATA_PATH_APR, unpack=True)
R_AU,  M_AU  = np.loadtxt(DATA_PATH_AU,  unpack=True)

def Radius_APR(M):
    return np.interp(M, M_APR, R_APR)

def Radius_AU(M):
    return np.interp(M, M_AU, R_AU)

# minimum and maximum allowed NS mass for different EoS:
M_APR_min, M_APR_max = M_APR.min(), M_APR.max()
M_AU_min,  M_AU_max  = M_AU.min(),  M_AU.max()

# ---------------------------------------------------------------------------
# Supramassive limit  -  Breu & Rezzolla (2016)
#   M_max,rot ~ 1.203 * M_TOV  (EOS-insensitive to ~2%)
# ---------------------------------------------------------------------------
BREU_REZZOLLA_FACTOR = 1.203

M_APR_TOV     = M_APR_max
M_APR_max_rot = BREU_REZZOLLA_FACTOR * M_APR_TOV

M_AU_TOV      = M_AU_max
M_AU_max_rot  = BREU_REZZOLLA_FACTOR * M_AU_TOV

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
        chi = np.clip(chi, 0.0, THORNE_LIMIT)
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
        chi  = np.clip(chi, 0.0, THORNE_LIMIT)
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
    v     = np.clip((Omega_S - sign * Nphi) * r_m / N, 0.0, 0.9999)
    gamma = 1.0 / np.sqrt(1.0 - v**2)

    E = (N + Nphi * v * r_m) * gamma                   # dimensionless (always > 0)
    J = sign * v * r_m * gamma / km                    # km  (positive prograde,
                                                        #      negative retrograde)
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
        chi  = np.clip(chi, 0.0, 1.0 - 1e-10)
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
        chi  = np.clip(chi, 0.0, 1.0 - 1e-10)
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

def evolve(Mi, Mf, NS, f=None, chi=0.0, dM=1e-3, eos='APR', prograde=True):
    """
    Evolve compact object spin via disk accretion.

    The object starts as a NS (NS=True) or BH (NS=False).  If NS=True and
    the mass exceeds M_max_rot the evolution automatically switches to the
    Bardeen (1970) BH regime.

    Parameters
    ----------
    Mi  : float - initial mass [M_sun]
    Mf  : float - final mass [M_sun]
    NS  : bool  - True if the initial object is a neutron star
    f   : float - initial spin frequency [Hz]  (provide f or chi, not both)
    chi : float - initial dimensionless spin    (default 0.0)
    dM       : float - mass step [M_sun] (sign set automatically from Mf-Mi)
    eos      : str   - 'APR' or 'AU', required if NS=True
    prograde : bool  - True for prograde accretion (spin-up), False for retrograde
                       (spin-down). When retrograde drives chi to 0 the orbit
                       automatically flips to prograde, matching Bardeen (1970).

    Returns
    -------
    dict with keys: NS, M, J, R, f, chi, fK, chiK
    """
    # --- input validation ---
    if f is not None and chi != 0.0:
        raise ValueError("Provide either f or chi as initial condition, not both.")
    if Mi == Mf:
        raise ValueError("Mi and Mf must differ.")
    if NS and eos is None:
        raise ValueError("NS=True requires an EOS: 'APR' or 'AU'.")

    # --- EOS selection ---
    if NS:
        if eos == 'APR':
            Radius    = Radius_APR
            M_TOV     = M_APR_TOV
            M_max_rot = M_APR_max_rot
        elif eos == 'AU':
            Radius    = Radius_AU
            M_TOV     = M_AU_TOV
            M_max_rot = M_AU_max_rot
        else:
            raise ValueError(f"Unknown EOS '{eos}'. Choose 'APR' or 'AU'.")
        if Mi > M_max_rot:
            raise ValueError(
                f"Mi={Mi:.3f} M_sun exceeds M_max_rot={M_max_rot:.3f} M_sun "
                f"for EOS={eos}. Use NS=False to treat as a BH."
            )
        R_TOV = Radius(M_TOV)   # freeze radius above M_TOV
    else:
        Radius    = None
        M_TOV     = None
        M_max_rot = -np.inf     # always in BH phase
        R_TOV     = None

    # --- dM sign follows direction of evolution ---
    dM = abs(dM) * np.sign(Mf - Mi)

    # --- initialise ---
    M  = float(Mi)
    Mg = geom_mass(M)

    if NS:
        R = Radius(min(M, M_TOV))
        if f is not None:
            chi = chi_from_frequency(M, R, f)
        chi = float(np.clip(chi, 0.0, THORNE_LIMIT))
    else:
        chi = float(np.clip(chi, 0.0, THORNE_LIMIT))
        R   = Mg * (1.0 + np.sqrt(1.0 - chi**2))   # outer horizon radius

    J = chi * Mg**2   # geometric angular momentum [km^2]

    # --- output storage ---
    NS_arr, M_arr, J_arr, R_arr      = [], [], [], []
    f_arr, chi_arr, fK_arr, chiK_arr = [], [], [], []

    def not_done(m):
        return m > Mf if dM < 0 else m < Mf

    # --- orbit direction (mutable: flips to prograde once chi reaches 0) ---
    # chi is always >= 0 (magnitude); prograde_now tracks direction.
    # When retrograde accretion drives chi to 0 the disk flips prograde,
    # matching Bardeen (1970) and the standalone evolve_spin_RK4.
    prograde_now = prograde

    # --- main loop ---
    while not_done(M):

        is_NS = NS and (M < M_max_rot)

        Mg  = geom_mass(M)
        chi = float(np.clip(J / Mg**2, 0.0, THORNE_LIMIT))

        if is_NS:
            # freeze radius in supramassive regime (M_TOV < M < M_max_rot)
            R    = Radius(min(M, M_TOV))
            f    = f_from_chi(M, chi, R=R, NS=True)
            fK   = f_kepler(M, NS=True, R=R)
            chiK = chi_kepler(M, NS=True, R=R)
        else:
            # BH: R = outer Kerr horizon
            R    = Mg * (1.0 + np.sqrt(1.0 - chi**2))
            f    = f_from_chi(M, chi, NS=False)
            fK   = f_kepler(M, NS=False)
            chiK = chi_kepler(M, NS=False)

        # --- store current state ---
        NS_arr.append(is_NS)
        M_arr.append(M)
        J_arr.append(J)
        R_arr.append(R)
        f_arr.append(f)
        chi_arr.append(chi)
        fK_arr.append(fK)
        chiK_arr.append(chiK)

        # --- advance ---
        if f < fK:
            if is_NS:
                e = E_isco(M, NS=True,  R=R, f=f, prograde=prograde_now)
                j = J_isco(M, NS=True,  R=R, f=f, prograde=prograde_now)
            else:
                e = E_isco(M, NS=False, chi=chi, prograde=prograde_now)
                j = J_isco(M, NS=False, chi=chi, prograde=prograde_now)

            # j < 0 for retrograde: J decreases. Clamp at 0 and flip to prograde.
            J_new = J + (j / e) * dM * Msun_to_km
            if J_new < 0.0 and not prograde_now:
                J_new        = 0.0
                prograde_now = True
            J  = J_new
            M += dM
        else:
            # at spin limit: clamp J to Kepler/Thorne value, advance M
            M  += dM
            Mg  = geom_mass(M)
            J   = chiK * Mg**2

    return {
        'NS'  : np.array(NS_arr),
        'M'   : np.array(M_arr),
        'J'   : np.array(J_arr),
        'R'   : np.array(R_arr),
        'f'   : np.array(f_arr),
        'chi' : np.array(chi_arr),
        'fK'  : np.array(fK_arr),
        'chiK': np.array(chiK_arr),
    }

# End of file.
