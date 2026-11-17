"""
Central massive BH (-eta; version 2.18.0).

Rapster describes the BH subsystem by a mean field (Breen & Heggie 2013): average BH mass,
half-mass and core radii, velocity dispersion, relaxation time. A BH much heavier than
the others (a -MBH seed, or a BH grown by mergers and accretion) breaks it: e.g. the BH
relaxation factor psi_BH = sum(m^5/2) / (N_BH <m>^5/2) is dominated by that one BH, and
the BH subsystem then relaxes 10^3-10^4 times too fast.

Phase 1 (promotion). At every step the heaviest BH (single, in a BBH, a BH-star pair or a
triple), of mass M, is promoted to a central massive BH once its influence radius
exceeds eta times its wandering radius in the core of the other BHs (promotion_ratio):
  r_infl = G M / sigma^2,
  r_wand^2 = <v^2> / omega^2 with <v^2> = 3 (m / M) sigma^2 (energy equipartition) and
  omega^2 = 4 pi G rho_c / 3 (the harmonic core) (Bahcall & Wolf 1976; Chatterjee,
  Hernquist & Loeb 2002),
with the mean field of the other BHs (sigma = v_BH / sqrt(3), rho_c = n_cBH <m_BH> +
n_star <m>, m = <m_BH>): r_infl >= eta r_wand <=> M^3 >= eta^2 9 m sigma^6 / (4 pi G^3
rho_c), and M >= 10 <m_BH> (MIN_MASS_RATIO: the estimate assumes M >> m). A promoted BH stays promoted while the heaviest BH is at least as heavy (until it
is ejected); a -MBH seed is promoted from the start. The mean field is then that of the
other BHs; the massive BH adds to the mass that sets the stellar velocity dispersion and
the relaxation time (not to M_cl), and is left out of the -fCO stopping condition. It
keeps forming binaries, merging and accreting.

Phase 2 (physics of the massive BH):
- escape velocity of the BHs and BBHs from the radius max(r_cBH, r_infl) in its
  potential, v_esc^2 = F_ESC^2 (v_star^2 + v_BH^2) + 2 G M / max(r_cBH, r_infl) (a merger
  remnant of the massive BH itself: without its own potential), functions.v_esc_cl;
- three-body binaries of two other BHs form only outside its influence radius: those of
  the mean-field rate are kept with probability f_3bb = 1 - min(1, 2 M / (N_cBH <m_BH>))
  (the fraction of the BH core outside r_infl, which holds 2 M); those with the massive BH
  are kept, as its binaries are the heat source of the core;
- BHs and stars reach it through its loss cone (loss_cone_rate), instead of the
  mean-field two-body captures and TDEs, from which it is left out: GW captures of BHs
  (IMRIs, merger channel 7) and TDEs of stars (tde_type 31) and white dwarfs (32).
"""

import numpy as np

from .constants import G_Newton, c_light, R_sun, white_dwarf_mass, tBH_form
from .functions import (T_GW, sample_angles, get_maxwell_sample, redshift_interp, lookback_interp,
                        v_esc_remnant, v_esc_cl, set_massive_bh)
from .remnant import merger_remnant
from .fastrandom import poisson_draw
from .bhpool import FAM
from .stellar_evolution import get_star, R_WhiteDwarf
from .tidal_disruptions import disrupt_star

# merger channel of the GW captures by the massive BH, TDE types of stars and WDs by it:
CHANNEL_IMRI = 7
TDE_TYPE_STAR = 31
TDE_TYPE_WD = 32

# a BH is promoted only if it is at least this many times heavier than the average of the
# other BHs: the wandering-radius estimate assumes M >> m, and in small BH cores (a few BHs)
# r_infl >= eta r_wand already holds for second- and third-generation remnants:
MIN_MASS_RATIO = 10.0

# keys appended to the evolution file with -eta > 0:
EVOLUTION_KEYS_MBH = ['M_MBH', 'chi_MBH', 'g_MBH', 'r_infl', 'f_3bb', 'N_IMRI', 'N_tdeMBH',
                      'Gamma_IMRI', 'Gamma_tdeMBH']

# prefactor of the focused TDE rate of Rapster's micro-TDEs (cluster_evolution.evolve_tdes):
_K_TDE = 2 * np.sqrt(2 * (3 * np.pi - 8))


def promotion_ratio(M, m, sigma2, rho_c):
    """r_infl / r_wand of a BH of mass M [Msun] at the center of a core of density rho_c
    [Msun/pc^3] and 1D velocity dispersion sigma2 [(km/s)^2], made of objects of mass m."""
    r_infl = G_Newton * M / sigma2
    r_wand = np.sqrt(9 * m * sigma2 / (4 * np.pi * G_Newton * rho_c * M))
    return r_infl / r_wand


def promote(state, config, i_aux1, N_BH, sum_m, sum_m25, mBH_max, n_star, m_avg, logLcl, mean_field):
    """(M, mean field of the other BHs) of the central massive BH, or (0, None) if there
    is none; mean_field(i_aux1, N, sum_m, sum_m25, Mcl, m_avg, rh, logLcl, c_xi) is the mean
    field of N BHs (cluster_evolution._bh_mean_field). Before the BH subsystem is active
    (i_aux1 != 1), a promoted BH (a -MBH seed) stays promoted, without a mean field."""
    eta = config.get('eta_MBH', 0.0)
    M_prev = state.get('M_MBH', 0.0)
    if eta <= 0:
        return 0.0, None
    if i_aux1 != 1:
        return M_prev, None
    if N_BH < 3:
        return 0.0, None
    triples = state['triples']
    M_tri = float(np.max(triples[:, 4:7])) if triples.size else 0.0
    M = float(max(mBH_max, M_tri))
    if M_tri > mBH_max:
        # (in a triple: not in the mean-field sums, only in N_BH)
        mf = mean_field(1, N_BH - 1, sum_m, sum_m25, state['Mcl'], m_avg, state['rh'], logLcl, config['c_xi'])
    else:
        mf = mean_field(1, N_BH - 1, sum_m - M, sum_m25 - M**(5/2), state['Mcl'], m_avg, state['rh'], logLcl, config['c_xi'])
    if M_prev > 0 and M >= M_prev:
        return M, mf
    sigma2 = mf['vBH']**2 / 3
    rho_c = mf['nc_BH'] * mf['mBH_avg'] + n_star * m_avg
    if M >= MIN_MASS_RATIO * mf['mBH_avg'] and sigma2 > 0 and rho_c > 0 \
            and promotion_ratio(M, mf['mBH_avg'], sigma2, rho_c) >= eta:
        return M, mf
    return 0.0, None


def phase2(config):
    """True if the physics of the central massive BH (Phase 2, -MBHph 1) is on."""
    return config.get('eta_MBH', 0.0) > 0 and config.get('MBH_physics', 0) == 1


def set_phase2_properties(state, config, M, mf):
    """Influence radius of the central massive BH of mass M (0: none) with the mean field
    mf of the other BHs, and with -MBHph 1 its escape-velocity term and f_3bb; stored in
    state."""
    if M > 0 and mf is not None and mf['vBH'] > 0:
        sigma2 = mf['vBH']**2 / 3
        r_infl = G_Newton * M / sigma2
        if phase2(config):
            set_massive_bh(M, 2 * G_Newton * M / max(mf['rc_BH'], r_infl))
            M_core = mf['Nc_BH'] * mf['mBH_avg']
            f_3bb = max(0.0, 1.0 - 2 * M / M_core) if M_core > 0 else 0.0
        else:
            set_massive_bh(0.0, 0.0)
            f_3bb = 1.0
    else:
        r_infl = 0.0
        set_massive_bh(0.0, 0.0)
        f_3bb = 1.0
    state['r_infl'] = r_infl
    state['f_3bb_MBH'] = f_3bb


def slot_in_pool(state):
    """Slot of the central massive BH in the pool of single BHs, or -1 (none, or in a
    binary, pair or triple)."""
    M = state.get('M_MBH', 0.0)
    BH = state['BH']
    if M > 0 and BH.size > 0 and BH.max_mass() >= M:
        return BH.D
    return -1


def take_out(state, config):
    """Take the central massive BH out of the pool (if it is there), e.g. for the
    mean-field GW captures and TDEs, from which it is left out; returns its (m, s, g, h),
    or None (also without -MBHph 1). put_back returns it to the pool."""
    if not phase2(config):
        return None
    k = slot_in_pool(state)
    if k < 0:
        return None
    BH = state['BH']
    rec = BH.get(k)
    BH.remove(k)
    return rec


def put_back(state, rec):
    if rec is not None:
        state['BH'].add(*rec)


def heaviest_object(state):
    """(m, s, g) of the heaviest compact object (single, in a BBH, a BH-star pair or a
    triple); (0, 0, 0) if there is none."""
    best = (0.0, 0.0, 0.0)
    BH = state['BH']
    if BH.size > 0:
        m, s, g, _ = BH.get(BH.D)
        best = (m, s, g)
    b = state['binaries']; p = state['pairs']; tr = state['triples']
    for arr, cm, cs, cg in ((b, 4, 6, 8), (b, 5, 7, 9), (p, 1, 2, 3)):
        if arr.shape[0] > 1:
            i = int(np.argmax(arr[:, cm]))
            if arr[i, cm] > best[0]:
                best = (float(arr[i, cm]), float(arr[i, cs]), float(arr[i, cg]))
    if tr.shape[0] > 1:
        # triples: [ind, ..., m0, m1, m2 (columns 4-6), s0, s1, s2 (7-9), g0, g1, g2 (10-12), ...]
        i, j = np.unravel_index(int(np.argmax(tr[:, 4:7])), (tr.shape[0], 3))
        if tr[i, 4 + j] > best[0]:
            best = (float(tr[i, 4 + j]), float(tr[i, 7 + j]), float(tr[i, 10 + j]))
    return best


# ----------------------------------------------------------------- loss-cone rates
def loss_cone_rate(rate_full, N_in, t_r, L):
    """Rate [1/Myr] at which objects reach the loss cone of the massive BH, interpolating
    between the full loss cone (rate_full: the gravitationally focused rate of encounters
    with pericenter within the loss-cone radius) and the empty one (diffusion by two-body
    relaxation, N_in / (t_r L), with N_in the objects within the influence radius, t_r the
    relaxation time there and L = ln(J_c^2 / J_lc^2)): 1 / (1 / rate_full + t_r L / N_in)
    (e.g. Merritt 2013, Sec. 6.1). With less than one object bound to the massive BH
    (N_in < 1), they arrive on unbound orbits, for which the loss cone is full: rate_full."""
    if rate_full <= 0 or t_r <= 0:
        return 0.0
    if N_in < 1:
        return rate_full
    return 1.0 / (1.0 / rate_full + t_r * L / N_in)


def gw_capture_pericenter(M, m, v):
    """Largest pericenter [pc] at which an encounter of m with M at relative velocity v
    [km/s] ends bound by GW emission (the energy radiated at pericenter, 85 pi / (12
    sqrt(2)) G^(7/2) mu^2 M_t^(5/2) / (c^5 r_p^(7/2)), equals mu v^2 / 2)."""
    Mt = M + m; mu = M * m / Mt
    return (85 * np.pi * mu * Mt**(5/2) * G_Newton**(7/2) / (6 * np.sqrt(2) * c_light**5 * v**2))**(2/7)


def imri_rate(M, m, n_BH, vBH, N_BH_other):
    """Rate [1/Myr] of GW captures of BHs (mass m, number density n_BH in the core, 3D
    velocity dispersion vBH) by the massive BH M (at rest): full loss cone, the capture
    cross section of Rapster's two-body captures (functions.Rate_cap) for masses M and m
    at relative velocity vBH; empty loss cone, the BHs within r_infl = G M / sigma^2 (mass
    2 M, singular isothermal sphere) diffusing in on the relaxation time there, t_r = 0.34
    sigma^3 / (G^2 m rho lnL) with rho = sigma^2 / (2 pi G r_infl^2)."""
    if M <= 0 or m <= 0 or n_BH <= 0 or vBH <= 0 or N_BH_other <= 0:
        return 0.0
    v = vBH
    Sigma_cap = 2 * np.pi * (85 * np.pi / 6 / np.sqrt(2))**(2/7) * G_Newton**2 * (M + m)**(10/7) * (M * m)**(2/7) \
        / c_light**(10/7) / v**(18/7)
    rate_full = n_BH * Sigma_cap * v
    sigma = vBH / np.sqrt(3)
    r_infl = G_Newton * M / sigma**2
    N_in = min(2 * M / m, N_BH_other)
    lnL = np.log(np.e + M / m)
    t_r = 0.68 * np.pi * G_Newton * M**2 / (m * sigma**3 * lnL)
    L = np.log(max(np.e, r_infl / (2 * gw_capture_pericenter(M, m, v))))
    return loss_cone_rate(rate_full, N_in, t_r, L)


def tde_rate(M, m_obj, R_obj, n_obj, v_obj, m_avg, n_star):
    """Rate [1/Myr] of TDEs of objects (mass m_obj, radius R_obj [pc], number density
    n_obj, 3D velocity dispersion v_obj) by the massive BH M (at rest): full loss cone,
    Rapster's focused micro-TDE rate for one BH of mass M; empty loss cone, the objects
    within r_infl = G M / sigma^2 (a Bahcall-Wolf cusp n ~ r^(-3/2) of the light objects
    around the heavy BH, (8 pi / 3) n r_infl^3) diffusing in on the stellar relaxation time
    there, t_r = 0.34 sigma^3 / (G^2 n_star <m>^2 lnL)."""
    if M <= 0 or n_obj <= 0 or v_obj <= 0 or n_star <= 0:
        return 0.0
    rate_full = _K_TDE * G_Newton * (M + m_obj) * n_obj * R_obj * (M / m_obj)**(1/3) / v_obj
    sigma = v_obj / np.sqrt(3)
    r_infl = G_Newton * M / sigma**2
    N_in = 8 * np.pi / 3 * n_obj * r_infl**3
    lnL = np.log(np.e + M / m_avg)
    t_r = 0.34 * sigma**3 / (G_Newton**2 * n_star * m_avg**2 * lnL)
    r_t = R_obj * (M / m_obj)**(1/3)
    L = np.log(max(np.e, r_infl / (2 * r_t)))
    return loss_cone_rate(rate_full, N_in, t_r, L)


# ----------------------------------------------------------------- events
def _imri(state, config, k_M):
    """GW capture of a BH by the massive BH (pool slot k_M): the BH drawn with weights
    m^(2/7) s(m) (as the captured BHs of Rapster's two-body captures), the orbit as in
    two_body_capture (relative velocity Maxwellian with 1D dispersion v_BH / sqrt(3), the
    massive BH at rest). The binary merges in this step if its GW time is shorter (remnant
    and kick from merger_remnant, at the actual mass ratio), else it joins the BBHs.
    Returns False if the massive BH left the pool (merged into a binary or ejected)."""
    BH = state['BH']
    t = state['t']; dt = state['dt']; z = state['z']
    vBH = state['vBH']; v_star = state['v_star']
    j = BH.drawer_without_D(FAM['c27'])()
    M, sM, gM, hM = BH.get(k_M)
    m, sm, gm, hm = BH.get(j)
    m12 = M + m; mu = M * m / m12
    v_rel = get_maxwell_sample(vBH / np.sqrt(3))
    b_max = (340 * np.pi / 3)**(1/7) * m12**(6/7) * mu**(1/7) / v_rel**(9/7) * G_Newton * c_light**(-5/7)
    for _ in range(1000):
        b = np.sqrt(np.random.rand() * b_max**2)
        rp = b**2 * v_rel**2 / 2 / G_Newton / m12
        E_gw = 85 * np.pi / 12 / np.sqrt(2) * mu**2 * m12**(5/2) / rp**(7/2) * G_Newton**(7/2) / c_light**5
        E_fin = mu * v_rel**2 / 2 - E_gw
        if 1 + 2 * E_fin * b**2 * v_rel**2 / m12**2 / mu / G_Newton**2 >= 0:
            break
    else:
        return True                     # (no bound orbit found: no capture)
    sma = - G_Newton * m12 * mu / 2 / E_fin
    eccen = np.sqrt(1 + 2 * E_fin * b**2 * v_rel**2 / m12**2 / mu / G_Newton**2)
    ind = np.random.randint(0, 999999999)
    state['N_IMRI'] = state.get('N_IMRI', 0) + 1
    T_gw = T_GW(M, m, sma, eccen)
    t_lb_form = lookback_interp(state['zCl_form'])
    t_left = t_lb_form - t
    if T_gw < min(dt, t_left):
        theta1, theta2, dPhi = sample_angles()
        m_rem, s_rem, vGW_kick = merger_remnant(M, m, sM, sm, theta1, theta2, dPhi)
        g_rem = max(gM, gm) + 1
        h_rem = hM + hm
        v_esc = v_esc_cl(v_star, vBH)
        retained = vGW_kick < v_esc_remnant(v_esc, v_star, vBH, M, m)
        BH.remove_two(k_M, j)
        if retained:
            BH.add(m_rem, s_rem, g_rem, h_rem)
            state['N_BH'] -= 1
            state['N_meRe'] += 1
            state['M_MBH'] = float(m_rem)
        else:
            state['N_BH'] -= 2
            state['N_meEj'] += 1
        state['N_me'] += 1
        s_eff = (M * sM * np.cos(theta1) + m * sm * np.cos(theta2)) / m12
        t_merge = t + T_gw
        t_lb = t_lb_form - t_merge
        z_merge = redshift_interp(t_lb) if t_lb >= 0.1 else 0.0
        state['mergers'].append([state['seed'], ind, CHANNEL_IMRI, sma, eccen, M, m, sM, sm, gM, gm, theta1, theta2, dPhi,
                                 t, z, t_merge, z_merge, m_rem, s_rem, g_rem, vGW_kick, s_eff, m / M, v_esc, hM, hm])
        return retained
    # a binary of the massive BH (hardened, and merging, with the other BBHs):
    BH.remove_two(k_M, j)
    state['binaries'] = np.append(state['binaries'], [[ind, CHANNEL_IMRI, sma, eccen, M, m, sM, sm, gM, gm, t, z, 0, hM, hm]], axis=0)
    state['N_BBH'] += 1
    return False


def _tde(state, config, k_M, tde_type, m_obj, R_obj, v_obj):
    """TDE of an object (a star drawn from the stellar mass function, or a WD) by the
    massive BH (pool slot k_M), as Rapster's TDEs (disrupt_star), the BH at rest."""
    BH = state['BH']
    M, s, g, h = BH.get(k_M)
    m_new, s_new, h_new = disrupt_star(state['seed'], state['t'], state['z'], tde_type, state['m_avg'], m_obj, R_obj,
                                       M, s, g, h, 0.0, v_obj, state['vBH'], state['tdes'], config['f_accreted'],
                                       config['EoS'])
    BH.set_mass(k_M, float(m_new)); BH.set_spin(k_M, float(s_new)); BH.set_tde_count(k_M, float(h_new))
    state['M_MBH'] = float(m_new)
    state['N_tdeMBH'] = state.get('N_tdeMBH', 0) + 1


def evolve_massive_bh(state, config):
    """GW captures of BHs (IMRIs) and TDEs of stars and WDs by the central massive BH in
    this step, Poisson numbers of the loss-cone rates, while it is a single BH (in a binary
    it takes part in the BBH evolution instead). Updates state in-place."""
    state['Gamma_IMRI'] = 0.0; state['Gamma_tdeMBH'] = 0.0
    if not phase2(config) or state.get('M_MBH', 0.0) <= 0 or state['i_aux1'] != 1:
        return
    k = slot_in_pool(state)
    BH = state['BH']
    if k < 0 or BH.size < 2:
        return
    dt = state['dt']
    M = BH.get(k)[0]
    N_other = BH.size - 1
    Gamma_IMRI = imri_rate(M, state['mBH_avg'], state['nc_BH'], state['vBH'], N_other)
    m_avg = state['m_avg']; n_star = state['n_star']; v_star = state['v_star']
    if config['with_tdes'] and n_star > 0:
        Gamma_star = tde_rate(M, m_avg, R_sun, n_star, v_star, m_avg, n_star)
        N_WD = state.get('N_WD', 0.0)
        if N_WD > 0:
            m_WD = white_dwarf_mass; R_WD = R_WhiteDwarf()
            v_WD = np.sqrt(m_avg / m_WD) * v_star
            nc_WD = 3 * N_WD / 4 / np.pi / (state['rh'] / 1.3)**3
            Gamma_WD = tde_rate(M, m_WD, R_WD, nc_WD, v_WD, m_avg, n_star)
        else:
            Gamma_WD = 0.0
    else:
        Gamma_star = Gamma_WD = 0.0
    state['Gamma_IMRI'] = Gamma_IMRI; state['Gamma_tdeMBH'] = Gamma_star + Gamma_WD

    k_imri = min(poisson_draw(min(dt * Gamma_IMRI, 1e8)), N_other)
    for _ in range(k_imri):
        k = slot_in_pool(state)
        if k < 0 or BH.size < 2 or not _imri(state, config, k):
            return
    k_star = poisson_draw(min(dt * Gamma_star, 1e8)) if Gamma_star > 0 else 0
    for _ in range(int(min(k_star, state['Mcl'] / m_avg))):
        k = slot_in_pool(state)
        if k < 0:
            return
        m_star, R_star = get_star(state['t'], tBH_form, config['m_min'], config['m_max'])
        _tde(state, config, k, TDE_TYPE_STAR, m_star, R_star, v_star)
    k_WD = poisson_draw(min(dt * Gamma_WD, 1e8)) if Gamma_WD > 0 else 0
    k_WD = int(min(k_WD, state.get('N_WD', 0.0)))
    for _ in range(k_WD):
        k = slot_in_pool(state)
        if k < 0:
            return
        _tde(state, config, k, TDE_TYPE_WD, white_dwarf_mass, R_WhiteDwarf(), np.sqrt(m_avg / white_dwarf_mass) * v_star)
        state['N_WD'] = state['N_WD'] - 1


def evolution_columns(state):
    """Values of EVOLUTION_KEYS_MBH for the current step."""
    M = state.get('M_MBH', 0.0)
    if M > 0:
        m, s, g = heaviest_object(state)
    else:
        m = s = g = 0.0
    return [M, s, g, state.get('r_infl', 0.0), state.get('f_3bb_MBH', 1.0), state.get('N_IMRI', 0),
            state.get('N_tdeMBH', 0), state.get('Gamma_IMRI', 0.0), state.get('Gamma_tdeMBH', 0.0)]
