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
from .functions import *
from .stellar_evolution import *
from .three_body_binary import three_body_binary
from .binary_evolution import evolve_BBHs
from .triples import evolve_triples
from .two_body_capture import two_body_capture
from .exchanges import StarStar_to_BHstar, BHstar_to_BBH
from .tidal_disruptions import BH_TidalDisruptions


def initialize_cluster(config):
    """Set up the initial cluster state from the configuration.

    Samples stellar masses from the Kroupa IMF, computes remnant masses using
    the selected prescription, applies supernova natal kicks, retains BHs below
    the escape velocity, initializes BH spins (uniform, monochromatic, or beta),
    optionally includes neutron stars and a massive BH seed, and allocates all
    data arrays (binaries, mergers, evolution, hardening, triples, pairs, tdes).

    Args:
        config (dict): Configuration dictionary from parse_args().

    Returns:
        dict: Simulation state dictionary containing all mutable variables
            (cluster properties, BH arrays, data arrays, counters, time variables).
    """
    # unpack configuration parameters into local variables:
    N = config['N']
    rh = config['rh']
    m_min = config['m_min']
    m_max = config['m_max']
    Z = config['Z']
    seed = config['seed']
    wSN_kick = config['wSN_kick']
    R_gal = config['R_gal']
    v_gal = config['v_gal']
    zCl_form = config['zCl_form']
    n_star = config['n_star']
    NKP = config['NKP']
    fb = config['fb']
    Bi = config['Bi']
    input_BH_file = config['input_BH_file']
    s1g_max = config['s1g_max']
    SD = config['SD']
    RP = config['RP']
    with_NSs = config['with_NSs']
    M_BH0 = config['M_BH0']
    s_BH0 = config['s_BH0']
    dt_min = config['dt_min']
    BMD = config['BMD']
    min_1g_bh_mass = config['min_1g_bh_mass']
    max_1g_bh_mass = config['max_1g_bh_mass']

    # initialize pseudo-random number generator:
    np.random.seed(seed)

    # initial half-mass radius:
    rh0 = rh

    # initial stellar number density:
    n_star0 = n_star

    # average stellar mass:
    Kroupa_norm = 1 / integrate.quad(IMF_kroupa, m_min, m_max)[0]
    m_avg = Kroupa_norm * integrate.quad(lambda x: x * IMF_kroupa(x), m_min, m_max)[0]

    # initial average stellar mass:
    m_avg0 = m_avg

    # initial cluster mass:
    Mcl = N * m_avg
    Mcl0 = Mcl

    # initial galactocentric radius:
    R_gal0 = R_gal

    # initial multimass relaxation factor:
    psi0 = integrate.quad(lambda x: x**(5/2) * IMF_kroupa(x), m_min, m_max)[0] \
        / integrate.quad(IMF_kroupa, m_min, m_max)[0] / m_avg**(5/2)

    # initial relaxation timescale:
    t_rlx0 = t_relax(Mcl0, rh0, m_avg, psi0, np.log(lc * N))

    # core collapse timescale:
    t_cc = k_cc * t_rlx0

    # number of binary stars:
    Nb = int(N * fb / 2)

    # minimum binary-star separation:
    ab_min = 3 * R_sun

    # maximum binary-star separation:
    ab_max = (3 / np.pi / 4 / n_star)**(1/3) / 10

    # initialize binary star population and filter hard binaries:
    if Nb>0:
        # semimajor axes of binaries (loguniform in [ab_min, ab_max]):
        ab = 10**(np.random.uniform(np.log10(ab_min), np.log10(ab_max), Nb))

        # stellar velocity dispersion:
        v_star = np.sqrt(0.4 * G_Newton * Mcl / rh)

        # binary star hard-soft boundary:
        ab_hard = G_Newton * m_avg / v_star

        # filter only hard binary stars:
        ab = ab[ab < ab_hard]
    else:
        ab = np.array([])
        v_star = np.sqrt(0.4 * G_Newton * Mcl / rh)

    # number of massive stars:
    N_massive = int(Mcl * integrate.quad(IMF_kroupa, mM_min, m_max)[0] \
                    / integrate.quad(lambda x: x * IMF_kroupa(x), m_min, m_max)[0])

    if BMD==0:
        # default: Kroupa IMF + stellar collapse + SN kicks

        # massive star masses drawn from Kroupa (2002) IMF:
        u_massive = np.random.rand(N_massive)
        m_massive = (u_massive * (m_max**(alphaIMF + 1) - mM_min**(alphaIMF + 1)) \
                     + mM_min**(alphaIMF + 1))**(1 / (alphaIMF + 1))

        # remnant masses:
        if   RP==0:
            m_rem = Mrem_SEVNdelayed(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==1:
            m_rem = Mrem_F12d(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==2:
            m_rem = Mrem_SEVNrapid(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==3:
            m_rem = Mrem_F12r(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)

        # Separate BHs from other remnants:
        mBH = m_rem[m_rem > mBH_min]

        # compute supernova kicks:
        if NKP==0:
            if   RP==0:
                vSN_kick = (1 - np.vectorize(f_fb_delayed)(m_massive[m_rem > mBH_min], np.vectorize(M_CO_SEVN)(m_massive[m_rem > mBH_min], Z))) * \
                           np.vectorize(get_SN_kick)(1.4 * np.ones(mBH.size), wSN_kick)
            elif RP==1:
                vSN_kick = (1 - np.vectorize(f_fb_delayed)(m_massive[m_rem > mBH_min], np.vectorize(M_CO_SSE)(m_massive[m_rem > mBH_min], Z))) * \
                           np.vectorize(get_SN_kick)(1.4 * np.ones(mBH.size), wSN_kick)
            elif RP==2:
                vSN_kick = (1 - np.vectorize(f_fb_rapid)(m_massive[m_rem > mBH_min], np.vectorize(M_CO_SEVN)(m_massive[m_rem > mBH_min], Z))) * \
                           np.vectorize(get_SN_kick)(1.4 * np.ones(mBH.size), wSN_kick)
            elif RP==3:
                vSN_kick = (1 - np.vectorize(f_fb_rapid)(m_massive[m_rem > mBH_min], np.vectorize(M_CO_SSE)(m_massive[m_rem > mBH_min], Z))) * \
                           np.vectorize(get_SN_kick)(1.4 * np.ones(mBH.size), wSN_kick)
        elif NKP==1:
            vSN_kick = np.vectorize(get_SN_kick)(mBH, wSN_kick)

        # retain BHs with SN kick < escape velocity:
        mBH = mBH[vSN_kick < v_esc(Mcl, rh)]

    else:
        # For BMD=1,2: run Kroupa path to determine N_BH, then replace masses.
        u_massive = np.random.rand(N_massive)
        m_massive = (u_massive * (m_max**(alphaIMF + 1) - mM_min**(alphaIMF + 1)) \
                     + mM_min**(alphaIMF + 1))**(1 / (alphaIMF + 1))
        if   RP==0:
            m_rem = Mrem_SEVNdelayed(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==1:
            m_rem = Mrem_F12d(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==2:
            m_rem = Mrem_SEVNrapid(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        elif RP==3:
            m_rem = Mrem_F12r(m_massive, Z) + 0.01 * np.random.rand(m_massive.size)
        N_BH_init = (m_rem > mBH_min).sum()

        if BMD==1:
            # uniform BH mass distribution in [min_1g_bh_mass, max_1g_bh_mass]:
            mBH = np.random.uniform(min_1g_bh_mass, max_1g_bh_mass, N_BH_init)

        elif BMD==2:
            # Salpeter power law (m^-2.35) in [min_1g_bh_mass, max_1g_bh_mass]:
            alpha_salpeter = -2.35
            u = np.random.rand(N_BH_init)
            mBH = (u * (max_1g_bh_mass**(alpha_salpeter + 1) - min_1g_bh_mass**(alpha_salpeter + 1)) \
                   + min_1g_bh_mass**(alpha_salpeter + 1))**(1 / (alpha_salpeter + 1))

        elif BMD==3:
            # log-uniform BH mass distribution in [min_1g_bh_mass, max_1g_bh_mass]:
            mBH = 10**(np.random.uniform(np.log10(min_1g_bh_mass), np.log10(max_1g_bh_mass), N_BH_init))

        # Apply momentum-conservation SN kicks and retain BHs below escape velocity.
        # Fallback kicks (NKP==0) are not available for BMD>0 because there is no
        # stellar progenitor or CO core mass to compute the fallback fraction from.
        # Momentum-conservation kicks only depend on the BH mass, so they apply regardless.
        vSN_kick = np.vectorize(get_SN_kick)(mBH, wSN_kick)
        mBH = mBH[vSN_kick < v_esc(Mcl, rh)]

    # optionally override BH masses from external file:
    if Bi==1:
        mBH = np.load(input_BH_file)['mBH_ini']
        mBH = mBH + 0.01 * np.random.rand(mBH.size)

    # BH spins:
    if SD==0:
        sBH = np.random.uniform(0, s1g_max, mBH.size)
    elif SD==1:
        sBH = s1g_max * np.ones(mBH.size)
    elif SD==2:
        rng = np.random.default_rng()
        sBH = rng.beta(1.4, 3.6, mBH.size) * s1g_max

    # optionally form and retain neutron stars based on natal kicks:
    if with_NSs==1:
        f_NS_form = Kroupa_norm*integrate.quad(IMF_kroupa, 8, 18)[0]
        N_NS_form = int(f_NS_form*N)
        v_NS_natal = maxwell.rvs(loc=0, scale=np.sqrt(3)*wSN_kick, size=N_NS_form)
        N_NS_ret = v_NS_natal[v_NS_natal < v_esc(Mcl, rh)].size

        if N_NS_ret > 0:
            mBH = np.concatenate((mBH, neutron_star_mass * np.ones(N_NS_ret)))
            sBH = np.concatenate((sBH, np.zeros(N_NS_ret)))

    # optionally add a massive BH seed to the initial conditions:
    if M_BH0>0:
        mBH = np.concatenate((mBH, [M_BH0]))
        sBH = np.concatenate((sBH, [s_BH0]))

    # BH generations:
    gBH = np.ones(mBH.size)
    hBH = np.zeros(mBH.size)

    # data arrays:
    binaries = np.zeros(shape=(1, 15))
    pairs = np.zeros(shape=(1, 5))
    triples = np.zeros(shape=(1, 24))
    mergers = np.zeros(shape=(1, 27))
    evolution = np.zeros(shape=(1, 69))
    hardening = np.zeros(shape=(1, 12))
    tdes = np.zeros(shape=(1, 18))

    # bundle all mutable simulation variables into the state dictionary:
    state = {
        # cluster properties:
        'Mcl': Mcl, 'Mcl0': Mcl0, 'rh': rh, 'rh0': rh0, 'R_gal': R_gal, 'R_gal0': R_gal0,
        'v_gal': v_gal, 'n_star': n_star, 'n_star0': n_star0,
        'm_avg': m_avg, 'm_avg0': m_avg0, 'v_star': v_star,
        'Kroupa_norm': Kroupa_norm, 'Nb': Nb, 'ab': ab,
        't_cc': t_cc, 't_rlx': 0.0,
        'N': N, 'Z': Z,
        # BH properties:
        'mBH': mBH, 'sBH': sBH, 'gBH': gBH, 'hBH': hBH, 
        'mBH_avg': 0.0, 'vBH': 0.0,
        # data arrays:
        'binaries': binaries, 'pairs': pairs, 'triples': triples,
        'mergers': mergers, 'evolution': evolution, 'hardening': hardening, 'tdes': tdes,
        # counters:
        'N_BH': 0, 'N_BBH': 0, 'N_3bb': 0, 'N_me': 0, 'N_meRe': 0, 'N_meEj': 0,
        'N_2cap': 0, 'N_3cap': 0, 'N_dis': 0, 'N_ex': 0, 'N_BHej': 0, 'N_BBHej': 0,
        'N_iter': 0, 'N_bb': 0, 'N_meFi': 0, 'N_me2b': 0,
        'N_ex1': 0, 'N_ex2': 0, 'N_BHstar': 0, 'N_pp': 0,
        'N_Triples': 0, 'N_ZLK': 0, 'N_WD': 0,
        'N_tdeBHWD': 0, 'N_tdeBHstar': 0, 'N_hardening': 0,
        # time:
        't': 0, 'z': zCl_form, 'dt': dt_min, 'zCl_form': zCl_form, 'seed': seed,
        # aux:
        'i_aux1': 0,
        # tracking lists:
        'simulation_times': [], 'black_hole_masses': [], 'black_hole_spins': [], 'black_hole_generations': [], 'black_hole_tdes': []
    }

    return state


def compute_cluster_properties(state, config):
    """Compute BH subsystem and cluster structural properties for the current timestep.

    Saves BH mass/spin/generation snapshots, computes the average and maximum BH
    mass, Spitzer mass-segregation parameter, BH half-mass and core radii, velocity
    dispersions, number densities, relaxation timescales, and the cluster Coulomb
    logarithm. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.

    Returns:
        bool: True if the simulation should continue, False if the BH subsystem
            has evaporated (N_BH <= 0).
    """

    # unpack current state into local variables:
    t = state['t']
    mBH = state['mBH']
    binaries = state['binaries']
    pairs = state['pairs']
    i_aux1 = state['i_aux1']
    N_BH = state['N_BH']
    N_BBH = state['N_BBH']
    N_BHstar = state['N_BHstar']
    N_Triples = state['N_Triples']
    Mcl = state['Mcl']
    rh = state['rh']
    m_avg = state['m_avg']
    Nb = state['Nb']
    ab = state['ab']
    n_star0 = state['n_star0']
    rh0 = state['rh0']
    Mcl0 = state['Mcl0']

    # save current BH properties:
    state['simulation_times'].append(t)
    state['black_hole_masses'].append(np.concatenate((mBH, binaries[:, 4], binaries[:, 5], pairs[:, 1])))
    state['black_hole_spins'].append(np.concatenate((state['sBH'], binaries[:, 6], binaries[:, 7], pairs[:, 2])))
    state['black_hole_generations'].append(np.concatenate((state['gBH'], binaries[:, 8], binaries[:, 9], pairs[:, 3])))

    # activate BH subsystem once BH formation time is reached:
    if t>tBH_form and i_aux1==0:
        N_BH = mBH.size
        i_aux1 = 1

    # check if BH subsystem has fully evaporated:
    if i_aux1==1 and N_BH<=0:
        print('BH SUBSYSTEM EVAPORATED')
        state['i_aux1'] = i_aux1
        state['N_BH'] = N_BH
        return False

    # average BH mass:
    if i_aux1==1:
        mBH_avg = (np.sum(mBH) + np.sum(binaries[:, 4]) + np.sum(binaries[:, 5]) + np.sum(pairs[:, 1]) ) / N_BH
    else:
        mBH_avg = 0.0

    # maximum BH mass:
    if i_aux1==1:
        try:
            mBHs_max = np.max(mBH)
        except Exception:
            mBHs_max = 0.0
        try:
            mBH1_max = np.max(binaries[:, 4])
        except Exception:
            mBH1_max = 0.0
        try:
            mBH2_max = np.max(binaries[:, 5])
        except Exception:
            mBH2_max = 0.0
        try:
            mBHp_max = np.max(pairs[:, 1])
        except Exception:
            mBHp_max = 0.0
        mBH_max = np.max([mBHs_max, mBH1_max, mBH2_max, mBHp_max])
    else:
        mBH_max = 0.0

    # total BH mass ratio:
    Q_BH = N_BH * mBH_avg / Mcl

    # individual BH mass ratio:
    q_BH = mBH_avg / m_avg

    # Cluster Coulomb logarithm:
    logLcl = 10.0

    # BH Coulomb logarithm:
    if i_aux1==1 and N_BH>0:
        logLBH = 1.0
    else:
        logLBH = 0.0

    # Spitzer parameter:
    S = Q_BH * q_BH**(3/2)

    # Spitzer instability condition:
    if S < 0.16 and S > 0:
        xi = 1
        S = Q_BH * q_BH
    else:
        if i_aux1==1:
            xi = q_BH**(3/5) * Q_BH**(2/5) * (logLBH / logLcl)**(-2/5)
        else:
            xi = 0.0

    # multimass relaxation factor:
    psi = 1 + S

    # BH relaxation factor:
    if i_aux1==1 and N_BH>0:
        psi_BH = (np.sum(mBH**(5/2)) + np.sum(binaries[:, 4]**(5/2)) + \
                  np.sum(binaries[:, 5]**(5/2)) + np.sum(pairs[:, 1]**(5/2)) ) \
                  / N_BH / mBH_avg**(5/2)
    else:
        psi_BH = 0.0

    # BH half-mass radius:
    if i_aux1==1:
        rh_BH = 1/xi * Q_BH * q_BH * rh
    else:
        rh_BH = 0.0

    # BH core radius:
    if i_aux1==1:
        Core_parameter = 1.76 / Hardening_constant / psi_BH * (12 * q_BH / xi - 1) * eta_min**(-11/2) * (1+3*eta_min) * (1 + 6*eta_min) * P_3bb
        rc_BH = rh_BH * N_BH**(-2/3) * (Core_parameter / zeta_BH / logLBH)**(1/3)
    else:
        rc_BH = 0.0

    # BH velocity dispersion:
    if i_aux1==1:
        vBH = np.sqrt(0.4 * G_Newton * mBH_avg * N_BH / rh_BH)
    else:
        vBH = 0.0

    # star velocity dispersion:
    v_star = np.sqrt(0.4 * G_Newton * Mcl / rh)

    # BH half-mass volume:
    Vh_BH = 4 / 3 * np.pi * rh_BH**3

    # BH half-mass density:
    if i_aux1==1:
        nh_BH = N_BH / 2 / Vh_BH
    else:
        nh_BH = 0.0

    # BH core density:
    if i_aux1==1:
        nc_BH = nh_BH * (rc_BH / rh_BH)**(-2)
    else:
        nc_BH = 0.0

    # BH core volume:
    Vc_BH = 4 / 3 * np.pi * rc_BH**3

    # mean BH density:
    if i_aux1==1:
        na_BH = nc_BH * 2/3
    else:
        na_BH = 0.0

    # stellar density updated:
    n_star = n_star0 * (rh0 / rh)**3 * (Mcl / Mcl0)

    # half-mass volume:
    Vh = 4 * np.pi / 3 * rh**3

    if Nb>0:
        nb = ab.size / Vh
    else:
        nb = 0.0

    # relaxation timescale:
    t_rlx = t_relax(Mcl, rh, m_avg, psi, logLcl)

    # BH relaxation timescale:
    if i_aux1==1:
        tBH_rlx = t_relax(N_BH*mBH_avg, rh_BH, mBH_avg, psi_BH, logLBH)
    else:
        tBH_rlx = 1e100

    # BH subsystem core collapse:
    if i_aux1==1:
        tBH_cc = k_cc * tBH_rlx
    else:
        tBH_cc = 1e100

    # Number of BHs in the core:
    Nc_BH = nc_BH * Vc_BH

    # update state:
    state['i_aux1'] = i_aux1
    state['N_BH'] = N_BH
    state['mBH_avg'] = mBH_avg
    state['mBH_max'] = mBH_max
    state['Q_BH'] = Q_BH
    state['q_BH'] = q_BH
    state['S'] = S
    state['xi'] = xi
    state['psi'] = psi
    state['psi_BH'] = psi_BH
    state['rh_BH'] = rh_BH
    state['rc_BH'] = rc_BH
    state['vBH'] = vBH
    state['v_star'] = v_star
    state['nh_BH'] = nh_BH
    state['nc_BH'] = nc_BH
    state['na_BH'] = na_BH
    state['Vc_BH'] = Vc_BH
    state['Nc_BH'] = Nc_BH
    state['n_star'] = n_star
    state['nb'] = nb
    state['t_rlx'] = t_rlx
    state['tBH_rlx'] = tBH_rlx
    state['logLcl'] = logLcl
    state['logLBH'] = logLBH

    return True


def compute_timescales(state, config):
    """Compute dynamical timescales and the adaptive timestep.

    Calculates the three-body binary formation (t_3bb), two-body gravitational
    capture (t_2cap), star-star to BH-star exchange (t_ex1), and BH-star to BBH
    exchange (t_ex2) timescales. Sets the adaptive timestep dt as the minimum of
    these, bounded by [dt_min, dt_max]. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """

    # unpack current state into local variables:
    t = state['t']
    t_cc = state['t_cc']
    N_BH = state['N_BH']
    N_BBH = state['N_BBH']
    N_BHstar = state['N_BHstar']
    N_Triples = state['N_Triples']
    mBH_avg = state['mBH_avg']
    nc_BH = state['nc_BH']
    vBH = state['vBH']
    Vc_BH = state['Vc_BH']
    Nc_BH = state['Nc_BH']
    nb = state['nb']
    v_star = state['v_star']
    m_avg = state['m_avg']
    binaries = state['binaries']
    pairs = state['pairs']
    ab = state['ab']
    dt_min = config['dt_min']
    dt_max = config['dt_max']

    # number of single (unbound) BHs available for interactions:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples

    # 3bb timescale:
    if t>t_cc and N_BHsin>2:
        try:
            t_3bb = 1 / Rate_3bb(mBH_avg, nc_BH, vBH) / Vc_BH
        except Exception:
            t_3bb = 1e100
    else:
        t_3bb = 1e100

    # 2-body capture timescale:
    if t>t_cc and N_BHsin>1:
        try:
            t_2cap = 1 / Rate_cap(mBH_avg, nc_BH, vBH) / Vc_BH
        except Exception:
            t_2cap = 1e100
    else:
        t_2cap = 1e100

    # star-star -> BH-star timescale:
    if t>t_cc and N_BHsin>0 and nb>0:
        try:
            t_ex1 = 1 / Rate_exc(m_avg, m_avg, mBH_avg, nb, np.sqrt(v_star**2 + vBH**2), np.mean(ab)) / Nc_BH / 2
        except Exception:
            t_ex1 = 1e100
    else:
        t_ex1 = 1e100

    # BH-star -> BH-BH timescale:
    if t>t_cc and N_BHstar>0 and N_BHsin>0:
        if N_BHstar>0:
            t_ex2 = 1 / Rate_exc(m_avg, np.mean(pairs[:, 1]), mBH_avg, nc_BH, vBH, np.mean(pairs[:, 0])) / N_BHstar
        else:
            t_ex2 = 1e100
    else:
        t_ex2 = 1e100

    # adaptive time-step:
    if t<t_cc:
        dt = dt_min
    else:
        dt = np.min([dt_max, np.max([dt_min, np.min([t_3bb, t_2cap, t_ex1, t_ex2])])])

    if dt > t and t>0:
        dt = t

    state['t_3bb'] = t_3bb
    state['t_2cap'] = t_2cap
    state['t_ex1'] = t_ex1
    state['t_ex2'] = t_ex2
    state['dt'] = dt


def form_binaries(state, config):
    """Form new binary BHs through three-body, two-body capture, and exchange channels.

    Draws Poisson-sampled event counts for each formation channel and executes:
      - Three-body binary (3bb) formation (channel 3)
      - Two-body gravitational wave capture (channel 2)
      - Star-star -> BH-star exchanges (channel 1, step 1)
      - BH-star -> BBH exchanges (channel 1, step 2)
    Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """
    
    # unpack current state into local variables:
    t = state['t']
    z = state['z']
    dt = state['dt']
    zCl_form = state['zCl_form']
    mBH_avg = state['mBH_avg']
    binaries = state['binaries']
    mBH = state['mBH']
    sBH = state['sBH']
    gBH = state['gBH']
    hBH = state['hBH']
    vBH = state['vBH']
    v_star = state['v_star']
    m_avg = state['m_avg']
    ab = state['ab']
    pairs = state['pairs']
    mergers = state['mergers']
    seed = state['seed']
    nc_BH = state['nc_BH']
    N_BH = state['N_BH']
    N_BBH = state['N_BBH']
    N_BHstar = state['N_BHstar']
    N_Triples = state['N_Triples']
    N_3bb = state['N_3bb']
    N_2cap = state['N_2cap']
    N_me = state['N_me']
    N_meRe = state['N_meRe']
    N_meEj = state['N_meEj']
    N_ex1 = state['N_ex1']
    N_ex2 = state['N_ex2']
    t_3bb = state['t_3bb']
    t_2cap = state['t_2cap']

    # number of single (unbound) BHs available for interactions:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples

    # number of 3bbs:
    k_3bb = np.min([poisson.rvs(mu=dt / t_3bb), int(N_BHsin / 3)])

    # 3bb formation:
    t, z, k_3bb, mBH_avg, binaries, mBH, sBH, gBH, hBH, vBH, N_3bb, N_BBH = three_body_binary(t, z, k_3bb, mBH_avg, binaries, mBH, sBH, gBH, hBH, vBH, N_3bb, N_BBH, random_pairing=config['random_pairing'])

    # number of 2-body captures:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples
    k_2cap = np.min([poisson.rvs(mu=dt / t_2cap), int(N_BHsin / 2)])

    # 2-body capture(s):
    seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, mBH, sBH, gBH, hBH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers = \
        two_body_capture(seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, mBH, sBH, gBH, hBH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers, random_pairing=config['random_pairing'])

    # number of star-star -> BH-star exchanges:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples
    k_ex1 = np.min([poisson.rvs(mu=dt / state['t_ex1']), int(N_BHsin)])

    # star-star -> BH-star exchange(s):
    if k_ex1 > 0:
        k_ex1, N_ex1, m_avg, mBH, sBH, gBH, hBH, ab, pairs, N_BHstar = StarStar_to_BHstar(k_ex1, N_ex1, m_avg, mBH, sBH, gBH, hBH, ab, pairs, N_BHstar)

    # number of BH-star -> BH-BH exchanges:
    N_BHsin = N_BH - 2*N_BBH - N_BHstar - 3*N_Triples
    k_ex2 = np.min([poisson.rvs(mu=dt / state['t_ex2']), int(N_BHsin), int(N_BHstar)])

    # BH-star -> BBH exchange(s):
    if k_ex2 > 0:
        t, z, k_ex2, N_ex2, m_avg, mBH, sBH, gBH, hBH, pairs, binaries, N_BBH, N_BHstar = BHstar_to_BBH(t, z, k_ex2, N_ex2, m_avg, mBH, sBH, gBH, hBH, pairs, binaries, N_BBH, N_BHstar)

    # write back:
    state['t'] = t; state['z'] = z; state['dt'] = dt; state['seed'] = seed
    state['mBH_avg'] = mBH_avg; state['binaries'] = binaries; state['mBH'] = mBH
    state['sBH'] = sBH; state['gBH'] = gBH; state['hBH'] = hBH; state['vBH'] = vBH; state['v_star'] = v_star
    state['m_avg'] = m_avg; state['ab'] = ab; state['pairs'] = pairs; state['mergers'] = mergers
    state['N_BH'] = N_BH; state['N_BBH'] = N_BBH; state['N_BHstar'] = N_BHstar
    state['N_3bb'] = N_3bb; state['N_2cap'] = N_2cap; state['N_me'] = N_me
    state['N_meRe'] = N_meRe; state['N_meEj'] = N_meEj
    state['N_ex1'] = N_ex1; state['N_ex2'] = N_ex2
    state['k_3bb'] = k_3bb; state['k_2cap'] = k_2cap
    state['k_ex1'] = k_ex1; state['k_ex2'] = k_ex2
    state['zCl_form'] = zCl_form


def evolve_interactions(state, config):
    """Evolve BBH binaries, pair-pair interactions, and hierarchical triples.

    Calls evolve_BBHs to harden/disrupt/merge existing BBHs through single-binary
    and binary-binary encounters, processes BH-star pair-pair interactions that can
    form new BBHs, and evolves hierarchical triples via von Zeipel-Lidov-Kozai
    (ZLK) oscillations. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """
    
    # unpack current state into local variables:
    seed = state['seed']
    t = state['t']; z = state['z']; dt = state['dt']; zCl_form = state['zCl_form']
    binaries = state['binaries']; hardening = state['hardening']; mergers = state['mergers']
    mBH = state['mBH']; sBH = state['sBH']; gBH = state['gBH']; hBH = state['hBH']
    n_star = state['n_star']; v_star = state['v_star']; vBH = state['vBH']
    t_rlx = state['t_rlx']; m_avg = state['m_avg']; mBH_avg = state['mBH_avg']
    na_BH = state['na_BH']; nc_BH = state['nc_BH']
    N_BH = state['N_BH']; N_BBH = state['N_BBH']; N_me = state['N_me']
    N_me2b = state['N_me2b']; N_3cap = state['N_3cap']; N_meFi = state['N_meFi']
    N_meRe = state['N_meRe']; N_meEj = state['N_meEj']
    N_dis = state['N_dis']; N_ex = state['N_ex']
    N_BHej = state['N_BHej']; N_BBHej = state['N_BBHej']
    N_hardening = state['N_hardening']; Vc_BH = state['Vc_BH']
    N_bb = state['N_bb']; triples = state['triples']; N_Triples = state['N_Triples']
    pairs = state['pairs']; N_BHstar = state['N_BHstar']
    N_pp = state['N_pp']; N_ZLK = state['N_ZLK']
    i_aux1 = state['i_aux1']
    t_cc = state['t_cc']

    # BBH evolution:
    seed, t, z, dt, zCl_form, binaries, hardening, mergers, mBH, sBH, gBH, hBH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples = evolve_BBHs(seed, t, z, dt, zCl_form, binaries, hardening, mergers, mBH, sBH, gBH, hBH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples)

    # average BBH number density:
    if i_aux1==1:
        n_BBH = N_BBH / Vc_BH
    else:
        n_BBH = 0.0

    # average BH-star number density:
    if i_aux1==1:
        n_BHstar = N_BHstar / Vc_BH
    else:
         n_BHstar = 0.0

    # Binary-binary interaction timescale:
    if t>t_cc and N_BBH>1:
        try:
            t_bb = 1 / Rate_int(np.mean(binaries[:, 4]+binaries[:, 5]), n_BBH, vBH, 2 * np.mean(binaries[:, 2])) / N_BBH
        except Exception:
            t_bb = 1e100
    else:
        t_bb = 1e100

    # Pair-pair interaction timescale:
    if t>t_cc and N_BHstar>1:
        try:
            t_pp = 1 / Rate_int(np.mean(pairs[:, 1]) + m_avg, n_BHstar, vBH, 2 * np.mean(pairs[:, 0])) / N_BHstar
        except Exception:
            t_pp = 1e100
    else:
        t_pp = 1e100

    # number of pair-pair interactions:
    dt = state['dt']
    k_pp = np.min([poisson.rvs(mu=dt / t_pp), int(N_BHstar / 2)])

    # execute pair-pair interactions: two BH-star pairs form a new BBH:
    if k_pp>0:
        for i in range(k_pp):
            N_pp+=1

            smas = pairs[:, 0]
            masses = pairs[:, 1]

            a1, a2 = np.random.choice(smas, size=2, replace=False, p=smas*masses / np.sum(smas*masses))

            k1 = np.squeeze(np.where(smas==a1))+0
            k2 = np.squeeze(np.where(smas==a2))+0

            if isinstance(k1, np.ndarray):
                k1=k1[0]
            if isinstance(k2, np.ndarray):
                k2=k2[0]

            m1 = pairs[k1][1]; s1 = pairs[k1][2]; g1 = pairs[k1][3]; h1 = pairs[k1][4]
            m2 = pairs[k2][1]; s2 = pairs[k2][2]; g2 = pairs[k2][3]; h2 = pairs[k2][4]

            if m2/a2 < m1/a1:
                sma = m2 / m_avg * a1
            else:
                sma = m1 / m_avg * a2

            eccen = np.sqrt(np.random.rand())

            binaries = np.append(binaries, [[np.random.randint(0, 999999999), 1, sma, eccen, m1, m2, s1, s2, g1, g2, t, z, 0, h1, h2]], axis=0)
            pairs = np.delete(pairs, [k1, k2], axis=0)
            N_BHstar = N_BHstar - 2
            N_BBH+=1

    # Triple evolution:
    seed, t, z, zCl_form, triples, binaries, mBH, sBH, gBH, hBH, mBH_avg, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, v_star, vBH, nc_BH, mergers = \
        evolve_triples(seed, t, z, zCl_form, triples, binaries, mBH, sBH, gBH, hBH, mBH_avg, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, v_star, vBH, nc_BH, mergers)

    # write back:
    state['seed'] = seed; state['t'] = t; state['z'] = z; state['dt'] = dt; state['zCl_form'] = zCl_form
    state['binaries'] = binaries; state['hardening'] = hardening; state['mergers'] = mergers
    state['mBH'] = mBH; state['sBH'] = sBH; state['gBH'] = gBH; state['hBH'] = hBH
    state['n_star'] = n_star; state['v_star'] = v_star; state['vBH'] = vBH
    state['t_rlx'] = t_rlx; state['m_avg'] = m_avg; state['mBH_avg'] = mBH_avg
    state['na_BH'] = na_BH; state['nc_BH'] = nc_BH
    state['N_BH'] = N_BH; state['N_BBH'] = N_BBH; state['N_me'] = N_me
    state['N_me2b'] = N_me2b; state['N_3cap'] = N_3cap; state['N_meFi'] = N_meFi
    state['N_meRe'] = N_meRe; state['N_meEj'] = N_meEj
    state['N_dis'] = N_dis; state['N_ex'] = N_ex
    state['N_BHej'] = N_BHej; state['N_BBHej'] = N_BBHej
    state['N_hardening'] = N_hardening; state['Vc_BH'] = Vc_BH
    state['N_bb'] = N_bb; state['triples'] = triples; state['N_Triples'] = N_Triples
    state['pairs'] = pairs; state['N_BHstar'] = N_BHstar
    state['N_pp'] = N_pp; state['N_ZLK'] = N_ZLK
    state['t_bb'] = t_bb; state['t_pp'] = t_pp; state['k_pp'] = k_pp


def evolve_tdes(state, config):
    """Evolve white dwarf population and process tidal disruption events.

    Computes WD formation, evaporation, and TDE rates. Draws Poisson-sampled
    BH-WD and BH-star TDE counts, executes each TDE (updating BH mass and spin
    via accretion), and records event parameters. Skipped if config['with_tdes']
    is False. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """
    
    # skip TDEs if disabled; zero out all TDE-related state:
    if not config['with_tdes']:
        state['N_WD'] = 0; state['v_WD'] = 0
        state['k_tdeBHWD'] = 0; state['N_tdeBHWD'] = state.get('N_tdeBHWD', 0)
        state['dN_WDformdt'] = 0; state['dN_WDevdt'] = 0; state['dN_tdeBHWDdt'] = 0
        state['k_tdeBHstar'] = 0; state['dN_tdeBHstardt'] = 0
        state['N_tdeBHstar'] = state.get('N_tdeBHstar', 0)
        return

    # unpack current state into local variables:
    seed = state['seed']; t = state['t']; z = state['z']; dt = state['dt']
    mBH = state['mBH']; sBH = state['sBH']; gBH = state['gBH']; hBH = state['hBH']
    vBH = state['vBH']; v_star = state['v_star']
    mBH_avg = state['mBH_avg']; m_avg = state['m_avg']
    N_BH = state['N_BH']; N_Triples = state['N_Triples']
    N_WD = state['N_WD']; N_tdeBHWD = state['N_tdeBHWD']; N_tdeBHstar = state['N_tdeBHstar']
    n_star = state['n_star']; t_rlx = state['t_rlx']
    Mcl = state['Mcl']; rh = state['rh']
    binaries = state['binaries']; pairs = state['pairs']; tdes = state['tdes']
    Kroupa_norm = state['Kroupa_norm']
    xi_e = state['xi_e']
    m_min = config['m_min']
    m_max = config['m_max']

    # white dwarf formation parameters (solar lifetime and max WD progenitor mass):
    solar_life = 1.0e4
    m_WN = 8.0
    t_WN = solar_life/m_WN**(2.5)
    N_strs = Mcl/m_avg

    # WD formation rate (/Myr):
    dN_WDformdt = N_strs/t/2.5 * Kroupa_norm*IMF_kroupa(np.array([(solar_life/t)**(1/2.5)]))[0] * (solar_life/t)**(1/2.5) if t>t_WN else 0.0

    # WD evaporation rate (/Myr):
    dN_WDevdt = xi_e*N_WD/t_rlx if N_WD>0 else 0.0

    # white dwarf properties and BH-WD TDE rate:
    m_WD = white_dwarf_mass
    R_WD = R_WhiteDwarf()
    v_WD = np.sqrt(m_avg/m_WD)*v_star
    rh_WD = rh
    nc_WD = 3*N_WD/4/np.pi/(rh_WD/1.3)**3
    v_BHWD = np.sqrt(vBH**2 + v_WD**2)
    dN_tdeBHWDdt = 2*np.sqrt(2*(3*np.pi-8))*G_Newton*(mBH_avg + m_WD)*N_BH*nc_WD*R_WD*(mBH_avg/m_WD)**(1/3)/v_BHWD if N_WD>0 else 0.0

    # effective WD evolution equation:
    dN_WDdt = dN_WDformdt - dN_WDevdt - dN_tdeBHWDdt
    N_WD = N_WD + dN_WDdt * dt

    # Poisson number of BH-WD TDEs:
    k_tdeBHWD = np.min([poisson.rvs(mu=dt*dN_tdeBHWDdt), int(N_BH-3*N_Triples), int(N_WD)]) if N_WD>0 else 0

    # execute BH-WD tidal disruption events:
    if k_tdeBHWD > 0:
        tde_type = 11
        seed, t, z, k_tdeBHWD, N_tdeBHWD, tde_type, m_avg, m_WD, R_WD, mBH, sBH, gBH, hBH, v_star, vBH, tdes, binaries, pairs = BH_TidalDisruptions(seed, t, z, k_tdeBHWD, N_tdeBHWD, tde_type, m_avg, m_WD, R_WD, mBH, sBH, gBH, hBH, v_WD, vBH, tdes, binaries, pairs)

    # micro-TDEs:
    v_BHstar = np.sqrt(v_star**2 + vBH**2)
    dN_tdeBHstardt = 2*np.sqrt(2*(3*np.pi-8))*G_Newton*(mBH_avg + m_avg)*n_star*N_BH*R_sun*(mBH_avg/m_avg)**(1/3)/v_BHstar if (N_BH>0)*(n_star>0) else 0.0

    # Poisson number of BH-star TDEs at this timestep:
    N_strs = Mcl/m_avg
    k_tdeBHstar = np.min([poisson.rvs(mu=dt*dN_tdeBHstardt), int(N_BH-3*N_Triples), int(N_strs)]) if (N_BH>0)*(n_star>0) else 0

    # execute BH-star tidal disruption events (micro-TDEs):
    if k_tdeBHstar > 0:
        tde_type = 1

        # Sample star mass from evolving mass function:
        m_star, R_star = get_star(t, tBH_form, m_min, m_max)

        seed, t, z, k_tdeBHstar, N_tdeBHstar, tde_type, m_avg, mstar, Rstar, mBH, sBH, gBH, hBH, v_star, vBH, tdes, binaries, pairs = BH_TidalDisruptions(seed, t, z, k_tdeBHstar, N_tdeBHstar, tde_type, m_avg, m_star, R_star, mBH, sBH, gBH, hBH, v_star, vBH, tdes, binaries, pairs)

    # write back:
    state['seed'] = seed; state['t'] = t; state['z'] = z
    state['mBH'] = mBH; state['sBH'] = sBH; state['gBH'] = gBH; state['hBH'] = hBH
    state['v_star'] = v_star; state['vBH'] = vBH
    state['binaries'] = binaries; state['pairs'] = pairs; state['tdes'] = tdes
    state['N_WD'] = N_WD; state['v_WD'] = v_WD
    state['N_tdeBHWD'] = N_tdeBHWD; state['N_tdeBHstar'] = N_tdeBHstar
    state['k_tdeBHWD'] = k_tdeBHWD; state['k_tdeBHstar'] = k_tdeBHstar
    state['dN_WDformdt'] = dN_WDformdt; state['dN_WDevdt'] = dN_WDevdt
    state['dN_tdeBHWDdt'] = dN_tdeBHWDdt; state['dN_tdeBHstardt'] = dN_tdeBHstardt


def record_evolution(state):
    """Append a row of 69 time-dependent quantities to the evolution array.

    Records the current cluster and BH subsystem state (masses, radii, densities,
    velocities, timescales, event counts, TDE statistics) as a single row in
    state['evolution']. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
    """
    
    # unpack current state into local variables:
    seed = state['seed']; t = state['t']; z = state['z']; dt = state['dt']
    m_avg = state['m_avg']; Mcl = state['Mcl']; rh = state['rh']
    R_gal = state['R_gal']; v_gal = state['v_gal']
    t_rlx = state['t_rlx']; tBH_rlx = state.get('tBH_rlx', 1e100)
    n_star = state['n_star']; N_BH = state['N_BH']
    mBH_avg = state['mBH_avg']; mBH_max = state.get('mBH_max', 0.0)
    rh_BH = state.get('rh_BH', 0.0); rc_BH = state.get('rc_BH', 0.0)
    S = state.get('S', 0.0); xi = state.get('xi', 0.0)
    psi = state.get('psi', 0.0); psi_BH = state.get('psi_BH', 0.0)
    t_3bb = state.get('t_3bb', 1e100); t_2cap = state.get('t_2cap', 1e100)
    k_3bb = state.get('k_3bb', 0); k_2cap = state.get('k_2cap', 0)
    N_me = state['N_me']; N_BBH = state['N_BBH']
    N_meRe = state['N_meRe']; N_meEj = state['N_meEj']
    v_star = state['v_star']; vBH = state['vBH']
    nh_BH = state.get('nh_BH', 0.0); nc_BH = state.get('nc_BH', 0.0)
    na_BH = state.get('na_BH', 0.0)
    N_3bb = state['N_3bb']; N_2cap = state['N_2cap']
    N_3cap = state['N_3cap']; N_BHej = state['N_BHej']
    N_BBHej = state['N_BBHej']; N_dis = state['N_dis']
    N_ex = state['N_ex']
    t_bb = state.get('t_bb', 1e100); N_bb = state['N_bb']
    N_meFi = state['N_meFi']; N_me2b = state['N_me2b']
    t_ex1 = state.get('t_ex1', 1e100); t_ex2 = state.get('t_ex2', 1e100)
    k_ex1 = state.get('k_ex1', 0); k_ex2 = state.get('k_ex2', 0)
    N_ex1 = state['N_ex1']; N_ex2 = state['N_ex2']
    N_BHstar = state['N_BHstar']
    t_pp = state.get('t_pp', 1e100); k_pp = state.get('k_pp', 0)
    N_pp = state['N_pp']
    N_Triples = state['N_Triples']; N_ZLK = state['N_ZLK']
    N_WD = state['N_WD']; v_WD = state.get('v_WD', 0.0)
    k_tdeBHWD = state.get('k_tdeBHWD', 0); N_tdeBHWD = state['N_tdeBHWD']
    dN_WDformdt = state.get('dN_WDformdt', 0.0); dN_WDevdt = state.get('dN_WDevdt', 0.0)
    dN_tdeBHWDdt = state.get('dN_tdeBHWDdt', 0.0)
    k_tdeBHstar = state.get('k_tdeBHstar', 0)
    dN_tdeBHstardt = state.get('dN_tdeBHstardt', 0.0)
    N_tdeBHstar = state['N_tdeBHstar']

    # append a row of 69 time-dependent quantities to the evolution array:
    state['evolution'] = np.append(state['evolution'], [[seed, t, z, dt, m_avg, Mcl, rh, R_gal, v_gal, t_rlx, tBH_rlx, n_star, N_BH, mBH_avg, mBH_max, rh_BH, rc_BH, S,
                                       xi, psi, psi_BH, t_3bb, t_2cap, k_3bb, k_2cap, N_me, N_BBH, N_meRe, N_meEj, v_star, vBH,
                                       nh_BH, nc_BH, na_BH, N_3bb, N_2cap, N_3cap, N_BHej, N_BBHej, N_dis, N_ex, t_bb, N_bb,
                                       N_meFi, N_me2b, t_ex1, t_ex2, k_ex1, k_ex2, N_ex1, N_ex2, N_BHstar, t_pp, k_pp, N_pp, 2*v_star,
                                       2*vBH, N_Triples, N_ZLK, N_WD, v_WD, k_tdeBHWD, N_tdeBHWD, dN_WDformdt, dN_WDevdt, dN_tdeBHWDdt, k_tdeBHstar,
                                       dN_tdeBHstardt, N_tdeBHstar]], axis=0)


def compute_external_params(state, config):
    """Compute external/environmental parameters for the cluster.

    Calculates the Jacobi (tidal) radius, the dimensionless escape rate xi_e,
    and the dynamical friction timescale t_df. These are needed before evolving
    TDEs (xi_e enters the WD evaporation rate) and before updating the cluster
    mass and radius. Updates state in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
    """
    
    # unpack current state into local variables:
    Mcl = state['Mcl']; rh = state['rh']; R_gal = state['R_gal']
    v_gal = state['v_gal']

    # Jacobi radius:
    rJ = (G_Newton * Mcl * R_gal**2 / 3 / v_gal**2)**(1/3)

    # dimensionless escape rate:
    xi_e = xi_e0 * np.exp(10 * rh / rJ)
    state['xi_e'] = xi_e

    # dynamical friction timescale:
    t_df = 0.45e3 * (R_gal / 1e3)**2 * v_gal / (Mcl / 1e5) / 2
    state['t_df'] = t_df


def update_cluster(state, config):
    """Update cluster mass, half-mass radius, galactocentric radius, and time.

    Applies stellar-evolution mass loss, relaxation-driven mass loss, adiabatic
    and relaxation-driven expansion, galactocentric inspiral via dynamical
    friction, and advances the simulation time and redshift. Updates state
    in-place.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.

    Returns:
        bool: True if the simulation should continue, False if the cluster has
            dissolved (Mcl < 0), reached the galaxy center (R_gal < 0), or
            the redshift has reached zero.
    """
    
    # unpack current state into local variables:
    t = state['t']; dt = state['dt']; z = state['z']
    Mcl = state['Mcl']; rh = state['rh']; R_gal = state['R_gal']
    m_avg = state['m_avg']; m_avg0 = state['m_avg0']
    t_rlx = state['t_rlx']; t_cc = state['t_cc']
    zCl_form = state['zCl_form']
    xi_e = state['xi_e']
    t_df = state['t_df']

    # average mass evolution:
    if t>t_sev:
        m_avg = m_avg0 * (t / t_sev)**nu_sev
    else:
        m_avg = m_avg0

    # stellar mass loss:
    if t>t_sev:
        dMcl_sev = - nu_sev * Mcl * dt / t
    else:
        dMcl_sev = 0.0

    # relaxation mass loss:
    dMcl_rlx = - xi_e * Mcl * dt / t_rlx

    # total mass loss:
    dMcl = dMcl_sev + dMcl_rlx

    # stellar evolution adiabatic expansion:
    drh_sev = - dMcl_sev * rh / Mcl

    # expansion due to relaxation:
    if t>t_cc:
        drh_rlx = zeta * rh * dt / t_rlx + 2 * dMcl_rlx * rh / Mcl
    else:
        drh_rlx = 0.0

    # total expansion:
    drh = drh_sev + drh_rlx

    # galactocentric radius step:
    dR_gal = - R_gal * dt / t_df

    # cluster mass update:
    Mcl = Mcl + dMcl

    if Mcl<0:
        print('CLUSTER DISSOLVED')
        state['Mcl'] = Mcl
        return False

    # half-mass radius update:
    rh = rh + drh

    # galactocentric radius update:
    R_gal = R_gal + dR_gal

    if R_gal<0:
        print('CLUSTER REACHED GALAXY CENTER')
        state['R_gal'] = R_gal
        return False

    # simulation time update:
    t+=dt

    # redshift update:
    try:
        z = redshift_interp(lookback_interp(zCl_form) - t)
    except Exception:
        print('REDSHIFT 0 REACHED')
        state['t'] = t
        return False

    state['N_iter'] = state['N_iter'] + 1
    state['t'] = t; state['z'] = z; state['dt'] = dt
    state['Mcl'] = Mcl; state['rh'] = rh; state['R_gal'] = R_gal
    state['m_avg'] = m_avg

    return True


def print_status(state, config, local_time_initial, simulation_time_initial):
    """Print a summary of the current iteration to stdout.

    Displays cluster properties (t, dt, z, Mcl, rh, R_gal), BH subsystem
    counts (N_BH, N_BBH, N_Triples, N_me, N_tdeBHWD, N_tdeBHstar), and
    timing information (step duration, cumulative runtime). Only prints if
    config['print_info'] == 1.

    Args:
        state (dict): Mutable simulation state.
        config (dict): Configuration dictionary.
        local_time_initial (float): Wall-clock time at the start of this iteration.
        simulation_time_initial (float): Wall-clock time at the start of the simulation.
    """
    
    # skip printing if disabled:
    if config['print_info']!=1:
        return

    t = state['t']; dt = state['dt']; z = state['z']
    Mcl = state['Mcl']; rh = state['rh']; R_gal = state['R_gal']
    N_BH = state['N_BH']; N_BBH = state['N_BBH']; N_Triples = state['N_Triples']
    N_me = state['N_me']; N_tdeBHWD = state['N_tdeBHWD']; N_tdeBHstar = state['N_tdeBHstar']
    N_iter = state['N_iter']

    local_time_final = time.time()

    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    print("Iteration #", N_iter)
    frmt_1 = '%.2f'
    frmt_2 = '%.0f'
    frmt_3 = '%.3f'
    frmt_4 = "%.1f"
    data_1 = {"t[Myr]": [frmt_1%t], "dt[Myr]": [frmt_1%dt], "z": [frmt_1%z], "Mcl[MMsun]": [frmt_1%(Mcl/1e6)], "rh[pc]": [frmt_1%rh], "R_gal[kpc]": [frmt_1%(R_gal/1e3)]}
    data_2 = {"N_BH": [frmt_2%N_BH], "N_BBH": [frmt_2%N_BBH], "N_Triples": [frmt_2%N_Triples], "N_me": [frmt_2%N_me], "N_tdeBHWD": [frmt_2%N_tdeBHWD], "N_tdeBHstar": [frmt_2%N_tdeBHstar]}
    data_3 = {"steptime[ms]": [frmt_4%(np.abs(local_time_final - local_time_initial)*1e3)], "runtime[s]": [frmt_3%np.abs(time.time() - simulation_time_initial)]}
    headers = [" "]
    df_1 = pd.DataFrame(data_1, headers)
    df_2 = pd.DataFrame(data_2, headers)
    df_3 = pd.DataFrame(data_3, headers)
    print(df_1)
    print(df_2)
    print(df_3)
    print("\n")


def write_output(state, config):
    """Export simulation results to the output directory.

    Writes the following files depending on config indicators:
      - outputBHs.pkl: BH masses, spins, and generations at every timestep.
      - tdes.txt: Tidal disruption event parameters.
      - mergers.txt: BBH merger source parameters (plus initial/final cluster state).
      - evolution.txt: Time-dependent cluster and BH subsystem quantities (69 columns).
      - hardening.txt: BBH hardening track details (12 columns).

    Args:
        state (dict): Simulation state at the end of the run.
        config (dict): Configuration dictionary.
    """
    mergers = np.delete(state['mergers'], 0, axis=0)
    evolution = np.delete(state['evolution'], 0, axis=0)
    hardening = np.delete(state['hardening'], 0, axis=0)
    tdes = np.delete(state['tdes'], 0, axis=0)

    CURRENT_WORKING_DIR = os.getcwd()
    RESULTS_DIR = os.path.join(CURRENT_WORKING_DIR, config['results_folder_name'])
    os.makedirs(RESULTS_DIR, exist_ok=True)

    if config['BOi']==1:
        data_to_save = {
            "t": state['simulation_times'],
            "mBH": state['black_hole_masses'],
            "sBH": state['black_hole_spins'],
            "gBH": state['black_hole_generations'],
            "hBH": state['black_hole_tdes']
        }
        BOF_path = os.path.join(RESULTS_DIR, config['BOF'] + '.pkl')
        with open(BOF_path, "wb") as f:
            pickle.dump(data_to_save, f)

    N_tdeBHWD = state['N_tdeBHWD']
    N_tdeBHstar = state['N_tdeBHstar']
    N_me = state['N_me']
    N_iter = state['N_iter']
    N_hardening = state['N_hardening']
    Mcl0 = state['Mcl0']; rh0 = state['rh0']; Z = state['Z']
    zCl_form = state['zCl_form']; R_gal0 = state['R_gal0']
    Mcl = state['Mcl']; rh = state['rh']; R_gal = state['R_gal']

    if config['Ti']==1:
        tdes_path = os.path.join(RESULTS_DIR, config['tdes_file'] + '.txt')
        with open(tdes_path, 'w') as f_tdes:
            f_tdes.write('# ' + ' '.join(tdes_keys) + '\n')
            for i in range(N_tdeBHWD+N_tdeBHstar):
                f_tdes.write(str(tdes[i][ 0])+' '+str(tdes[i][ 1])+' '+str(tdes[i][ 2])+' '+str(tdes[i][ 3])+' '+\
                             str(tdes[i][ 4])+' '+str(tdes[i][ 5])+' '+str(tdes[i][ 6])+' '+str(tdes[i][ 7])+' '+\
                             str(tdes[i][ 8])+' '+str(tdes[i][ 9])+' '+str(tdes[i][10])+' '+str(tdes[i][11])+' '+\
                             str(tdes[i][12])+' '+str(tdes[i][13])+' '+str(tdes[i][14])+' '+str(tdes[i][15])+' '+\
                             str(tdes[i][16])+' '+str(tdes[i][17]))
                f_tdes.write('\n')

    if config['Mi']==1:
        mergers_path = os.path.join(RESULTS_DIR, config['mergers_file'] + '.txt')
        with open(mergers_path, 'w') as f_mergers:
            f_mergers.write('# ' + ' '.join(merger_keys) + '\n')
            for i in range(N_me):
                f_mergers.write(str(mergers[i][0 ])+' '+str(mergers[i][1 ])+' '+str(mergers[i][2 ])+' '+str(mergers[i][3 ])+' '+str(mergers[i][4 ])+' '+str(mergers[i][5 ])+' '+str(mergers[i][6 ])+' '+\
                                str(mergers[i][7 ])+' '+str(mergers[i][8 ])+' '+str(mergers[i][9 ])+' '+str(mergers[i][10])+' '+str(mergers[i][11])+' '+str(mergers[i][12])+' '+str(mergers[i][13])+' '+\
                                str(mergers[i][14])+' '+str(mergers[i][15])+' '+str(mergers[i][16])+' '+str(mergers[i][17])+' '+str(mergers[i][18])+' '+str(mergers[i][19])+' '+str(mergers[i][20])+' '+\
                                str(mergers[i][21])+' '+str(mergers[i][22])+' '+str(mergers[i][23])+' '+str(mergers[i][24])+' '+str(mergers[i][25])+' '+str(mergers[i][26])+' '+str(Mcl0)+' '+str(rh0)+' '+str(Z)+' '+str(zCl_form)+' '+str(R_gal0)+' '+str(Mcl)+' '+str(rh)+' '+str(R_gal))
                f_mergers.write('\n')

    if config['Ei']==1:
        evolution_path = os.path.join(RESULTS_DIR, config['evolution_file'] + '.txt')
        with open(evolution_path, 'w') as f_evolution:
            f_evolution.write('# ' + ' '.join(evolution_keys) + '\n')
            for i in range(N_iter-1):
                f_evolution.write(str(evolution[i][0 ])+' '+str(evolution[i][1 ])+' '+str(evolution[i][2 ])+' '+str(evolution[i][3 ])+' '+str(evolution[i][4 ])+' '+str(evolution[i][5 ])+' '+\
                                  str(evolution[i][6 ])+' '+str(evolution[i][7 ])+' '+str(evolution[i][8 ])+' '+str(evolution[i][9 ])+' '+str(evolution[i][10])+' '+str(evolution[i][11])+' '+\
                                  str(evolution[i][12])+' '+str(evolution[i][13])+' '+str(evolution[i][14])+' '+str(evolution[i][15])+' '+str(evolution[i][16])+' '+str(evolution[i][17])+' '+\
                                  str(evolution[i][18])+' '+str(evolution[i][19])+' '+str(evolution[i][20])+' '+str(evolution[i][21])+' '+str(evolution[i][22])+' '+str(evolution[i][23])+' '+\
                                  str(evolution[i][24])+' '+str(evolution[i][25])+' '+str(evolution[i][26])+' '+str(evolution[i][27])+' '+str(evolution[i][28])+' '+str(evolution[i][29])+' '+\
                                  str(evolution[i][30])+' '+str(evolution[i][31])+' '+str(evolution[i][32])+' '+str(evolution[i][33])+' '+str(evolution[i][34])+' '+str(evolution[i][35])+' '+\
                                  str(evolution[i][36])+' '+str(evolution[i][37])+' '+str(evolution[i][38])+' '+str(evolution[i][39])+' '+str(evolution[i][40])+' '+str(evolution[i][41])+' '+\
                                  str(evolution[i][42])+' '+str(evolution[i][43])+' '+str(evolution[i][44])+' '+str(evolution[i][45])+' '+str(evolution[i][46])+' '+str(evolution[i][47])+' '+\
                                  str(evolution[i][48])+' '+str(evolution[i][49])+' '+str(evolution[i][50])+' '+str(evolution[i][51])+' '+str(evolution[i][52])+' '+str(evolution[i][53])+' '+\
                                  str(evolution[i][54])+' '+str(evolution[i][55])+' '+str(evolution[i][56])+' '+str(evolution[i][57])+' '+str(evolution[i][58])+' '+str(evolution[i][59])+' '+\
                                  str(evolution[i][60])+' '+str(evolution[i][61])+' '+str(evolution[i][62])+' '+str(evolution[i][63])+' '+str(evolution[i][64])+' '+str(evolution[i][65])+' '+\
                                  str(evolution[i][66])+' '+str(evolution[i][67])+' '+str(evolution[i][68]))
                f_evolution.write('\n')

    if config['Hi']==1:
        hardening_path = os.path.join(RESULTS_DIR, config['hardening_file'] + '.txt')
        with open(hardening_path, 'w') as f_hardening:
            f_hardening.write('# ' + ' '.join(hardening_keys) + '\n')
            for i in range(N_hardening):
                f_hardening.write(str(hardening[i][0 ])+' '+str(hardening[i][1 ])+' '+str(hardening[i][2 ])+' '+str(hardening[i][3 ])+' '+str(hardening[i][4 ])+' '+\
                                  str(hardening[i][5 ])+' '+str(hardening[i][6 ])+' '+str(hardening[i][7 ])+' '+str(hardening[i][8 ])+' '+str(hardening[i][9 ])+' '+\
                                  str(hardening[i][10])+' '+str(hardening[i][11]))
                f_hardening.write('\n')

# End of file.
