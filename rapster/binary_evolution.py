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
import sys
from .constants import *
from .functions import *
from .remnant import *
from .compact_accretion import *  # provides the evolve() spin-accretion routine
from .stellar_evolution import *  # provides the get_star() field-star sampling routine
from .tidal_disruptions import *  # provides the try_BBH_star_disruption() helper
from .fastrandom import choice_weighted
from .bhpool import FAM

# constant factors of T_GW (functions.py): 5 c^5 and G^3
_C5 = 5 * c_light**5
_G3 = G_Newton**3

# constant factor (85 pi / 3 / sqrt(2))^(2/7) of the critical pericenter of GW capture in IMS:
_RP3BM_CONST = (85 * np.pi / 3 / np.sqrt(2))**(2/7)

# indices of the pool weight families of the partner draw:
_I_P15, _I_P25, _I_P17, _I_P27 = FAM['p15'], FAM['p25'], FAM['p17'], FAM['p27']
assert (_I_P25, _I_P17, _I_P27) == (_I_P15 + 1, _I_P15 + 2, _I_P15 + 3)

def _inv(x):
    """1 / x for a Python float x, as NumPy would compute it for a float64: inf (or
    -inf, nan) instead of ZeroDivisionError for x = 0."""
    try:
        return 1 / x
    except ZeroDivisionError:
        return float(np.float64(1) / np.float64(x))

def _pow(x, y):
    """x ** y for Python floats, as NumPy would compute it (inf instead of
    OverflowError)."""
    try:
        return x ** y
    except OverflowError:
        return float(np.float64(x) ** np.float64(y))

def _div(x, y):
    """x / y for Python floats, as NumPy would compute it (inf or nan for y = 0)."""
    try:
        return x / y
    except ZeroDivisionError:
        return float(np.float64(x) / np.float64(y))

def sample_BBH_partner(m12, BH):
    """
    Draw the single BH that interacts with a BBH of total mass m12, with weights
        p3 = (m12 + m) / sqrt(m12**(-2/5) + m**(-2/5)) * m**(3/2)
    over the single-BH masses m (as np.random.choice(mBH, p=p3/np.sum(p3)) did).

    Rejection sampling from the pool BH (bhpool.BHPool): since
    (m12**(-2/5) + m**(-2/5))**(-1/2) <= min(m12, m)**(1/5), p3 is bounded by
        B1 = m12**(1/5) (m12 m^(3/2) + m^(5/2))   or   B2 = m12 m^(17/10) + m^(27/10),
    both mixtures of two pool weight families. The bound with the smaller total is
    used; a candidate is accepted with probability p3 / bound, which is at least
    1/sqrt(2) for B1 when m >= m12 and for B2 when m <= m12, and ~0.4-0.9 otherwise.

    @in m12: total BBH mass m1 + m2 [Msun]
    @in BH: pool of single BHs (not empty)

    @out m3: mass of the sampled single BH [Msun]
    @out k3: its slot in the pool
    """
    M = float(m12)
    BH.prepare(_I_P15, 4)      # (cumulative sums of p15, p25, p17, p27 in one call)
    S15 = BH.total(_I_P15); S25 = BH.total(_I_P25)
    S17 = BH.total(_I_P17); S27 = BH.total(_I_P27)
    M02 = M**0.2
    M_25 = M**(-2/5)
    T1 = M02 * (M * S15 + S25)
    T2 = M * S17 + S27
    m = BH.m
    draw = BH.draw
    # (the candidate's weight and bound are evaluated with Python floats: the same
    # IEEE operations as with NumPy scalars, with less overhead)
    while True:
        if T1 <= T2:
            k = draw(_I_P15) if np.random.rand() * (M * S15 + S25) < M * S15 else draw(_I_P25)
            mk = m[k]
            x = float(mk)
            bound = M02 * (M + x) * x**1.5
        else:
            k = draw(_I_P17) if np.random.rand() * T2 < M * S17 else draw(_I_P27)
            mk = m[k]
            x = float(mk)
            bound = (M + x) * x**1.7
        p3 = (M + x) / math.sqrt(M_25 + x**(-2/5)) * x**(3/2)
        if np.random.rand() * bound < p3:
            return mk, k

def evolve_BBHs(seed, t, z, dt, zCl_form, binaries, hardening, mergers, BH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples, tdes, N_tdeBBHstar, m_min, m_max, f_accreted, EoS):
    """
    @in seed: simulation seed number
    @in t: simulation time
    @in z: redshift
    @in dt: simulation timestep
    @in zCl_form: cluster formation redshift
    @in binaries: array of BBHs
    @in hardening: array of BBH evolution
    @in mergers: array of BBH mergers
    @in BH: pool of single BHs (bhpool.BHPool: masses, spins, generations, TDE counts)
    @in m_avg: average mass
    @in n_star: central stellar density
    @in v_star: 3D star velocity dispersion
    @in vBH: 3D BH velocity dispersion
    @in t_rlx: relaxation timescale
    @in m_avg: average mass
    @in mBH_avg: average BH mass
    @in nc_BH: central BH number density
    @in na_BH: average BH number density
    @in N_BH: number of BHs
    @in N_BBH: number of BBHs
    @in N_me: number of mergeres
    @in N_me2b: number of 2-body mergers
    @in N_3cap: number of 3-captures
    @in N_meFi: number of field mergers
    @in N_meRe: number of retained merger remnants
    @in N_meEj: number of ejected merger remnants
    @in N_dis: number of disruptions
    @in N_ex: number of exchanges
    @in N_BHej: number of single BH ejections
    @in N_BBHej: number of BBH ejections
    @in N_hardening: number of hardening interactions
    @in Vc_BH: BH core volume
    @in N_bb: number of BBH-BBH interactions
    @in triples: array of hierarchical triple-BH systems
    @in N_Triples: number of hierarchical triple-BHs
    @in tdes: tdes array [seed, t, z, type, mstar, Rstar, m_BH, s_BH, g_BH, r_t, r_p, beta, iota, r_mb, dm, s_new, v_rel, h_BH]
    @in N_tdeBBHstar: cumulative number of BBH-star TDEs
    @in m_min: minimum stellar mass
    @in m_max: maximum stellar mass
    @in f_accreted: fraction of disrupted star mass accreted
    @in EoS: NS equation of state; a name registered in compact_accretion.EOS_TABLES

    @out: all inputs
    """
    
    new_triples = []   # triples formed in BBH-BBH interactions (appended at the end)

    # BBH evolution:
    if N_BBH>0:
        
        # shuffle binaries to avoid biases:
        np.random.shuffle(binaries[1:])  # shuffle all BBH rows excluding the placeholder first row

        # mean total BBH mass np.mean(binaries[:, 4]+binaries[:, 5]), recomputed only
        # when the binaries array is replaced or its masses are modified in place:
        mtot_mean = None
        mtot_mean_of = None   # the binaries array mtot_mean was computed for

        # Quantities that stay constant during this call, computed once instead of at
        # every interaction (same operations in the same order, so the values are
        # identical to evaluating them in place):
        v_esc = float(v_esc_cl(v_star, vBH))        # cluster escape velocity
        t_lb_form = float(lookback_interp(zCl_form))  # lookback time at cluster formation
        # interaction rates Rate_int(m, n, v, rp) = n * (2 pi G m rp / v_rel^2) * v_rel,
        # v_rel = sqrt(2) v, written out with the v-dependent factors precomputed:
        two_pi_G = 2 * np.pi * G_Newton
        vrel_star = float(np.sqrt(2) * np.sqrt(v_star**2 + vBH**2))   # (star-BBH relative velocity)
        vrel_star2 = vrel_star**2
        vrel_BH = float(np.sqrt(2) * vBH)
        vrel_BH2 = vrel_BH**2
        # leading factors of the Maxwellian dispersions sqrt(m_avg v_star^2 / 3 / m) and
        # sqrt(mBH_avg vBH^2 / 3 / m):
        disp_star = float(m_avg * v_star**2 / 3)
        disp_BH = float(mBH_avg * vBH**2 / 3)

        # The interactions below are computed with Python floats: the binary's row is
        # read with tolist(), and the inputs used in the loop are converted once here.
        # Python floats and NumPy scalars give the same IEEE results for +, -, *, /,
        # ** and sqrt (math.sqrt is correctly rounded, as np.sqrt), at a fraction of the
        # overhead; the NumPy functions (cos, arccos, ...) are still called where they
        # were. (The function returns its inputs unchanged.) Square roots of arguments
        # that rounding could make slightly negative keep np.sqrt (NaN, not an error).
        n_star_f = float(n_star); nc_BH_f = float(nc_BH); Vc_BH_f = float(Vc_BH)
        m_avg_f = float(m_avg); mBH_avg_f = float(mBH_avg)
        t_f = float(t); dt_f = float(dt); t_rlx_f = float(t_rlx); vBH_f = float(vBH)
        rand = np.random.rand
        uniform = np.random.uniform
        sqrt = math.sqrt
        h_append = hardening.data.extend   # (= hardening.append: rows of 12 values)
        if hardening.ncols != 12:
            h_append = hardening.append
        h_set_last = hardening.set_last
        two_vBH = 2 * vBH_f
        # (degenerate inputs keep NumPy's inf/nan semantics where they are divisors:)
        if not vrel_star2 > 0:
            vrel_star = np.float64(vrel_star); vrel_star2 = np.float64(vrel_star2)
        if not vrel_BH2 > 0:
            vrel_BH = np.float64(vrel_BH); vrel_BH2 = np.float64(vrel_BH2)
        if not Vc_BH_f > 0:
            Vc_BH_f = np.float64(Vc_BH_f)

        # initialize iteration index:
        i=1
        
        # iterate over all BBHs:
        while i<N_BBH+1:
            
            # binary status:
            condition=0
            
            # initialize local timescale:
            t_local=0
            
            # binary convection timescale:
            t_conv = 0
            
            # while binary is available, evolve it:
            while N_BBH>0:
                
                if i==0:
                    i=1
                
                # unwrap parameters of current binary:
                ind, channel, a, e, m1, m2, s1, s2, g1, g2, t_form, z_form, Nex, h1, h2 = binaries[i].tolist()
                q = (m2 / m1) if m2 < m1 else (m1 / m2)    # = min_of(m1, m2)/max_of(m1, m2)
                
                # BBH-star interaction timescale (encounter cross-section set by the binary separation a);
                # = 1 / Rate_int(m1 + m2 + m_avg, n_star, sqrt(v_star**2 + vBH**2), kp_max * a):
                if n_star_f>0:
                    t_BBH_star = _inv(n_star_f * (two_pi_G * (m1 + m2 + m_avg_f) * (kp_max * a) / vrel_star2) * vrel_star)
                else:
                    t_BBH_star = 1e100

                # BBH-BH interaction timescale, = 1 / Rate_int(m1+m2+mBH_avg, nc_BH, vBH, kp_max * a) + t_conv:
                if BH.n>0:
                    t_BBH_BH = _inv(nc_BH_f * (two_pi_G * (m1+m2+mBH_avg_f) * (kp_max * a) / vrel_BH2) * vrel_BH) + t_conv
                else:
                    t_BBH_BH = 1e100
                
                # BBH-BBH interaction timescale:
                if N_BBH>1:
                    # BBH number density:
                    n_BBH = (N_BBH - 1) / Vc_BH_f
                    
                    if mtot_mean_of is not binaries:
                        mtot_tot = binaries[:, 4]+binaries[:, 5]
                        mtot_mean = float(np.add.reduce(mtot_tot) / mtot_tot.shape[0])   # = np.mean(mtot_tot)
                        mtot_mean_of = binaries
                    t_BBH_BBH = _inv(n_BBH * (two_pi_G * (m1 + m2 + mtot_mean) * (kp_max * a) / vrel_BH2) * vrel_BH) + t_conv
                else:
                    t_BBH_BBH = 1e100
                    
                # competing rates for the three possible encounter types:
                rate_star = _inv(t_BBH_star)
                rate_BH = _inv(t_BBH_BH)
                rate_BBH = _inv(t_BBH_BBH)
                rate_tot = rate_star + rate_BH + rate_BBH

                # probabilities proportional to each channel's rate:
                p_BBH_star = _div(rate_star, rate_tot)
                p_BBH_BH = _div(rate_BH, rate_tot)
                p_BBH_BBH = _div(rate_BBH, rate_tot)

                u_int = rand()

                if u_int < p_BBH_star: # do BBH-star interaction
                    dt_local = t_BBH_star
                    type_int = 1
                elif u_int < p_BBH_star + p_BBH_BH: # do BBH-BH interaction
                    dt_local = t_BBH_BH
                    type_int = 2
                else: # do BBH-BBH interaction
                    dt_local = t_BBH_BBH
                    type_int = 3
                    
                h_append([t, dt, t_local, dt_local, ind, a, e, m1, m2, q, condition, Nex])
                N_hardening+=1
                
                if condition>0:
                    break
                
                # GW coalescence time T_GW(m1, m2, a, e), with the expressions of T_GW for
                # ordinary values written out here (same operations, called ~1e6 times):
                if 0.0 < m1 and 0.0 < m2 and 0.0 <= e < 1.0 and 0.0 <= a < 1e60:
                    Tc_ = _C5 * a**4 / 256 / _G3 / m1 / m2 / (m1 + m2)
                    T_gw = Tc_ * ((1 + 0.27 * e**10 + 0.33 * e**20 + 0.2 * e**1000) * (1 - e**2)**(7/2))
                else:
                    T_gw = T_GW(m1, m2, a, e)
                if T_gw < min_of(dt_local, dt_f, t_lb_form - t_f - t_local): # 2-body in-cluster merger during the current time-step
                    
                    theta1, theta2, dPhi = sample_angles()
                    
                    m_rem, s_rem, vGW_kick = merger_remnant(m1, m2, s1, s2, theta1, theta2, dPhi)
                    g_rem = np.max([g1, g2]) + 1
                    h_rem = h1 + h2

                    if vGW_kick < v_esc: # merger remnant retained in cluster
                        
                        BH.add(m_rem, s_rem, g_rem, h_rem)
                        
                        N_BH = N_BH - 1
                        
                        N_meRe+=1
                        
                    else: # merger remnant ejected
                        
                        N_BH = N_BH - 2
                        
                        N_meEj+=1
                        
                    N_me+=1
                    N_me2b+=1
                    
                    # order BHs by mass:
                    mA = m1; sA = s1; gA = g1; hA = h1; thetaA = theta1
                    mB = m2; sB = s2; gB = g2; hB = h2; thetaB = theta2
                    if mA>mB:
                        m1 = mA; s1 = sA; g1 = gA; h1 = hA; theta1 = thetaA
                        m2 = mB; s2 = sB; g2 = gB; h2 = hB; theta2 = thetaB
                    else:
                        m1 = mB; s1 = sB; g1 = gB; h1 = hB; theta1 = thetaB
                        m2 = mA; s2 = sA; g2 = gA; h2 = hA; theta2 = thetaA
                        
                    # mass ratio:
                    q = m2 / m1
                    
                    # effective spin parameter:
                    s_eff = (m1 * s1 * np.cos(theta1) + m2 * s2 * np.cos(theta2)) / (m1 + m2)
                    
                    # merger time and redshift:
                    t_merge = t_f + t_local + T_GW(m1, m2, a, e)     # absolute merger time [Myr]
                    t_lb = t_lb_form - t_merge     # lookback time remaining at merger
                    z_merge = redshift_interp(t_lb) if t_lb >= 0.1 else 0.0     # merger redshift (0 if within 0.1 Myr of today)

                    # append merger:
                    mergers.append([seed, ind, channel, a, e, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form,
                                    t_merge, z_merge, m_rem, s_rem, g_rem, vGW_kick, s_eff, q, v_esc, h1, h2])
                    
                    binaries = delete_row_inplace(binaries, i)
                    
                    i = i - 1
                    
                    N_BBH = N_BBH - 1
                    
                    condition=2
                    h_set_last(10, condition)  # the row appended for this interaction
                    break
                
                # probability for interaction:
                p_int = _div(dt_f, dt_f + dt_local)
                
                # update local time:
                t_local+=dt_local
                
                random_number = rand()

                if random_number > p_int:
                    # interaction will not happen
                    condition = 1
                    t_local = dt
                    break
                else:
                    # interaction will happen
                    if t_local >= dt_f:
                        # local time exceeds the parent iteration time-step
                        t_local = dt
                        condition=1
                    else:
                        # local time does not exceed the parent iteration time-step
                        pass
                    
                if type_int==3: # BBH-BBH interaction occurs
                    
                    N_bb+=1
                    
                    # semimajor axis of current BBH:
                    a1 = a

                    # sample the other BBH with probability ~ (m1 + m2) a, among the BBHs
                    # other than the current one (row i, weight set to 0; the placeholder
                    # row 0 has weight 0), located directly in the cumulative weights:
                    w_bb = (binaries[:, 4] + binaries[:, 5]) * binaries[:, 2]
                    w_bb[i] = 0.0
                    cw_bb = np.add.accumulate(w_bb)   # (= np.cumsum(w_bb), without its Python wrapper)
                    total_bb = cw_bb[-1]
                    if not total_bb > 0:
                        continue
                    k2 = int(cw_bb.searchsorted(rand() * total_bb, side='right'))
                    while k2 >= cw_bb.shape[0]:      # (rounding at the end: draw again)
                        k2 = int(cw_bb.searchsorted(rand() * total_bb, side='right'))
                    a2 = binaries.item(k2, 2)

                    # index of current BBH:
                    k1 = i
                        
                    # order based on hardness:
                    if a1 < a2:
                        k_soft = k2; a_soft = a2
                        k_hard = k1; a_hard = a1
                    else:
                        k_soft = k1; a_soft = a1
                        k_hard = k2; a_hard = a2
                        
                    # BHs of tighter BBH:
                    rh_ = binaries[k_hard].tolist()
                    m0 = rh_[4]; s0 = rh_[6]; g0 = rh_[8]; h0 = rh_[13]
                    m1 = rh_[5]; s1 = rh_[7]; g1 = rh_[9]; h1 = rh_[14]
                
                    # BHs of softer BBH:
                    rs_ = binaries[k_soft].tolist()
                    if rs_[4] < rs_[5]:
                        m_freed = rs_[4]; s_freed = rs_[6]; g_freed = rs_[8]; h_freed = rs_[13]
                        m2 = rs_[5]; s2 = rs_[7]; g2 = rs_[9]; h2 = rs_[14]
                    else:
                        m_freed = rs_[5]; s_freed = rs_[7]; g_freed = rs_[9]; h_freed = rs_[14]
                        m2 = rs_[4]; s2 = rs_[6]; g2 = rs_[8]; h2 = rs_[13]
                    
                    BH.add(m_freed, s_freed, g_freed, h_freed)
                    
                    # delete softer BBH:
                    binaries = delete_row_inplace(binaries, k_soft)
                    
                    N_dis+=1
                    
                    N_BBH = N_BBH - 1
                    
                    if k_hard > k_soft:
                        k_hard = k_hard - 1
                        
                    a_inner = a_hard
                    a_outer = a_soft * (m0 + m1) / m_freed
                    e_inner = sqrt(rand())
                    e_outer = sqrt(rand())
                    inclination1 = np.arccos(uniform(-1, 1))
                    inclination2 = np.arccos(uniform(-1, 1))
                    inclination = inclination1 + inclination2
                    
                    i_old = i
                    
                    # check Mardling & Aarseth (2001) criterion:
                    if a_outer/a_inner > 2.8 * (1 + m2 / (m1 + m0))**(2/5) * (1 + e_outer)**(2/5) / (1 - e_outer)**(6/5) * (1 - 0.3 * inclination / np.pi): # hierarchical triple
                        
                        # append triple:
                        # (appended to triples after the loop, which does not read it: the
                        # same rows as appending them one by one)
                        new_triples.append([a_inner, a_outer, e_inner, e_outer, m0, m1, m2, s0, s1, s2, g0, g1, g2, inclination1, inclination2, t, z,
                                            rh_[0], rh_[1], rh_[10], rh_[11], h0, h1, h2])
                        
                        N_Triples+=1
                        
                        binaries = delete_row_inplace(binaries, k_hard)
                        
                        N_BBH = N_BBH - 1
                        
                        i = i - 1
                        
                    else: # `inner' binary is freed:
                        
                        BH.add(m2, s2, g2, h2)
                        
                        # append binary:
                        Delta = 0.38 # Zevin et al., ApJ 871 (2019), 91.
                        binaries[k_hard, 2] = a_inner/(1 + Delta * m_freed * m2 / m1 / m0 * a_inner / a_outer)
                        binaries[k_hard, 3] = e_inner
                        
                    if k_soft < i_old:
                        i = i - 1
                        continue
                    elif k_soft == i_old:
                        i = i - 1
                        condition=4
                        h_set_last(10, condition)  # the row appended for this interaction
                        break
                    else:
                        continue
                    
                # at this point BBH-single interaction occurs (signle= star or BH)

                if type_int==1: # BBH-star interaction occurs: a field star is gravitationally
                                 # focused toward the binary (encounter cross-section ~ a).
                                 # Disruption is checked later, once the pericenter rp is
                                 # sampled (shared with the resonance check below), so this
                                 # block only samples the star and sets up m3 for the
                                 # generic 3-body dynamics that follow if no disruption occurs.

                    # sample a field star (mass and radius) from the current stellar mass function:
                    m_star, R_star = get_star(t, tBH_form, m_min, m_max)

                    # treat the star as a generic 3rd body, using its actual sampled mass
                    # rather than the population average:
                    m3 = m_star

                    # single velocity before interaction, sampled from a Maxwellian whose
                    # dispersion reflects energy equipartition with the star's own mass,
                    # anchored to the stellar population's characteristic energy scale
                    # (m_avg * v_star^2), mirroring the analogous BH formula below:
                    vS_before = float(get_maxwell_sample(sqrt(disp_star / m3)))

                else: # BBH-BH occurs

                    if BH.size == 0:  # no single BHs left to interact with, skip
                        continue
                    # sample the single BH, with weight
                    # p3 = (m1 + m2 + m) / sqrt((m1 + m2)**(-2/5) + m**(-2/5)) * m**(3/2),
                    # and its slot k3 in the pool:
                    m3, k3 = sample_BBH_partner(m1 + m2, BH)
                    m3, s3, g3, h3 = BH.get(k3)        # (Python floats; m3 as drawn)
                    
                    # sample single velocity before interaction:
                    vS_before = float(get_maxwell_sample(sqrt(disp_BH / m3)))
                    
                # sample binary velocity before interaction:
                vB_before = float(get_maxwell_sample(sqrt(disp_BH / (m1 + m2))))

                # sample cosine angle before interaction:
                cos_theta_before = uniform(-1, 1)
                
                # relative velocity before interaction:
                v2_ = vB_before**2 + vS_before**2 - 2 * vB_before * vS_before * cos_theta_before
                v_rel_before = sqrt(v2_) if v2_ >= 0.0 else float(np.sqrt(v2_))
                
                # reduced mass before interaction:
                mu_before = (m1 + m2) * m3 / (m1 + m2 + m3)
                
                if 1/2 * mu_before * v_rel_before**2 > G_Newton * m1 * m2 / 2 / a: # binary is ionized
                    
                    N_dis+=1
                    
                    binaries = delete_row_inplace(binaries, i)
                    
                    i = i - 1
                    
                    BH.add(m1, s1, g1, h1)
                    BH.add(m2, s2, g2, h2)
                    
                    N_BBH = N_BBH - 1
                    
                    condition=3
                    h_set_last(10, condition)  # the row appended for this interaction
                    break

                # sample pericenter of interaction (shared between TDE and dynamical checks,
                # since both refer to the same physical passage); uniform in rp itself,
                # consistent with the linear (gravitational-focusing) cross-section scaling:
                rp = uniform(0, kp_max * a)

                # critical pericenter for resonant interaction:
                rp_c = max_of(m1, m2) / (m1 + m2) * a

                if type_int==1: # check for tidal disruption on this passage
                    disrupted, m1, m2, s1, s2, h1, h2, tdes = try_BBH_star_disruption(
                        rp, a, m1, m2, s1, s2, g1, g2, h1, h2, m_star, R_star, f_accreted, EoS,
                        seed, t, z, tdes, m_avg, vBH, v_star)

                    if disrupted:
                        N_tdeBBHstar+=1
                        binaries[i, 4] = m1; binaries[i, 6] = s1; binaries[i, 13] = h1
                        binaries[i, 5] = m2; binaries[i, 7] = s2; binaries[i, 14] = h2
                        mtot_mean_of = None  # masses changed
                        continue # star consumed; move to next interaction draw for this binary

                    if rp < rp_c: # resonant encounter: multiple close passages possible
                        for _ in range(N_IMS - 1): # passage 1 already checked above
                            rp = uniform(0, kp_max * a)
                            disrupted, m1, m2, s1, s2, h1, h2, tdes = try_BBH_star_disruption(
                                rp, a, m1, m2, s1, s2, g1, g2, h1, h2, m_star, R_star, f_accreted, EoS,
                                seed, t, z, tdes, m_avg, vBH, v_star)
                            if disrupted:
                                N_tdeBBHstar+=1
                                binaries[i, 4] = m1; binaries[i, 6] = s1; binaries[i, 13] = h1
                                binaries[i, 5] = m2; binaries[i, 7] = s2; binaries[i, 14] = h2
                                mtot_mean_of = None  # masses changed
                                break
                        if disrupted:
                            continue # star consumed during one of the resonant passages

                        # no disruption across all N_IMS passages: fall through to ordinary
                        # hardening below (no capture/merger physics -- stars never merge
                        # with BHs, enforced by the type_int==2 gates further down)

                if rp < rp_c and type_int==2: # interaction is resonant
                    
                    u_3bm = np.random.rand(3).tolist()

                    # For each of the three pairs (m1, m2), (m2, m3), (m3, m1) of the
                    # intermediate states (IMS): semimajor axis assuming energy conservation,
                    # critical pericenter for a GW capture merger, critical eccentricity and
                    # probability of a 3-body merger over N_IMS IMS; the pair with the
                    # largest ratio p_3bm / u_3bm is selected (written with Python floats,
                    # pair by pair, instead of three-element arrays):
                    mm12 = m1 * m2
                    rp_pre = _RP3BM_CONST * a**(2/7) / 2
                    mm12_25 = mm12**(2/5)
                    pair_m = ((m1, m2), (m2, m3), (m3, m1))
                    j_3bm = 0
                    ratio_max = None
                    a_3bm = [0.0, 0.0, 0.0]
                    e_3bm_c = [0.0, 0.0, 0.0]
                    for j in range(3):
                        x, y = pair_m[j]
                        a_3bm[j] = a / mm12 * (x * y)
                        rp_c3 = rp_pre * (2 * G_Newton * (x * y)**(4/5) * (x + y)**(1/5) / mm12_25 / c_light**2)**(5/7)
                        e_3bm_c[j] = 1 - rp_c3 / a_3bm[j]
                        p_3bm = 1 - _pow(e_3bm_c[j], 2 * N_IMS)
                        u = u_3bm[j]
                        ratio = p_3bm / u if u > 0.0 else (math.inf if p_3bm > 0 else -math.inf)
                        if ratio_max is None or ratio > ratio_max:
                            ratio_max = ratio
                            j_3bm = j

                    # eccentricities of the IMS (drawn for all three as before; NaN, as before,
                    # if the critical eccentricity is negative and the draw exceeds 1/(1 - e_c)):
                    v_3bm = np.random.rand(3)

                    m_3bm_1 = (m1, m2, m3); s_3bm_1 = (s1, s2, s3); g_3bm_1 = (g1, g2, g3); h_3bm_1 = (h1, h2, h3)
                    m_3bm_2 = (m2, m3, m1); s_3bm_2 = (s2, s3, s1); g_3bm_2 = (g2, g3, g1); h_3bm_2 = (h2, h3, h1)
                    m_3bm_3 = (m3, m1, m2); s_3bm_3 = (s3, s1, s2); g_3bm_3 = (g3, g1, g2); h_3bm_3 = (h3, h1, h2)
                    
                    if ratio_max > 1: # binary-single GW capture merger occurs
                        
                        BH.remove(k3)
                        
                        mA = m_3bm_1[j_3bm]; sA = s_3bm_1[j_3bm]; gA = g_3bm_1[j_3bm]; hA = h_3bm_1[j_3bm]
                        mB = m_3bm_2[j_3bm]; sB = s_3bm_2[j_3bm]; gB = g_3bm_2[j_3bm]; hB = h_3bm_2[j_3bm]
                        mC = m_3bm_3[j_3bm]; sC = s_3bm_3[j_3bm]; gC = g_3bm_3[j_3bm]; hC = h_3bm_3[j_3bm]
                            
                        thetaA, thetaB, dPhi = sample_angles()
                            
                        m_rem, s_rem, vGW_kick = merger_remnant(mA, mB, sA, sB, thetaA, thetaB, dPhi)
                        g_rem = np.max([gA, gB]) + 1
                        h_rem = h1 + h2
                        
                        sma = a_3bm[j_3bm]
                        eccen = np.sqrt(1 - (1 - e_3bm_c[j_3bm]) * v_3bm[j_3bm])
                        
                        if vGW_kick < v_esc: # merger remnant retained in cluster
                            
                            BH.add(m_rem, s_rem, g_rem, h_rem)
                            
                            N_BH = N_BH - 1
                            N_meRe+=1
                            
                        else:
                            
                            N_BH = N_BH - 2
                            N_meEj+=1
                            
                        N_3cap+=1
                        N_me+=1
                        
                        # order BHs my mass:
                        if mA>mB:
                            m1 = mA; s1 = sA; g1 = gA; h1 = hA; theta1 = thetaA
                            m2 = mB; s2 = sB; g2 = gB; h2 = hB; theta2 = thetaB
                        else:
                            m1 = mB; s1 = sB; g1 = gB; h1 = hB; theta1 = thetaB
                            m2 = mA; s2 = sA; g2 = gA; h2 = hA; theta2 = thetaA
                            
                        # mass ratio:
                        q = m2 / m1
                        
                        # effective spin parameter:
                        s_eff = (m1 * s1 * np.cos(theta1) + m2 * s2 * np.cos(theta2)) / (m1 + m2)
                        
                        # merger time and redshift:
                        t_merge = t_f + t_local + T_GW(mA, mB, sma, eccen)  # absolute merger time [Myr]
                        t_lb = t_lb_form - t_merge          # lookback time remaining at merger
                        z_merge = redshift_interp(t_lb) if t_lb >= 0.1 else 0.0  # merger redshift (0 if within 0.1 Myr of today)

                        # append merger:
                        mergers.append([seed, ind, 6, sma, eccen, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form, 
                                        t_merge, z_merge, m_rem, s_rem, g_rem, vGW_kick, s_eff, q, v_esc, h1, h2])
                        
                        # delete BBH:
                        binaries = delete_row_inplace(binaries, i)
                        
                        N_BBH = N_BBH - 1
                        
                        i = i - 1
                        
                        BH.add(mC, sC, gC, hC)
                        
                        condition=6
                        h_set_last(10, condition)  # the row appended for this interaction
                        break
                        
                if m3>min_of(m1, m2) and rp<rp_c and type_int==2: # exchange occurs

                    N_ex+=1
                    binaries[i, 12]+=1
                    
                    if m1<m2:

                        ms = m1; ss = s1; gs = g1; hs = h1
                        mr = m2; sr = s2; gr = g2; hr = h2
                        
                    else:
                        
                        ms = m2; ss = s2; gs = g2; hs = h2
                        mr = m1; sr = s1; gr = g1; hr = h1
                        
                    # impose binding energy conservation after substitution:
                    a = m3/ms * a                  # (= binaries[i][2], which is a)
                    binaries[i, 2] = a

                    # the binary member replaced by the single BH takes its place in the pool:
                    BH.replace(k3, ms, ss, gs, hs)

                    # update binary masses:
                    binaries[i, 4] = mr
                    binaries[i, 5] = m3
                    mtot_mean_of = None  # masses changed
                    
                    # update binary spins:
                    binaries[i, 6] = sr
                    binaries[i, 7] = s3
                    
                    # update binary generations:
                    binaries[i, 8] = gr
                    binaries[i, 9] = g3
                    
                    # update binary tde counters:
                    binaries[i, 13] = hr
                    binaries[i, 14] = h3
                    
                    m1 = mr; s1 = sr; g1 = gr; h1 = hr
                    m2 = m3; s2 = s3; g2 = g3; h2 = h3
                    m3 = ms; s3 = ss; g3 = gs; h3 = hs
                    
                    # (e = binaries[i][3], unchanged)
                    
                # Binary hardening:
                a_new = a / (1 + Hardening_constant * m3 / (m1 + m2))
                e = sqrt(rand())
                binaries[i, 2] = a_new
                binaries[i, 3] = e
                
                # energy extracted:
                dE_b = G_Newton * m1 * m2 / 2 * (1/a_new - 1/a)
                
                a = a_new
                
                # reduced mass after interaction:
                mu_after = (m1 + m2) * m3 / (m1 + m2 + m3)
                
                # relative velocity after interaction:
                v2_ = mu_before / mu_after * v_rel_before**2 + 2 / mu_after * dE_b
                v_rel_after = sqrt(v2_) if v2_ >= 0.0 else float(np.sqrt(v2_))
                
                # velocity of single after interaction:
                v3_after = (m1 + m2) / (m1 + m2 + m3) * v_rel_after
                
                # velocity of binary after interaction:
                v12_after = m3 / (m1 + m2 + m3) * v_rel_after
                
                # check if single is ejected:
                if v3_after > v_esc and type_int==2:
                    
                    BH.remove(k3)

                    N_BH = N_BH - 1
                    
                    N_BHej+=1
                    
                if v12_after > two_vBH \
                   and v12_after < v_esc: # binary convection
                    
                    t_conv = m_avg_f / (m1 + m2) * t_rlx_f
                    
                # check if binary is ejected:
                if v12_after > v_esc:
                    
                    # check if BBH mergers in the field:
                    if t_f + t_local + T_GW(m1, m2, a, e) < t_lb_form: # BBH merges
                        
                        theta1, theta2, dPhi = sample_angles()
                        
                        m_rem, s_rem, vGW_kick = merger_remnant(m1, m2, s1, s2, theta1, theta2, dPhi)
                        g_rem = np.max([g1, g2]) + 1
                        h_rem = h1 + h2

                        N_me+=1
                        N_meFi+=1
                        
                        # merger time:
                        t_merge = t_f + t_local + T_GW(m1, m2, a, e)
                        
                        # order BHs by mass:
                        mA = m1; sA = s1; gA = g1; hA = h1; thetaA = theta1
                        mB = m2; sB = s2; gB = g2; hB = h2; thetaB = theta2
                        if mA>mB:
                            m1 = mA; s1 = sA; g1 = gA; h1 = hA; theta1 = thetaA
                            m2 = mB; s2 = sB; g2 = gB; h2 = hB; theta2 = thetaB
                        else:
                            m1 = mB; s1 = sB; g1 = gB; h1 = hB; theta1 = thetaB
                            m2 = mA; s2 = sA; g2 = gA; h2 = hA; theta2 = thetaA
                            
                        # mass ratio:
                        q = m2 / m1
                        
                        # effective spin parameter:
                        s_eff = (m1 * s1 * np.cos(theta1) + m2 * s2 * np.cos(theta2)) / (m1 + m2)
                        
                        # lookback time and redshift at field merger:
                        t_lb = t_lb_form - t_merge  # lookback time remaining at merger
                        z_merge = redshift_interp(t_lb) if t_lb >= 0.1 else 0.0  # merger redshift (0 if within 0.1 Myr of today)

                        # append merger:
                        mergers.append([seed, ind, -channel, a, e, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form, t_merge,
                                        z_merge, m_rem, s_rem, g_rem, vGW_kick, s_eff, q, v_esc, h1, h2])
                        
                    binaries = delete_row_inplace(binaries, i)
                    
                    i = i - 1
                    
                    N_BBH = N_BBH - 1
                    
                    N_BH = N_BH - 2
                    
                    N_BBHej+=1
                    
                    condition=5
                    h_set_last(10, condition)  # the row appended for this interaction
                    break

            h_append([t, dt, t_local, dt_local, ind, a, e, m1, m2, q, condition, Nex])
            N_hardening+=1
            i+=1
            
    if new_triples:
        triples = np.concatenate((triples, np.array(new_triples)), axis=0)

    return seed, t, z, dt, zCl_form, binaries, hardening, mergers, BH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples, tdes, N_tdeBBHstar, m_min, m_max, f_accreted, EoS

# End of file.
