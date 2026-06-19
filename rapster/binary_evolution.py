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

inport sys
from .constants import *
from .functions import *
from .remnant import *
from .compact_accretion import *  # provides the evolve() spin-accretion routine
from .stellar_evolution import *  # provides the get_star() field-star sampling routine

def evolve_BBHs(seed, t, z, dt, zCl_form, binaries, hardening, mergers, mBH, sBH, gBH, hBH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples, tdes, N_tdeBBHstar, m_min, m_max, f_accreted, EoS):
    """
    @in seed: simulation seed number
    @in t: simulation time
    @in z: redshift
    @in dt: simulation timestep
    @in zCl_form: cluster formation redshift
    @in binaries: array of BBHs
    @in hardening: array of BBH evolution
    @in mergers: array of BBH mergers
    @in mBH: array of single BH masses
    @in sBH: array of single BH spins
    @in gBH: array of single BH generations
    @in hBH: array of BH tdes count
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
    @in EoS: equation of state for neutron stars; either 'APR' or 'AU'

    @out: all inputs
    """
    
    # BBH evolution:
    if N_BBH>0:
        
        # shuffle binaries to avoid biases:
        np.random.shuffle(binaries[1:])  # shuffle all BBH rows excluding the placeholder first row
        
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
                ind = binaries[i][0]
                channel = binaries[i][1]
                a = binaries[i][2]
                e = binaries[i][3]
                m1 = binaries[i][4]
                m2 = binaries[i][5]
                s1 = binaries[i][6]
                s2 = binaries[i][7]
                g1 = binaries[i][8]
                g2 = binaries[i][9]
                q = np.min([m1, m2])/np.max([m1, m2])
                t_form = binaries[i][10]
                z_form = binaries[i][11]
                Nex = binaries[i][12]
                h1 = binaries[i][13]
                h2 = binaries[i][14]
                
                # BBH-star interaction timescale (encounter cross-section set by the binary separation a):
                if n_star>0:
                    v_rel = np.sqrt(v_star**2 + vBH**2)  # relative velocity between star and BBH center of mass
                    t_BBH_star = 1 / Rate_int(m1 + m2 + m_avg, n_star, v_rel, kp_max * a)
                else:
                    t_BBH_star = 1e100

                # BBH-BH interaction timescale:
                if mBH.size>0:
                    t_BBH_BH = 1 / Rate_int(m1+m2+mBH_avg, nc_BH, vBH, kp_max * a) + t_conv
                else:
                    t_BBH_BH = 1e100
                
                # BBH-BBH interaction timescale:
                if N_BBH>1:
                    # BBH number density:
                    n_BBH = (N_BBH - 1) / Vc_BH
                    
                    t_BBH_BBH = 1 / Rate_int(m1 + m2 + np.mean(binaries[:, 4]+binaries[:, 5]), n_BBH, vBH, kp_max * a) + t_conv
                else:
                    t_BBH_BBH = 1e100
                    
                # competing rates for the three possible encounter types:
                rate_star = 1 / t_BBH_star
                rate_BH = 1 / t_BBH_BH
                rate_BBH = 1 / t_BBH_BBH
                rate_tot = rate_star + rate_BH + rate_BBH

                # probabilities proportional to each channel's rate:
                p_BBH_star = rate_star / rate_tot
                p_BBH_BH = rate_BH / rate_tot
                p_BBH_BBH = rate_BBH / rate_tot

                u_int = np.random.rand()

                if u_int < p_BBH_star: # do BBH-star interaction
                    dt_local = t_BBH_star
                    type_int = 1
                elif u_int < p_BBH_star + p_BBH_BH: # do BBH-BH interaction
                    dt_local = t_BBH_BH
                    type_int = 2
                else: # do BBH-BBH interaction
                    dt_local = t_BBH_BBH
                    type_int = 3
                    
                hardening = np.append(hardening, [[t, dt, t_local, dt_local, ind, a, e, m1, m2, q, condition, Nex]], axis=0)
                N_hardening+=1
                
                if condition>0:
                    break
                
                if T_GW(m1, m2, a, e) < np.min([dt_local, dt, lookback_interp(zCl_form) - t - t_local]): # 2-body in-cluster merger during the current time-step
                    
                    theta1, theta2, dPhi = sample_angles()
                    
                    m_rem, s_rem, vGW_kick = merger_remnant(m1, m2, s1, s2, theta1, theta2, dPhi)
                    g_rem = np.max([g1, g2]) + 1
                    h_rem = h1 + h2

                    if vGW_kick < 2 * np.sqrt(v_star**2 + vBH**2): # merger remnant retained in cluster
                        
                        mBH = np.append(mBH, m_rem)
                        sBH = np.append(sBH, s_rem)
                        gBH = np.append(gBH, g_rem)
                        hBH = np.append(hBH, h_rem)
                        
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
                    
                    # append merger:
                    mergers = np.append(mergers, [[seed, ind, channel, a, e, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form,
                                                   t + t_local + T_GW(m1, m2, a, e), redshift_interp(lookback_interp(zCl_form) - t - t_local - T_GW(m1, m2, a, e)), m_rem, s_rem, g_rem, vGW_kick, s_eff, q, 2*v_star, h1, h2]], axis=0)
                    
                    binaries = np.delete(binaries, i, axis=0)
                    
                    i = i - 1
                    
                    N_BBH = N_BBH - 1
                    
                    condition=2
                    hardening[i][10]=condition
                    break
                
                # probability for interaction:
                p_int = dt / (dt + dt_local)
                
                # update local time:
                t_local+=dt_local
                
                random_number = np.random.rand()

                if random_number > p_int:
                    # interaction will not happen
                    condition = 1
                    t_local = dt
                    break
                else:
                    # interaction will happen
                    if t_local >= dt:
                        # local time exceeds the parent iteration time-step
                        t_local = dt
                        condition=1
                    else:
                        # local time does not exceed the parent iteration time-step
                        pass
                    
                if type_int==1: # BBH-star interaction occurs: a field star is gravitationally
                                 # focused toward the binary (encounter cross-section ~ a),
                                 # and may or may not end up tidally disrupted

                    # sample a field star (mass and radius) from the current stellar mass function:
                    m_star, R_star = get_star(t, tBH_form, m_min, m_max)

                    # whole-binary (monopole) tidal radius: treats the BBH as a single point mass
                    # of total mass (m1+m2), valid when the star's pericenter stays well outside
                    # the binary's own separation a (see physics discussion: monopole vs quadrupole):
                    r_t_bin = R_star * ((m1+m2)/m_star)**(1/3)

                    # pericenter of the star's encounter with the binary, sampled uniformly in
                    # impact-parameter-squared out to the binary's own encounter cross-section a
                    # (NOT the tidal radius -- the star is gravitationally focused on the scale
                    # of the binary's separation, since a >> r_t_bin is typical for NSC BBHs):
                    r_p = np.sqrt(np.random.rand()) * a

                    if r_p > r_t_bin:
                        # star's pericenter is outside the binary's tidal radius: it is merely
                        # deflected by the BBH's gravity but NOT disrupted. No mass transfer,
                        # no counter increment -- this is not a TDE.
                        pass

                    elif r_p > a:
                        # CHANNEL A: "whole-binary" disruption.
                        # The star's pericenter is inside r_t_bin but still outside the binary's
                        # own separation a, so to leading order (monopole approximation) the star
                        # is disrupted by the combined tidal field of the binary as a whole, not
                        # by either BH individually. tde_type = 22.

                        N_tdeBBHstar+=1 # update cumulative count of BBH-star TDEs (Channel A disruption)

                        # debris mass fraction accreted by each BH, weighted by its own mass
                        # (derived from combining geometric debris-proximity ~1/m_i and
                        # gravitational-focusing efficiency ~m_i^2, giving a net m_i scaling):
                        f1 = m1 / (m1+m2)
                        f2 = m2 / (m1+m2)

                        # mass increment accreted by each BH (only a fraction f_accreted of the
                        # disrupted star's mass is actually captured; the rest is ejected/lost):
                        dm1 = f1 * f_accreted * m_star
                        dm2 = f2 * f_accreted * m_star

                        # determine the neutron-star mass range for the chosen equation of state,
                        # needed by the spin-evolution routine to decide BH vs NS accretion physics:
                        if EoS=='APR':
                            M_NS_min = M_APR_min; M_NS_max = M_APR_max
                        elif EoS=='AU':
                            M_NS_min = M_AU_min; M_NS_max = M_AU_max
                        else: # invalid EoS string supplied by the user
                            sys.exit("Invalid EoS; please use 'APR' or 'AU'")

                        # classify each component as NS or BH based on its current mass:
                        NS1 = True if m1<M_NS_max else False
                        NS2 = True if m2<M_NS_max else False

                        # evolve the spin of each BH/NS as it accretes its share of the debris
                        # (prograde accretion assumed, consistent with the disk-formation picture
                        # used for the existing single-BH TDE channel):
                        evo1 = evolve(Mi=m1, Mf=m1+dm1, NS=NS1, f=None, chi=s1, dM=dm1/100, eos=EoS, prograde=True)
                        evo2 = evolve(Mi=m2, Mf=m2+dm2, NS=NS2, f=None, chi=s2, dM=dm2/100, eos=EoS, prograde=True)

                        # final spins after accretion:
                        s1_new = evo1['chi'][-1]
                        s2_new = evo2['chi'][-1]

                        # relative velocity between the star and the binary's center of mass,
                        # same functional form as used in the single-BH TDE rate calculation:
                        v_rel_tde = np.sqrt(v_star**2*m_avg/m_star + (m1+m2)/m_star*vBH**2)

                        # record one TDE entry per accreting BH (tde_type=22 flags this as a
                        # whole-binary/Channel-A disruption event); beta = r_t_bin/r_p is the
                        # penetration parameter, angles set to 0 since this isn't a single-BH
                        # marginally-bound-orbit calculation:
                        tdes = np.append(tdes, [[seed, t, z, 22, m_star, R_star, m1, s1, g1, r_t_bin, r_p, r_t_bin/r_p, 0.0, 0.0, dm1, s1_new, v_rel_tde, h1]], axis=0)
                        tdes = np.append(tdes, [[seed, t, z, 22, m_star, R_star, m2, s2, g2, r_t_bin, r_p, r_t_bin/r_p, 0.0, 0.0, dm2, s2_new, v_rel_tde, h2]], axis=0)

                        # update both BH masses, spins, and TDE counters after accretion:
                        m1 = m1 + dm1; s1 = s1_new; h1 = h1 + 1
                        m2 = m2 + dm2; s2 = s2_new; h2 = h2 + 1

                        # write the updated values back into the binaries array in place:
                        binaries[i][4] = m1; binaries[i][6] = s1; binaries[i][13] = h1
                        binaries[i][5] = m2; binaries[i][7] = s2; binaries[i][14] = h2

                    else:
                        # CHANNEL B: "single-component" disruption.
                        # The star's pericenter is small enough (r_p <= a) that it penetrates
                        # inside the binary's orbit and effectively has a close encounter with
                        # just one of the two BHs -- same physics as the existing single-BH TDE
                        # channel (BH_TidalDisruptions), just triggered during a BBH encounter.
                        # tde_type = 21.

                        N_tdeBBHstar+=1 # update cumulative count of BBH-star TDEs (Channel B disruption)

                        # select which BH is the disruptor, weighted by mass (same probability
                        # form used in the existing single-BH TDE function):
                        p1 = (m_star + m1) * m1**(1/3)
                        p2 = (m_star + m2) * m2**(1/3)
                        disrupt_1 = np.random.rand() < p1/(p1+p2)

                        # pull out the mass/spin/generation/tde-count of the chosen disruptor:
                        m_d = m1 if disrupt_1 else m2
                        s_d = s1 if disrupt_1 else s2
                        g_d = g1 if disrupt_1 else g2
                        h_d = h1 if disrupt_1 else h2

                        # standard single-BH tidal radius for the chosen disruptor:
                        r_t = R_star * (m_d/m_star)**(1/3)

                        # mass increment accreted by the disruptor:
                        dm = f_accreted * m_star

                        # determine NS mass range for the chosen EoS, as above:
                        if EoS=='APR':
                            M_NS_min = M_APR_min; M_NS_max = M_APR_max
                        elif EoS=='AU':
                            M_NS_min = M_AU_min; M_NS_max = M_AU_max
                        else:
                            sys.exit("Invalid EoS; please use 'APR' or 'AU'")

                        NS = True if m_d<M_NS_max else False

                        # evolve the disruptor's spin during accretion:
                        evo = evolve(Mi=m_d, Mf=m_d+dm, NS=NS, f=None, chi=s_d, dM=dm/100, eos=EoS, prograde=True)
                        s_new = evo['chi'][-1]

                        # relative velocity between star and the mean BH mass of the pair:
                        v_rel_tde = np.sqrt(v_star**2*m_avg/m_star + np.mean([m1, m2])/m_star*vBH**2)

                        # record the TDE entry (tde_type=21 flags this as a single-component/
                        # Channel-B disruption event), beta = r_t/r_p is the penetration parameter:
                        tdes = np.append(tdes, [[seed, t, z, 21, m_star, R_star, m_d, s_d, g_d, r_t, r_p, r_t/r_p, 0.0, 0.0, dm, s_new, v_rel_tde, h_d]], axis=0)

                        # write the updated mass/spin/tde-count back to whichever component
                        # was the disruptor, leaving the other BH untouched:
                        if disrupt_1:
                            m1 = m_d + dm; s1 = s_new; h1 = h_d + 1
                            binaries[i][4] = m1; binaries[i][6] = s1; binaries[i][13] = h1
                        else:
                            m2 = m_d + dm; s2 = s_new; h2 = h_d + 1
                            binaries[i][5] = m2; binaries[i][7] = s2; binaries[i][14] = h2


                if type_int==3: # BBH-BBH interaction occurs
                    
                    N_bb+=1
                    
                    # BBH semimajor axes (excluding the current one):
                    other_binaries = np.delete(binaries, i, axis=0)
                    smas = other_binaries[:, 2]

                    # BBH masses (excluding the current one):
                    masses = other_binaries[:, 4] + other_binaries[:, 5]
                    
                    # semimajor axis of current BBH:
                    a1 = a
                    
                    # sample binary:
                    valid = np.where(masses * smas > 0)[0]
                    if len(valid) < 1:
                        continue
                    a2 = np.random.choice(smas[valid], p=(masses * smas)[valid] / np.sum((masses * smas)[valid]))

                    # index of current BBH:
                    k1 = i
                    
                    # find index location of binary:
                    k2 = np.squeeze(np.where(binaries[:, 2]==a2))+0
                    
                    k2 = int(np.atleast_1d(k2)[0])
                        
                    # order based on hardness:
                    if a1 < a2:
                        k_soft = k2; a_soft = a2
                        k_hard = k1; a_hard = a1
                    else:
                        k_soft = k1; a_soft = a1
                        k_hard = k2; a_hard = a2
                        
                    # BHs of tighter BBH:
                    m0 = binaries[k_hard][4]; s0 = binaries[k_hard][6]; g0 = binaries[k_hard][8]; h0 = binaries[k_hard][13]
                    m1 = binaries[k_hard][5]; s1 = binaries[k_hard][7]; g1 = binaries[k_hard][9]; h1 = binaries[k_hard][14]
                
                    # BHs of softer BBH:
                    if binaries[k_soft][4] < binaries[k_soft][5]:
                        m_freed = binaries[k_soft][4]; s_freed = binaries[k_soft][6]; g_freed = binaries[k_soft][8]; h_freed = binaries[k_soft][13]
                        m2 = binaries[k_soft][5]; s2 = binaries[k_soft][7]; g2 = binaries[k_soft][9]; h2 = binaries[k_soft][14]
                    else:
                        m_freed = binaries[k_soft][5]; s_freed = binaries[k_soft][7]; g_freed = binaries[k_soft][9]; h_freed = binaries[k_soft][14]
                        m2 = binaries[k_soft][4]; s2 = binaries[k_soft][6]; g2 = binaries[k_soft][8]; h2 = binaries[k_soft][13]
                    
                    mBH = np.append(mBH, m_freed)
                    sBH = np.append(sBH, s_freed)
                    gBH = np.append(gBH, g_freed)
                    hBH = np.append(hBH, h_freed)
                    
                    # delete softer BBH:
                    binaries = np.delete(binaries, k_soft, axis=0)
                    
                    N_dis+=1
                    
                    N_BBH = N_BBH - 1
                    
                    if k_hard > k_soft:
                        k_hard = k_hard - 1
                        
                    a_inner = a_hard
                    a_outer = a_soft * (m0 + m1) / m_freed
                    e_inner = np.sqrt(np.random.rand())
                    e_outer = np.sqrt(np.random.rand())
                    inclination1 = np.arccos(np.random.uniform(-1, 1))
                    inclination2 = np.arccos(np.random.uniform(-1, 1))
                    inclination = inclination1 + inclination2
                    
                    i_old = i
                    
                    # check Mardling & Aarseth (2001) criterion:
                    if a_outer/a_inner > 2.8 * (1 + m2 / (m1 + m0))**(2/5) * (1 + e_outer)**(2/5) / (1 - e_outer)**(6/5) * (1 - 0.3 * inclination / np.pi): # hierarchical triple
                        
                        # append triple:
                        triples = np.append(triples, [[a_inner, a_outer, e_inner, e_outer, m0, m1, m2, s0, s1, s2, g0, g1, g2, inclination1, inclination2, t, z,
                                                       binaries[k_hard][0], binaries[k_hard][1], binaries[k_hard][10], binaries[k_hard][11], h0, h1, h2]], axis=0)
                        
                        N_Triples+=1
                        
                        binaries = np.delete(binaries, k_hard, axis=0)
                        
                        N_BBH = N_BBH - 1
                        
                        i = i - 1
                        
                    else: # `inner' binary is freed:
                        
                        mBH = np.append(mBH, m2)
                        sBH = np.append(sBH, s2)
                        gBH = np.append(gBH, g2)
                        hBH = np.append(hBH, h2)
                        
                        # append binary:
                        Delta = 0.38 # Zevin et al., ApJ 871 (2019), 91.
                        binaries[k_hard][2] = a_inner/(1 + Delta * m_freed * m2 / m1 / m0 * a_inner / a_outer)
                        binaries[k_hard][3] = e_inner
                        
                    if k_soft < i_old:
                        i = i - 1
                        continue
                    elif k_soft == i_old:
                        i = i - 1
                        condition=4
                        hardening[i][10]=condition
                        break
                    else:
                        continue
                    
                # at this point BBH-BH interaction occurs
                
                if type_int==1:
                    
                    m3 = m_avg
                    
                    # single velocity before interaction:
                    vS_before = v_star

                else:
                    if mBH.size == 0:  # no single BHs left to interact with, skip
                        continue
                    p3 = (m1 + m2 + mBH) / np.sqrt((m1 + m2)**(-2/5) + mBH**(-2/5)) * mBH**(3/2)
                    m3 = np.random.choice(mBH, replace=False, p=p3/np.sum(p3))
                    
                    k3 = np.squeeze(np.where(mBH==m3))+0
                    
                    k3 = int(np.atleast_1d(k3)[0])
                    
                    s3 = sBH[k3]
                    g3 = gBH[k3]
                    h3 = hBH[k3]
                    
                    # sample single velocity before interaction:
                    vS_before = get_maxwell_sample(np.sqrt(mBH_avg * vBH**2 / 3 / m3))
                    
                # sample binary velocity before interaction:
                vB_before = get_maxwell_sample(np.sqrt(mBH_avg * vBH**2 / 3 / (m1 + m2)))

                # sample cosine angle before interaction:
                cos_theta_before = np.random.uniform(-1, 1)
                
                # relative velocity before interaction:
                v_rel_before = np.sqrt(vB_before**2 + vS_before**2 - 2 * vB_before * vS_before * cos_theta_before)
                
                # reduced mass before interaction:
                mu_before = (m1 + m2) * m3 / (m1 + m2 + m3)
                
                if 1/2 * mu_before * v_rel_before**2 > G_Newton * m1 * m2 / 2 / a: # binary is ionized
                    
                    N_dis+=1
                    
                    binaries = np.delete(binaries, i, axis=0)
                    
                    i = i - 1
                    
                    mBH = np.append(mBH, m1); mBH = np.append(mBH, m2)
                    sBH = np.append(sBH, s1); sBH = np.append(sBH, s2)
                    gBH = np.append(gBH, g1); gBH = np.append(gBH, g2)
                    hBH = np.append(hBH, h1); hBH = np.append(hBH, h2)
                    
                    N_BBH = N_BBH - 1
                    
                    condition=3
                    hardening[i][10]=condition
                    break
                
                # sample pericenter of interaction:
                rp = np.random.uniform(0, kp_max * a)
                
                # critical pericenter for resonant interaction:
                rp_c = np.max([m1, m2]) / (m1 + m2) * a
                
                if rp < rp_c and type_int==2: # interaction is resonant
                    
                    u_3bm = np.random.rand(3)
                    
                    # semimajor axes of IMS assuming energy conservation:
                    a_3bm = a / (m1*m2) * np.array([m1*m2, m2*m3, m3*m1])
                    
                    # critical pericenter distances for GW capture merger during binary-single interaction:
                    rp_3bm_c = (85 * np.pi / 3 / np.sqrt(2))**(2/7) * a**(2/7) / 2 \
                        * np.array([(2 * G_Newton * (m1 * m2)**(4/5) *(m1 + m2)**(1/5) /(m1*m2)**(2/5) /c_light**2)**(5/7),
                                    (2 * G_Newton * (m2 * m3)**(4/5) *(m2 + m3)**(1/5) /(m1*m2)**(2/5) /c_light**2)**(5/7),
                                    (2 * G_Newton * (m3 * m1)**(4/5) *(m3 + m1)**(1/5) /(m1*m2)**(2/5) /c_light**2)**(5/7)])
                    
                    # critical eccentricities for 3-body merger:
                    e_3bm_c = 1 - rp_3bm_c / a_3bm
                    
                    # probabilities for 3-body merger:
                    p_3bm = (1 - e_3bm_c**(2 * N_IMS))
                    
                    e_3bm = np.sqrt(1 - (1 - e_3bm_c) * np.random.rand(3))
                    
                    ratios_3bm = p_3bm / u_3bm
                    j_3bm = np.argmax(ratios_3bm)
                    
                    m_3bm_1 = np.array([m1, m2, m3]); s_3bm_1 = np.array([s1, s2, s3]); g_3bm_1 = np.array([g1, g2, g3]); h_3bm_1 = np.array([h1, h2, h3])
                    m_3bm_2 = np.array([m2, m3, m1]); s_3bm_2 = np.array([s2, s3, s1]); g_3bm_2 = np.array([g2, g3, g1]); h_3bm_2 = np.array([h2, h3, h1])
                    m_3bm_3 = np.array([m3, m1, m2]); s_3bm_3 = np.array([s3, s1, s2]); g_3bm_3 = np.array([g3, g1, g2]); h_3bm_3 = np.array([h3, h1, h2])
                    
                    if ratios_3bm[j_3bm] > 1: # binary-single GW capture merger occurs
                        
                        mBH = np.delete(mBH, k3)
                        sBH = np.delete(sBH, k3)
                        gBH = np.delete(gBH, k3)
                        hBH = np.delete(hBH, k3)
                        
                        mA = m_3bm_1[j_3bm]; sA = s_3bm_1[j_3bm]; gA = g_3bm_1[j_3bm]; hA = h_3bm_1[j_3bm]
                        mB = m_3bm_2[j_3bm]; sB = s_3bm_2[j_3bm]; gB = g_3bm_2[j_3bm]; hB = h_3bm_2[j_3bm]
                        mC = m_3bm_3[j_3bm]; sC = s_3bm_3[j_3bm]; gC = g_3bm_3[j_3bm]; hC = h_3bm_3[j_3bm]
                            
                        thetaA, thetaB, dPhi = sample_angles()
                            
                        m_rem, s_rem, vGW_kick = merger_remnant(mA, mB, sA, sB, thetaA, thetaB, dPhi)
                        g_rem = np.max([gA, gB]) + 1
                        h_rem = h1 + h2
                        
                        sma = a_3bm[j_3bm]
                        eccen = e_3bm[j_3bm]
                        
                        if vGW_kick < 2 * np.sqrt(v_star**2 + vBH**2): # merger remnant retained in cluster
                            
                            mBH = np.append(mBH, m_rem)
                            sBH = np.append(sBH, s_rem)
                            gBH = np.append(gBH, g_rem)
                            hBH = np.append(hBH, h_rem)
                            
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
                        
                        # append merger:
                        mergers = np.append(mergers,
                                            [[seed, ind, 6, sma, eccen, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form, 
                                              t + t_local + T_GW(mA, mB, sma, eccen), redshift_interp(lookback_interp(zCl_form) - t - t_local - T_GW(mA, mB, sma, eccen)), m_rem, s_rem, g_rem, vGW_kick, s_eff, q, 2*v_star, h1, h2]], axis=0)
                        
                        # delete BBH:
                        binaries = np.delete(binaries, i, axis=0)
                        
                        N_BBH = N_BBH - 1
                        
                        i = i - 1
                        
                        mBH = np.append(mBH, mC)
                        sBH = np.append(sBH, sC)
                        gBH = np.append(gBH, gC)
                        hBH = np.append(hBH, hC)
                        
                        condition=6
                        hardening[i][10]=condition
                        break
                        
                if m3>np.min([m1, m2]) and rp<rp_c and type_int==2: # exchange occurs

                    N_ex+=1
                    binaries[i][12]+=1
                    
                    if m1<m2:

                        ms = m1; ss = s1; gs = g1; hs = h1
                        mr = m2; sr = s2; gr = g2; hr = h2
                        
                    else:
                        
                        ms = m2; ss = s2; gs = g2; hs = h2
                        mr = m1; sr = s1; gr = g1; hr = h1
                        
                    # impose binding energy conservation after substitution:
                    binaries[i][2] = m3/ms * binaries[i][2]

                    mBH = np.delete(mBH, k3)
                    sBH = np.delete(sBH, k3)
                    gBH = np.delete(gBH, k3)
                    hBH = np.delete(hBH, k3)

                    mBH = np.append(mBH, ms)
                    sBH = np.append(sBH, ss)
                    gBH = np.append(gBH, gs)
                    hBH = np.append(hBH, hs)

                    # update binary masses:
                    binaries[i][4] = mr
                    binaries[i][5] = m3
                    
                    # update binary spins:
                    binaries[i][6] = sr
                    binaries[i][7] = s3
                    
                    # update binary generations:
                    binaries[i][8] = gr
                    binaries[i][9] = g3
                    
                    # update binary tde counters:
                    binaries[i][13] = hr
                    binaries[i][14] = h3
                    
                    m1 = mr; s1 = sr; g1 = gr; h1 = hr
                    m2 = m3; s2 = s3; g2 = g3; h2 = h3
                    m3 = ms; s3 = ss; g3 = gs; h3 = hs
                    
                    k3 = mBH.size - 1
                    
                    a = binaries[i][2]
                    e = binaries[i][3]
                    
                # Binary hardening:
                binaries[i][2] = a / (1 + Hardening_constant * m3 / (m1 + m2))
                binaries[i][3] = np.sqrt(np.random.rand())
                
                # energy extracted:
                dE_b = G_Newton * m1 * m2 / 2 * (1/binaries[i][2] - 1/a)
                
                a = binaries[i][2]
                e = binaries[i][3]
                
                # reduced mass after interaction:
                mu_after = (m1 + m2) * m3 / (m1 + m2 + m3)
                
                # relative velocity after interaction:
                v_rel_after = np.sqrt(mu_before / mu_after * v_rel_before**2 + 2 / mu_after * dE_b)
                
                # velocity of single after interaction:
                v3_after = (m1 + m2) / (m1 + m2 + m3) * v_rel_after
                
                # velocity of binary after interaction:
                v12_after = m3 / (m1 + m2 + m3) * v_rel_after
                
                # check if single is ejected:
                if v3_after > 2 * np.sqrt(v_star**2 + vBH**2) and type_int==2:
                    
                    mBH = np.delete(mBH, k3)
                    sBH = np.delete(sBH, k3)
                    gBH = np.delete(gBH, k3)
                    hBH = np.delete(hBH, k3)

                    N_BH = N_BH - 1
                    
                    N_BHej+=1
                    
                if v12_after > 2 * vBH \
                   and v12_after < 2 * np.sqrt(v_star**2 + vBH**2): # binary convection
                    
                    t_conv = m_avg / (m1 + m2) * t_rlx
                    
                # check if binary is ejected:
                if v12_after > 2 * np.sqrt(v_star**2 + vBH**2):
                    
                    # check if BBH mergers in the field:
                    if t + t_local + T_GW(m1, m2, a, e) < lookback_interp(zCl_form): # BBH merges
                        
                        theta1, theta2, dPhi = sample_angles()
                        
                        m_rem, s_rem, vGW_kick = merger_remnant(m1, m2, s1, s2, theta1, theta2, dPhi)
                        g_rem = np.max([g1, g2]) + 1
                        h_rem = h1 + h2

                        N_me+=1
                        N_meFi+=1
                        
                        # merger time:
                        t_merge = t + t_local + T_GW(m1, m2, a, e)
                        
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
                        
                        # append merger:
                        mergers = np.append(mergers, [[seed, ind, -channel, a, e, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form, t_merge,
                                                       redshift_interp(lookback_interp(zCl_form) - t_merge), m_rem, s_rem, g_rem, vGW_kick, s_eff, q, 2*v_star, h1, h2]], axis=0)
                        
                    binaries = np.delete(binaries, i, axis=0)
                    
                    i = i - 1
                    
                    N_BBH = N_BBH - 1
                    
                    N_BH = N_BH - 2
                    
                    N_BBHej+=1
                    
                    condition=5
                    hardening[i][10]=condition
                    break

            hardening = np.append(hardening, [[t, dt, t_local, dt_local, ind, a, e, m1, m2, q, condition, Nex]], axis=0)
            N_hardening+=1
            i+=1
            
    return seed, t, z, dt, zCl_form, binaries, hardening, mergers, mBH, sBH, gBH, hBH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples, tdes, N_tdeBBHstar, m_min, m_max, f_accreted, EoS

# End of file.
