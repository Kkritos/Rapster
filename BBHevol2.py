'''
 Copyright (C) 2023  Konstantinos Kritos <kkritos1@jhu.edu>

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

from constants2 import *
from functions2 import *

#@njit
def evolve_BBHs(seed, t, z, dt, zCl_form, binaries, hardening, mergers, mBH, sBH, gBH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples):
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
    
    @out: all inputs
    """
    
    # BBH evolution:
    if N_BBH>0:
        
        # shuffle binaries to avoid biases:
        np.random.shuffle(binaries[1:binaries.size])
        
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
            while N_BH>0:
                
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
                
                # BBH-star interaction timescale:
                t_BBH_star = 1e100 #1 / Rate_int(m1+m2+m_avg, n_star, v_star, kp_max * a)
                
                # BBH-BH interaction timescale:
                if mBH.size>0:
                    t_BBH_BH = 1 / Rate_int(m1+m2+mBH_avg, nc_BH, vBH, kp_max * a) + t_conv
                else:
                    t_BBH_BH = 1e100
                
                # BBH-BBH interaction timescale:
                if N_BBH>1:
                    # BBH number density:
                    n_BBH = (N_BBH - 1) / Vc_BH
                    
                    t_BBH_BBH = 1 / Rate_int(m1 + m2 + np.mean(np.transpose(binaries)[:][4]+np.transpose(binaries)[:][5]), n_BBH, vBH, kp_max * a) + t_conv
                else:
                    t_BBH_BBH = 1e100
                    
                # probability for BBH-star encounter:
                p_BBH_star = t_BBH_BH * t_BBH_BBH / (t_BBH_star * t_BBH_BH + t_BBH_BH * t_BBH_BBH + t_BBH_BBH * t_BBH_star)
                
                # probability for BBH-BH encounter:
                p_BBH_BH = t_BBH_BBH / (t_BBH_BBH + t_BBH_BH) #t_BBH_star * t_BBH_BBH / (t_BBH_star * t_BBH_BH + t_BBH_BH * t_BBH_BBH + t_BBH_BBH * t_BBH_star)
                
                # probability for BBH-BBH encounter:
                p_BBH_BBH = 1 - p_BBH_BH #- p_BBH_star
                
                u_int = np.random.rand()
                
                if u_int < p_BBH_BH: # do BBH-BH interaction
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
                    
                    if vGW_kick < 2 * np.sqrt(v_star**2 + vBH**2): # merger remnant retained in cluster
                        
                        mBH = np.append(mBH, m_rem)
                        sBH = np.append(sBH, s_rem)
                        gBH = np.append(gBH, g_rem)
                        
                        N_BH = N_BH - 1
                        
                        N_meRe+=1
                        
                    else: # merger remnant ejected
                        
                        N_BH = N_BH - 2
                        
                        N_meEj+=1
                        
                    N_me+=1
                    N_me2b+=1
                    
                    # order BHs by mass:
                    mA = m1; sA = s1; gA = g1; thetaA = theta1
                    mB = m2; sB = s2; gB = g2; thetaB = theta2
                    if mA>mB:
                        m1 = mA; s1 = sA; g1 = gA; theta1 = thetaA
                        m2 = mB; s2 = sB; g2 = gB; theta2 = thetaB
                    else:
                        m1 = mB; s1 = sB; g1 = gB; theta1 = thetaB
                        m2 = mA; s2 = sA; g2 = gA; theta2 = thetaA
                        
                    # mass ratio:
                    q = m2 / m1
                    
                    # effective spin parameter:
                    s_eff = (m1 * s1 * np.cos(theta1) + m2 * s2 * np.cos(theta2)) / (m1 + m2)
                    
                    # append merger:
                    mergers = np.append(mergers, [[seed, ind, channel, a, e, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form,
                                                   t + t_local + T_GW(m1, m2, a, e), redshift_interp(lookback_interp(zCl_form) - t - t_local - T_GW(m1, m2, a, e)), m_rem, s_rem, g_rem, vGW_kick, s_eff, q]], axis=0)
                    
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
                    
                if type_int==3: # BBH-BBH interaction occurs
                    
                    N_bb+=1
                    
                    # BBH semimajor axes (excluding the current one):
                    smas = np.transpose(np.delete(binaries, i, axis=0))[:][2]
                    
                    # BBH masses (excluding the current one):
                    masses = np.transpose(np.delete(binaries, i, axis=0))[:][4] + np.transpose(np.delete(binaries, i, axis=0))[:][5]
                    
                    # semimajor axis of current BBH:
                    a1 = a
                    
                    # sample binary:
                    a2 = np.random.choice(smas, p=masses * smas / np.sum(masses * smas))

                    # index of current BBH:
                    k1 = i
                    
                    # find index location of binary:
                    k2 = np.squeeze(np.where(np.transpose(binaries)[:][2]==a2))+0
                    
                    if isinstance(k2, np.ndarray):
                        k2 = k2[0]
                        
                    # order based on hardness:
                    if a1 < a2:
                        k_soft = k2; a_soft = a2
                        k_hard = k1; a_hard = a1
                    else:
                        k_soft = k1; a_soft = a1
                        k_hard = k2; a_hard = a2
                        
                    # BHs of tighter BBH:
                    m0 = binaries[k_hard][4]; s0 = binaries[k_hard][6]; g0 = binaries[k_hard][8]
                    m1 = binaries[k_hard][5]; s1 = binaries[k_hard][7]; g1 = binaries[k_hard][9]
                
                    # BHs of softer BBH:
                    if binaries[k_soft][4] < binaries[k_soft][5]:
                        m_freed = binaries[k_soft][4]; s_freed = binaries[k_soft][6]; g_freed = binaries[k_soft][8]
                        m2 = binaries[k_soft][5]; s2 = binaries[k_soft][7]; g2 = binaries[k_soft][9]
                    else:
                        m_freed = binaries[k_soft][5]; s_freed = binaries[k_soft][7]; g_freed = binaries[k_soft][9]
                        m2 = binaries[k_soft][4]; s2 = binaries[k_soft][6]; g2 = binaries[k_soft][8]
                    
                    mBH = np.append(mBH, m_freed)
                    sBH = np.append(sBH, s_freed)
                    gBH = np.append(gBH, g_freed)
                    
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
                    if a_outer/a_inner > 2.8 * (1 + m2 / (m1 + m0))**(2/5) * (1 + e_outer)**(2/5) / (1 - e_outer)**(6/5) * (1 - 0.3 * inclination / np.pi) and 1<2: # hierarchical triple
                        
                        # append triple:
                        triples = np.append(triples, [[a_inner, a_outer, e_inner, e_outer, m0, m1, m2, s0, s1, s2, g0, g1, g2, inclination1, inclination2, t, z,
                                                       binaries[k_hard][0], binaries[k_hard][1], binaries[k_hard][10], binaries[k_hard][11]]], axis=0)
                        
                        N_Triples+=1
                        
                        binaries = np.delete(binaries, k_hard, axis=0)
                        
                        N_BBH = N_BBH - 1
                        
                        i = i - 1
                        
                    else: # `inner' binary is freed:
                        
                        mBH = np.append(mBH, m2)
                        sBH = np.append(sBH, s2)
                        gBH = np.append(gBH, g2)
                        
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
                    p3 = (m1 + m2 + mBH) / np.sqrt((m1 + m2)**(-2/5) + mBH**(-2/5)) * mBH**(3/2)
                    m3 = np.random.choice(mBH, replace=False, p=p3/np.sum(p3))
                    
                    k3 = np.squeeze(np.where(mBH==m3))+0
                    
                    if isinstance(k3, np.ndarray):
                        k3 = k3[0]
                    
                    s3 = sBH[k3]
                    g3 = gBH[k3]
                    
                    # sample single velocity before interaction:
                    vS_before = maxwell.rvs(loc=0, scale=np.sqrt(mBH_avg * vBH**2 / 3 / m3))
                    
                # sample binary velocity before interaction:
                vB_before = maxwell.rvs(loc=0, scale=np.sqrt(mBH_avg * vBH**2 / 3 / (m1 + m2)))

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
                    
                    m_3bm_1 = np.array([m1, m2, m3]); s_3bm_1 = np.array([s1, s2, s3]); g_3bm_1 = np.array([g1, g2, g3])
                    m_3bm_2 = np.array([m2, m3, m1]); s_3bm_2 = np.array([s2, s3, s1]); g_3bm_2 = np.array([g2, g3, g1])
                    m_3bm_3 = np.array([m3, m1, m2]); s_3bm_3 = np.array([s3, s1, s2]); g_3bm_3 = np.array([g3, g1, g2])
                    
                    if ratios_3bm[j_3bm] > 1: # binary-single GW capture merger occurs
                        
                        mBH = np.delete(mBH, k3)
                        sBH = np.delete(sBH, k3)
                        gBH = np.delete(gBH, k3)
                        
                        mA = m_3bm_1[j_3bm]; sA = s_3bm_1[j_3bm]; gA = g_3bm_1[j_3bm]
                        mB = m_3bm_2[j_3bm]; sB = s_3bm_2[j_3bm]; gB = g_3bm_2[j_3bm]
                        mC = m_3bm_3[j_3bm]; sC = s_3bm_3[j_3bm]; gC = g_3bm_3[j_3bm]
                            
                        thetaA, thetaB, dPhi = sample_angles()
                            
                        m_rem, s_rem, vGW_kick = merger_remnant(mA, mB, sA, sB, thetaA, thetaB, dPhi)
                        g_rem = np.max([gA, gB]) + 1
                        
                        sma = a_3bm[j_3bm]
                        eccen = e_3bm[j_3bm]
                        
                        if vGW_kick < 2 * np.sqrt(v_star**2 + vBH**2): # merger remnant retained in cluster
                            
                            mBH = np.append(mBH, m_rem)
                            sBH = np.append(sBH, s_rem)
                            gBH = np.append(gBH, g_rem)
                            
                            N_BH = N_BH - 1
                            N_meRe+=1
                            
                        else:
                            
                            N_BH = N_BH - 2
                            N_meEj+=1
                            
                        N_3cap+=1
                        N_me+=1
                        
                        # order BHs my mass:
                        if mA>mB:
                            m1 = mA; s1 = sA; g1 = gA; theta1 = thetaA
                            m2 = mB; s2 = sB; g2 = gB; theta2 = thetaB
                        else:
                            m1 = mB; s1 = sB; g1 = gB; theta1 = thetaB
                            m2 = mA; s2 = sA; g2 = gA; theta2 = thetaA
                            
                        # mass ratio:
                        q = m2 / m1
                        
                        # effective spin parameter:
                        s_eff = (m1 * s1 * np.cos(theta1) + m2 * s2 * np.cos(theta2)) / (m1 + m2)
                        
                        # append merger:
                        mergers = np.append(mergers,
                                            [[seed, ind, 6, sma, eccen, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form, 
                                              t + t_local + T_GW(mA, mB, sma, eccen), redshift_interp(lookback_interp(zCl_form) - t - t_local - T_GW(mA, mB, sma, eccen)), m_rem, s_rem, g_rem, vGW_kick, s_eff, q]], axis=0)
                        
                        # delete BBH:
                        binaries = np.delete(binaries, i, axis=0)
                        
                        N_BBH = N_BBH - 1
                        
                        i = i - 1
                        
                        mBH = np.append(mBH, mC)
                        sBH = np.append(sBH, sC)
                        gBH = np.append(gBH, gC)
                        
                        condition=6
                        hardening[i][10]=condition
                        break
                        
                if m3>np.min([m1, m2]) and rp<rp_c and type_int==2: # exchange occurs

                    N_ex+=1
                    binaries[i][12]+=1
                    
                    if m1<m2:

                        ms = m1; ss = s1; gs = g1
                        mr = m2; sr = s2; gr = g2
                        
                    else:
                        
                        ms = m2; ss = s2; gs = g2
                        mr = m1; sr = s1; gr = g1
                        
                    # impose binding energy conservation after substitution:
                    binaries[i][2] = m3/ms * binaries[i][2]

                    mBH = np.delete(mBH, k3)
                    sBH = np.delete(sBH, k3)
                    gBH = np.delete(gBH, k3)
                    
                    mBH = np.append(mBH, ms)
                    sBH = np.append(sBH, ss)
                    gBH = np.append(gBH, gs)
                    
                    # update binary masses:
                    binaries[i][4] = mr
                    binaries[i][5] = m3
                    
                    # update binary spins:
                    binaries[i][6] = sr
                    binaries[i][7] = s3
                    
                    # update binary generations:
                    binaries[i][8] = gr
                    binaries[i][9] = g3
                    
                    m1 = mr; s1 = sr; g1 = gr
                    m2 = m3; s2 = s3; g2 = g3
                    m3 = ms; s3 = ss; g3 = gs
                    
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
                        
                        N_me+=1
                        N_meFi+=1
                        
                        # merger time:
                        t_merge = t + t_local + T_GW(m1, m2, a, e)
                        
                        # order BHs by mass:
                        mA = m1; sA = s1; gA = g1; thetaA = theta1
                        mB = m2; sB = s2; gB = g2; thetaB = theta2
                        if mA>mB:
                            m1 = mA; s1 = sA; g1 = gA; theta1 = thetaA
                            m2 = mB; s2 = sB; g2 = gB; theta2 = thetaB
                        else:
                            m1 = mB; s1 = sB; g1 = gB; theta1 = thetaB
                            m2 = mA; s2 = sA; g2 = gA; theta2 = thetaA
                            
                        # mass ratio:
                        q = m2 / m1
                        
                        # effective spin parameter:
                        s_eff = (m1 * s1 * np.cos(theta1) + m2 * s2 * np.cos(theta2)) / (m1 + m2)
                        
                        # append merger:
                        mergers = np.append(mergers, [[seed, ind, -channel, a, e, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form, t_merge,
                                                       redshift_interp(lookback_interp(zCl_form) - t_merge), m_rem, s_rem, g_rem, vGW_kick, s_eff, q]], axis=0)
                        
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
            
    return seed, t, z, dt, zCl_form, binaries, hardening, mergers, mBH, sBH, gBH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples

# end of file
