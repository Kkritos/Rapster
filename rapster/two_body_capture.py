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
from .remnant import *
from .bhpool import FAM

def sample_2capture(BH):
    """
    Sample the two single BHs that form a binary by GW capture, from the pool BH
    (bhpool.BHPool), with the distribution of the sampler of versions <= 2.13.1
    (fast_sample_2capture): the pair (i, j), i != j, drawn as
    np.random.choice(size=2, replace=False) with proposal probabilities p ~ m^(2/7)
    (i first, then j among the others) and accepted with probability
    ~ (m_i + m_j)^(10/7), i.e.
        P(i, j) ~ p_i p_j / (1 - p_i) (m_i + m_j)^(10/7).

    Rejection sampling: with w = m^(2/7) and (m_i + m_j)^(10/7) <=
    2^(3/7) (m_i^(10/7) + m_j^(10/7)) (convexity), candidates are drawn from
        w_i w_j (m_i^(10/7) + m_j^(10/7)) = c127_i c27_j + c27_i c127_j
    (pool families c27 = m^(2/7), c127 = m^(12/7); equal pairs i = j rejected) and
    accepted with probability
        (m_i + m_j)^(10/7) / (2^(3/7) (m_i^(10/7) + m_j^(10/7))) * (1 - p_max) / (1 - p_i),
    with p_max the largest p. The first factor is at least 2^(-3/7) ~ 0.74, the second
    is close to 1 unless one BH dominates the m^(2/7) weights.

    @in BH: pool of single BHs (at least 2)

    @out: m1, m2, slot of m1, slot of m2
    """
    f27 = FAM['c27']; f127 = FAM['c127']
    S27 = BH.total(f27)
    q_max = 1 - BH.weight(f27, BH.D) / S27      # 1 - p_max (the heaviest BH has the largest weight)
    m = BH.m
    c = 2**(3/7)
    while True:
        if np.random.rand() < 0.5:
            i = BH.draw(f127); j = BH.draw(f27)
        else:
            i = BH.draw(f27); j = BH.draw(f127)
        if i == j:
            continue
        mi = m[i]; mj = m[j]
        accept = (mi + mj)**(10/7) / (c * (mi**(10/7) + mj**(10/7))) * q_max / (1 - BH.weight(f27, i) / S27)
        if np.random.rand() < accept:
            return mi, mj, i, j

def two_body_capture(seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, BH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers, random_pairing=False):
    """
    @in seed: simulation seed number
    @in t: simulation time
    @in dt: simulation time step
    @in z: simulation redshift
    @in zCl_form: cluster formation redshift
    @in k_2cap: number of 2-captures in current step
    @in mBH_avg: average BH mass
    @in binaries: array of BBHs
    @in BH: pool of single BHs (bhpool.BHPool: masses, spins, generations, TDE counts)
    @in vBH: 3D BH velocity dispersion
    @in v_star: 3D star velocity dispersion
    @in N_2cap: number of 2-captures
    @in N_BH: number of BHs
    @in N_BBH: number of BBHs
    @in N_me: number of mergers
    @in N_meRe: number of retained mergers
    @in N_meEj: number of ejected mergers
    @in mergers: array of mergers: [seed, ind, channel, a, e, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form, t_merge, z_merge, m_rem, s_rem, g_rem, vGW_kick, s_eff, q, v_esc, h1, h2]
    @in random_pairing: if True, use uniform random pairing instead of mass-weighted (m^2)

    @out: all inputs
    """
    
    if k_2cap > 0:

        retained = []   # merger remnants retained in the cluster (added after the loop)

        for i in range(k_2cap):
            if BH.size < 2:
                break # avoid infinite loop in the sampler
            
            # sample the BHs that form the captured binary (slots k1 != k2):
            if random_pairing:
                k1 = BH.draw_uniform()
                k2 = BH.draw_uniform()
                while k2 == k1:
                    k2 = BH.draw_uniform()
                m1, m2 = BH.m[k1], BH.m[k2]
            else:
                m1, m2, k1, k2 = sample_2capture(BH)
                
            s1 = BH.s[k1]; g1 = BH.g[k1]; h1 = BH.h[k1]
            s2 = BH.s[k2]; g2 = BH.g[k2]; h2 = BH.h[k2]
            
            ind = np.random.randint(0, 999999999)
            
            theta1, theta2, dPhi = sample_angles()
            
            m_rem, s_rem, vGW_kick = merger_remnant(m1, m2, s1, s2, theta1, theta2, dPhi)
            g_rem = np.max([g1, g2]) + 1
            h_rem = h1 + h2

            # relative velocity:
            v_rel = get_maxwell_sample(np.sqrt(2/3) * vBH)
            
            # total mass:
            m12 = m1 + m2
            
            # reduced mass:
            mu = m1*m2/m12
            
            # maximum impact parameter for capture:
            b_max = (340 * np.pi / 3)**(1/7) * m12**(6/7) * mu**(1/7) / v_rel**(9/7) * G_Newton * c_light**(-5/7)

            # impact parameter sampled from uniform in b^2 distribution:
            b = np.sqrt(np.random.rand() * b_max**2)
            
            # pericenter distance:
            rp = b**2 * v_rel**2 / 2 / G_Newton / m12
            
            # GW energy released:
            E_gw = 85 * np.pi / 12 / np.sqrt(2) * mu**2 * m12**(5/2) / rp**(7/2) * G_Newton**(7/2) / c_light**5
            
            # final energy:
            E_fin = mu * v_rel**2 / 2 - E_gw
            
            max_iter = 1000  # safety cap to prevent infinite loop for extreme mass/velocity combinations
            n_iter = 0
            # make sure eccentricity is strictly smaller than unity:
            while 1 + 2 * E_fin * b**2 * v_rel**2 / m12**2 / mu / G_Newton**2 < 0:
                n_iter += 1
                if n_iter > max_iter:  # could not find valid eccentricity, skip this capture event
                    break
                
                # impact parameter sampled from uniform in b^2 distribution:
                b = np.sqrt(np.random.rand() * b_max**2)
                
                # pericenter distance:
                rp = b**2 * v_rel**2 / 2 / G_Newton / m12
                
                # GW energy released:
                E_gw = 85 * np.pi / 12 / np.sqrt(2) * mu**2 * m12**(5/2) / rp**(7/2) * G_Newton**(7/2) / c_light**5
                
                # final energy:
                E_fin = mu * v_rel**2 / 2 - E_gw
                
            if n_iter > max_iter:  # capture event skipped due to no valid eccentricity found
                continue
                
            # semimajor axis at formation:
            sma = - G_Newton * m12 * mu / 2 / E_fin
            
            # eccentricity at formation:
            eccen = np.sqrt(1 + 2 * E_fin * b**2 * v_rel**2 / m12**2 / mu / G_Newton**2)
            
            # delete captured BHs:
            BH.remove_two(k1, k2)
            
            N_2cap+=1
            
            # check if binary merges within the current step:
            if T_GW(m1, m2, sma, eccen) < np.min([dt, lookback_interp(zCl_form) - t]):
                
                if vGW_kick < v_esc_cl(v_star, vBH): # merger remnant retained in cluster
                    
                    retained.append((m_rem, s_rem, g_rem, h_rem))
                    
                    N_BH = N_BH - 1
                    
                    N_meRe+=1
                        
                else: # merger remnant ejected from cluster
                        
                    N_BH = N_BH - 2
                        
                    N_meEj+=1
                        
                N_me+=1
                    
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
                
                # merger time and redshift for 2-body GW capture:
                t_merge = t + T_GW(m1, m2, sma, eccen)  # absolute merger time [Myr]
                t_lb = lookback_interp(zCl_form) - t_merge  # lookback time remaining at merger
                z_merge = redshift_interp(t_lb) if t_lb >= 0.1 else 0.0  # merger redshift (0 if within 0.1 Myr of today)

                # append merger:
                mergers.append([seed, ind, 2, sma, eccen, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t, z, t_merge,
                                z_merge, m_rem, s_rem, g_rem, vGW_kick, s_eff, q, v_esc_cl(v_star, vBH), h1, h2])

            else:
                
                # append binary:
                binaries = np.append(binaries, [[ind, 2, sma, eccen, m1, m2, s1, s2, g1, g2, t, z, 0, h1, h2]], axis=0)
                
                N_BBH+=1
                
        for remnant in retained:
            BH.add(*remnant)
        
    return seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, BH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers

# End of file.
