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

def p_2capture(m1, m2):
    """
    Joint probability density function for masses m1 and m2 to be captured.
    The result is not normalized.
    """
    
    return (m1 + m2)**(10/7) * m1**(2/7) * m2**(2/7)
p_2capture = np.vectorize(p_2capture)

def two_body_capture(seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, mBH, sBH, gBH, hBH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers, random_pairing=False):
    """
    @in seed: simulation seed number
    @in t: simulation time
    @in dt: simulation time step
    @in z: simulation redshift
    @in zCl_form: cluster formation redshift
    @in k_2cap: number of 2-captures in current step
    @in mBH_avg: average BH mass
    @in binaries: array of BBHs
    @in mBH: array of single BH masses
    @in sBH: array of single BH spins
    @in gBH: array of single BH generations
    @in hBH: array of BH tdes count
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

        mBH_temp = []
        sBH_temp = []
        gBH_temp = []
        hBH_temp = []

        for i in range(k_2cap):
            
            # sample the masses that form the captured binary:
            if random_pairing:
                m1, m2 = np.random.choice(mBH, size=2, replace=False)
            else:
                p1 = np.array([p_2capture(m_1, mBH[mBH!=m_1]).sum() for m_1 in mBH]) # marginalized probability
                p1 /= p1.sum() # normalize margninalized probability
                m1 = np.random.choice(mBH, size=1, replace=False, p=p1)[0] # sample first mass
                p2 = p_2capture(m1, mBH[mBH!=m1]) # conditional probability
                p2 /= p2.sum() # normalize conditional probability
                m2 = np.random.choice(mBH[mBH!=m1], size=1, replace=False, p=p2)[0] # sample second mass
            
            # find index locations of the sampled BHs:
            k1 = np.squeeze(np.where(mBH==m1))+0
            k2 = np.squeeze(np.where(mBH==m2))+0
            
            if isinstance(k1, np.ndarray):
                k1=k1[0]
            if isinstance(k2, np.ndarray):
                k2=k2[0]
                
            s1 = sBH[k1]; g1 = gBH[k1]; h1 = hBH[k1]
            s2 = sBH[k2]; g2 = gBH[k2]; h2 = hBH[k2]
            
            ind = np.random.randint(0, 999999999)
            
            theta1, theta2, dPhi = sample_angles()
            
            m_rem, s_rem, vGW_kick = merger_remnant(m1, m2, sBH[k1], sBH[k2], theta1, theta2, dPhi)
            g_rem = np.max([gBH[k1], gBH[k2]]) + 1
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
            
            # make sure eccentricity is strictly smaller than unity:
            while 1 + 2 * E_fin * b**2 * v_rel**2 / m12**2 / mu / G_Newton**2 < 0:
                
                # impact parameter sampled from uniform in b^2 distribution:
                b = np.sqrt(np.random.rand() * b_max**2)
                
                # pericenter distance:
                rp = b**2 * v_rel**2 / 2 / G_Newton / m12
                
                # GW energy released:
                E_gw = 85 * np.pi / 12 / np.sqrt(2) * mu**2 * m12**(5/2) / rp**(7/2) * G_Newton**(7/2) / c_light**5
                
                # final energy:
                E_fin = mu * v_rel**2 / 2 - E_gw
                
            # semimajor axis at formation:
            sma = - G_Newton * m12 * mu / 2 / E_fin
            
            # eccentricity at formation:
            eccen = np.sqrt(1 + 2 * E_fin * b**2 * v_rel**2 / m12**2 / mu / G_Newton**2)
            
            # delete captured BHs:
            mBH = np.delete(mBH, [k1, k2])
            sBH = np.delete(sBH, [k1, k2])
            gBH = np.delete(gBH, [k1, k2])
            hBH = np.delete(hBH, [k1, k2])
            
            N_2cap+=1
            
            # check if binary merges within the current step:
            if T_GW(m1, m2, sma, eccen) < np.min([dt, lookback(zCl_form) - t]):
                
                if vGW_kick < 2 * np.sqrt(v_star**2 + vBH**2): # merger remnant retained in cluster
                    
                    mBH_temp.append(m_rem)
                    sBH_temp.append(s_rem)
                    gBH_temp.append(g_rem)
                    hBH_temp.append(h_rem)
                    
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
                
                # append merger:
                mergers = np.append(mergers, [[seed, ind, 2, sma, eccen, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t, z, t + T_GW(m1, m2, sma, eccen),
                                               redshift(lookback(zCl_form) - t + T_GW(m1, m2, sma, eccen)), m_rem, s_rem, g_rem, vGW_kick, s_eff, q, 2*v_star, h1, h2]], axis=0)

            else:
                
                # append binary:
                binaries = np.append(binaries, [[ind, 2, sma, eccen, m1, m2, s1, s2, g1, g2, t, z, 0, h1, h2]], axis=0)
                
                N_BBH+=1
                
        mBH_temp = np.array(mBH_temp)
        sBH_temp = np.array(sBH_temp)
        gBH_temp = np.array(gBH_temp)
        hBH_temp = np.array(hBH_temp)

        mBH = np.append(mBH, mBH_temp)
        sBH = np.append(sBH, sBH_temp)
        gBH = np.append(gBH, gBH_temp)
        hBH = np.append(hBH, hBH_temp)
        
    return seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, mBH, sBH, gBH, hBH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers

# End of file.
