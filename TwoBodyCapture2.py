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
def two_body_capture(seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, mBH, sBH, gBH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers):
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
    @in vBH: 3D BH velocity dispersion
    @in v_star: 3D star velocity dispersion
    @in N_2cap: number of 2-captures
    @in N_BH: number of BHs
    @in N_BBH: number of BBHs
    @in N_me: number of mergers
    @in N_meRe: number of retained mergers
    @in N_meEj: number of ejected mergers
    @in mergers: array of mergers
    
    @out: all inputs
    """
    
    if k_2cap > 0:

        mBH_temp = []
        sBH_temp = []
        gBH_temp = []

        for i in range(k_2cap):
            
            m1, m2 = np.random.choice(mBH, size=2, replace=False, p=mBH**(2)/np.sum(mBH**(2)))
            
            # find index locations of the sampled BHs:
            k1 = np.squeeze(np.where(mBH==m1))+0
            k2 = np.squeeze(np.where(mBH==m2))+0
            
            if isinstance(k1, np.ndarray):
                k1=k1[0]
            if isinstance(k2, np.ndarray):
                k2=k2[0]
                
            s1 = sBH[k1]; g1 = gBH[k1]
            s2 = sBH[k2]; g2 = gBH[k2]
            
            ind = np.random.randint(0, 999999999)
            
            theta1, theta2, dPhi = sample_angles()
            
            m_rem, s_rem, vGW_kick = merger_remnant(m1, m2, sBH[k1], sBH[k2], theta1, theta2, dPhi)
            g_rem = np.max([gBH[k1], gBH[k2]]) + 1
            
            # relative velocity:
            v_rel = get_maxwell_sample(np.sqrt(2/3) * vBH)
            
            # total mass:
            m12 = m1 + m2
            
            # reduced mass:
            mu = m1*m2/m12
            
            # maximum impact parameter for capture:
            b_max = (340 * np.pi / 3)**(1/7) * m12 * mu**(1/7) / v_rel**(9/7) * G_Newton * c_light**(-5/7)

            # impact parameter sampled from uniform in b^2 distribution:
            b = np.sqrt(np.random.rand() * b_max**2)
            
            # pericenter distance:
            rp = b**2 * v_rel**2 / 2 / G_Newton / m12
            
            # GW energy released:
            E_gw = 85 * np.pi / 12 / np.sqrt(2) * mu**2 * m12**(9/2) / rp**(7/2) * G_Newton**(7/2) / c_light**5
            
            # final energy:
            E_fin = mu * v_rel**2 / 2 - E_gw
            
            # make sure eccentricity is strictly smaller than unity:
            while 1 + 2 * E_fin * b**2 * v_rel**2 / m12**3 / mu / G_Newton**2 < 0:
                
                # impact parameter sampled from uniform in b^2 distribution:
                b = np.sqrt(np.random.rand() * b_max**2)
                
                # pericenter distance:
                rp = b**2 * v_rel**2 / 2 / G_Newton / m12
                
                # GW energy released:
                E_gw = 85 * np.pi / 12 / np.sqrt(2) * mu**2 * m12**(9/2) / rp**(7/2) * G_Newton**(7/2) / c_light**5
                
                # final energy:
                E_fin = mu * v_rel**2 / 2 - E_gw
                
            # semimajor axis at formation:
            sma = - G_Newton * m12 * mu / 2 / E_fin
            
            # eccentricity at formation:
            eccen = np.sqrt(1 + 2 * E_fin * b**2 * v_rel**2 / m12**3 / mu / G_Newton**2)
            
            # delete captured BHs:
            mBH = np.delete(mBH, [k1, k2])
            sBH = np.delete(sBH, [k1, k2])
            gBH = np.delete(gBH, [k1, k2])
            
            N_2cap+=1
            
            # check if binary merges within the current step:
            if T_GW(m1, m2, sma, eccen) < np.min([dt, lookback(zCl_form) - t]):
                
                if vGW_kick < 2 * np.sqrt(v_star**2 + vBH**2): # merger remnant retained in cluster
                    
                    mBH_temp.append(m_rem)
                    sBH_temp.append(s_rem)
                    gBH_temp.append(g_rem)
                    
                    N_BH = N_BH - 1
                    
                    N_meRe+=1
                        
                else: # merger remnant ejected from cluster
                        
                    N_BH = N_BH - 2
                        
                    N_meEj+=1
                        
                N_me+=1
                    
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
                mergers = np.append(mergers, [[seed, ind, 2, sma, eccen, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t, z, t + T_GW(m1, m2, sma, eccen),
                                               redshift(lookback(zCl_form) - t + T_GW(m1, m2, sma, eccen)), m_rem, s_rem, g_rem, vGW_kick, s_eff, q]], axis=0)

            else:
                
                # append binary:
                binaries = np.append(binaries, [[ind, 2, sma, eccen, m1, m2, s1, s2, g1, g2, t, z, 0]], axis=0)
                
                N_BBH+=1
                
        mBH_temp = np.array(mBH_temp)
        sBH_temp = np.array(sBH_temp)
        gBH_temp = np.array(gBH_temp)
        
        mBH = np.append(mBH, mBH_temp)
        sBH = np.append(sBH, sBH_temp)
        gBH = np.append(gBH, gBH_temp)
        
    return seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, mBH, sBH, gBH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers

# end of file
