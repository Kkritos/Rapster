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
def evolve_triples(seed, t, z, zCl_form, triples, binaries, mBH, sBH, gBH, mBH_avg, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, v_star, vBH, nc_BH, mergers):
    """
    @in seed: simulation seed number
    @in t: simulation time
    @in z: redshift
    @in zCl_form: cluster formation redshift
    @in triples: array of triples
    @in binaries: array of BBHs
    @in mBH: array of single BH masses
    @in sBH: array of single BH spins
    @in gBH: array of single BH generations
    @in mBH_avg: average BH mass
    @in N_Triples: number of triples
    @in N_BBH: number of BBHs
    @in N_BH: number of BHs
    @in N_me: number of mergers
    @in N_meRe: number of retained mergers
    @in N_meEj: number of ejected mergers
    @in N_ZLK: number of ZLK mergers
    @in v_star: 3D star velocity dispersion
    @in vBH: 3D BH velocity dispersion
    @in nc_BH: central BH density
    @in mergers: array of mergers
    
    @out: all inputs
    """
    
    if N_Triples>0:
        
        # initialize iteration index:
        i=1
        
        # iterate over all available triples:
        while i < N_Triples+1:
            
            # Unwrap current triples's parameters:
            a_in = triples[i][0]
            a_out = triples[i][1]
            e_in = triples[i][2]
            e_out = triples[i][3]
            m0 = triples[i][4]
            m1 = triples[i][5]
            m2 = triples[i][6]
            s0 = triples[i][7]
            s1 = triples[i][8]
            s2 = triples[i][9]
            g0 = triples[i][10]
            g1 = triples[i][11]
            g2 = triples[i][12]
            inclination1 = triples[i][13]
            inclination2 = triples[i][14]
            ind_in = triples[i][17]
            channel_in = triples[i][18]
            t_form_in = triples[i][19]
            z_form_in = triples[i][20]
            
            # inner pair mass:
            m01 = m0 + m1
            
            # total triple mass:
            m012 = m01 + m2
            
            # epsilon parameter of inner pair:
            epsilon = 1 - e_in**2
            
            # reduced mass of inner pair:
            mu_in = m0 * m1 / m01
            
            # reduced mass of outer pair:
            mu_out = m01 * m2 / m012
            
            # auxiliary parameters:
            beta = mu_out * np.sqrt(m012) / mu_in / np.sqrt(m01) * np.sqrt(a_out / a_in * (1 - e_out**2))
            alpha = np.sqrt(epsilon) * np.cos(inclination1) + beta * np.cos(inclination2)
            
            # total inclination:
            inclination = inclination1 + inclination2
            
            # cosine of total inclination:
            cos_inclination = (alpha**2 - beta**2 - epsilon) / 2 / beta / np.sqrt(epsilon)
            
            # proxy for maximum eccentricity of inner pair:
            epsilon_min = 5/3 * np.cos(inclination)**2
            
            # ZLK merger timescale:
            t_ZLK = 9e32 / m01**2 / mu_in * a_in**4 * epsilon_min**3
            
            # triple - single interaction timescale:
            t_3 = 1 / Rate_int(m012 + mBH_avg, nc_BH, vBH, kp_max * a_out)
            
            # condition for GR precession:
            ZLK_not_Destroyed_By_GR_precession = a_out / a_in < 34 * (a_in / 1e-2)**(1/3) * (m01 / 2e6)**(-1/3) * (2 * m2 / m01)**(1/3) * ((1 - e_in**2) / (1 - e_out**2))**(1/2)
            
            # check whether inner pair merges before next interaction with another object and GR precession does not destroy ZLK:
            if t_ZLK < t_3 and ZLK_not_Destroyed_By_GR_precession: # ZLK merger occurs
                
                # sample spin orientations:
                theta0, theta1, dPhi = sample_angles()
                
                # merger remnant properties:
                m_rem, s_rem, vGW_kick = merger_remnant(m0, m1, s0, s1, theta0, theta1, dPhi)
                g_rem = np.max([g1, g2]) + 1
                
                N_me+=1
                N_ZLK+=1
                
                # merger time:
                t_merge = t + T_GW(m0, m1, a_in, e_in)
                
                # order BHs by mass:
                mA = m0; sA = s0; gA = g0; thetaA = theta0
                mB = m1; sB = s1; gB = g1; thetaB = theta1
                if mA>mB:
                    m0 = mA; s0 = sA; g0 = gA; theta0 = thetaA
                    m1 = mB; s1 = sB; g1 = gB; theta1 = thetaB
                else:
                    m0 = mB; s0 = sB; g0 = gB; theta0 = thetaB
                    m1 = mA; s1 = sA; g1 = gA; theta1 = thetaA
                    
                # mass ratio:
                q = m1 / m0
                
                # effective spin parameter:
                s_eff = (m0 * s0 * np.cos(theta0) + m1 * s1 * np.cos(theta1)) / (m0 + m1)
                
                # append merger:
                mergers = np.append(mergers, [[seed, ind_in, 4, e_in, e_in, m0, m1, s0, s1, g0, g1, theta0, theta1, dPhi, t_form_in, z_form_in, t_merge,
                                               redshift(lookback(zCl_form) - t_merge), m_rem, s_rem, g_rem, vGW_kick, s_eff, q]], axis=0)
                
                triples = np.delete(triples, i, axis=0)

                i = i - 1
                
                N_Triples = N_Triples - 1
                
                N_BH = N_BH - 1
                
                if vGW_kick > np.sqrt(G_Newton * (m_rem + m2) / a_out): # new outer binary disrupts due to merger kick
                    
                    mBH = np.append(mBH, m2)
                    sBH = np.append(sBH, s2)
                    gBH = np.append(gBH, g2)
                    
                    if vGW_kick > 2 * np.sqrt(v_star**2 + vBH**2): # merger remnant is ejected

                        N_BH = N_BH - 1

                        N_meEj+=1

                    else: # merger remnant is retained
                        
                        N_meRe+=1
                        
                        mBH = np.append(mBH, m_rem)
                        sBH = np.append(sBH, s_rem)
                        gBH = np.append(gBH, g_rem)
                        
                else: # new outer binary survives
                    
                    # append binary:
                    binaries = np.append(binaries, [[np.random.randint(0, 999999999), 5, a_out, np.sqrt(np.random.rand()), m2, m_rem, s2, s_rem, g2, g_rem, t + t_ZLK, redshift(lookback(zCl_form) - t - t_ZLK), 0]], axis=0)
                    
                    N_BBH+=1
                    N_meRe+=1

            else: # triple breaks and inner binary is released
                
                triples = np.delete(triples, i, axis=0)
                
                i = i - 1
                
                N_Triples = N_Triples - 1
                
                N_BBH+=1
                
                # release outer member into the cluster:
                mBH = np.append(mBH, m2)
                sBH = np.append(sBH, s2)
                gBH = np.append(gBH, g2)
                
                binaries = np.append(binaries, [[ind_in, channel_in, a_in, e_in, m0, m1, s0, s1, g0, g1, t, z, 0]], axis=0)
                
            i+=1
            
    return seed, t, z, zCl_form, triples, binaries, mBH, sBH, gBH, mBH_avg, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, v_star, vBH, nc_BH, mergers

# end of file
