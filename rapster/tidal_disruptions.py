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

def BH_TidalDisruptions(seed, t, z, k_tde, N_tde, tde_type, m_avg, m_star, R_star, mBH, sBH, gBH, hBH, vSTAR, vBH, tdes, binaries, pairs):
    """
    @in seed: seed number of the main simulation
    @in t: current time (Myr)
    @in z: current redshift
    @in k_tde: number of BH TDE occurances in the current step
    @in N_tde: cumulative number of BH TDEs
    @in tde_type: stellar type (integer)
    @in m_avg: average stellar mass (Msun)
    @in m_star: stellar mass (Msun)
    @in R_star: stellar radius (pc)
    @in mBH: single BH masses (Msun)
    @in sBH: single BH dimensionless spins
    @in gBH: single BH generations
    @in hBH: array of BH tdes count
    @in vSTAR: velocity dispersion of stars (km/s)
    @in vBH: velocity dispersion of BHs (km/s)
    @in tdes: tdes array [seed, t, z, type, mstar, Rstar, m_BH, s_BH, g_BH, r_t, r_p, beta, iota, r_mb, dm, s_new, v_rel, h_BH]
    @in binaries: [ind, channel, a, e, m1, m2, s1, s2, g1, g2, t_form, z_form, Nex, h1, h2]
    @in pairs: [a, m, s, g, h]

    @out: all inputs
    """

    if k_tde>0: # perform BH TDE(s)
        
       	for i in range(k_tde):
            
            N_tde+=1 # update number of BH TDEs
            
            # sample single BH mass:
            p = (m_star + mBH)*mBH**(1/3) # weighting probability
            m = np.random.choice(mBH, p=p/np.sum(p))
            
            # find index location of that BH mass:
            k = np.squeeze(np.where(mBH==m))+0
            
            if isinstance(k, np.ndarray):
                k=k[0]
                
            s = sBH[k] # get BH's (dimensionless) spin
            g = gBH[k] # get BH's generation
            h = hBH[k] # get BH's tdes count
            
            # orbital inclination:
            iota = np.arccos(np.random.uniform(-1, 1)) # radians
            
            # marginally bound orbit:
            r_mb = R_mb(m, s, iota)
            
            # tidal radius:
            r_t = R_star * (m/m_star)**(1/3)
            
            # pericenter radius:
            r_p = np.sqrt(np.random.rand())*r_t
            
            # penetration parameter:
            beta = r_t/r_p
            
            f_accreted = 0.5 # fraction of stellar mass accreted
            # mass increment:
            dm = f_accreted*m_star
            
            prograde = +1 if np.cos(iota)>0 else -1
            # final spin:
            s_new = evolve_spin_RungeKutta(m, m+dm, s, prograde, dM=dm/100)
            
            # relative velocity:
            v_rel = np.sqrt(vSTAR**2*m_avg/m_star + np.mean(mBH)/m*vBH**2)
            
            # append tde:
            tdes = np.append(tdes, [[seed, t, z, tde_type, m_star, R_star, m, s, g, r_t, r_p, beta, iota, r_mb, dm, s_new, v_rel, h]], axis=0)
            
            # update BH mass:
            mBH[k] = m + dm
            
            # update BH spin:
            sBH[k] = s_new

            # update BH tdes count:
            hBH[k] += 1

    return seed, t, z, k_tde, N_tde, tde_type, m_avg, m_star, R_star, mBH, sBH, gBH, hBH, vSTAR, vBH, tdes, binaries, pairs

# End of file.
