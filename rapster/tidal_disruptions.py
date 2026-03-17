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

from constants import *
from functions import *

# Marginally bound orbit:
# =======================

def x_mb(a_BH, iota):
    """
    Marginally bound orbit normalized to gravitational radius.
    """
    def fun(x):
        return x**4 - 4*x**3 - a_BH**2 * (1 - 3*np.cos(iota)**2)*x**2 + a_BH**4 * np.cos(iota)**2 + 4 * a_BH * np.sin(iota) * np.sqrt(x**5 - a_BH**2 * x**3 * np.cos(iota)**2)
        
    x = fsolve(fun, 100000)
    return x
x_mb = np.vectorize(x_mb)

# bulk of domain:
N_points_bulk = 10**3

a_BH_vec = np.random.uniform(-1, 1, N_points_bulk)
iota_vec = np.random.uniform(0, np.pi/2, N_points_bulk)

# boundary of parameter space:
N_points_boundary = 10**3

a_BH_vec = np.concatenate([a_BH_vec, np.ones(N_points_boundary)])
iota_vec = np.concatenate([iota_vec, np.random.uniform(0, np.pi/2, N_points_boundary)])

a_BH_vec = np.concatenate([a_BH_vec, -np.ones(N_points_boundary)])
iota_vec = np.concatenate([iota_vec, np.random.uniform(0, np.pi/2, N_points_boundary)])

a_BH_vec = np.concatenate([a_BH_vec, np.random.uniform(-1, 1, N_points_boundary)])
iota_vec = np.concatenate([iota_vec, np.pi/2*np.ones(N_points_boundary)])

a_BH_vec = np.concatenate([a_BH_vec, np.random.uniform(-1, 1, N_points_boundary)])
iota_vec = np.concatenate([iota_vec, np.zeros(N_points_boundary)])

# corners of domain:
a_BH_vec = np.concatenate([a_BH_vec, np.array([-1, -1, 1, 1])])
iota_vec = np.concatenate([iota_vec, np.array([0, np.pi/2, 0, np.pi/2])])

# compute values in the scatter grid:
x_mb_vec = x_mb(a_BH_vec, iota_vec)

# interpolate auxiliary function:
x_mb_interp = CloughTocher2DInterpolator(np.c_[a_BH_vec, iota_vec], x_mb_vec)

def R_mb(M_BH, a_BH, iota):
    """
    Marginally bound orbit, from interpolated function above.
    """
    
    rg = G_Newton*M_BH/c_light**2 # gravitational radius
    rmb = x_mb_interp(a_BH, iota) * rg
    return rmb
R_mb = np.vectorize(R_mb)

def BH_TidalDisruptions(seed, t, z, k_tde, N_tde, type, m_star, R_star, mBH, sBH, gBH, vSTAR, vBH, tdes, binaries, pairs):
    """
    @in seed: seed number of the main simulation
    @in t: current time (Myr)
    @in z: current redshift
    @in k_tde: number of BH TDE occurances in the current step
    @in N_tde: cumulative number of BH TDEs
    @in type: stellar type (integer)
    @in m_star: stellar mass (Msun)
    @in R_star: stellar radius (pc)
    @in mBH: single BH masses (Msun)
    @in sBH: single BH dimensionless spins
    @in gBH: single BH generations
    @in vSTAR: velocity dispersion of stars (km/s)
    @in vBH: velocity dispersion of BHs (km/s)
    @in tdes: tdes array [seed, t, z, type, m_star, R_star, m_BH, s_BH, g_BH, r_t, r_p, beta, iota, r_mb, dm, ds, t_fb, eta_R, L_pk, v_rel]
    @in binaries: [ind, channel, a, e, m1, m2, s1, s2, g1, g2, t_form, z_form, Nex]
    @in pairs: [a, m, s, g]

    @out: all inputs
    """

    if k_tde>0: # perform BH TDE(s)
        
       	for i in range(k_tde):
            
            N_tde+=1 # update number of BH TDEs
            
            # sample single BH mass:
            p = mBH**(4/3) # weighting probability
            m = np.random.choice(mBH, p=p/np.sum(p))
            
            # find index location of that BH mass:
            k = np.squeeze(np.where(mBH==m))+0
            
            if isinstance(k, np.ndarray):
                k=k[0]
                
            s = sBH[k] # get BH's (dimensionless) spin
            g = gBH[k] # get BH's generation
            
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
            # spin change:
            ds = s - evolve_spin_RungeKutta(m, m+dm, s, prograde, dM=dm/100)
            
            ENERGY = G_Newton*m*R_star/r_p**2 # energy of most-bound debris
            t_fb = 2*np.pi*G_Newton*(ENERGY)**(-3/2) # fallback time
            eta_R = 2*ENERGY/c_light**2 # radiative efficiency
            L_pk = eta_R*dm*c_light**2/t_fb # peak luminosity
            
            # relative velocity:
            v_rel = np.sqrt(vSTAR**2 + np.mean(mBH)/m*vBH**2)
            
            # append tde:
            tdes = np.append(tdes, [[seed, t, z, type, m_star, R_star, m, s, g, r_t, r_p, beta, iota, r_mb, dm, ds, t_fb, eta_R, L_pk, v_rel]], axis=0)
            
            # update BH mass:
            mBH[k] = m + dm
            
            # update BH spin:
            sBH[k] = s + ds

    return seed, t, z, k_tde, N_tde, type, m_star, R_star, mBH, sBH, gBH, vSTAR, vBH, tdes, binaries, pairs

# End of file.
