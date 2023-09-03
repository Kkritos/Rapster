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
from ThreeBodyBinary2 import three_body_binary
from BBHevol2 import evolve_BBHs
from triples2 import evolve_triples
from TwoBodyCapture2 import two_body_capture
from Exchanges2 import StarStar_to_BHstar, BHstar_to_BBH

parser = argparse.ArgumentParser(description="Rapster2 input parameters")

parser.add_argument('-N', '--number', type=float, metavar=' ', default=1e6, help='Initial number of stars')

parser.add_argument('-r', '--half_mass_radius', type=float, metavar=' ', default=1, help='Initial half-mass radius [pc]')

parser.add_argument('-mm', '--minimum_star_mass', type=float, metavar=' ', default=0.08, help='Smallest ZAMS mass [Msun]')

parser.add_argument('-mM', '--maximum_star_mass', type=float, metavar=' ', default=150, help='Largest ZAMS mass [Msun]')

parser.add_argument('-Z', '--metallicity', type=float, metavar=' ', default=0.001, help='Absolute metallicity')

parser.add_argument('-z', '--cluster_formation_redshift', type=float, metavar=' ', default=3, help='Redshift of cluster formation')

parser.add_argument('-n', '--central_stellar_density', type=float, metavar=' ', default=1e6, help='Central stellar number density [pc^-3]')

parser.add_argument('-fb', '--binary_fraction', type=float, metavar=' ', default=0.1, help='Initial binary star fraction')

parser.add_argument('-S', '--seed', type=int, metavar=' ', default=1234567890, help='Seed number')

parser.add_argument('-dtm', '--minimum_time_step', type=float, metavar=' ', default=0.1, help='Minimum simulation time-step [Myr]')

parser.add_argument('-dtM', '--maximum_time_step', type=float, metavar=' ', default=50.0, help='Maximum simulation time-step [Myr]')

parser.add_argument('-tM', '--maximum_time', type=float, metavar=' ', default=14000, help='Maximum simulation time [Myr]')

parser.add_argument('-wK', '--supernova_kick_parameter', type=float, metavar=' ', default=265, help='One-dimensional supernova kick parameter [km/s]')

parser.add_argument('-K', '--natal_kick_prescription', type=int, metavar=' ', default=1, help='Natal kick prescription (0 for fallback, 1 for momentum conservation)')

parser.add_argument('-R', '--galactocentric_radius', type=float, metavar=' ', default=8, help='Initial galactocentric radius [kpc]')

parser.add_argument('-vg', '--galactocentric_velocity', type=float, metavar=' ', default=220.0, help='Galactocentric circular velocity [km/s]')

parser.add_argument('-s', '--spin_parameter', type=float, metavar=' ', default=0.0, help='Natal spin parameter of first generation (1g) BHs')

parser.add_argument('-SD', '--spin_distribution', type=int, metavar=' ', default=0, help='Natal spin distribution model (0 for uniform, 1 for monochromatic)')

parser.add_argument('-P', '--print_information', type=int, metavar=' ', default=1, help='Print runtime information (0 for no, 1 for yes)')

parser.add_argument('-Mi', '--mergers_file_indicator', type=int, metavar=' ', default=1, help='Export mergers file (0 for no, 1 for yes)')

parser.add_argument('-MF', '--mergers_file_name', type=str, metavar=' ', default='mergers', help='Name of .txt output file with BBH merger source parameters')

parser.add_argument('-Ei', '--evolution_file_indicator', type=int, metavar=' ', default=1, help='Export evolution file (0 for no, 1 for yes)')

parser.add_argument('-EF', '--evolution_file_name', type=str, metavar=' ', default='evolution', help='Name of .txt output file with time-dependent quantities')

parser.add_argument('-Hi', '--hardening_file_indicator', type=int, metavar=' ', default=1, help='Export hardening file (0 for no, 1 for yes)')

parser.add_argument('-HF', '--hardening_file_name', type=str, metavar=' ', default='hardening', help='Name of .txt output file with BBH time evolution information')

parser.add_argument('-BIi', '--blackholes_in_file_indicator', type=int, metavar=' ', default=0, help='Use external BH file (0 for no, 1 for yes)')

parser.add_argument('-BIF', '--blackholes_in_file_name', type=str, metavar=' ', default='input_BHs.npz', help='Name of .npz input file with initial BH masses')

parser.add_argument('-BOi', '--blackholes_out_file_indicator', type=int, metavar=' ', default=1, help='Export BH masses file (0 for no, 1 for yes)')

parser.add_argument('-BOF', '--blackholes_out_file_name', type=str, metavar=' ', default='output_BHs.npz', help='Name of .npz file with the masses of all BHs in solar masses')

parser.add_argument('-RP', '--remnant_mass_prescription', type=int, metavar=' ', default=1, help='Remnant mass prescription (0 for SEVN delayed, 1 for Fryer+2012 delayed, 2 for Fryer+2012 rapid)')

args = parser.parse_args()

N = args.number
rh = args.half_mass_radius
m_min = args.minimum_star_mass
m_max = args.maximum_star_mass
Z = args.metallicity
seed = args.seed
dt_min = args.minimum_time_step
dt_max = args.maximum_time_step
t_max = args.maximum_time
wSN_kick = args.supernova_kick_parameter
R_gal = args.galactocentric_radius * 1e3
v_gal = args.galactocentric_velocity
zCl_form = args.cluster_formation_redshift
n_star = args.central_stellar_density
NKP = args.natal_kick_prescription
mergers_file = args.mergers_file_name
evolution_file = args.evolution_file_name
hardening_file = args.hardening_file_name
fb = args.binary_fraction
print_info = args.print_information
Bi = args.blackholes_in_file_indicator
input_BH_file = args.blackholes_in_file_name
s1g_max = args.spin_parameter
SD = args.spin_distribution
Mi = args.mergers_file_indicator
Ei = args.evolution_file_indicator
Hi = args.hardening_file_indicator
BOi = args.blackholes_out_file_indicator
BOF = args.blackholes_out_file_name
RP = args.remnant_mass_prescription

# initial conditions:

if __name__ == "__main__":
    
    # initialize pseudo-random number generator:
    np.random.seed(seed)
    
    # initial half-mass radius:
    rh0 = rh
    
    # initial stellar number density:
    n_star0 = n_star
    
    # average stellar mass:
    m_avg = integrate.quad(lambda x: x * IMF_kroupa(np.array([x])), m_min, m_max)[0] \
    / integrate.quad(lambda x: IMF_kroupa(np.array([x])), m_min, m_max)[0]
    
    # initial average stellar mass:
    m_avg0 = m_avg
    
    # initial cluster mass:
    Mcl = N * m_avg
    Mcl0 = Mcl
    
    # initial multimass relaxation factor:
    psi0 = integrate.quad(lambda x: x**(5/2) * IMF_kroupa(np.array([x])), m_min, m_max)[0] \
        / integrate.quad(lambda x: IMF_kroupa(np.array([x])), m_min, m_max)[0] / m_avg**(5/2)
    
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
    
    if Nb>0: # there are binary stars
        
        # semimajor axes of binaries (loguniform in [ab_min, ab_max]):
        ab = 10**(np.random.uniform(np.log10(ab_min), np.log10(ab_max), Nb))
        
        # stellar velocity dispersion:
        v_star = np.sqrt(0.4 * G_Newton * Mcl / rh)
        
        # binary star hard-soft boundary:
        ab_hard = G_Newton * m_avg / v_star
        
        # filter only hard binary stars:
        ab = ab[ab < ab_hard]
        
    # number of massive stars:
    N_massive = int(Mcl * integrate.quad(lambda x: IMF_kroupa(np.array([x])), mM_min, m_max)[0] \
                    / integrate.quad(lambda x: x * IMF_kroupa(np.array([x])), m_min, m_max)[0])
    
    # massive star masses drawn from Kroupa (2002) IMF:
    u_massive = np.random.rand(N_massive)
    m_massive = (u_massive * (m_max**(alphaIMF + 1) - mM_min**(alphaIMF + 1)) \
                 + mM_min**(alphaIMF + 1))**(1 / (alphaIMF + 1))
    
    # remnant masses:
    m_rem = np.zeros(m_massive.size)
    for i in range(m_massive.size):
        if RP==0:
            m_rem[i] = Mrem_SEVN(m_massive[i], Z) + 0.01 * np.random.rand()
        if RP==1:
            m_rem[i] = Mrem_F12d(m_massive[i], Z) + 0.01 * np.random.rand()
        if RP==2:
            m_rem[i] = Mrem_F12r(m_massive[i], Z) + 0.01 * np.random.rand()
            
    # separate BHs from other remnants:
    mBH = m_rem[m_rem > mBH_min]
    
    # compute supernova kicks:
    if NKP==0:
        vSN_kick = (1 - np.vectorize(f_fb)(m_massive[m_rem > mBH_min])) * maxwell.rvs(loc=0, scale=np.sqrt(3) * wSN_kick, size=mBH.size)
    elif NKP==1:
        vSN_kick = 1.4 / mBH * maxwell.rvs(loc=0, scale=np.sqrt(3) * wSN_kick, size=mBH.size)
        
    # retain BHs with SN kick < escape velocity:
    mBH = mBH[vSN_kick < v_esc(Mcl, rh)]
    
    if Bi==1:
        # load BH masses from file:
        mBH = np.load(input_BH_file)['mBH_ini']
        
        # perturb BH masses to avoid same masses:
        mBH = mBH + 0.01 * np.random.rand(mBH.size)
        
    # BH spins:
    if SD==0:
        # draw spins from uniform distribution:
        sBH = np.random.uniform(0, s1g_max, mBH.size)
    elif SD==1:
        # draw spins from monochromatic distribution:
        sBH = s1g_max * np.ones(mBH.size)
        
    # BH generations:
    gBH = np.ones(mBH.size)
    
    N_BH = 0 # initial number of BHs
    N_BBH = 0 # number of BBHs
    N_3bb = 0 # number of 3bbs formed
    N_me = 0 # total number of mergers
    N_meRe = 0 # number of mergers retained
    N_meEj = 0 # number of mergers ejected
    N_2cap = 0 # number of 2-body captures
    N_3cap = 0 # number of 3-body captures
    N_dis = 0 # number of BBH disruptions
    N_ex = 0 # number of BBH-BH exchanges
    N_BHej = 0 # number of post-interaction BH ejections
    N_BBHej = 0 # number of post-interaction BBH ejections
    N_iter = 0 # number of iterations
    N_bb = 0 # number of BBH-BBH interactions
    N_meFi = 0 # number of field BBH mergers
    N_me2b = 0 # number of 2-body in-cluster mergers
    N_ex1 = 0 # number of star-star -> BH-star exchanges
    N_ex2 = 0 # number of BH-star -> BH-BH exchanges
    N_BHstar = 0 # number of BH-star pairs
    N_pp = 0 # number of BHstar - BHstar interactions
    N_Triples = 0 # number of BHBHBH hierarchical triples
    N_ZLK = 0 # number of ZLK mergers
    
    # binaries [ind, channel, a, e, m1, m2, s1, s2, g1, g2, t_form, z_form, Nex]:
    binaries = np.zeros(shape=(1, 13))
    
    # pairs [a, m, s, g]:
    pairs = np.zeros(shape=(1, 4))
    
    # triples [a_inner, a_outer, e_inner, e_outer, m0, m1, m2, s0, s1, s2, g0, g1, g2, i1, i2, t_form, z_form, ind_inner, channel_inner, t_form_inner, z_form_inner]:
    triples = np.zeros(shape=(1, 21))
    
    # mergers [seed, ind, channel, a, e, m1, m2, s1, s2, g1, g2, theta1, theta2, dPhi, t_form, z_form, t_merge, z_merge, m_rem, s_rem, g_rem, vGW_kick, s_eff, q]:
    mergers = np.zeros(shape=(1, 24))
    
    # evolution [t, z, dt, m_avg, Mcl, rh, R_gal, v_gal, t_rlx, tBH_rlx, n_star, N_BH, mBH_avg, mBH_max, rh_BH, rc_BH, S, xi, psi, psi_BH, t_3bb, t_2cap, k_3bb, k_2cap, N_me, N_BBH, N_meRe, N_meEj, v_star, vBH, nh_BH, nc_BH, na_BH, N_3bb, N_2cap, N_3cap, N_BHej, N_BBHej, N_dis, N_ex, t_bb, N_bb, N_meFi, N_me2b, t_ex1, t_ex2, k_ex1, k_ex2, N_ex1, N_ex2, N_BHstar, t_pp, k_pp, N_pp, v_esc, vBH_esc, N_Triples, N_ZLK]:
    evolution = np.zeros(shape=(1, 58))
    
    # hardening [t, dt, t_local, dt_local, ind, a, e, m1, m2, q, condition, Nex]:
    hardening = np.zeros(shape=(1, 12))
    N_hardening = 0
    
    # initialize time:
    t = 0
    
    # initialize redshift:
    z = zCl_form
    
    # initialize time-step:
    dt = dt_min
    
    i_aux1 = 0
    
    # start global clock:
    global_time_initial = time.time()
    
    mBH_ini = mBH
    
    # Simulation:
    while t<t_max and R_gal>0 and Mcl>0 and mBH.sum()<fBH_max*Mcl:
        
        # start local clock:
        local_time_initial = time.time()
        
        if t>tBH_form and i_aux1==0:
            N_BH = mBH.size
            i_aux1 = 1
            
        if i_aux1==1 and N_BH==0:
            print('BH SUBSYSTEM EVAPORATED')
            break
        
        # average BH mass:
        if i_aux1==1:
            mBH_avg = (np.sum(mBH) + np.sum(np.transpose(binaries)[:][4]) + np.sum(np.transpose(binaries)[:][5]) + np.sum(np.transpose(pairs)[:][1]) ) / N_BH
        else:
            mBH_avg = 0.0
            
        # maximum BH mass:
        if i_aux1==1:
            try:
                mBHs_max = np.max(mBH)
            except:
                mBHs_max = 0.0
            try:
                mBH1_max = np.max(np.transpose(binaries)[:][4])
            except:
                mBH1_max = 0.0
            try:
                mBH2_max = np.max(np.transpose(binaries)[:][5])
            except:
                mBH2_max = 0.0
            try:
                mBHp_max = np.max(np.transpose(pairs)[:][1])
            except:
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
        if S < 0.16 and S > 0: # equipartition is achieved
            xi = 1
            S = Q_BH * q_BH
        else: # equipartition is not achieved
            # equipartition parameter:
            if i_aux1==1:
                xi = q_BH**(3/5) * Q_BH**(2/5) * (logLBH / logLcl)**(-2/5)
            else:
                xi = 0.0
                
        # multimass relaxation factor:
        psi = 1 + S
        
        # BH relaxation factor:
        if i_aux1==1 and N_BH>0:
            psi_BH = (np.sum(mBH**(5/2)) + np.sum(np.transpose(binaries)[:][4]**(5/2)) + \
                      np.sum(np.transpose(binaries)[:][5]**(5/2)) + np.sum(np.transpose(pairs)[:][1]**(5/2)) ) \
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
        
        # half-mass radius:
        Vh = 4 * np.pi / 3 * rh**3
        
        if Nb>0:
            # binary star number density:
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
            
        # 3bb timescale:
        if t>t_cc and N_BH-2*N_BBH-N_BHstar-3*N_Triples>2:
            try:
                t_3bb = 1 / Rate_3bb(mBH_avg, nc_BH, vBH) / Vc_BH
            except:
                t_3bb = 1e100
        else:
            t_3bb = 1e100
            
        # 2-body capture timescale:
        if t>t_cc and N_BH-2*N_BBH-N_BHstar-3*N_Triples>1:
            try:
                t_2cap = 1 / Rate_cap(mBH_avg, nc_BH, vBH) / Vc_BH
            except:
                t_2cap = 1e100
        else:
            t_2cap = 1e100
            
        # star-star -> BH-star timescale:
        if t>t_cc and N_BH-2*N_BBH-N_BHstar-3*N_Triples>0 and nb>0:
            try:
                t_ex1 = 1 / Rate_exc(m_avg, m_avg, mBH_avg, nb, np.sqrt(v_star**2 + vBH**2), np.mean(ab)) / Nc_BH / 2
            except:
                t_ex1 = 1e100
        else:
            t_ex1 = 1e100
            
        # BH-star -> BH-BH timescale:
        if t>t_cc and N_BHstar>0 and N_BH-2*N_BBH-N_BHstar-3*N_Triples>0:
            if N_BHstar>0:
                t_ex2 = 1 / Rate_exc(m_avg, np.mean(np.transpose(pairs)[:][1]), mBH_avg, nc_BH, vBH, np.mean(np.transpose(pairs)[:][0])) / N_BHstar
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
            
        # Binary BH formation:
        
        # number of 3bbs (cannot exceed 3*N_BHsin):
        k_3bb = np.min([poisson.rvs(mu=dt / t_3bb), int((N_BH-2*N_BBH-N_BHstar-3*N_Triples) / 3)])
        
        # 3bb formation:
        t, z, k_3bb, mBH_avg, binaries, mBH, sBH, gBH, vBH, N_3bb, N_BBH = three_body_binary(t, z, k_3bb, mBH_avg, binaries, mBH, sBH, gBH, vBH, N_3bb, N_BBH)
        
        # number of 2-body captures (cannot exceed 2*N_BHsin):
        k_2cap = np.min([poisson.rvs(mu=dt / t_2cap), int((N_BH-2*N_BBH-N_BHstar-3*N_Triples) / 2)])
        
        # 2-body capture(s):
        seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, mBH, sBH, gBH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers = \
            two_body_capture(seed, t, dt, z, zCl_form, k_2cap, mBH_avg, binaries, mBH, sBH, gBH, vBH, v_star, N_2cap, N_BH, N_BBH, N_me, N_meRe, N_meEj, mergers)
        
        # number of star-star -> BH-star exchanges (cannot exceed N_BHsin)
        k_ex1 = np.min([poisson.rvs(mu=dt / t_ex1), int(N_BH-2*N_BBH-N_BHstar-3*N_Triples)])
        
        # star-star -> BH-star exchange(s):
        if k_ex1 > 0:
            k_ex1, N_ex1, m_avg, mBH, sBH, gBH, ab, pairs, N_BHstar = StarStar_to_BHstar(k_ex1, N_ex1, m_avg, mBH, sBH, gBH, ab, pairs, N_BHstar)
            
        # number of BH-star -> BH-BH exchanges (cannot exceed N_BHsin or N_BHstar)
        k_ex2 = np.min([poisson.rvs(mu=dt / t_ex2), int(N_BH-2*N_BBH-N_BHstar-3*N_Triples), int(N_BHstar)])
        
        # BH-star -> BBH exchange(s):
        if k_ex2 > 0:
            t, z, k_ex2, N_ex2, m_avg, mBH, sBH, gBH, pairs, binaries, N_BBH, N_BHstar = BHstar_to_BBH(t, z, k_ex2, N_ex2, m_avg, mBH, sBH, gBH, pairs, binaries, N_BBH, N_BHstar)
            
        # BBH evolution:
        seed, t, z, dt, zCl_form, binaries, hardening, mergers, mBH, sBH, gBH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples = evolve_BBHs(seed, t, z, dt, zCl_form, binaries, hardening, mergers, mBH, sBH, gBH, n_star, v_star, vBH, t_rlx, m_avg, mBH_avg, na_BH, nc_BH, N_BH, N_BBH, N_me, N_me2b, N_3cap, N_meFi, N_meRe, N_meEj, N_dis, N_ex, N_BHej, N_BBHej, N_hardening, Vc_BH, N_bb, triples, N_Triples)
        
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
        if t>t_cc and N_BBH>1: # need at least two BBHs
            try:
                t_bb = 1 / Rate_int(np.mean(np.transpose(binaries)[:][4]+np.transpose(binaries)[:][5]), n_BBH, vBH, 2 * np.mean(np.transpose(binaries)[:][2])) / N_BBH
            except:
                t_bb = 1e100
        else:
            t_bb = 1e100
        
        # Pair-pair interaction timescale:
        if t>t_cc and N_BHstar>1: # need at least two BH-stars
            try:
                t_pp = 1 / Rate_int(np.mean(np.transpose(pairs)[:][1]) + m_avg, n_BHstar, vBH, 2 * np.mean(np.transpose(pairs)[:][0])) / N_BHstar
            except:
                t_pp = 1e100
        else:
            t_pp = 1e100

        # number of pair-pair interactions (cannot exceed 2*N_BHstar):
        k_pp = np.min([poisson.rvs(mu=dt / t_pp), int(N_BHstar / 2)])
        
        if k_pp>0:
            
            for i in range(k_pp):
                
                N_pp+=1
                
                smas = np.transpose(pairs)[:][0]
                masses = np.transpose(pairs)[:][1]
                
                a1, a2 = np.random.choice(smas, size=2, replace=False, p=smas*masses / np.sum(smas*masses))
                
                k1 = np.squeeze(np.where(smas==a1))+0
                k2 = np.squeeze(np.where(smas==a2))+0
                
                if isinstance(k1, np.ndarray):
                    k1=k1[0]
                if isinstance(k2, np.ndarray):
                    k2=k2[0]
                    
                m1 = pairs[k1][1]; s1 = pairs[k1][2]; g1 = pairs[k1][3]
                m2 = pairs[k2][1]; s2 = pairs[k2][2]; g2 = pairs[k2][3]
                
                if m2/a2 < m1/a1: # pair 1 is softer and controls the energy of the new BBH
                    sma = m2 / m_avg * a1
                else:
                    sma = m1 / m_avg * a2
                    
                eccen = np.sqrt(np.random.rand())
                
                # append binary:
                binaries = np.append(binaries, [[np.random.randint(0, 999999999), 1, sma, eccen, m1, m2, s1, s2, g1, g2, t, z, 0]], axis=0)
                
                # destroy pairs:
                pairs = np.delete(pairs, [k1, k2], axis=0)
                
                N_BHstar = N_BHstar - 2
                
                N_BBH+=1
                
        # Triple evolution:
        seed, t, z, zCl_form, triples, binaries, mBH, sBH, gBH, mBH_avg, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, v_star, vBH, nc_BH, mergers = \
            evolve_triples(seed, t, z, zCl_form, triples, binaries, mBH, sBH, gBH, mBH_avg, N_Triples, N_BBH, N_BH, N_me, N_meRe, N_meEj, N_ZLK, v_star, vBH, nc_BH, mergers)
        
        # external parameters:
        
        # Jacobi radius:
        rJ = (G_Newton * Mcl * R_gal**2 / 3 / v_gal**2)**(1/3)
        
        # dimensionless escape rate:
        xi_e = xi_e0 * np.exp(10 * rh / rJ)
        
        # dynamical friction timescale:
        t_df = 0.45e3 * (R_gal / 1e3)**2 * v_gal / (Mcl / 1e5) / 2
        
        # append evolution:
        evolution = np.append(evolution, [[t, z, dt, m_avg, Mcl, rh, R_gal, v_gal, t_rlx, tBH_rlx, n_star, N_BH, mBH_avg, mBH_max, rh_BH, rc_BH, S, 
                                           xi, psi, psi_BH, t_3bb, t_2cap, k_3bb, k_2cap, N_me, N_BBH, N_meRe, N_meEj, v_star, vBH, 
                                           nh_BH, nc_BH, na_BH, N_3bb, N_2cap, N_3cap, N_BHej, N_BBHej, N_dis, N_ex, t_bb, N_bb, 
                                           N_meFi, N_me2b, t_ex1, t_ex2, k_ex1, k_ex2, N_ex1, N_ex2, N_BHstar, t_pp, k_pp, N_pp, 2*v_star, 
                                           2*vBH, N_Triples, N_ZLK]], axis=0)
        
        # Cluster evolution:
        
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
            break
        
        # half-mass radius update:
        rh = rh + drh
        
        # galactocentric radius update:
        R_gal = R_gal + dR_gal
        
        if R_gal<0:
            print('CLUSTER REACHED GALAXY CENTER')
            break
        
        # simulation time update:
        t+=dt
        
        # redshift update:
        try:
            z = redshift_interp(lookback_interp(zCl_form) - t)
        except:
            print('REDSHIFT 0 REACHED')
            break
        
        N_iter+=1
        
        local_time_final = time.time()
        
        if print_info==1:
            
            # Print runtime information:
            print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
            print("Iteration #", N_iter)
            frmt_1 = '%.2f'
            frmt_2 = '%.0f'
            frmt_3 = '%.3f'
            frmt_4 = "%.1f"
            data_1 = {"t[Myr]": [frmt_1%t], "dt[Myr]": [frmt_1%dt], "z": [frmt_1%z], "Mcl[MMsun]": [frmt_1%(Mcl/1e6)], "rh[pc]": [frmt_1%rh], "R_gal[kpc]": [frmt_1%(R_gal/1e3)]}
            data_2 = {"N_BH": [frmt_2%N_BH], "N_BBH": [frmt_2%N_BBH], "N_Triples": [frmt_2%N_Triples], "N_me": [frmt_2%N_me]}
            data_3 = {"steptime[ms]": [frmt_4%(np.abs(local_time_final - local_time_initial)*1e3)], "runtime[s]": [frmt_3%np.abs(time.time() - global_time_initial)]}
            headers = [" "]
            df_1 = pd.DataFrame(data_1, headers)
            df_2 = pd.DataFrame(data_2, headers)
            df_3 = pd.DataFrame(data_3, headers)
            print(df_1)
            print(df_2)
            print(df_3)
            print("\n")
            
    # ----------------------------------------------------------------------------------------------------------------------
    
    mergers = np.delete(mergers, 0, axis=0)
    evolution = np.delete(evolution, 0, axis=0)
    hardening = np.delete(hardening, 0, axis=0)
    
    # exporting output files:
    
    mBH_fin = mBH
    
    if BOi==1:
        np.savez(BOF, mBH_ini=mBH_ini, mBH_fin=mBH_fin)
    
    if Mi==1: # export mergers file
        with open(mergers_file+'.txt', 'w') as f_mergers:
            for i in range(N_me):
                f_mergers.write(str(mergers[i][0 ])+' '+str(mergers[i][1 ])+' '+str(mergers[i][2 ])+' '+str(mergers[i][3 ])+' '+str(mergers[i][4 ])+' '+str(mergers[i][5 ])+' '+str(mergers[i][6 ])+' '+\
                                str(mergers[i][7 ])+' '+str(mergers[i][8 ])+' '+str(mergers[i][9 ])+' '+str(mergers[i][10])+' '+str(mergers[i][11])+' '+str(mergers[i][12])+' '+str(mergers[i][13])+' '+\
                                str(mergers[i][14])+' '+str(mergers[i][15])+' '+str(mergers[i][16])+' '+str(mergers[i][17])+' '+str(mergers[i][18])+' '+str(mergers[i][19])+' '+str(mergers[i][20])+' '+\
                                str(mergers[i][21])+' '+str(mergers[i][22])+' '+str(mergers[i][23]))
                f_mergers.write('\n')

    if Ei==1: # export evolution file
        with open(evolution_file+'.txt', 'w') as f_evolution:
            for i in range(N_iter):
                f_evolution.write(str(evolution[i][0 ])+' '+str(evolution[i][1 ])+' '+str(evolution[i][2 ])+' '+str(evolution[i][3 ])+' '+str(evolution[i][4 ])+' '+str(evolution[i][5 ])+' '+\
                                  str(evolution[i][6 ])+' '+str(evolution[i][7 ])+' '+str(evolution[i][8 ])+' '+str(evolution[i][9 ])+' '+str(evolution[i][10])+' '+str(evolution[i][11])+' '+\
                                  str(evolution[i][12])+' '+str(evolution[i][13])+' '+str(evolution[i][14])+' '+str(evolution[i][15])+' '+str(evolution[i][16])+' '+str(evolution[i][17])+' '+\
                                  str(evolution[i][18])+' '+str(evolution[i][19])+' '+str(evolution[i][20])+' '+str(evolution[i][21])+' '+str(evolution[i][22])+' '+str(evolution[i][23])+' '+\
                                  str(evolution[i][24])+' '+str(evolution[i][25])+' '+str(evolution[i][26])+' '+str(evolution[i][27])+' '+str(evolution[i][28])+' '+str(evolution[i][29])+' '+\
                                  str(evolution[i][30])+' '+str(evolution[i][31])+' '+str(evolution[i][32])+' '+str(evolution[i][33])+' '+str(evolution[i][34])+' '+str(evolution[i][35])+' '+\
                                  str(evolution[i][36])+' '+str(evolution[i][37])+' '+str(evolution[i][38])+' '+str(evolution[i][39])+' '+str(evolution[i][40])+' '+str(evolution[i][41])+' '+\
                                  str(evolution[i][42])+' '+str(evolution[i][43])+' '+str(evolution[i][44])+' '+str(evolution[i][45])+' '+str(evolution[i][46])+' '+str(evolution[i][47])+' '+\
                                  str(evolution[i][48])+' '+str(evolution[i][49])+' '+str(evolution[i][50])+' '+str(evolution[i][51])+' '+str(evolution[i][52])+' '+str(evolution[i][53])+' '+\
                                  str(evolution[i][54])+' '+str(evolution[i][55])+' '+str(evolution[i][56])+' '+str(evolution[i][57]))
                f_evolution.write('\n')
            
    if Hi==1: # export hardening file
        with open(hardening_file+'.txt', 'w') as f_hardening:
            for i in range(N_hardening):
                f_hardening.write(str(hardening[i][0 ])+' '+str(hardening[i][1 ])+' '+str(hardening[i][2 ])+' '+str(hardening[i][3 ])+' '+str(hardening[i][4 ])+' '+\
                                  str(hardening[i][5 ])+' '+str(hardening[i][6 ])+' '+str(hardening[i][7 ])+' '+str(hardening[i][8 ])+' '+str(hardening[i][9 ])+' '+\
                                  str(hardening[i][10])+' '+str(hardening[i][11]))
                f_hardening.write('\n')
    
    print('END OF SIMULATION. RUNTIME:', np.abs(time.time() - global_time_initial), 's')
    print('\n')
    
# end of file
