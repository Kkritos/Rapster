'''
 Copyright (C) 2022  Konstantinos Kritos <kkritos1@jhu.edu>
 
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

# Imports and global constants
# ----------------------------------------------------------------------------------------------------------------------------

from functions import *

# Hydrogen mass (kg):
mH = 1.67e-27

# centimeter (m):
cm = 1e-2

# Boltzmann constant (SI):
kB = 1.38e-23

# proton mass (kg):
m_proton = mH

# Thomson scattering cross section (m^2):
sThomson = 6.65e-29

# User inputs
# ----------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Cluster code input parameters')

parser.add_argument('-Mcl' ,'--ClusterMass'          ,type=float,metavar=' ',default=1e6       ,\
    help='Initial cluster mass in solar masses (Msun)')

parser.add_argument('-rh'  ,'--HalfMassRadius'       ,type=float,metavar=' ',default=1         ,\
    help='Initial half-mass radius in parsec (pc)')

parser.add_argument('-rhoC','--CentralDensity'       ,type=float,metavar=' ',default=4e5       ,\
    help='Initial central star density in Msun/pc^3')

parser.add_argument('-Rgal','--GalactocentricRadius' ,type=float,metavar=' ',default=8         ,\
    help='Initial galactocentric radius in kpc')

parser.add_argument('-Z'   ,'--Metallicity'          ,type=float,metavar=' ',default=.1        ,\
    help='Cluster metallicity in units of solar metallicity')

parser.add_argument('-fB'  ,'--BinaryFraction'       ,type=float,metavar=' ',default=.1        ,\
    help='Initial binary star fraction')

parser.add_argument('-w'   ,'--NatalKickParameter'   ,type=float,metavar=' ',default=265       ,\
    help='Natal velocity parameter of black holes (BH) in km/s')

parser.add_argument('-chi' ,'--NatalSpinParameter'   ,type=float,metavar=' ',default=0         ,\
    help='Natal spin parameter of first genetation BHs')

parser.add_argument('-SM'  ,'--NatalSpinDistribution',type=int  ,metavar=' ',default=0         ,\
    help='Natal spin distribution (=1 for monochromatic, =0 for uniform)')

parser.add_argument('-tMax','--SimulationTime'       ,type=float,metavar=' ',default=13800     ,\
    help='Maximum simulation time in mega years (Myr)')

parser.add_argument('-dtMin','--MinimumTimeStep'     ,type=float,metavar=' ',default=.1        ,\
    help='Minimum simulation timestep in Myr')

parser.add_argument('-dtMax','--MaximumTimeStep'     ,type=float,metavar=' ',default=50        ,\
    help='Maximum simulation timestep in Myr')

parser.add_argument('-z'    ,'--FormationRedshift'   ,type=float,metavar=' ',default=3         ,\
    help='Redshfit of cluster formation')

parser.add_argument('-aIMF' ,'--HighMassIMFslope'    ,type=float,metavar=' ',default=-2.3      ,\
    help='High mass initial star mass function slope')

parser.add_argument('-ZAMSmax'   ,'--MaximumZAMSmass',type=float,metavar=' ',default=150       ,\
    help='Maximum ZAMS star mass in Msun (less than 340 Msun)')

parser.add_argument('-c'    ,'--ConcentrationNFW'    ,type=float,metavar=' ',default=10        ,\
    help='Concentration parameter of the Navarro-Frenk-White profile')

parser.add_argument('-Rs'   ,'--ScaleRadiusNFW'      ,type=float,metavar=' ',default=50        ,\
    help='Scale radius of the Navarro-Frenk-White profile in kpc')

parser.add_argument('-Mh'   ,'--DarkMatterHaloMass'  ,type=float,metavar=' ',default=1e12      ,\
    help='Total dark matter halo mass of the host galaxy in Msun')

parser.add_argument('-s'    ,'--Seed'                ,type=int  ,metavar=' ',default=123456789 ,\
    help='Random number generator seed')

parser.add_argument('-MF'   ,'--MergersFile'         ,type=str  ,metavar=' ',default='mergers'    ,\
    help='Name of output .txt file containing binary BH merger source parameters')

parser.add_argument('-EF'   ,'--EvolutionFile'       ,type=str  ,metavar=' ',default='evolution'    ,\
    help='Name of output .txt file containing time-dependent quantities of interest')

parser.add_argument('-BF'   ,'--BlackHoleFile'      ,type=str  ,metavar=' ',default='blackholes'    ,\
    help='Name of output .npz file containing all the masses of black holes at the initial and final state')

parser.add_argument('-CF'   ,'--CollisionsFile'     ,type=str  ,metavar=' ',default='collisions' ,\
    help='Name of output .npz file containing information regarding initial phase of stellar collisions giving rise to runaway')

#parser.add_argument('-rhog', '--gasDensity', type=float, metavar=' ', default=0, help='Ambient gas density of hydrogen in cm^-3')
#parser.add_argument('-cs', '--soundSpeed', type=float, metavar=' ', default=100, help='Sound speed of ISM in km/s')
parser.add_argument('-ep', '--radiativeEfficiency', type=float, metavar=' ', default=0.1, help='Radiative efficiency in black hole accretion')
parser.add_argument('-eSF', '--starFormationEfficiency', type=float, metavar=' ', default=0.3, help='Star formation efficiency')
parser.add_argument('-tM', '--expulsionTimescale', type=float, metavar=' ', default=-1, help='Expulsion type in Myr (<0 for default residual gas removal)')

args = parser.parse_args()

Mcl           = args.ClusterMass * Msun
rh            = args.HalfMassRadius * pc
rhoC          = args.CentralDensity * Msun/pc**3
Rgal          = args.GalactocentricRadius * kpc
Z             = args.Metallicity * Zsun
fB            = args.BinaryFraction
wNatalKick    = args.NatalKickParameter * 1e3
chiNatal      = args.NatalSpinParameter
spinModel     = args.NatalSpinDistribution
tMax          = args.SimulationTime * Myr
dtMin         = args.MinimumTimeStep * Myr
dtMax         = args.MaximumTimeStep * Myr
zClForm       = args.FormationRedshift
alphaIMF      = args.HighMassIMFslope
Mstar_max     = args.MaximumZAMSmass
conc          = args.ConcentrationNFW
Rscale        = args.ScaleRadiusNFW * kpc
Mhalo         = args.DarkMatterHaloMass * Msun
SEED          = args.Seed
mergersFile   = args.MergersFile
evolutionFile = args.EvolutionFile
blackholeFile = args.BlackHoleFile

collisionFile = args.CollisionsFile

#rho_gas     = args.gasDensity * mH / cm**3
#c_sound     = args.soundSpeed * 1e3
epsilon_acc = args.radiativeEfficiency
epsilon_SF  = args.starFormationEfficiency
tM          = args.expulsionTimescale * Myr

Mcl0  = Mcl  # initial cluster mass (stars+gas)
rh0   = rh   # initial half-mass radius
rhoC0 = rhoC # initial central density

# initial gas half-mass radius:
rh_gas0 = rh0

# initial residual gas not converted into stars:
M_gas0 = (1-epsilon_SF) * Mcl0

# mass converted in stars:
Mcl_stars = epsilon_SF * Mcl0

# initial central gas density:
rho_gas0 = 3 * M_gas0 / (4 * np.pi * (rh0 / 1.3)**3)

# current gas mass:
M_gas = M_gas0

# current gas radius:
rh_gas = rh_gas0

# Initialize pseudo-number generator
# ----------------------------------------------------------------------------------------------------------------------------

np.random.seed(SEED)

if __name__=="__main__":

    # Cluster set-up
    # ---------------------------------------------------------------------------------------------------------------------------

    # lightest star mass (solar masses):
    Mstar_min = 0.08

    # average stellar mass in solar masses:
    mAvg = integrate.quad(lambda x: x*IMF_kroupa(np.array([x])),Mstar_min,Mstar_max)[0]\
          /integrate.quad(lambda x:   IMF_kroupa(np.array([x])),Mstar_min,Mstar_max)[0] * Msun

    # total initial number of stars in cluster:
    Nstar = Mcl_stars/mAvg

    # central number density of stars:
    nStar = rhoC/mAvg

    # total number of original binary stars:
    N_StarStar = int(fB*Nstar/2)

    # central number density of binary stars:
    n_StarStar = fB*nStar/2

    # velocity dispersion of stars:
    vStar = veloDisp(mAvg,1,mAvg,Mcl,rh)
    
    # initial half-mass relaxation timescale:
    tRel0 = tRelax(Mcl0, epsilon_SF * Mcl0 / mAvg,rh0,mAvg)

    # lightest massive star (solar masses) DO NOT CHANGE THIS VALUE should be set at `20`:
    MstarMassive_min = 20
    
    # Core collapse:
    # ---------------------------------------------------------------------------------------------------------------------------
    
    # core-collapse time:
    t_cc = 0.20 * tRel0
    
    # initialize simulation time:
    t = t_cc
    
    # initialize redshift:
    z = redd(t_lbb(zClForm)-t)
    
    # early phase of stellar mergers:
    # ---------------------------------------------------------------------------------------------------------------------------
    
    # runaway star mass:
    m_r = 0.0
    
    if t_cc < 3*Myr: # stellar mergers dominate the initial evolution of the cluster
     
        mZAMS = np.linspace(Mstar_min, Mstar_max, 10**6)
        pdf_mZAMS = IMF_kroupa(mZAMS)
        cdf_mZAMS = integrate.cumulative_trapezoid(pdf_mZAMS, mZAMS, initial=0)
        inv_cdf_mZAMS = interpolate.interp1d(cdf_mZAMS/cdf_mZAMS[-1], mZAMS)

        # sample all cluster stars:
        mZAMS = inv_cdf_mZAMS(np.random.rand(Nstar)) * Msun

        # Compute stellar ages:
        tLives = np.zeros(mZAMS.size)
        for k in range(tLives.size):
            tLives[k] = tMS(mZAMS[k], Z)
     
        # seed runaway star:
        m_seed = mZAMS.max()
       
        # index position of seed runaway:
        jHeaviest = int(np.where(mZAMS==m_seed)[0])
        
        # remove mass of seed from star list:
        mZAMS = np.delete(mZAMS, jHeaviest)
        
        # remove age of seed from star list:
        tLives = np.delete(tLives, jHeaviest)
        
        # runaway mass:
        m_r = m_seed
        
        # rejuvenation efficiency:
        frejuv = 1.0
        
        # collision timescale:
        t_coll = 1/(2.2e-4 * Nstar / tRel0)
        
        # age of runaway star when it first collides (at moment of core collapse):
        t_r = t_cc
        
        i = 0
        
        redshift__ = []
        time__ = []
        m_runaway__ = []
        m_increase__ = []
        t_collision__ = []
        t_age__ = []
        
        mZAMS_collapsed = []
        
        while tMS(m_r, Z)-t_r > t_coll and M_cl > 0: # evolve the runaway star
        
            # Remove stars that evolve beyond MS to their death (leaving behind a remnant NS or BH)
            jEvolved = np.where(tLives < t)
            mZAMS_evolved = mZAMS[jEvolved]
            mZAMS_collapsed = mZAMS_collapsed + list(mZAMS_evolved)
            mZAMS = np.delete(mZAMS, jEvolved)
            tLives = np.delete(tLives, jEvolved)
            
            # collision timescale:
            t_coll = 1 / (2.2e-4 * Nstar / tRel0)
            
            # local simulation timestep:
            dt = t_coll
            
            # smallest star segregated by time t:
            m_f = 1.9*Msun * Myr/t * (rh0/pc)**(3/2) * (Mcl/Msun)**(1/2) / np.log(0.1*Nstar)

            # segregated stars by time t:
            mZAMS_candidates = mZAMS[mZAMS > m_f]

            # sampling probability:
            power = 1.5
            p = mZAMS_candidates**power / np.sum(mZAMS_candidates**power)
            
            # mass increase in the collision:
            dm = np.random.choice(mZAMS_candidates, size=1, p=p)[0]

            jDelete = int(np.where(mZAMS==dm)[0][0])
            mZAMS = np.delete(mZAMS, jDelete)
            tLives = np.delete(tLives, jDelete)

            # mass ratio:
            q = np.min([m_r, dm])/np.max([m_r, dm])

            redshift__.append(z)
            time__.append(t/Myr)
            m_runaway__.append(m_r/Msun)
            m_increase__.append(dm/Msun)
            t_collision__.append(t_coll/Myr)
            t_age__.append(t_r/Myr)

            # update number of stars:
            Nstar = Nstar - 1

            # collision remnant age:
            t_r = frejuv * tMS(m_r+dm, Z)/(m_r+dm)*(m_r*t_r/tMS(m_r, Z) + dm*t/tMS(dm, Z))

            # update runaway mass (some fraction is lost):
            m_r = (m_r + dm)*(1 - 0.3*q/(1+q)**2)

            # time update:
            t = t + dt

            # redshift update:
            z = redd(t_lbb(zClForm)-t)
            
            i = i + 1
            
        redshift__ = np.array(redshift__)
        time__ = np.array(time__)
        m_runaway__ = np.array(m_runaway__)
        m_increase__ = np.array(m_increase__)
        t_collision__ = np.array(t_collision__)
        t_age__ = np.array(t_age__)
        
        np.savez(collisionFile, z=redshift__, t=time__, mr=m_runaway__, dm=m_increase__, tc=t_collision__, ta=t_age__)
    
        # remaining massive stars that collapse to compact objects:
        starMasses = np.array(mZAMS_collapsed + list(mZAMS[mZAMS/Msun>MstarMassive_min]))
    
    else: # no runaway stellar collisions immediately after core collapse
     
        # initial number of stars with mass above MstarMassive_min = 20 solar masses:
        Nstar_massive = int(Mcl_stars/Msun*integrate.quad(lambda x:   IMF_kroupa(np.array([x])), MstarMassive_min, Mstar_max)[0]\
                                          /integrate.quad(lambda x: x*IMF_kroupa(np.array([x])), Mstar_min       , Mstar_max)[0])

        # star masses drawn from Kroupa (2002) IMF with inverse sampling (we allow for possibly different high-mass slope):
        u_star = np.random.rand(Nstar_massive)
        starMasses = (u_star*(Mstar_max**(alphaIMF+1)-MstarMassive_min**(alphaIMF+1))\
                      +MstarMassive_min**(alphaIMF+1))**(1/(alphaIMF+1))
    
    # Black holes in cluster
    # ---------------------------------------------------------------------------------------------------------------------------
    
    # remnant masses from SEVN interpolant:
    remnantMasses = Mrem(starMasses,Z) * Msun

    # minimum natal black hole mass:
    Mbh_min = 3*Msun

    # discard zero masses or smaller than some minimum value (black holes masses in kg):
    bhMasses = remnantMasses[remnantMasses>Mbh_min]

    # ZAMS star masses that give rise to compact remnants >3 solar masses (BHs):
    star_masses = starMasses[remnantMasses>Mbh_min] * Msun

    # total number of black holes initially to form:
    NBH = bhMasses.size

    # mean (initial) BH mass Produced:
    mBH_mean = np.mean(bhMasses)

    aBin_min = 3*Rsun                      # minimum original binary separation
    aBin_max = (3/np.pi/4/nStar)**(1/3)/10 # maximum original binary separation

    # probability that a star that gives birth to a BH (with mass above 3 solar masses):
    proba_BH = NBH/Nstar

    # As a fraction of the total number of original binary stars:
    ###########################################################################################################################
    # In this version of the code, the neglect original BBHs and BH-star pairs.
    ###########################################################################################################################
    N_BHstar_original = 0 #N_StarStar * proba_BH * (1-proba_BH)       # number of original binaries of the type BH-star
    N_BBH_original    = 0 #N_StarStar * proba_BH**2                   # number of original binaries of the type BH-BH
    N_BH_iso          = NBH - N_BHstar_original - 2*N_BBH_original # number of isolated BHs

    aBHstar_har    = G_Newt * mAvg    /4/vStar**2 # BHstar    semimajor axis at the hard-soft boundary
    aBBH_har       = G_Newt * mBH_mean/4/vStar**2 # BBH       semimajor axis at the hard-soft boundary
    aStarStar_har  = G_Newt * mAvg    /4/vStar**2 # star-star semimajor axis at the hard-soft boundary

    # Number of original star-star binaries that do not evolve into BHs (re-defined):
    N_StarStar_original = N_StarStar - N_BHstar_original - N_BBH_original
    N_StarStar = int(N_StarStar_original) # re-name variable

    # star-star semimajor axes (log-flat distributed):
    a_StarStar_original = 10**(np.random.uniform(np.log10(aBin_min),np.log10(aBin_max),int(N_StarStar)))
    aStarStar_original = a_StarStar_original[a_StarStar_original < aStarStar_har] # store only the hard binary semimajor axes
    aStarStar = aStarStar_original # re-name array

    # BH-star semimajor axes (log-flat distributed):
    a_BHstar_original = 10**(np.random.uniform(np.log10(aBin_min),np.log10(aBin_max),int(N_BHstar_original)))
    aBHstar_original = a_BHstar_original[a_BHstar_original < aBHstar_har] # store only the hard binary semimajor axes

    # BBH semimajor axes (log-flat distributed):
    a_BBH_original = 10**(np.random.uniform(np.log10(aBin_min),np.log10(aBin_max),int(N_BBH_original)))
    aBBH_original = a_BBH_original[a_BBH_original < aBBH_har] # store only the hard binary semimajor axes

    N_BHstar_original_hard = aBHstar_original.size                 # number of original hard binaries of the type BH-star
    N_BBH_original_hard    = aBBH_original.size                    # number of original hard binaries of the type BH-BH
    N_BH_iso               = NBH - N_BHstar_original_hard - 2*N_BBH_original_hard # number of isolated BHs (re-calculated)

    # Redefine variable names to unclutter notation (original binary here means hard binary):
    N_BHstar_original = N_BHstar_original_hard
    N_BBH_original    = N_BBH_original_hard

    # indices of ZAMS masses giving birth to BHs:
    indices_array = np.arange(NBH)

    # chose indices of stars to participate in binaries: (draw randomly to imitate a uniform mass ratio):
    indices_BHstar_original = np.random.choice(indices_array,size=N_BHstar_original,replace=False)
    indices_BBH_original    = np.random.choice(indices_array,size=2*N_BBH_original ,replace=False)

    mask_BHstar_original = np.zeros(indices_array.shape,dtype=bool)
    mask_BHstar_original[indices_BHstar_original] = True

    mask_BBH_original = np.zeros(indices_array.shape,dtype=bool)
    mask_BBH_original[indices_BBH_original] = True

    # Array with BH masses born in original BHstar binaries:
    bhMasses_BHstar_original = bhMasses[mask_BHstar_original]
    
    # Array with BH progenitor ZAMS masses in original BHstar binaries:
    star_masses_BHstar_original = star_masses[mask_BHstar_original]

    # Array with BH masses born in original BBH binaries:
    bhMasses_BBH_original = bhMasses[mask_BBH_original]
    
    # Array with BH progenitor ZAMS masses in original BBH binaries:
    star_masses_BBH_original = bhMasses[mask_BBH_original]

    mask_iso = np.zeros(indices_array.shape,dtype=bool)
    mask_iso[np.append(indices_BHstar_original,indices_BBH_original)] = True

    # Array with BH masses born from stars in isolation (within the cluster):
    bhMasses_iso = bhMasses[~mask_iso]
    
    # Array with BH progenitor ZAMS masses of stars in isolation (within the cluster):
    star_masses_iso = star_masses[~mask_iso]

    # Natal kicks of BHs:
    # ---------------------------------------------------------------------------------------------------------------------------

    # Final natal kick imprinted onto isolated 1g BHs after the initial fallback:
    vNatalKickFinal_iso = np.zeros(star_masses_iso.size)
    for i in range(0,star_masses_iso.size):
        # sample BH natal kick in m/s:
        vNatalKick = maxwell.rvs(loc=0,scale=wNatalKick/np.sqrt(3),size=1)
        vNatalKickFinal_iso[i] = (1-f_fb(star_masses_iso[i]/Msun))*vNatalKick

    # Final natal kick imprinted onto 1g BHs in original BHstar binaries after the initial fallback:
    vNatalKickFinal_BHstar_original = np.zeros(star_masses_BHstar_original.size)
    
    for i in range(0,star_masses_BHstar_original.size):
        
        # sample BH natal kick in m/s:
        vNatalKick = maxwell.rvs(loc=0,scale=wNatalKick/np.sqrt(3),size=1)
        vNatalKickFinal_BHstar_original[i] = (1-f_fb(star_masses_BHstar_original[i]/Msun))*vNatalKick

    # Final natal kick imprinted onto 1g BHs in original BBH binaries after the initial fallback:
    vNatalKickFinal_BBH_original = np.zeros(star_masses_BBH_original.size)
    
    for i in range(0,star_masses_BBH_original.size):
        
        # sample BH natal kick in m/s:
        vNatalKick = maxwell.rvs(loc=0,scale=wNatalKick/np.sqrt(3),size=1)
        vNatalKickFinal_BBH_original[i] = (1-f_fb(star_masses_BBH_original[i]/Msun))*vNatalKick

    # BH retention:
    # ---------------------------------------------------------------------------------------------------------------------------

    # Isolated BHs initially retained in the cluster:
    bhMasses_iso_ret = bhMasses_iso[vNatalKickFinal_iso < vEscape(Mcl,rh)] # check if black hole is retained
    N_BH_iso_ret = bhMasses_iso_ret.size # number of isolated BHs retained

    # Original BBHs initially retained in the cluster:
    bhMasses_BBH_original_ret = bhMasses_BBH_original
    aBBH_original_ret = aBBH_original
    
    # number of original retained BBHs:
    N_BBH_original_ret = int(bhMasses_BBH_original_ret.size/2)

    # Original BHstar initially retained in the cluster:
    bhMasses_BHstar_original_ret = bhMasses_BHstar_original[vNatalKickFinal_BHstar_original < vEscape(Mcl,rh)]
    aBHstar_original_ret = aBHstar_original[vNatalKickFinal_BHstar_original < vEscape(Mcl,rh)]
    
    # number of original retained BHstar binaries:
    N_BHstar_original_ret = bhMasses_BHstar_original_ret.size

    # Total number of retained BHs:
    N_BH_ret = N_BH_iso_ret + 2*N_BBH_original_ret + N_BHstar_original_ret
    N_BH = N_BH_ret # change of variable name

    # Set initial conditions
    # ---------------------------------------------------------------------------------------------------------------------------

    mBH = bhMasses_iso_ret # initialize dynamical mass array of isolated retained BHs

    if spinModel==1: # natal spin distribution chosen monochromatic:
        
        # dynamical array of spins of initially isolated BHs:
        sBH = chiNatal * np.ones(N_BH_iso_ret)
    
    else: # natal spin distribution chosen uniform:
        
        sBH = np.random.uniform(0,chiNatal,size=N_BH_iso_ret)

    # dynamical array of generations (initially all black holes are first generation by definition):
    gBH = np.ones(N_BH_iso_ret)

    # Collapse of the runaway star into BH:
    if m_r > 0.0:
        
        # mass of runaway remnant:
        m_r_BH = Mrem_Fryer2012(m_r, Z) * Msun
        
        # check if remnant is a BH:
        if m_r_BH > Mbh_min:
        
            # append runaway remnant into BH masses array:
            mBH = np.append(mBH, m_r_BH)

            if spinModel==1: # natal spin distribution chosen monochromatic:
                s_r_BH = chiNatal
            else: # natal spin distribution chosen uniform:
                s_r_BH = np.random.uniform(0,chiNatal,size=1)

            # append runaway spin into BH spins array:
            sBH = np.append(sBH, s_r_BH)

            # append generation array:
            gBH = np.append(gBH, 1)

            # increase number of BHs in cluster by 1, corresponding to the runaway star collapsing into a BH:
            N_BH += 1
    
    # define `empty` allocatable arrays of bound systems in cluster:

    # binaries = [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,zForm,Nhar,Nsub]
    binaries = np.zeros(shape=(1,13))

    # pairs = [a,m,chi,g]
    pairs = np.zeros(shape=(1,4))

    # triples = [aInner,aOuter,eInner,eInner,m0,m1,m2,s0,s1,s2,g0,g1,g2,inclination1,inclination2,tForm,zForm]
    triples = np.zeros(shape=(1,17))

    # mergers = [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,tMerge,zForm,zMerge,Nhar,Nsub,q,chiEff,theta1,theta2,dPhi,
    #            mRem,sRem,gRem,vGW,j1,j2]
    mergers = np.zeros(shape=(1,26))

    # evolut = [t,z,Mcl,rh,Rgal,N_BH,N_BH_sin,N_BHstar,N_BBH,N_Triples,N_me,N_meRe,N_meEj,
    #           N_meOut,N_ZLK,N_cap,N_ej,vEsc,nStarStar,meanVseg,xi,vStar,vBH,tCap,tEx1,tEx2,t3bb,tBB,tPP, N_ex1, N_ex2, N_3bb, N_bb, N_pp, mBHmean]
    evolut = np.zeros(shape=(1,35))

    # BBH assembly channel: (- sign means BBH was ejected from the cluster)
    #      0: original binary
    #   (-)1: exchange processes
    #   (-)2: two-body capture
    #   (-)3: three-BH binary induced
    #   (-)4: von Zeipel-Lidov-Kozai merger

    # Append original BBHs retained in the cluster:
    for i in range(0,N_BBH_original_ret):

        if spinModel==1: # natal spin distribution chosen monochromatic:
            
            sBH_1 = chiNatal
            sBH_2 = chiNatal
        
        else: # natal spin distribution chosen uniform:
            
            sBH_1 = np.random.uniform(0,chiNatal,size=1)[0]
            sBH_2 = np.random.uniform(0,chiNatal,size=1)[0]

        mBH_1 = bhMasses_BBH_original_ret[i]                    # first  BH member mass
        mBH_2 = bhMasses_BBH_original_ret[i+N_BBH_original_ret] # second BH member mass

        # append binary [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,zForm,Nhar,Nsub]:
        binaries = np.append(binaries,[[0,aBBH_original_ret[i],0,mBH_1,mBH_2,sBH_1,sBH_2,1,1,0,zClForm,0,0]],axis=0)

    # Append original BHstar pairs retained in the cluster:
    for i in range(0,N_BHstar_original_ret):

        if spinModel==1: # natal spin distribution chosen monochromatic:
            
            sBH_ = chiNatal
        
        else: # natal spin distribution chosen uniform:
            
            sBH_ = np.random.uniform(0,chiNatal,size=1)[0]

        # append pair [a,m,chi,g]:
        pairs = np.append(pairs,[[aBHstar_original_ret[i],bhMasses_BHstar_original_ret[i],sBH_,1]],axis=0)

    N_BBH     = N_BBH_original_ret    # current number of BBHs     in cluster
    N_BHstar  = N_BHstar_original_ret # current number of BH-stars in cluster
    N_Triples = 0                     # current number of triples  in cluster

    N_me       = 0 # number of binary black hole mergers
    N_cap      = 0 # number of BHBH captures
    N_meRe     = 0 # number of merger remnants retained
    N_meEj     = 0 # number of merger remnants ejected
    N_ZLK      = 0 # number of von Zeipel-Lidov-Kozai mergers
    N_dis      = 0 # number of binary black hole disruptions
    N_ej       = 0 # number of binary black holes ejected
    N_meOut    = 0 # number of ejected binaries that merge in the field
    N_exch     = 0 # number of BBH-BH dynamical exchanges
    N_sinBHeje = 0 # number of single BHs ejected from binary-single interactions

    N_bb = 0 # cumulative number of BBH - BBH interactions
    N_pp = 0 # cumulative number of BHstar - BHstar interactions
    N_3bb = 0 # cumulative number of 3bb interactions
    N_ex1 = 0 # cumulative number of star-star -> BH-star interactions
    N_ex2 = 0 # cumulative number of BH-star -> BBH interactions
    
    N_me_1_1 = 0 # number of 1g+1g mergers
    N_me_1_2 = 0 # number of 1g+2g mergers
    N_me_2_2 = 0 # number of 2g+2g mergers
    N_me_1_3 = 0 # number of 1g+3g mergers
    N_me_2_3 = 0 # number of 2g+3g mergers
    N_me_3_3 = 0 # number of 3g+3g mergers
    
    # initialize number of iterations:
    Niter = 0
    
    # Cumulative number of associations formed:
    N_bhstars_cumul = N_BHstar_original_ret # cumul. num. of BH-stars
    N_bbhs_cumul    = N_BBH_original_ret    # cumul. num. of BBHs
    N_triples_cumul = 0                     # cumul. num. of triples

    # half-mass relaxation timescale, updated:
    tRel = tRelax(Mcl, Mcl_stars/mAvg,rh,mAvg)

    # velocity dispersion of stars (updated):
    vStar = veloDisp(mAvg,1,mAvg,Mcl,rh)

    # update central density of stars (self-similarly):
    rhoC = (Mcl_stars / epsilon_SF / Mcl0) * (rh0/rh)**3 * rhoC0

    # central number density of stars (updated):
    nStar = rhoC/mAvg

    # current number of single black holes:
    N_BH_sin = mBH.size

    # mean BH mass:
    meanBHmass = (np.sum(mBH) + np.sum(np.transpose(binaries)[:][3]+np.transpose(binaries)[:][4]) + np.sum(np.transpose(pairs)[:][1])) / N_BH
    
    # individual mass ratio:
    q__ = meanBHmass / mAvg
    
    # total mass ratio:
    Q__ = meanBHmass*N_BH / Mcl
    
    # temperature ratio:
    xi  = q__*Q__**(2/3)*(27/4)**(1/3)*(1+3/2*Q__/q__)**(1/3)*(1+5/2*Q__)**(-1)

    if xi<1:
        
        # equipartition can be established:
        xi = 1

    # mean segregation volume:
    meanVseg = Vseg(meanBHmass,meanBHmass*N_BH,xi,mAvg,Mcl,rh)

    # number density of single BHs:
    nBHsin = N_BH_sin / meanVseg
    
    # binary star semimajor axis at the hard-soft boundary, updated:
    aStarStar_har = G_Newt*mAvg/4/vStar**2

    # (as cluster evolves and velocity decreases some previously hard binary stars become soft)
    # regard only hard binary binaries:
    aStarStar = aStarStar[aStarStar < aStarStar_har]

    # central number density of hard binary stars (within BH segregation volume):
    nStarStar = ((3*meanVseg/4/np.pi)**(1/3)/rh)**3*aStarStar.size / meanVseg

    # Calculation of timescales
    # -----------------------------------------------------------------------------------------------------------------------

    # Single-single capture:
    if N_BH_sin < 2: # no captures possible if there are less than two single black holes available in the system
    # (ignore captures from multiple body interactions)
            
        tauCap = 1e100*yr
        
    else:
            
        tauCap = 1/RateCapture(meanBHmass,meanBHmass,N_BH_sin,nBHsin,xi,mAvg,Mcl,rh)

    # BHBHBH induced BBH:
    if N_BH_sin < 3: # need at least three single BHs for this channel
            
        tau3bBH = 1e100*yr
        
    else:
            
        # minimum binary hardness ratio:
        etaMin = 5
            
        # (division by the mean segregation volume allows to convert rate density into total rate):
        tau3bBH = 1/Rate3body(meanBHmass,nBHsin,veloDisp(meanBHmass,xi,mAvg,Mcl,rh),etaMin)/meanVseg

    # Binary star -> BH-star exchange:
    if N_BH_sin < 1 or len(aStarStar)==0\
        or nStarStar*meanVseg<1: # need at least two single BHs to perform the exchange and at least one hard star binary
            
        tauPair = 1e100*yr
        
    else:
            
        # (factor of `2` included to account for the possibility of either star being exchanged):
        tauPair = 1/RateExchange(mAvg,mAvg,meanBHmass,nStarStar,N_BH_sin,np.mean(aStarStar),xi,mAvg,Mcl,rh)/2

    # BH-star -> BBH exchange:
    if N_BHstar < 1 or N_BH_sin<1: # need at least one BH-star pair and a single BH for this channel
            
        tauEx = 1e100*yr
        
    else:
            
        # mean BH mass which participates in BH-star pairs:
        meanBHpairMass = np.mean(np.transpose(pairs)[:][1])
            
        # BH-star -> BBH exchange timescale:
        tauEx = 1/RateExchange(mAvg,meanBHpairMass,meanBHmass,N_BHstar/meanVseg,\
            N_BH_sin,np.mean(np.transpose(pairs)[:][0]),xi,mAvg,Mcl,rh)

    # BBH-BBH strong encounter:
    if N_BBH < 2: # need at least two BBHs to have a binary-binary encounter
            
        tauBB = 1e100*yr
        
    else:
            
        # mean BBH mass:
        meanBBHmass = np.mean(np.transpose(binaries)[:][3]+np.transpose(binaries)[:][4])

        # mean BBH sma:
        meanBBHsma = np.mean(np.transpose(binaries)[:][1])

        # number density of BBHs:
        nBBH = (N_BBH-1)/meanVseg

        # BBH-BBH strong encounter timescale:
        tauBB = 1/RateInter(meanBBHmass,meanBBHmass,N_BBH,nBBH,7*meanBBHsma,xi,mAvg,Mcl,rh)

    # BHstar - BHstar strong encounter:
    if N_BHstar < 2: # need at least two BH-star pairs to have a pair-pair encounter
            
        tauPP = 1e100*yr
        
    else:
            
        # mean BH-star mass:
        meanBHstarMass = np.mean(np.transpose(pairs)[:][1])+mAvg

        # mean BH-star sma:
        meanBHstarSma = np.mean(np.transpose(pairs)[:][0])

        # number density of BH-stars:
        nBHstar = (N_BHstar-1)/meanVseg

        # BHstar - BHstar mean interaction timescale:
        tauPP = 1/RateInter(meanBHstarMass,meanBHstarMass,N_BHstar,nBHstar,7*meanBHstarSma,xi,mAvg,Mcl,rh)
    
    # append evolution params [t,z,Mcl,rh,Rgal,N_BH,N_BH_sin,N_BHstar,N_BBH,N_Triples,N_me,
    #                          N_meRe,N_meEj,N_meOut,N_ZLK,N_cap,N_ej,vEsc,nStarStar,meanVseg,xi,vStar,vBH,tCap,tEx1,tEx2,t3bb,tBB,tPP,
    #                          N_ex1, N_ex2, N_3bb, N_bb, N_pp, mBHmean]:
    evolut = np.append(evolut,[[t/Myr,z,Mcl/Msun,rh/pc,Rgal/kpc,\
        N_BH,N_BH_sin,N_BHstar,N_BBH,N_Triples,N_me,N_meRe,N_meEj,N_meOut,N_ZLK,N_cap,N_ej,vEscape(Mcl,rh)/1e3,\
            nStarStar*pc**3,meanVseg/pc**3,xi,np.sqrt(0.4*G_Newt*Mcl/rh)/1e3,veloDisp(meanBHmass,xi,mAvg,Mcl,rh)/1e3,\
                               tauCap/Myr,tauPair/Myr,tauEx/Myr,tau3bBH/Myr,tauBB/Myr,tauPP/Myr, N_ex1, N_ex2, N_3bb, N_bb, N_pp, \
                                meanBHmass/Msun]],axis=0)

    # save masses of 1g black holes:
    mBH_1g = mBH
    
    # Simulation block
    # ---------------------------------------------------------------------------------------------------------------------------

    # start global clock:
    startTimeGlobal = time.time()

    # simulate while termination criteria is not met:
    while t < tMax and N_BH > 1 and z > 0 and Rgal > 0 and Mcl>0:

        # start local clock:
        startTimeLocal = time.time()

        # half-mass relaxation timescale, updated:
        tRel = tRelax(Mcl, Mcl_stars/mAvg,rh,mAvg)

        # velocity dispersion of stars (updated):
        vStar = veloDisp(mAvg,1,mAvg,Mcl,rh)

        # update central density of stars (self-similarly):
        rhoC = (Mcl_stars / epsilon_SF / Mcl0) * (rh0/rh)**3 * rhoC0

        # central number density of stars (updated):
        nStar = rhoC/mAvg
        
        # current number of stars:
        Nstar = Mcl_stars / mAvg

        # current number of single black holes:
        N_BH_sin = mBH.size
        
        # mean BH mass:
        meanBHmass = (np.sum(mBH) + np.sum(np.transpose(binaries)[:][3]+np.transpose(binaries)[:][4]) + np.sum(np.transpose(pairs)[:][1])) / N_BH

        # individual mass ratio:
        q__ = meanBHmass / mAvg
        
        # total mass ratio:
        Q__ = meanBHmass*N_BH / Mcl
        
        # temperature ratio:
        xi  = q__*Q__**(2/3)*(27/4)**(1/3)*(1+3/2*Q__/q__)**(1/3)*(1+5/2*Q__)**(-1)

        if xi<1:
            
            # equipartition can be established:
            xi = 1
        
        # mean segregation volume:
        meanVseg = Vseg(meanBHmass,meanBHmass*N_BH,xi,mAvg,Mcl,rh)

        # number density of single BHs:
        nBHsin = N_BH_sin / meanVseg

        # binary star semimajor axis at the hard-soft boundary, updated:
        aStarStar_har = G_Newt*mAvg/4/vStar**2

        # (as cluster evolves and velocity decreases some previously hard binary stars become soft)
        # regard only hard binary binaries:
        aStarStar = aStarStar[aStarStar < aStarStar_har]

        # central number density of hard binary stars (within BH segregation volume):
        nStarStar = ((3*meanVseg/4/np.pi)**(1/3)/rh)**3*aStarStar.size / meanVseg

        # Calculation of timescales
        # -----------------------------------------------------------------------------------------------------------------------

        # Single-single capture:
        if N_BH_sin < 2: # no captures possible if there are less than two single black holes available in the system
        # (ignore captures from multiple body interactions)
            
            tauCap = 1e100*yr
        
        else:
            
            tauCap = 1/RateCapture(meanBHmass,meanBHmass,N_BH_sin,nBHsin,xi,mAvg,Mcl,rh)

        # BHBHBH induced BBH:
        if N_BH_sin < 3: # need at least three single BHs for this channel
            
            tau3bBH = 1e100*yr
        
        else:
            
            # minimum binary hardness ratio:
            etaMin = 5
            
            # (division by the mean segregation volume allows to convert rate density into total rate):
            tau3bBH = 1/Rate3body(meanBHmass,nBHsin,veloDisp(meanBHmass,xi,mAvg,Mcl,rh),etaMin)/meanVseg

        # Binary star -> BH-star exchange:
        if N_BH_sin < 1 or len(aStarStar)==0\
            or nStarStar*meanVseg<1: # need at least two single BHs to perform the exchange and at least one hard star binary
            
            tauPair = 1e100*yr
        
        else:
            
            # (factor of `2` included to account for the possibility of either star being exchanged):
            tauPair = 1/RateExchange(mAvg,mAvg,meanBHmass,nStarStar,N_BH_sin,np.mean(aStarStar),xi,mAvg,Mcl,rh)/2

        # BH-star -> BBH exchange:
        if N_BHstar < 1 or N_BH_sin<1: # need at least one BH-star pair and a single BH for this channel
            
            tauEx = 1e100*yr
        
        else:
            
            # mean BH mass which participates in BH-star pairs:
            meanBHpairMass = np.mean(np.transpose(pairs)[:][1])
            
            # BH-star -> BBH exchange timescale:
            tauEx = 1/RateExchange(mAvg,meanBHpairMass,meanBHmass,N_BHstar/meanVseg,\
                N_BH_sin,np.mean(np.transpose(pairs)[:][0]),xi,mAvg,Mcl,rh)

        # BBH-BBH strong encounter:
        if N_BBH < 2: # need at least two BBHs to have a binary-binary encounter
            
            tauBB = 1e100*yr
        
        else:
            
            # mean BBH mass:
            meanBBHmass = np.mean(np.transpose(binaries)[:][3]+np.transpose(binaries)[:][4])

            # mean BBH sma:
            meanBBHsma = np.mean(np.transpose(binaries)[:][1])

            # number density of BBHs:
            nBBH = (N_BBH-1)/meanVseg

            # BBH-BBH strong encounter timescale:
            tauBB = 1/RateInter(meanBBHmass,meanBBHmass,N_BBH,nBBH,7*meanBBHsma,xi,mAvg,Mcl,rh)

        # BHstar - BHstar strong encounter:
        if N_BHstar < 2: # need at least two BH-star pairs to have a pair-pair encounter
            
            tauPP = 1e100*yr
        
        else:
            
            # mean BH-star mass:
            meanBHstarMass = np.mean(np.transpose(pairs)[:][1])+mAvg

            # mean BH-star sma:
            meanBHstarSma = np.mean(np.transpose(pairs)[:][0])

            # number density of BH-stars:
            nBHstar = (N_BHstar-1)/meanVseg

            # BHstar - BHstar mean interaction timescale:
            tauPP = 1/RateInter(meanBHstarMass,meanBHstarMass,N_BHstar,nBHstar,7*meanBHstarSma,xi,mAvg,Mcl,rh)

        # Calculate time step (should not be below minimum or above some maximum threshold)
        # -----------------------------------------------------------------------------------------------------------------------

        dt1 = np.min([dtMax , np.max([ dtMin , np.min([tauCap,tauPair,tauEx,tau3bBH,tauBB,tauPP]) ]) ])

        # Draw number of incidents for each channel in the given time step
        # -----------------------------------------------------------------------------------------------------------------------

        kCap  = poisson.rvs(mu=dt1/tauCap , size=1)[0] # num. of BH+BH     -> BH-BH capture     events in this step
        kPair = poisson.rvs(mu=dt1/tauPair, size=1)[0] # num. of star-star -> BH-star exchange  events in this step
        kEx   = poisson.rvs(mu=dt1/tauEx  , size=1)[0] # num. of BH-star   -> BH-BH exchange    events in this step
        k3bBH = poisson.rvs(mu=dt1/tau3bBH, size=1)[0] # num. of BH+BH+BH  -> BH-BH + BH 3bb    events in this step
        kBB   = poisson.rvs(mu=dt1/tauBB  , size=1)[0] # num. of BH-BH + BH-BH            interactions in this step
        kPP   = poisson.rvs(mu=dt1/tauPP  , size=1)[0] # num. of BH-star + BH-star        interactions in this step
        
        # simulation time update:
        t = t + dt1
        
        # check if maximum time already surpassed:
        if t > tMax:
            
            print('Simulation time surpassed within a time step.')
            
            break

        # redshfit update:
        z = redd(t_lbb(zClForm)-t)

        # Cluster evolution
        # -----------------------------------------------------------------------------------------------------------------------
        # Ref. [O.Y.Gnedin et al., ApJ 785 (2014), 71.]

        # central density scale:
        rhoScale = Mhalo / ( 4*np.pi*Rscale**3 ) / (np.log(1+conc)-conc/(1+conc))

        # circular orbital velocity in m/s:
        Vcirc = veloNFW(Rgal,Rscale,rhoScale)

        # dynamical friction timescale (factor of 0.5 accounts for eccentric orbits):
        tDF = 0.45*Gyr*(Rgal/kpc)**2*(Vcirc/1e3)/(Mcl/(1e5*Msun))*0.5;

        # evolve cluster's galactrocentric radius:
        Rgal = Rgal-dt1*Rgal/2/tDF
        
        if(Rgal < 0):
            
            print('Galactocentric radius negative during time step.')
            
            break

        # auxiliary factor:
        Pfactor = 41.4*(Rgal/kpc)/(Vcirc/1e3)

        # tidal timescale:
        tTid = 10*Gyr*(Mcl/(2e5*Msun))**(2/3)*Pfactor

        # internal timescale:
        tIso = 17*Gyr*(Mcl/(2e5*Msun))
        
        # residual gas mass removed by galactic tides:
        dM_gas = M_gas * dt1 / tTid
        
        # stellar mass removed by galactic tides and internal processes:
        dMcl_stars = Mcl_stars * dt1 / np.min([tTid,tIso])
        
        # evolve gas mass:
        M_gas = M_gas - dM_gas
        
        # evolve stellar mass:
        Mcl_stars = Mcl_stars - dMcl_stars
        
        # evolve cluster mass:
        Mcl = Mcl - dMcl_stars - dM_gas
        
        if Mcl < 0 or Mcl_stars < 0:
            
            print('Cluster mass negative during time step.')
            
            break

        # moment black hole burning initiates; assume `immediately`:
        tBurn = t_cc
        
        # burning efficinency parameter (Henon, 1965):
        zetaBurn = 0.0926
        
        # evolve cluster's half-mass radius (due to BH heating):
        #rh = rh0*(1+3/2*zetaBurn*(t-tBurn)/tRel0)**(2/3)
        rh = rh * (1 + zetaBurn * dt1 / tRelax(Mcl, Mcl_stars/mAvg, rh, mAvg))

        # Ubdate number density of star-star binaries (creation by triple-star interaction)
        # -----------------------------------------------------------------------------------------------------------------------

        try:
            
            # triple-star timescale:
            tauTripleStar = 1/Rate3body(mAvg,nStar,veloDisp(mAvg,1,mAvg,Mcl,rh),5)/meanVseg

            # num. of star+star+star -> star-star events in this step:
            kTripleStar = poisson.rvs(mu=dt1/tauTripleStar,size=1)[0]

            if kTripleStar>0 and aStarStar_min < aStarStar_har: # new hard binary stars form;
            # also make sure the hardness sma is above the minimum allowed binary star separation

                for i in range(0,kTripleStar):

                    # sample hardness ratio:
                    etaValue = hardness_sampler(np.random.rand(), mAvg, mAvg, mAvg, mAvg, 0.4, etamin=5)
                    
                    # new sma:
                    smaNEW = G_Newt*mAvg/etaValue/veloDisp(mAvg,1,mAvg,Mcl,rh)**2

                    # make sure stars in binary are not very close together; if so, draw hardness ratio again:
                    while smaNEW < aStarStar_min:
                        
                        # sample hardness ratio:
                        etaValue = hardness_sampler(np.random.rand(), mAvg, mAvg, mAvg, mAvg, 0.4, etamin=5)

                        # new sma:
                        smaNEW = G_Newt*mAvg/etaValue/veloDisp(mAvg,1,mAvg,Mcl,rh)**2

                    # append new sma:
                    aStarStar = np.append(aStarStar,smaNEW)

        except:
            
            print('An exception occured while evolving the number of hard binary stars.')

        # Binary assembly via CAPTURE
        # -----------------------------------------------------------------------------------------------------------------------

        # make sure there are available BHs as they evolve in the current step:
        kCap = np.min([kCap,int(N_BH_sin/2)])
        
        if kCap > 0:

            # temporary lists, holding merger remnants and avoid sampling a BH remnant during the same step:
            mBH_temporary = []
            sBH_temporary = []
            gBH_temporary = []

            for i in range(0,kCap):

                # a binary forms:
                N_bbhs_cumul+=1

                # mass-bias array for capture; p_cap(m) ~ <Σ_cap*v_rel> ~ m^(5/2):
                p_ = (mBH/Msun)**(5/2)/((mBH/Msun)**(5/2)).sum()
                
                # biased selection of BHs captured:
                masses = np.random.choice(mBH,size=2,replace=False,p=p_)
                
                m1 = masses[0] # first  BH mass
                m2 = masses[1] # second BH mass
                
                k_1 = np.squeeze(np.where(mBH==m1))+0 # index position of first  BH
                k_2 = np.squeeze(np.where(mBH==m2))+0 # index position of second BH

                if isinstance(k_1,np.ndarray):
                    k1 = k_1[0]
                else:
                    k1 = k_1
                
                if isinstance(k_2,np.ndarray):
                    k2 = k_2[0]
                else:
                    k2 = k_2
                
                s1 = sBH[k1]  # spin parameter of first  BH
                s2 = sBH[k2]  # spin parameter of second BH
                
                g1 = int(gBH[k1]) # generation of first  BH
                g2 = int(gBH[k2]) # generation of second BH

                # Set up spin orientation parameters for PRECESSION:
                
                costheta1 = np.random.uniform(-1,1) # random value in range [-1,1]
                costheta2 = np.random.uniform(-1,1) # random value in range [-1,1]
                
                theta1 = np.arccos(costheta1) # convert angle in radians (angle between spin 1 and ang. mom.)
                theta2 = np.arccos(costheta2) # convert angle in radians (angle between spin 2 and ang. mom.)
                
                dPhi = np.random.uniform(0,2*np.pi) # random value in range [0,2*pi]

                # get remnant mass, remnant spin and GW kick velocity:
                mRem,sRem,vGW = mergerRemnant(m1,m2,s1,s2,theta1,theta2,dPhi)

                # remnant generation definition:
                gRem = np.max([g1,g2])+1

                if(g1 == 1 and g2 == 1):
                    N_me_1_1 += 1
                if((g1 == 2 and g2 == 1) or (g1 == 1 and g2 == 2)):
                    N_me_1_2 += 1
                if(g1 == 2 and g2 == 2):
                    N_me_2_2 += 1
                if((g1 == 3 and g2 == 1) or (g1 == 1 and g2 == 3)):
                    N_me_1_3 += 1
                if((g1 == 3 and g2 == 2) or (g1 == 2 and g2 == 3)):
                    N_me_2_3 += 1
                if(g1 == 3 and g2 == 3):
                    N_me_3_3 += 1

                # Relations below from arXiv:0807.2638[astro-ph]
                
                # rel. velo. sampled from Maxwellian:
                wBhBh = maxwell.rvs(loc=0,scale=veloDispRel(m1,m2,xi,mAvg,Mcl,rh)/np.sqrt(3),size=1)[0]

                # total mass:
                Mtot = m1+m2
                
                # reduced mass:
                mu = m1*m2/Mtot**2
                
                # max. impact param. for cap.:
                bMax = (340*np.pi/3)**(1/7)*Mtot*mu**(1/7)/wBhBh**(9/7) * G_Newt**(1)*c_light**(-5/7)

                # impact parameter sampled from uniform in b^2 distribution:
                b = np.sqrt(np.random.rand()*bMax**2)
                
                # pericenter distance:
                rp = b**2*wBhBh**2/2/G_Newt/Mtot
                
                # GW energy released:
                dE_gw = 85*np.pi/12/np.sqrt(2)*mu**2*Mtot**(9/2)/rp**(7/2) * G_Newt**(7/2)/c_light**5

                # final two-body energy:
                Efinal = Mtot*mu*wBhBh**2/2 - dE_gw

                # make sure eccentricity is strictly smaller than unity:
                while 1 + 2*Efinal*b**2*wBhBh**2/Mtot**3/mu/G_Newt**2<0:
                    
                    # impact parameter sampled from uniform in b^2 distribution:
                    b = np.sqrt(np.random.rand()*bMax**2)
                    
                    # pericenter distance:
                    rp = b**2*wBhBh**2/2/G_Newt/Mtot
                    
                    # GW energy released:
                    dE_gw = 85*np.pi/12/np.sqrt(2)*mu**2*Mtot**(9/2)/rp**(7/2) * G_Newt**(7/2)/c_light**5
                    
                    # final two-body energy:
                    Efinal = Mtot*mu*wBhBh**2/2 - dE_gw
                    
                # semimajor axis at formation of captured pair:
                sma = -G_Newt*Mtot**2*mu/2/Efinal
                
                # eccentricity at formation of captured pair:
                eccen = np.sqrt(1 + 2*Efinal*b**2*wBhBh**2/Mtot**3/mu/G_Newt**2)

                tCoal = T_coal(m1,m2,sma,eccen) # GW coalescence timescale
                chiEff = (m1*s1*np.cos(theta1)+m2*s2*np.cos(theta2))/(m1+m2) # effective spin parameter

                mBH = np.delete(mBH,[k1,k2]) # delete captured masses
                sBH = np.delete(sBH,[k1,k2]) # delete captured spins
                gBH = np.delete(gBH,[k1,k2]) # delete captured generations

                if vGW < vEscape(Mcl,rh): # merger remnant retained in cluster
                    
                    mBH_temporary.append(mRem) # append remnant mass
                    sBH_temporary.append(sRem) # append remnant spin
                    gBH_temporary.append(gRem) # append remnant generation
                    
                    # update number of BHs:
                    N_BH = N_BH - 1
                    
                    # increase number of retained merger remnants:
                    N_meRe+=1
                    
                else: # merger remnant ejected from cluster
                    
                    # update number of BHs:
                    N_BH = N_BH - 2
                    
                    # increase number of ejected merger remnants:
                    N_meEj+=1

                N_cap+=1 # update number of captures
                N_me +=1 # update number of mergers

                # Find location of BH progentitors in the merger tree:
                jBH_1 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m1))+0
                if jBH_1.size==0: # then this is a 1g BH
                    jBH_1 = -1
                jBH_2 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m2))+0
                if jBH_2.size==0: # then this is 1g BH
                    jBH_2 = -1

                # Append merger [channel,a,e,m1,m2,s1,s2,g1,g2,tForm,tMerge,zForm,zMerge,Nhar,Nsub,q,chiEff,theta1,theta2,dPhi,
                #                mRem,sRem,gRem,vGW,jBH_1,jBH_2]
                mergers = np.append(mergers,[[2,sma,eccen,m1,m2,s1,s2,g1,g2,t,t+tCoal,z,redd(t_lbb(zClForm)-t-tCoal),0,0,\
                    np.min([m1,m2])/np.max([m1,m2]),chiEff,theta1,theta2,dPhi,\
                                             mRem,sRem,gRem,vGW,jBH_1,jBH_2]],axis=0)

            # convert temporary lists into numpy arrays:
            mBH_temporary = np.array(mBH_temporary)
            sBH_temporary = np.array(sBH_temporary)
            gBH_temporary = np.array(gBH_temporary)

            # append merger products into bank of free BHs:
            mBH = np.append(mBH, mBH_temporary) # append remnant masses
            sBH = np.append(sBH, sBH_temporary) # append remnant spins
            gBH = np.append(gBH, gBH_temporary) # append remnant generations

        # all BHs could have depleted at this point:
        if mBH.size==0:
            break
            
        # Binary assembly via BHBHBH -> BBH
        # -----------------------------------------------------------------------------------------------------------------------

        N_BH_sin = mBH.size
        
        # make sure there are available BHs as they evolve in the current step:
        k3bBH = np.min([k3bBH,int(N_BH_sin/3)])
        
        if k3bBH > 0:

            for i in range(0,k3bBH):
           
                N_3bb += 1

                # a binary forms:
                N_bbhs_cumul+=1

                # mass-bias array for 3bb; p_3bb(m) ~ <Σ_3bb*v_rel> ~ m^(19/2):
                p_ = (mBH/Msun)**(19/2)/((mBH/Msun)**(19/2)).sum()
                
                # biased selection of BH masses:
                masses = np.random.choice(mBH,size=2,replace=False,p=p_)
                
                m1 = masses[0] # first  BH mass
                m2 = masses[1] # second BH mass
                
                k_1 = np.squeeze(np.where(mBH==m1))+0 # index position of first  BH
                k_2 = np.squeeze(np.where(mBH==m2))+0 # index position of second BH
                
                if isinstance(k_1,np.ndarray):
                    k1 = k_1[0]
                else:
                    k1 = k_1
                
                if isinstance(k_2,np.ndarray):
                    k2 = k_2[0]
                else:
                    k2 = k_2
                
                s1 = sBH[k1] # spin parameter of first  BH
                s2 = sBH[k2] # spin parameter of second BH
                
                g1 = gBH[k1] # generation of first  BH
                g2 = gBH[k2] # generation of second BH

                # sample initial hardness of newly formed BBH:
                try:
                    eta = hardness_sampler(np.random.rand(),m1,m2,meanBHmass,mAvg,xi,etamin=5)
                except:
                    eta = 5
                
                # BBH semimajor axis:
                sma = G_Newt*m1*m2/eta/mAvg/veloDisp(mAvg,1,mAvg,Mcl,rh)**2
                
                #eccen = np.sqrt(1-np.random.rand()**2) # super-thermal eccentricity
                eccen = np.sqrt(np.random.rand())      #       thermal eccentricity

                # append binary [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,zForm,Nhar,Nsub]:
                binaries = np.append(binaries,[[3,sma,eccen,m1,m2,s1,s2,g1,g2,t,z,0,0]],axis=0)
                
                # unpdate number of available BBHs in cluster:
                N_BBH+=1

                mBH = np.delete(mBH,[k1,k2]) # delete captured masses
                sBH = np.delete(sBH,[k1,k2]) # delete captured spins
                gBH = np.delete(gBH,[k1,k2]) # delete captured generations

        # Pair assembly via STAR-STAR -> BH-STAR
        # -----------------------------------------------------------------------------------------------------------------------

        # make sure there are available BHs as they evolve in the current step:
        kPair = np.min([kPair,N_BH_sin])
        
        if kPair > 0:

            for i in range(0,kPair):
           
                N_ex1 += 1

                # make sure there are single BHs in the system:
                if mBH.size<1:

                    break

                # a binary forms:
                N_bhstars_cumul+=1

                # mass-bias array for exchange; p_ex(m) ~ <Σ_ex*v_rel> ~ m^(3/2):
                p_ = (mBH/Msun)**(3/2)/((mBH/Msun)**(3/2)).sum()
                
                # biased selection of BH mass:
                m1 = np.random.choice(mBH,size=1,replace=False,p=p_)[0]
                
                # index position of the BH drawn:
                k_1 = np.squeeze(np.where(mBH==m1))+0
                
                if isinstance(k_1,np.ndarray):
                    k1 = k_1[0]
                else:
                    k1 = k_1
                
                # spin parameter of the BH:
                s1 = sBH[k1]
                
                # generation of the BH:
                g1 = gBH[k1]

                # sma-bias array for exchange; p_ex(a) ~ a:
                p_ = aStarStar/aStarStar.sum()
                
                # sampled semimajor axis value:
                smaStarBinary = np.random.choice(aStarStar,size=1,replace=False,p=p_)[0]
                
                # index position of drawn binary star:
                ksmaStarBinary = np.squeeze(np.where(aStarStar==smaStarBinary))+0
                
                # delete sampled binary-star sma:
                aStarStar = np.delete(aStarStar,ksmaStarBinary)

                # BH-star semimajor axis (from binding energy conservation during direct exchanges):
                sma = m1/mAvg * smaStarBinary

                # append pair [a,m,chi,g]:
                pairs = np.append(pairs,[[sma,m1,s1,g1]],axis=0)

                # update number of available BH-stars in cluster:
                N_BHstar+=1

                mBH = np.delete(mBH,k1) # delete mass
                sBH = np.delete(sBH,k1) # delete spin
                gBH = np.delete(gBH,k1) # delete generation
                
                # update cumulative number of BH-star pairs:
                N_bhstars_cumul+=1

        # Binary assembly via BH-STAR -> BBH
        # -----------------------------------------------------------------------------------------------------------------------

        # update number of single BHs (may have changed in this step):
        N_BH_sin = mBH.size

        # make sure there are available BH-stars or BHs as they evolve in the current step:
        kEx = np.min([kEx,N_BHstar,N_BH_sin])
        
        if kEx > 0 and mBH.size>0:

            for i in range(0,kEx):
           
                N_ex2 += 1

                # a binary forms:
                N_bbhs_cumul+=1

                # mass-bias array for exchange; p_ex(m) ~ <Σ_ex*v_rel> ~ m^(3/2):
                p_ = (mBH/Msun)**(3/2)/((mBH/Msun)**(3/2)).sum()
                
                # biased selection of BH mass:
                m2 = np.random.choice(mBH,size=1,replace=False,p=p_)[0]
                
                # index position of the BH drawn:
                k_2 = np.squeeze(np.where(mBH==m2))+0
                
                if isinstance(k_2,np.ndarray):
                    k2 = k_2[0]
                else:
                    k2 = k_2
                
                # spin parameter of the BH:
                s2 = sBH[k2]
                
                # generation of the BH:
                g2 = gBH[k2]

                # draw an available BH-star pair; index of pair (excluding j=0):
                j = int(np.random.rand()*(N_BHstar-1)+1)

                # BBH semimajor axis (from binding energy conservation during direct exchanges):
                sma = pairs[j][0]*m2/mAvg

                eccen = np.sqrt(np.random.rand()) # thermal eccentricity

                m1 = pairs[j][1] # BH mass in pair
                s1 = pairs[j][2] # BH spin in pair
                g1 = int(pairs[j][3]) # BH generation in pair

                # append binary [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,zForm,Nhar,Nsub]:
                binaries = np.append(binaries,[[1,sma,eccen,m1,m2,s1,s2,g1,g2,t,z,0,0]],axis=0)
                
                # unpdate number of available BBHs in cluster:
                N_BBH+=1
                
                # update number of available BH-stars in cluster:
                N_BHstar = N_BHstar-1

                # delete BH-star pair that exchanged into BBH:
                pairs = np.delete(pairs,j,axis=0)

                mBH = np.delete(mBH,k2) # delete mass
                sBH = np.delete(sBH,k2) # delete spin
                gBH = np.delete(gBH,k2) # delete generation


        # (some binaries formed in this step might undergo binary-binary interaction before they harden significantly)
        # Re-evaluate BBH-BBH and BHstar-BHstar interaction timescales:
        # -----------------------------------------------------------------------------------------------------------------------

        # BBH-BBH strong encounter:
        if N_BBH < 2: # need at least two BBHs to have a binary-binary encounter
            
            tauBB = 1e100*yr
        
        else:
            
            # mean BBH mass:
            meanBBHmass = np.mean(np.transpose(binaries)[:][3]+np.transpose(binaries)[:][4])
            
            # mean BBH sma:
            meanBBHsma = np.mean(np.transpose(binaries)[:][1])
            
            # number density of BBHs:
            nBBH = (N_BBH-1)/meanVseg
            
            # BBH-BBH strong encounter timescale:
            tauBB = 1/RateInter(meanBBHmass,meanBBHmass,N_BBH,nBBH,7*meanBBHsma,xi,mAvg,Mcl,rh)

        # BHstar - BHstar strong encounter:
        if N_BHstar < 2: # need at least two BH-star pairs to have a pair-pair encounter
            
            tauPP = 1e100*yr
        
        else:
            
            # mean BH-star mass:
            meanBHstarMass = np.mean(np.transpose(pairs)[:][1])+mAvg
            
            # mean BH-star sma:
            meanBHstarSma = np.mean(np.transpose(pairs)[:][0])
            
            # number density of BH-stars:
            nBHstar = (N_BHstar-1)/meanVseg
            
            # BHstar - BHstar mean binary-binary interaction timescale:
            tauPP = 1/RateInter(meanBHstarMass,meanBHstarMass,N_BHstar,nBHstar,7*meanBHstarSma,xi,mAvg,Mcl,rh)

        kBB = poisson.rvs(mu=dt1/tauBB, size=1)[0] # num. of BH-BH + BH-BH     interactions in this step
        kPP = poisson.rvs(mu=dt1/tauPP, size=1)[0] # num. of BH-star + BH-star interactions in this step

        # Pair - pair interactions and BBH assembly
        # -----------------------------------------------------------------------------------------------------------------------

        # make sure there are at least two BH-stars for a pair-pair enounter to occur:
        kPP = np.min([kPP,int(N_BHstar/2)])
        
        if kPP > 0:
            
            for i in range(0,kPP):

                N_pp += 1
            
                # array of BH-star semimajor axes (smas):
                smas = np.transpose(pairs)[:][0]
                
                # array of BH-star masses:
                BhStarMasses = np.transpose(pairs)[:][1]

                # mass- & sma-biased array:
                p_ = (BhStarMasses/Msun)**(3/2)*smas/AU/((BhStarMasses/Msun)**(3/2)*smas/AU).sum()
                
                # biased sampling of two available binaries based on their size and mass:
                a = np.random.choice(smas,size=2,replace=False,p=p_)
                
                a1 = a[0] # first  sma
                a2 = a[1] # second sma
                
                i_1 = np.squeeze(np.where(smas==a1))+0 # find index of first  BBH
                i_2 = np.squeeze(np.where(smas==a2))+0 # find index of second BBH
                
                if isinstance(i_1,np.ndarray):
                    i1 = i_1[0]
                else:
                    i1 = i_1
                
                if isinstance(i_2,np.ndarray):
                    i2 = i_2[0]
                else:
                    i2 = i_2

                m1 = pairs[i1][1]; s1 = pairs[i1][2]; g1 = pairs[i1][3] # properties of first  BH
                m2 = pairs[i2][1]; s2 = pairs[i2][2]; g2 = pairs[i2][3] # properties of second BH
                
                # sma of new BBH (from energy conservation considerations):
                if m2/a2 < m1/a1: # pair 1 is softer and controls the energy of the BBH
                    
                    aBBH = m2/mAvg*a1
                
                else: # pair 2 is softer
                    
                    aBBH = m1/mAvg*a2
                
                # eccentricity of new BBH (thermal):
                eBBH = np.sqrt(np.random.rand())

                # append binary [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,zForm,Nhar,Nsub]:
                binaries = np.append(binaries,[[1,aBBH,eBBH,m1,m2,s1,s2,g1,g2,t,z,0,0]],axis=0)

                # delete BH-star pairs:
                pairs = np.delete(pairs,[i1,i2],axis=0)

                # update number of BH-stars:
                N_BHstar = N_BHstar - 2
                
                # update number of BBHs:
                N_BBH+=1

        # Binary - binary interactions and hierarchical triple assembly
        # -----------------------------------------------------------------------------------------------------------------------
        # Adopt the `Theory of Binary-Binary Interactions` by Lyman Spitzer, Jr. and Robert D.Mathieu, ApJ 241(1980),618

        # make sure there are at least two BBHs for a binary-binary encounter to occur:
        kBB = np.min([kBB,int(N_BBH/2)])
        
        if kBB > 0:
            
            for i in range(0,kBB):

                N_bb += 1
            
                # array of BBH semimajor axes (smas):
                smas = np.transpose(binaries)[:][1]
                
                # array of BBH masses:
                BBHmasses = np.transpose(binaries)[:][3] + np.transpose(binaries)[:][4]

                # mass- & sma-biased array:
                p_ = (BBHmasses/Msun)**(3/2)*smas/AU/((BBHmasses/Msun)**(3/2)*smas/AU).sum()

                # biased sampling of two available binaries based on their size and masses:
                a = np.random.choice(smas,size=2,replace=False,p=p_)
                
                a1 = a[0] # first  sma
                a2 = a[1] # second sma
                
                i_1 = np.squeeze(np.where(smas==a1))+0 # find index of first  BBH
                i_2 = np.squeeze(np.where(smas==a2))+0 # find index of second BBH
                
                if isinstance(i_1,np.ndarray):
                    i1 = i_1[0]
                else:
                    i1 = i_1
                
                if isinstance(i_2,np.ndarray):
                    i2 = i_2[0]
                else:
                    i2 = i_2

                # find out which binary is harder and which softer of the two:
                if a1<a2:
                
                    iHard=i1
                    iSoft=i2
                
                else:
                    
                    iHard=i2
                    iSoft=i1

                m0 = binaries[iHard][3] # `inner` pair mass 0
                m1 = binaries[iHard][4] # `inner` pair mass 1

                s0 = binaries[iHard][5] # `inner` pair spin 0
                s1 = binaries[iHard][6] # `inner` pair spin 1

                g0 = int(binaries[iHard][7]) # `inner` pair generation 0
                g1 = int(binaries[iHard][8]) # `inner` pair generation 1

                # Determine params. of BH that is `freed` from bound systems:
                if binaries[iSoft][3] < binaries[iSoft][4]:
                    
                    mFreed =     binaries[iSoft][3]  ; m2 =     binaries[iSoft][4]
                    sFreed =     binaries[iSoft][5]  ; s2 =     binaries[iSoft][6]
                    gFreed = int(binaries[iSoft][7]) ; g2 = int(binaries[iSoft][8])
                
                else:
                
                    mFreed =     binaries[iSoft][4]  ; m2 =     binaries[iSoft][3]
                    sFreed =     binaries[iSoft][6]  ; s2 =     binaries[iSoft][5]
                    gFreed = int(binaries[iSoft][8]) ; g2 = int(binaries[iSoft][7])

                # semimajor axis of `inner` BBH:
                aInner = a1
                
                # outer semimajor axis determined from energy conservation during the direct exchange:
                aOuter = a2*(m0+m1)/mFreed

                # append lighter component of `soft` BBH into the arrays of single BHs in the cluster:
                mBH = np.append(mBH,mFreed)
                sBH = np.append(sBH,sFreed)
                gBH = np.append(gBH,gFreed)

                # delete `softer` BBH:
                binaries = np.delete(binaries,iSoft,axis=0)
                
                # effectively, the softer BBH is disrupted and the harder BBH may tightens a bit more:
                N_dis += 1
                
                # number of BBHs updated (dissociation of softer binary):
                N_BBH = N_BBH - 1
                
                if iHard > iSoft:
                    iHard = iHard - 1

                eInner = np.sqrt(np.random.rand()) # Inner eccen. modified here
                eOuter = np.sqrt(np.random.rand()) # Outer eccen. modified as well

                cosI1 = np.random.uniform(low=-1,high=1,size=None)
                cosI2 = np.random.uniform(low=-1,high=1,size=None)

                inclination1 = np.arccos(cosI1)
                inclination2 = np.arccos(cosI2)
                
                # angle between inner and outer ang. mom.:
                inclination = inclination1 + inclination2

                # Check if hierarchical triple is stable;
                # if not then a harder binary forms [Mardling & Aarseth (2001) criterion]:
                if aOuter/aInner > 2.8*(1+m2/(m1+m0))**(2/5)*(1+eOuter)**(2/5)/(1-eOuter)**(6/5)*(1-0.3*inclination/np.pi):

                    N_triples_cumul+=1

                    # append triple [aInner,aOuter,eInner,eInner,m0,m1,m2,s0,s1,s2,g0,g1,g2,
                    #                inclination1,inclination2,tForm,zForm]:
                    triples = np.append(triples,[[aInner,aOuter,eInner,eInner,m0,m1,m2,s0,s1,s2,g0,g1,g2,\
                        inclination1,inclination2,t,z]],axis=0)
                    
                    # update number of triples in cluster:
                    N_Triples+=1

                    # delete hard binary as well:
                    binaries = np.delete(binaries,iHard,axis=0)
                    
                    # number of BBHs updated (inner binary participates as inner pair of triple):
                    N_BBH = N_BBH - 1

                else: # inner binary is `freed`

                    # wider binary breaks
                    mBH = np.append(mBH,m2)
                    sBH = np.append(sBH,s2)
                    gBH = np.append(gBH,g2)

                    # harder binary tightens even more, up to 50% here: Zevin et al., ApJ 871 (2019), 91.
                    Delta=0.50 # fraction of energy that gets trasferred into the inner pair
                    binaries[iHard][1] = aInner/(1+Delta*mFreed*m2/m1/m0*aInner/aOuter)
                    binaries[iHard][2] = eInner

        # append evolution params [t,z,Mcl,rh,Rgal,N_BH,N_BH_sin,N_BHstar,N_BBH,N_Triples,N_me,
        #                          N_meRe,N_meEj,N_meOut,N_ZLK,N_cap,N_ej,vEsc,nStarStar,meanVseg,xi,vStar,vBH,tCap,tEx1,tEx2,t3bb,tBB,tPP,
        #                          N_ex1,N_ex2,N_3bb,N_bb,N_pp,mBHmean]:
        evolut = np.append(evolut,[[t/Myr,z,Mcl/Msun,rh/pc,Rgal/kpc,\
            N_BH,N_BH_sin,N_BHstar,N_BBH,N_Triples,N_me,N_meRe,N_meEj,N_meOut,N_ZLK,N_cap,N_ej,vEscape(Mcl,rh)/1e3,\
                nStarStar*pc**3,meanVseg/pc**3,xi,np.sqrt(0.4*G_Newt*Mcl/rh)/1e3,veloDisp(meanBHmass,xi,mAvg,Mcl,rh)/1e3,\
                                   tauCap/Myr,tauPair/Myr,tauEx/Myr,tau3bBH/Myr,tauBB/Myr,tauPP/Myr,N_ex1,N_ex2,N_3bb,N_bb,N_pp,\
                                   meanBHmass/Msun]],axis=0)

        # BBH evolution (hardening, exchanges, dissociation, ejection)
        # -----------------------------------------------------------------------------------------------------------------------

        if N_BBH > 0:

            # Shuffle binaries randomizing the nonzero rows to avoid any biases:
            np.random.shuffle(binaries[1:binaries.size])

            # initialize iteration index:
            i = 1
            
            # iterate over all available binaries (start from 1 as the the first row is a buffer with zeros):
            while i < N_BBH+1:

                # indicator for binary status (=0 binary available to evolve):
                condition = 0
                
                # initialize local timescale:
                tLocal = 0
                
                # while binary is available, evolve it:
                while condition==0:

                    # Unwrap current binary's parameters:
                    channel = int(binaries[i][0]) # binary's assembly channel
                    a       =     binaries[i][1]  # binary's semimajor axis
                    e       =     binaries[i][2]  # binary's eccentricity
                    m1      =     binaries[i][3]  # first  member's mass
                    m2      =     binaries[i][4]  # second member's mass
                    s1      =     binaries[i][5]  # first  member's spin
                    s2      =     binaries[i][6]  # second member's spin
                    g1      = int(binaries[i][7]) # first  member's generation
                    g2      = int(binaries[i][8]) # second member's generation

                    # update number of single BHs in cluster:
                    N_BH_sin = mBH.size

                    # BBH-BH interaction timescale:
                    dt2_BBHBH   = 1/RateInter(meanBHmass ,m1+m2,1,N_BH_sin/meanVseg,7*a,xi,mAvg,Mcl,rh)
                    
                    # BBH-star interaction timescale:
                    dt2_BBHstar = 1/RateInter(mAvg       ,m1+m2,1,rhoC/mAvg        ,7*a,xi,mAvg,Mcl,rh)

                    pBBHBH   = dt2_BBHstar/(dt2_BBHBH+dt2_BBHstar) # probability BBH-BH   interaction occurs
                    pBBHstar = dt2_BBHBH  /(dt2_BBHBH+dt2_BBHstar) # probability BBH-star interaction occurs

                    # draw number from [0,1]:
                    u = np.random.rand()

                    if u < pBBHstar or mBH.size==0: # BBH-star interaction happens

                        # set time step to BBH-star collision timescale:
                        dt2 = dt2_BBHstar

                        # available time window (until redshfit z=0):
                        tAvail = t_lbb(zClForm)-t-tLocal
                        
                        # merger time < BBH-single interaction time:
                        if T_coal(m1,m2,a,e) < np.min([dt1,dt2,tAvail]) and condition==0:

                            # Set up spin orientation parameters for PRECESSION:
                            
                            costheta1 = np.random.uniform(-1,1) # random value in range [-1,1]
                            costheta2 = np.random.uniform(-1,1) # random value in range [-1,1]
                            
                            theta1 = np.arccos(costheta1) # convert angle in radians (angle between spin 1 and ang. mom.)
                            theta2 = np.arccos(costheta2) # convert angle in radians (angle between spin 2 and ang. mom.)
                            
                            # random value in range [0,2*pi] (angle between projected spins on orbital plane):
                            dPhi = np.random.uniform(0,2*np.pi)

                            # get remnant mass, remnant spin and GW kick velocity:
                            mRem,sRem,vGW = mergerRemnant(m1,m2,s1,s2,theta1,theta2,dPhi)

                            # remnant generation definition:
                            gRem = int(np.max([g1,g2]))+1

                            if(g1 == 1 and g2 == 1):
                                N_me_1_1 += 1
                            if((g1 == 2 and g2 == 1) or (g1 == 1 and g2 == 2)):
                                N_me_1_2 += 1
                            if(g1 == 2 and g2 == 2):
                                N_me_2_2 += 1
                            if((g1 == 3 and g2 == 1) or (g1 == 1 and g2 == 3)):
                                N_me_1_3 += 1
                            if((g1 == 3 and g2 == 2) or (g1 == 2 and g2 == 3)):
                                N_me_2_3 += 1
                            if(g1 == 3 and g2 == 3):
                                N_me_3_3 += 1

                            # GW coalescence timescale:
                            tCoal = T_coal(m1,m2,a,e)
                            
                            # effective spin parameter:
                            chiEff = (m1*s1*np.cos(theta1)+m2*s2*np.cos(theta2))/(m1+m2)

                            if vGW < vEscape(Mcl,rh): # merger remnant retained in cluster
                                
                                mBH = np.append(mBH,mRem) # append remnant mass
                                sBH = np.append(sBH,sRem) # append remnant spin
                                gBH = np.append(gBH,gRem) # append remnant generation
                                
                                # update number of BHs:
                                N_BH = N_BH - 1
                                
                                # increase number of retained merger remnants:
                                N_meRe+=1
                                
                            else: # merger remnant ejected from cluster
                                
                                # update number of BHs:
                                N_BH = N_BH - 2
                                
                                # increase number of ejected merger remnants:
                                N_meEj+=1

                            # update number of mergers:
                            N_me +=1

                            # Find location of BH progentitors in the merger tree:
                            jBH_1 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m1))+0
                            if jBH_1.size==0: # then this is a 1g BH
                                jBH_1 = -1
                            jBH_2 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m2))+0
                            if jBH_2.size==0: # then this is 1g BH
                                jBH_2 = -1

                            # Append merger [channel,a,e,m1,m2,s1,s2,g1,g1,tForm,tMerge,zForm,
                            #                zMerge,Nhar,Nsub,q,chiEff,theta1,theta2,dPhi,mRem,sRem,gRem,vGW,j1,j2]:
                            mergers = np.append(mergers,[[channel,a,e,m1,m2,s1,s2,g1,g2,\
                                binaries[i][9],t+tLocal+T_coal(m1,m2,a,e),binaries[i][10],\
                                    redd(t_lbb(zClForm)-t-tLocal-T_coal(m1,m2,a,e)),binaries[i][11],binaries[i][12],\
                                        np.min([m1,m2])/np.max([m1,m2]),chiEff,theta1,theta2,dPhi,\
                                                         mRem,sRem,gRem,vGW,jBH_1,jBH_2]],axis=0)

                            binaries = np.delete(binaries,i,axis=0) # delete binary from BBH list
                            i = i - 1 # update running iteration index
                            N_BBH = N_BBH - 1 # update number of BBHs
                            condition=2 # binary merges
                            break

                        else: # binary-single collision occurs

                            # sample velocity of binary before the BBH-star interaction:
                            vB = maxwell.rvs(loc=0,scale=veloDisp(m1+m2,xi,mAvg,Mcl,rh)/np.sqrt(3),size=1)
                            
                            # sample velocity of single before the BBH-star interaction:
                            vS = maxwell.rvs(loc=0,scale=veloDisp(mAvg ,1 ,mAvg,Mcl,rh)/np.sqrt(3),size=1)

                            thetaBefore = np.random.uniform(0,np.pi) # sample relative angle before the collision
                            thetaAfter  = np.random.uniform(0,np.pi) # sample relative angle after  the collision

                            # reduced mass in system before the encounter:
                            mu = (m1+m2)*mAvg/(m1+m2+mAvg)
                            
                            # relative velocity squared before the interaction:
                            vRel2 = vB**2 + vS**2 - 2*vB*vS*np.cos(thetaBefore)

                            # check if single is energetic enough to ionize the binary:
                            if mu*vRel2 > G_Newt*m1*m2/a: # binary ionizes
                                
                                # update number of BBH disruptions:
                                N_dis+=1
                                
                                # delete binary from BBH list:
                                binaries = np.delete(binaries,i,axis=0)
                                
                                # update running iteration index:
                                i = i - 1
                                
                                mBH = np.append(mBH,m1); mBH = np.append(mBH,m2) # append masses
                                sBH = np.append(sBH,s1); sBH = np.append(sBH,s2) # append spins
                                gBH = np.append(gBH,g1); gBH = np.append(gBH,g2) # append generations
                                
                                # update number of BBHs:
                                N_BBH = N_BBH - 1
                                
                                # binary disrupts:
                                condition=4
                                
                                break

                            if tLocal > dt1: # if total time exceeds parent iteration time step
                                
                                # end local iteration:
                                condition=1
                                
                                break

                            # reduced mass of binary-single system after the collision (here it is the same):
                            muNEW = mu

                            # binary hardness semimajor axis:
                            aHard = G_Newt*np.min([m1,m2])/4/veloDisp(mAvg,1,mAvg,Mcl,rh)**2

                            # hardening rate:
                            Hrate = H_rate(a/aHard,  np.min([m1,m2])/np.max([m1,m2]))
                            
                            # eccentricity growth rate:
                            Krate = K_rate(a/aHard,e,np.min([m1,m2])/np.max([m1,m2]))

                            # Hardening evolution lines:
                            aNEW = a / (1+Hrate/2/np.pi*mAvg/(m1+m2))
                            eNEW = e + Krate/(1+2*np.pi/Hrate*(m1+m2)/mAvg)
                            binaries[i][1] = aNEW # update semimajor axis
                            binaries[i][2] = eNEW # update eccentricity
                            
                            # number of binary-single hardening interactions updated:
                            binaries[i][11]+=1

                            # update local time:
                            tLocal+=dt2

                            # check is binary is catapulted out of the cluster after the collision:
                            if a<aEj(m1,m2,mAvg,Hrate,vEscape(Mcl,rh)): # binary kicked out of cluster

                                # Ejected BBH may merge outside cluster in the field
                                
                                # check if ej. bin. merges in the available time
                                if t+tLocal+T_coal(m1,m2,a,e) < t_lbb(zClForm):

                                    costheta1 = np.random.uniform(-1,1) # random value in range [-1,1]
                                    costheta2 = np.random.uniform(-1,1) # random value in range [-1,1]
                                    
                                    theta1 = np.arccos(costheta1) # convert angle in radians (angle between spin 1 and ang. mom.)
                                    theta2 = np.arccos(costheta2) # convert angle in radians (angle between spin 2 and ang. mom.)
                                    
                                    dPhi = np.random.uniform(0,2*np.pi) # random value in range [0,2*pi]

                                    # get remnant mass, remnant spin and GW kick velocity:
                                    mRem,sRem,vGW = mergerRemnant(m1,m2,s1,s2,theta1,theta2,dPhi)

                                    # effective spin parameter:
                                    chiEff = (s1*m1*np.cos(theta1)+s2*m2*np.cos(theta2))/(m1+m2)

                                    # remnant generation:
                                    gRem = np.max([g1, g2]) + 1

                                    # update tot. number of binaries that merge:
                                    N_me += 1
                                    
                                    # update number of BBHs that merge outside of cluster:
                                    N_meOut += 1
                                    
                                    if(g1 == 1 and g2 == 1):
                                        N_me_1_1 += 1
                                    if((g1 == 2 and g2 == 1) or (g1 == 1 and g2 == 2)):
                                        N_me_1_2 += 1
                                    if(g1 == 2 and g2 == 2):
                                        N_me_2_2 += 1
                                    if((g1 == 3 and g2 == 1) or (g1 == 1 and g2 == 3)):
                                        N_me_1_3 += 1
                                    if((g1 == 3 and g2 == 2) or (g1 == 2 and g2 == 3)):
                                        N_me_2_3 += 1
                                    if(g1 == 3 and g2 == 3):
                                        N_me_3_3 += 1

                                    # Find location of BH progentitors in the merger tree:
                                    jBH_1 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m1))+0
                                    if jBH_1.size==0: # then this is a 1g BH
                                        jBH_1 = -1
                                    jBH_2 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m2))+0
                                    if jBH_2.size==0: # then this is 1g BH
                                        jBH_2 = -1

                                    # Append merger [channel,a,e,m1,m2,s1,s2,g1,g1,tForm,tMerge,zForm,zMerge,
                                    #                Nhar,Nsub,q,chiEff,theta1,theta2,dPhi,mRem,sRem,gRem,vGW,j1,j2],
                                    # (negative `channel` means BBH merges outside outside):
                                    mergers = np.append(mergers,[[-channel,a,e,m1,m2,s1,s2,g1,g2,\
                                        binaries[i][9],t+tLocal+T_coal(m1,m2,a,e),binaries[i][10],\
                                            redd(t_lbb(zClForm)-t-tLocal-T_coal(m1,m2,a,e)),\
                                                binaries[i][11],binaries[i][12],np.min([m1,m2])/np.max([m1,m2]),chiEff,\
                                                    theta1,theta2,dPhi,\
                                                                 mRem,sRem,gRem,vGW,jBH_1,jBH_2]],axis=0)

                                # delete binary from BBH list:
                                binaries = np.delete(binaries,i,axis=0)
                                
                                # update running iteration index:
                                i = i - 1

                                # update number of BBHs in cluster:
                                N_BBH = N_BBH - 1
                                
                                # update total number of BHs in cluster:
                                N_BH = N_BH - 2
                                
                                # update number of BBH ejections:
                                N_ej+=1
                                
                                # BBH-BH ejection:
                                condition=5

                                break

                    else: # BBH-BH interaction happens

                        # set time step to BBH-BH collision timescale:
                        dt2 = dt2_BBHBH

                        # available time window (until redshfit z=0):
                        tAvail = t_lbb(zClForm)-t-tLocal
                        
                        # merger time < BBH-single interaction time:
                        if T_coal(m1,m2,a,e) < np.min([dt1,dt2,tAvail]) and condition==0:

                            # Set up spin orientation parameters for PRECESSION:
                            
                            costheta1 = np.random.uniform(-1,1) # random value in range [-1,1]
                            costheta2 = np.random.uniform(-1,1) # random value in range [-1,1]
                            
                            theta1 = np.arccos(costheta1) # convert angle in radians (angle between spin 1 and ang. mom.)
                            theta2 = np.arccos(costheta2) # convert angle in radians (angle between spin 2 and ang. mom.)
                            
                            dPhi = np.random.uniform(0,2*np.pi) # random value in range [0,2*pi]

                            # get remnant mass, remnant spin and GW kick velocity:
                            mRem,sRem,vGW = mergerRemnant(m1,m2,s1,s2,theta1,theta2,dPhi)

                            # remnant generation definition:
                            gRem = int(np.max([g1,g2]))+1

                            if(g1 == 1 and g2 == 1):
                                N_me_1_1 += 1
                            if((g1 == 2 and g2 == 1) or (g1 == 1 and g2 == 2)):
                                N_me_1_2 += 1
                            if(g1 == 2 and g2 == 2):
                                N_me_2_2 += 1
                            if((g1 == 3 and g2 == 1) or (g1 == 1 and g2 == 3)):
                                N_me_1_3 += 1
                            if((g1 == 3 and g2 == 2) or (g1 == 2 and g2 == 3)):
                                N_me_2_3 += 1
                            if(g1 == 3 and g2 == 3):
                                N_me_3_3 += 1

                            # GW coalescence timescale:
                            tCoal = T_coal(m1,m2,a,e)
                            
                            # effective spin parameter:
                            chiEff = (m1*s1*np.cos(theta1)+m2*s2*np.cos(theta2))/(m1+m2)

                            if vGW < vEscape(Mcl,rh): # merger remnant retained in cluster
                                
                                mBH = np.append(mBH,mRem) # append remnant mass
                                sBH = np.append(sBH,sRem) # append remnant spin
                                gBH = np.append(gBH,gRem) # append remnant generation
                                
                                # update number of BHs:
                                N_BH = N_BH - 1
                                
                                # increase number of retained merger remnants:
                                N_meRe+=1
                                
                            else: # merger remnant ejected from cluster
                                
                                # update number of BHs:
                                N_BH = N_BH - 2
                                
                                # increase number of ejected merger remnants:
                                N_meEj+=1

                            # update number of mergers:
                            N_me +=1

                            # Find location of BH progentitors in the merger tree:
                            jBH_1 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m1))+0
                            
                            if jBH_1.size==0: # then this is a 1g BH
                                jBH_1 = -1
                            jBH_2 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m2))+0
                            if jBH_2.size==0: # then this is 1g BH
                                jBH_2 = -1

                            # Append merger [channel,a,e,m1,m2,s1,s2,g1,g1,tForm,
                            #                tMerge,zForm,zMerge,Nhar,Nsub,q,chiEff,theta1,theta2,dPhi,mRem,sRem,gRem,vGW,j1,j2]:
                            mergers = np.append(mergers,[[channel,a,e,m1,m2,s1,s2,g1,g2,\
                                binaries[i][9],t+tLocal+T_coal(m1,m2,a,e),binaries[i][10],\
                                    redd(t_lbb(zClForm)-t-tLocal-T_coal(m1,m2,a,e)),binaries[i][11],binaries[i][12],\
                                        np.min([m1,m2])/np.max([m1,m2]),chiEff,theta1,theta2,dPhi,\
                                                         mRem,sRem,gRem,vGW,jBH_1,jBH_2]],axis=0)

                            # delete binary from BBH list:
                            binaries = np.delete(binaries,i,axis=0)
                            
                            # update running iteration index:
                            i = i - 1
                            
                            # update number of BBHs:
                            N_BBH = N_BBH - 1
                            
                            # binary merges:
                            condition=2
                            
                            break

                        else: # binary-single collision occurs

                            # mass-bias array for interaction; p_int(m) ~ <Σ_int*v_rel> ~ m^(3/2):
                            p_ = (mBH/Msun)**(3/2)/((mBH/Msun)**(3/2)).sum()
                            
                            # biased sampling of available single BH mass:
                            m3 = np.random.choice(mBH,size=1,replace=False,p=p_)
                            
                            # index of sampled mass:
                            k3 = np.squeeze(np.where(mBH==m3))+0
                            
                            if isinstance(k3,np.ndarray):
                                k3 = k3[0]
                            
                            s3 =     sBH[k3]  # spin parameter of third BH
                            g3 = int(gBH[k3]) # generation     of third BH

                            # sample velocity of binary before the BBH-BH interaction:
                            vB = maxwell.rvs(loc=0,scale=veloDisp(m1+m2,xi,mAvg,Mcl,rh)/np.sqrt(3),size=1)
                            
                            # sample velocity of single before the BBH-BH interaction:
                            vS = maxwell.rvs(loc=0,scale=veloDisp(m3   ,xi,mAvg,Mcl,rh)/np.sqrt(3),size=1)

                            thetaBefore = np.random.uniform(0,np.pi) # sample relative angle before the collision
                            thetaAfter  = np.random.uniform(0,np.pi) # sample relative angle after  the collision

                            # reduced mass in system before the encounter:
                            mu = (m1+m2)*m3/(m1+m2+m3)
                            
                            # relative velocity squared before the interaction:
                            vRel2 = vB**2 + vS**2 - 2*vB*vS*np.cos(thetaBefore)

                            # check if single is energetic enough to ionize the binary:
                            if mu*vRel2 > G_Newt*m1*m2/a: # binary ionizes
                                
                                # update number of BBH disruptions:
                                N_dis+=1
                                
                                # delete binary from BBH list:
                                binaries = np.delete(binaries,i,axis=0)
                                
                                # update running iteration index:
                                i = i - 1
                                
                                mBH = np.append(mBH,m1); mBH = np.append(mBH,m2) # append masses
                                sBH = np.append(sBH,s1); sBH = np.append(sBH,s2) # append spins
                                gBH = np.append(gBH,g1); gBH = np.append(gBH,g2) # append generations
                                
                                # update number of BBHs:
                                N_BBH = N_BBH - 1
                                
                                # binary disrupts:
                                condition=4
                                
                                break

                            if tLocal > dt1: # if total time exceeds parent iteration time step
                                
                                # end local iteration:
                                condition=1
                                
                                break

                            if m3<np.min([m1,m2]): # fly-by encounter (based on J.G.Hills and L.W.Fullerton, AJ 85, 1281 (1980))

                                # reduced mass of binary-single system after the collision (here it is the same):
                                muNEW = mu

                                # binary hardness semimajor axis:
                                aHard = G_Newt*np.min([m1,m2])/4/veloDisp(meanBHmass,xi,mAvg,Mcl,rh)**2

                                # hardening rate:
                                Hrate = H_rate(a/aHard,  np.min([m1,m2])/np.max([m1,m2]))
                                
                                # eccentricity growth rate:
                                Krate = K_rate(a/aHard,e,np.min([m1,m2])/np.max([m1,m2]))

                                # Hardening evolution lines:
                                aNEW = a / (1+Hrate/2/np.pi*m3/(m1+m2))
                                eNEW = e + Krate/(1+2*np.pi/Hrate*(m1+m2)/m3)
                                binaries[i][1] = aNEW # update semimajor axis
                                binaries[i][2] = eNEW # update eccentricity
                                
                                # number of binary-single hardening interactions updated:
                                binaries[i][11]+=1

                                # update local time:
                                tLocal+=dt2

                                # energy extracted in the reduced body system:
                                y = G_Newt*m1*m2*(1/aNEW-1/a)/2
                                
                                # new single velocity, after the fly-by:
                                vSnew = np.sqrt(vS**2+2*y/m3)

                                # check if single is catapulted out of the cluster after the collision:
                                if vSnew > vEscape(Mcl,rh):

                                    mBH = np.delete(mBH,k3) # delete third BH mass
                                    sBH = np.delete(sBH,k3) # delete third BH spin
                                    gBH = np.delete(gBH,k3) # delete third BH generation
                                    
                                    # one BH ejected away:
                                    N_BH = N_BH - 1
                                    
                                    # count number of ejected single BHs from BBHs:
                                    N_sinBHeje += 1

                                # check is binary is catapulted out of the cluster after the collision:
                                if a<aEj(m1,m2,m3,Hrate,vEscape(Mcl,rh)): # binary kicked out of cluster

                                    # Ejected BBH may merge outside cluster in the field.
                                    
                                    # check if ej. bin. merges in the available time:
                                    if t+tLocal+T_coal(m1,m2,a,e) < t_lbb(zClForm):

                                        costheta1 = np.random.uniform(-1,1) # random value in range [-1,1]
                                        costheta2 = np.random.uniform(-1,1) # random value in range [-1,1]
                                        
                                        # convert angle in radians (angle between spin 1 and ang. mom.):
                                        theta1 = np.arccos(costheta1)
                                        
                                        # convert angle in radians (angle between spin 2 and ang. mom.):
                                        theta2 = np.arccos(costheta2)
                                        
                                        # random value in range [0,2*pi]:
                                        dPhi = np.random.uniform(0,2*np.pi)

                                        # effective spin parameter:
                                        chiEff = (s1*m1*np.cos(theta1)+s2*m2*np.cos(theta2))/(m1+m2)

                                        # update tot. number of binaries that merge:
                                        N_me += 1
                                        
                                        # update number of BBHs that merge outside of cluster:
                                        N_meOut += 1
                                        
                                        if(g1 == 1 and g2 == 1):
                                            N_me_1_1 += 1
                                        if((g1 == 2 and g2 == 1) or (g1 == 1 and g2 == 2)):
                                            N_me_1_2 += 1
                                        if(g1 == 2 and g2 == 2):
                                            N_me_2_2 += 1
                                        if((g1 == 3 and g2 == 1) or (g1 == 1 and g2 == 3)):
                                            N_me_1_3 += 1
                                        if((g1 == 3 and g2 == 2) or (g1 == 2 and g2 == 3)):
                                            N_me_2_3 += 1
                                        if(g1 == 3 and g2 == 3):
                                            N_me_3_3 += 1

                                        # Find location of BH progentitors in the merger tree:
                                        jBH_1 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m1))+0
                                        if jBH_1.size==0: # then this is a 1g BH
                                            jBH_1 = -1
                                        jBH_2 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m2))+0
                                        if jBH_2.size==0: # then this is 1g BH
                                            jBH_2 = -1

                                        # get remnant mass, remnant spin and GW kick velocity:
                                        mRem,sRem,vGW = mergerRemnant(m1,m2,s1,s2,theta1,theta2,dPhi)

                                        # merger remnant generation:
                                        gRem = np.max([g1, g2]) + 1
                                        
                                        # Append merger [channel,a,e,m1,m2,s1,s2,g1,g1,tForm,tMerge,zForm,
                                        #                zMerge,Nhar,Nsub,q,chiEff,theta1,theta2,dPhi,mRem,sRem,gRem,vGW,j1,j2],
                                        # (negative `channel` means BBH merges outside outside):
                                        mergers = np.append(mergers,[[-channel,a,e,m1,m2,s1,s2,g1,g2,\
                                            binaries[i][9],t+tLocal+T_coal(m1,m2,a,e),binaries[i][10],\
                                                redd(t_lbb(zClForm)-t-tLocal-T_coal(m1,m2,a,e)),\
                                                    binaries[i][11],binaries[i][12],np.min([m1,m2])/np.max([m1,m2]),chiEff,\
                                                        theta1,theta2,dPhi,\
                                                                     mRem,sRem,gRem,vGW,jBH_1,jBH_2]],axis=0)

                                    # delete binary from BBH list:
                                    binaries = np.delete(binaries,i,axis=0)
                                    
                                    # update running iteration index:
                                    i = i - 1

                                    # update number of BBHs in cluster:
                                    N_BBH = N_BBH - 1
                                    
                                    # update total number of BHs in cluster:
                                    N_BH = N_BH - 2
                                    
                                    # update number of BBH ejections:
                                    N_ej+=1
                                    
                                    # BBH-BH ejection:
                                    condition=5

                                    break

                            else: # exchange occurs
                                
                                # update cumulative number of exchanges occured in the cluster:
                                N_exch+=1

                                if(m1<m2):
                                    
                                    ms = m1; ss = s1; gs = g1 # exchanged for m3
                                    mr = m2; sr = s2; gr = g2 # retained in binary
                                
                                else:
                                    
                                    ms = m2; ss = s2; gs = g2 # exchanged for m3
                                    mr = m1; sr = s1; gr = g1 # retained in binary

                                # impose binding energy conservation after substitution:
                                binaries[i][1] = m3/ms*binaries[i][1]
                                
                                # eccentricity re-sampled:
                                binaries[i][2] = np.sqrt(np.random.rand())

                                mBH = np.delete(mBH,k3) # delete exchanged mass
                                sBH = np.delete(sBH,k3) # delete exchanged spin
                                gBH = np.delete(gBH,k3) # delete exchanged generation

                                mBH = np.append(mBH,ms) # append released mass
                                sBH = np.append(sBH,ss) # append released spin
                                gBH = np.append(gBH,gs) # append released generation

                                # update binary masses:
                                binaries[i][3] = mr
                                binaries[i][4] = m3

                                # update binary spins:
                                binaries[i][5] = sr
                                binaries[i][6] = s3

                                # update binaries generations:
                                binaries[i][7] = int(gr)
                                binaries[i][8] = int(g3)

                                # number of exchanges increased:
                                binaries[i][12]+=1

                # update iteration index:
                i+=1

        # Triple evolution
        # -----------------------------------------------------------------------------------------------------------------------
        # Ref. [Miller, M. C., & Hamilton, D. P. 2002, ApJ, 576, 894].
        # Triples don't survive for long in dynamical systems, however there's a slim chance they merge before the next interaction

        if N_Triples > 0:

            # Shuffle triples randomizing the nonzero rows to avoid any biases:
            np.random.shuffle(triples[1:triples.size])

            # initialize iteration index:
            i = 1
            
            # iterate over all available triples (start from 1 as the the first row is a buffer with zeros):
            while i < N_Triples+1:

                # indicator for triple status (=0 triple available to evolve):
                condition = 0
                
                # initialize local timescale:
                tLocal = 0
                
                # while triple is available, evolve it:
                while condition==0:

                    # Unwrap current triples's parameters:
                    aIn          = triples[i][0 ]
                    aOut         = triples[i][1 ]
                    eIn          = triples[i][2 ]
                    eOut         = triples[i][3 ]
                    m0           = triples[i][4 ]
                    m1           = triples[i][5 ]
                    m2           = triples[i][6 ]
                    s0           = triples[i][7 ]
                    s1           = triples[i][8 ]
                    s2           = triples[i][9 ]
                    g0           = triples[i][10]
                    g1           = triples[i][11]
                    g2           = triples[i][12]
                    inclination1 = triples[i][13]
                    inclination2 = triples[i][14]

                    m01  = m0  + m1 # inner pair mass
                    m012 = m01 + m2 # total triple mass

                    # epsilon paramerer of inner pair:
                    epsilon = 1 - eIn**2

                    muIn  = m0*m1 /m01  # reduced mass of inner pair
                    muOut = m2*m01/m012 # reduced mass of outer pair

                    # auxiliary params:
                    beta  = muOut*np.sqrt(m012)/muIn/np.sqrt(m01)*np.sqrt(aOut/aIn*(1-eOut**2))
                    alpha = np.sqrt(epsilon)*np.cos(inclination1)+beta*np.cos(inclination2)

                    # total inclination:
                    inclination = inclination1 + inclination2
                    
                    # cosine of tot. inclination:
                    cosInclination = (alpha**2-beta**2-epsilon)/2/beta/np.sqrt(epsilon)

                    # proxy for max. eccentricity of inner pair, at the Newtonian order:
                    epsilonMIN = 5/3*np.cos(inclination)**2

                    # ZLK merger timescale:
                    tauZLKmerger = 5e17*yr*(Msun**3/m01**2/muIn)*(aIn/AU)**4*epsilonMIN**3
                    
                    # Hill's approx.: regard inner and outer orbits as two nested binaries
                    dt3 = 1/(RateInter(m012,mAvg,1,nStar,7*aOut,xi,mAvg,Mcl,rh)+\
                        RateInter(m012,meanBHmass,1,N_BH/meanVseg,7*aOut,xi,mAvg,Mcl,rh))
                    aHard = G_Newt*np.min(np.array([m2,m01]))/4/veloDisp(mAvg,1,mAvg,Mcl,rh)**2
                    Hrate  = H_rate(aOut/aHard,     np.min([m01,m2])/np.max([m01,m2]))
                    Krate  = K_rate(aOut/aHard,eOut,np.min([m01,m2])/np.max([m01,m2]))

                    # check whether inner pair merges before next interaction with another object
                    # and GR precession does not destroy ZLK:
                    ZLK_not_Destroyed_By_GR_precession = aOut/aIn<34*(aIn/1e-2/pc)**(1/3)*((m0+m1)/2e6/Msun)**(-1/3)\
                        *(2*m2/(m0+m1))**(1/3)*((1-eIn**2)/(1-eOut**2))**(1/2)
                    
                    if tauZLKmerger < dt3 and condition == 0 and ZLK_not_Destroyed_By_GR_precession:

                        # Set up spin orientation parameters for PRECESSION:
                        
                        costheta0 = np.random.uniform(-1,1) # random value in range [-1,1] (angle between spin 1 and ang. mom.)
                        costheta1 = np.random.uniform(-1,1) # random value in range [-1,1]
                        
                        theta0 = np.arccos(costheta0) # convert angle in radians
                        theta1 = np.arccos(costheta1) # convert angle in radians
                        
                        dPhi = np.random.uniform(0,2*np.pi) # random value in range [0,2*pi]

                        # get remnant mass, remnant spin and GW kick velocity:
                        mRem,sRem,vGW = mergerRemnant(m0,m1,s0,s1,theta0,theta1,dPhi)

                        # remnant generation definition:
                        gRem = int(np.max([g0,g1]))+1

                        if(g0 == 1 and g1 == 1):
                            N_me_1_1 += 1
                        if((g0 == 2 and g1 == 1) or (g0 == 1 and g1 == 2)):
                            N_me_1_2 += 1
                        if(g0 == 2 and g1 == 2):
                            N_me_2_2 += 1
                        if((g0 == 3 and g1 == 1) or (g0 == 1 and g1 == 3)):
                            N_me_1_3 += 1
                        if((g0 == 3 and g1 == 2) or (g0 == 2 and g1 == 3)):
                            N_me_2_3 += 1
                        if(g0 == 3 and g1 == 3):
                            N_me_3_3 += 1

                        # GW coalescence timescale:
                        tCoal = T_coal(m0,m1,aIn,eIn)
                        
                        # effective spin parameter:
                        chiEff = (m0*s0*np.cos(theta0)+m1*s1*np.cos(theta1))/(m0+m1)

                        # update number of tot. mergers:
                        N_me+=1
                        
                        # update number of ZLK mergers:
                        N_ZLK+=1

                        # delete triple with merging inner pair:
                        triples = np.delete(triples,i,axis=0)
                        
                        # update running iteration index:
                        i = i - 1
                        
                        # number of triples updated:
                        N_Triples = N_Triples - 1
                        
                        # (when a merger occurs, the number of BHs reduces by 1)
                        # Update number of BHs in cluster:
                        N_BH = N_BH - 1

                        # Find location of BH progentitors in the merger tree:
                        jBH_0 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m0))+0
                        if jBH_0.size==0: # then this is a 1g BH
                            jBH_0 = -1
                        jBH_1 = np.squeeze(np.where(np.transpose(mergers)[:][20]==m1))+0
                        if jBH_1.size==0: # then this is 1g BH
                            jBH_1 = -1

                        # append mergers [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,tMerge,zForm,zMerge,
                        #                 Nhar,Nsub,q,chiEff,theta1,theta2,dPhi,mRem,sRem,gRem,vGW,j1,j2]:
                        mergers = np.append(mergers,[[4,aIn,np.sqrt(1-epsilonMIN),m0,m1,s0,s1,g0,g1,\
                            triples[i][15],t+tLocal+tauZLKmerger,triples[i][16],redd(t_lbb(zClForm)-t-tLocal-tauZLKmerger),\
                                0,0,np.min([m0,m1])/np.max([m0,m1]),chiEff,theta0,theta1,dPhi,\
                                                     mRem,sRem,gRem,vGW,jBH_0,jBH_1]],axis=0)

                        # new outer binary disrupts due to merger kick:
                        if vGW > np.sqrt(1*G_Newt*(m2+mRem)/aOut):
                            
                            # update number of disruptions:
                            N_dis+=1

                            # append outer BH params. into single list:
                            mBH = np.append(mBH,m2)
                            sBH = np.append(sBH,s2)
                            gBH = np.append(gBH,g2)

                            if vGW > vEscape(Mcl,rh): # BH remnant ejected from cluster
                                
                                # further reduce the number of BHs in the cluster due to ejection:
                                N_BH = N_BH - 1
                                
                                # update number of mergers ejected:
                                N_meEj+=1
                                
                            else:
                                
                                # update number of mergers retained:
                                N_meRe+=1

                                # append inner merger remnant BH params. into single list:
                                mBH = np.append(mBH,mRem)
                                sBH = np.append(sBH,sRem)
                                gBH = np.append(gBH,gRem)

                        else: # new outer binary survives
                            
                            # append binary [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,zForm,Nhar,Nsub]:
                            binaries = np.append(binaries,[[5,aOut,np.sqrt(np.random.rand()),mRem,m2,sRem,s2,gRem,g2,\
                                t+tLocal+tauZLKmerger,redd(t_lbb(zClForm)-t-tLocal-tauZLKmerger),0,0]],axis=0)
                            
                            # update number of BBHs:
                            N_BBH+=1
                            
                            # update number of mergers retained:
                            N_meRe+=1

                        # inner pair of triple merges:
                        condition=1
                        
                        break

                    if tLocal > dt1 and condition==0:
                        
                        # local time range depleted:
                        condition=2
                        
                        break

                    # Hardening lines
                    triples[i][1] = triples[i][1] / (1+Hrate/2/np.pi*mAvg/m012) # outer semimajor axis evolution
                    triples[i][3] = triples[i][3] + Krate/(1+2*np.pi/Hrate*m012/mAvg) # outer eccentricity evolution
                    
                    # and the triple reorients:
                    triples[i][13] = np.arccos(np.random.uniform(low=-1, high=1, size=None))
                    triples[i][14] = np.arccos(np.random.uniform(low=-1, high=1, size=None))
                    
                    # update local time:
                    tLocal+=dt3

                    # Triple ionized:
                    
                    # critical velocity for breakup:
                    vCritical = np.sqrt(G_Newt*m01*m2*(m012+mAvg)/aOut/m012/mAvg)
                    
                    # star velo. from Maxwellian:
                    v_star = maxwell.rvs(loc=0,scale=veloDisp(mAvg,1,mAvg,Mcl,rh)/np.sqrt(3),size=1)[0]
                    
                    # check if star velocity drawn exceeds the critical value:
                    if v_star > vCritical and condition == 0:

                        # delete triple:
                        triples = np.delete(triples,i,axis=0)
                        
                        # update running iteration index:
                        i = i - 1
                        
                        # update number of triples:
                        N_Triples = N_Triples - 1

                        # append BBH: [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,zForm,Nhar,Nsub]
                        binaries = np.append(binaries,[[5,aIn,eIn,m0,m1,s0,s1,g0,g1,t+tLocal,\
                            redd(t_lbb(zClForm)-t-tLocal),0,0]],axis=0) # inner binary is `freed`
                        
                        # update number of BBHs:
                        N_BBH+=1

                        # append tertiary parms. into list of single BHs:
                        # (outer BH is `freed`)
                        mBH = np.append(mBH,m2)
                        sBH = np.append(sBH,s2)
                        gBH = np.append(gBH,g2)

                        # update number of disruptions:
                        N_dis+=1
                        
                        # triple destroyed by energetic encounter:
                        condition=3
                        
                        break

                    # Triple is not hierarchical anymore and breaks into binary+single
                    # Triple's outer orbit hardens and triple may become non-hierarchical
                    if aOut<aIn*3.3/(1-eOut)*(2/3*(1+m2/(m01))*(1+eOut)/np.sqrt(1-eOut))**(2/5)\
                        *(1-0.3*inclination/np.pi) and condition == 0:

                        # delete triple:
                        triples = np.delete(triples,i,axis=0)
                        
                        # update running iteration index:
                        i = i - 1
                        
                        # update number of triples:
                        N_Triples = N_Triples - 1

                        # append BBH: [channel,a,e,m1,m2,chi1,chi2,g1,g2,tForm,zForm,Nhar,Nsub]
                        binaries = np.append(binaries,[[5,aIn,eIn,m0,m1,s0,s1,g0,g1,t+tLocal,\
                            redd(t_lbb(zClForm)-t-tLocal),0,0]],axis=0) # inner binary is `freed`

                        # update number of BBHs:
                        N_BBH+=1

                        # append tertiary parms. into list of single BHs:
                        # (outer BH is `freed`)
                        mBH = np.append(mBH,m2)
                        sBH = np.append(sBH,s2)
                        gBH = np.append(gBH,g2)

                        # update number of disruptions:
                        N_dis+=1
                        
                        # triple destroyed by energetic encounter:
                        condition=4
                        
                        break

                # update iteration index:
                i+=1

        # update number of single BHs in cluster (before printing information):
        N_BH_sin = mBH.size

        # Print information on each time step:
        print('Iteration step: #',int(Niter))
        print('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
        print('t/Myr','dt/Myr','z','Mcl/Msun','rh/pc','Rgal/kpc','Rseg/pc',sep='    ')
        print(format(t/Myr,'.2f'),'   ',format(dt1/Myr,'.2f'),'   ',format(z,'.2f'),'  ',format(Mcl/Msun,'.0f'),'  ',\
            format(rh/pc,'.2f'),'   ',format(Rgal/kpc,'.2f'),\
                '   ',format(Rseg(meanBHmass,meanBHmass*N_BH,xi,mAvg,Mcl,rh)/pc,'.3f'))
        print()
        print('N_BH','N_BH_sin','N_BHstar','N_BBH','N_Triples','N_me','N_cap','N_ZLK',sep='   ')
        print(N_BH,'      ',N_BH_sin,'      ',N_BHstar,'      ',N_BBH,'      ',N_Triples,'      ',N_me,'      ',N_cap,'      ',N_ZLK)
        print()
        print('N_ej','N_dis','N_exch','N_bhstars_cumul','N_bbhs_cumul','N_triples_cumul',sep='   ')
        print(N_ej,'      ',N_dis,'      ',N_exch,'      ',N_bhstars_cumul,'      ',N_bbhs_cumul,'      ',N_triples_cumul)
        print()
        print('kCap','kPair','kEx','k3bBH','kBB','kPP',sep='    ')
        print(kCap,'      ',kPair,'      ',kEx,'     ',k3bBH,'     ',kBB,'     ',kPP)
        print()
        print('N_sinBHeje =',N_sinBHeje)
        print()
        print('Iteration  time:',format(time.time()-startTimeLocal ,'.3f'),'sec')
        print('Simulation time:',format(time.time()-startTimeGlobal,'.2f'),\
            'sec (',format((time.time()-startTimeGlobal)/60,'.2f'),'min )')
        print('\n')
        
        # Bondi accretion onto BHs (limited by the Eddington limit) and intracluster gas evolution:
        # -----------------------------------------------------------------------------------------------------------------------
        
        # sound speed (virial temperature):
        c_sound = np.sqrt(0.4 * G_Newt * Mcl / rh)
        
        # compute gas density in the core of the cluster:
        rho_gas = 3 * M_gas / (4 * np.pi * (rh / 1.3)**3)
        
        # Eddington accretion limit:
        dMdt_Edd = 4 * np.pi * G_Newt * mBH * m_proton / epsilon_acc / sThomson / c_light
        
        # Bondi accretion rate:
        dMdt_Bondi = 4 * np.pi * rho_gas * (G_Newt * mBH)**2 / c_sound**3
        
        # Bondi rate should not exceed the Eddington limit:
        def accretion_rate(r_bon_i, r_edd_i):
            return np.piecewise(r_bon_i, [(r_bon_i<r_edd_i), (r_bon_i>=r_edd_i)], [r_bon_i, r_edd_i]) + 0
        accretion_rate_vector = np.vectorize(accretion_rate)
        
        # accretion rate:
        dMdt = accretion_rate_vector(dMdt_Bondi, dMdt_Edd)
        
        # accreted matter by all BHs:
        dM_gas_accreted = np.sum(dt1 * dMdt)
        
        # adjust BH masses due to gas accretion:
        mBH = mBH + dt1 * dMdt
        
        # gas expulsion timescale:
        if tM<0:
            t_expulsion = 7.1e-3 * Myr * (1 - epsilon_SF) * Mcl / epsilon_SF / 1e5 / Msun * pc / rh
        else:
            t_expulsion = tM
        
        # update half-mass radius of the gas:
        rh_gas = rh_gas0 * (1 + np.sqrt(t / t_expulsion))
        
        # update gas reservoir in cluster:
        M_gas1 = M_gas0 * np.exp(-t / t_expulsion) # expontial expulsion
        
        # amoung of gas ejected due to stellar feedback (UV, winds, SNe):
        dM_gas_removed = M_gas - M_gas1
        
        dM_gas = dM_gas_removed + dM_gas_accreted
        
        # subtract matter accreted:
        M_gas = M_gas - dM_gas
        
        # check if M_gas is negative value:
        if M_gas<0:
            M_gas = 0.0
        
        # update cluster mass:
        Mcl = Mcl - dM_gas
        
        M_gas = M_gas1
        
        # -----------------------------------------------------------------------------------------------------------------------
        
        # update number of iterations performed:
        Niter+=1

    # Primary->1, Secondary->2
    # ---------------------------------------------------------------------------------------------------------------------------

    # discard first lines of `junk` zeroes:
    mergers = np.delete(mergers,0,axis=0)
    evolut = np.delete(evolut,0,axis=0)

    for i in range(0,N_me):

        channel = mergers[i][0 ]
        sma     = mergers[i][1 ]
        ecc     = mergers[i][2 ]
        m1      = mergers[i][3 ]
        m2      = mergers[i][4 ]
        chi1    = mergers[i][5 ]
        chi2    = mergers[i][6 ]
        g1      = mergers[i][7 ]
        g2      = mergers[i][8 ]
        tForm   = mergers[i][9 ]
        tMerge  = mergers[i][10]
        zForm   = mergers[i][11]
        zMerge  = mergers[i][12]
        Nhar    = mergers[i][13]
        Nsub    = mergers[i][14]
        q       = mergers[i][15]
        chiEff  = mergers[i][16]
        theta1  = mergers[i][17]
        theta2  = mergers[i][18]
        dPhi    = mergers[i][19]
        mRem    = mergers[i][20]
        sRem    = mergers[i][21]
        gRem    = mergers[i][22]
        vGW     = mergers[i][23]
        jBH_1   = mergers[i][24]
        jBH_2   = mergers[i][25]

        if m2 > m1:
            mergers[i][3 ] = m2
            mergers[i][4 ] = m1
            mergers[i][5 ] = chi2
            mergers[i][6 ] = chi1
            mergers[i][7 ] = int(g2)
            mergers[i][8 ] = int(g1)
            mergers[i][24] = int(jBH_2)
            mergers[i][25] = int(jBH_1)

    # Exporting data
    # ---------------------------------------------------------------------------------------------------------------------------

    with open(str(mergersFile)+'.txt','w') as f_mergers:
        for i in range(0,N_me):
            f_mergers.write(str(mergers[i][0])+' '+str(mergers[i][1]/AU)+' '+\
                str(mergers[i][2])+' '+str(mergers[i][3]/Msun)+' '+str(mergers[i][4]/Msun)+' '+\
                    str(mergers[i][5])+' '+str(mergers[i][6])+' '+str(mergers[i][7])+' '+\
                        str(mergers[i][8])+' '+str(mergers[i][9]/Myr)+' '+\
                            str(mergers[i][10]/Myr)+' '+str(mergers[i][11])+' '+\
                                str(mergers[i][12])+' '+str(mergers[i][13])\
                                    +' '+str(mergers[i][14])+' '+str(mergers[i][15])\
                                        +' '+str(mergers[i][16])+' '+str(mergers[i][17])+' '+str(mergers[i][18])+' '+\
                  str(mergers[i][19])+' '+str(mergers[i][20]/Msun)+' '+str(mergers[i][21])+' '+\
                  str(mergers[i][22])+' '+str(mergers[i][23]/1e3)+' '+str(mergers[i][24])+' '+str(mergers[i][25])+' '+\
                            str(Mcl0/Msun)+' '+str(zClForm))
            f_mergers.write('\n')

    with open(str(evolutionFile)+'.txt','w') as f_evolution:
        for i in range(0,Niter):
            f_evolution.write( str(evolut[i][0 ])+' '+str(evolut[i][1 ])+' '\
                              +str(evolut[i][2 ])+' '+str(evolut[i][3 ])+' '\
                              +str(evolut[i][4 ])+' '+str(evolut[i][5 ])+' '\
                              +str(evolut[i][6 ])+' '+str(evolut[i][7 ])+' '\
                              +str(evolut[i][8 ])+' '+str(evolut[i][9 ])+' '\
                              +str(evolut[i][10])+' '+str(evolut[i][11])+' '\
                              +str(evolut[i][12])+' '+str(evolut[i][13])+' '\
                              +str(evolut[i][14])+' '+str(evolut[i][15])+' '\
                              +str(evolut[i][16])+' '+str(evolut[i][17])+' '\
                              +str(evolut[i][18])+' '+str(evolut[i][19])+' '\
                              +str(evolut[i][20])+' '+str(evolut[i][21])+' '\
                              +str(evolut[i][22])+' '+str(evolut[i][23])+' '\
                              +str(evolut[i][24])+' '+str(evolut[i][25])+' '\
                              +str(evolut[i][26])+' '+str(evolut[i][27])+' '\
                              +str(evolut[i][28])+' '+str(evolut[i][29])+' '\
                              +str(evolut[i][30])+' '+str(evolut[i][31])+' '\
                              +str(evolut[i][32])+' '+str(evolut[i][33])+' '\
                              +str(evolut[i][34]))
            f_evolution.write('\n')

    np.savez(blackholeFile, mBH_ini=mBH_1g/Msun, mBH_fin=mBH/Msun)
    
# End of source code
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
