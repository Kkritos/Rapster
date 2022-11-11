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

# Imports
# ----------------------------------------------------------------------------------------------------------------------------

import numpy as np
import precession as pre
import scipy.integrate as integrate
import scipy.special as scs
import time
import math
import astropy.units as u
import argparse
from scipy.stats import poisson
from scipy.stats import maxwell
from cmath import isnan
from scipy import interpolate
from astropy.cosmology import FlatLambdaCDM
from scipy.optimize import fsolve
from subprocess import Popen, PIPE

# Global constants in S.I.
# ----------------------------------------------------------------------------------------------------------------------------

yr       = 365*24*60*60  # year
kyr      = 1e3*yr        # kilo year
Myr      = 1e6*yr        # mega year
Gyr      = 1e9*yr        # giga year
pc       = 3.01e16       # parsec
kpc      = 1e3*pc        # kilo parsec
Mpc      = 1e6*pc        # mega parsec
Gpc      = 1e9*pc        # giga parsec
Msun     = 1.99e30       # solar mass
Zsun     = 0.014         # solar metallicity
Rsun     = 6.957e8       # solar radius
AU       = 1.5e11        # astronomical unit
c_light  = 3e8           # speed of light
G_Newt   = 6.7e-11       # Newton constant
OmegaK0  = 0             # curvature density parameter
OmegaR0  = 5.38e-5       # radiation density parameter
OmegaV   = 0.685         # vacuum density parameter
OmegaM0  = 0.315         # matter density paramter
h        = 0.674         # Hubble reduced parameter
H0       = 100*h*1e3/Mpc # Hubble constant
zEq      = 3402          # radiation-matter redshift
cosmo    = FlatLambdaCDM(H0=100*h * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=OmegaM0) # define cosmology
T_Hubble = cosmo.age(0).value * Gyr # Hubble time

def E_Hubble(z):
    '''
    Auxiliary function for cosmology (Hogg, Distance measures in cosmology, 2000).

    @in z: redshift
    '''

    return np.sqrt( (1+z)**3*OmegaM0 + (1+z)**2*OmegaK0 + (1+z)**4*OmegaR0 + OmegaV )

def t_lbb(z):
    '''
    Lookback time given the redshift.
    
    @in z: redshift, scalar
    
    @out : look-back time in seconds
    '''

    out = integrate.quad(lambda zz: 1/H0/(1+zz)/E_Hubble(zz),0,z)[0]
    
    return out

def redd(look):
    '''
    Redshift given the lookback time.
    
    @in look: cosmic time in seconds, scalar
    
    @out    : redshift at time `look`
    '''

    s=0
    dz=0.0001
    z=0
    
    while(s<look):
        s=s+dz/H0/(1+z)/E_Hubble(z)
        z=z+dz
    
    return z

def D_lumin(z):
    '''
    Luminosity distance.

    @in z: redshift
    
    @out : luminosity distance in meters
    '''

    ret = np.zeros(z.size)
    
    for i in range(0,z.size):
        ret[i] = integrate.quad(lambda zz: 1/E_Hubble(zz),0,z[i])[0]
    
    return ret*(1+z)*c_light/H0

def veloNFW(r,Rscale,rhoScale):
    '''
    Navaro-Frenk-White circular velocity profile.
    
    @in r       : observation radius
    @in Rscale  : effective radius scale
    @in rhoScale: central density scale
    
    @out        : circular velocity in m/s
    '''

    Mr = 4*np.pi*rhoScale*Rscale**3*(np.log(1+r/Rscale)-1/(1+Rscale/r)) # mass within radius r
    
    return np.sqrt(2*G_Newt*Mr/r)

def IMF_kroupa(m):
    '''
    Kroupa (2002) initial mass function.
    
    @in m: stellar mass in solar masses, array input
    
    @out : number dN of stars in mass bin [m,m+dm] in units of 1/Msun
    '''

    # mass boundaries (in solar masses):
    m1=0.08
    m2=0.50
    m3=1.00
    
    # spectral indices (broken power law; central values):
    a0=-0.3
    a1=-1.3
    a2=-2.3
    a3=-2.3
    
    # normalization constants:
    c1=m1**a0/m1**a1
    c2=c1*m2**a1/m2**a2
    c3=c2*m3**a2/m3**a3
    
    out=np.zeros(m.size)
    
    for i in range(0,m.size):
        
        if  (m[i]<=m1):
        
            out[i]=m[i]**a0
        
        elif(m[i]<=m2 and m[i]>m1):
            
            out[i]=c1*m[i]**a1
       
        elif(m[i]<=m3 and m[i]>m2):
         
            out[i]=c2*m[i]**a2
        
        elif(m[i]>=m3):
            
            out[i]=c3*m[i]**a3
    
    return out

def T_coal(m1,m2,a0,e0):
    '''
    GW coalescence timescale (I. Mandel 2021 fit to Peters timescale),
    including 1st order pN correction effects [Zwick et al., MNRAS 495, 2321 (2020)]
    
    @in m1: primary   mass component in kg
    @in m2: secondary mass component in kg
    @in a0: initial semimajor axis in meters
    @in e0: initial eccentricity in range [0,1)
    
    @out  : coalescence timescale in seconds    
    '''

    # coalescence timescale for circular orbit:
    Tc = 5*c_light**5*a0**4/(256*G_Newt**3*m1*m2*(m1+m2))
    
    # 1st order pN correction:
    S = 8**(1-np.sqrt(1-e0)) * np.exp( 5*G_Newt*(m1+m2)/c_light**2/a0/(1-e0) )
    
    return Tc*(1+0.27*e0**10+0.33*e0**20+0.2*e0**1000)*(1-e0**2)**(7/2) * S

def fGW(a,e,M):
    '''
    Get GW frequency, given Keplerian parameters [Wen (2003)]

    @in a: semimajor axis of binary in meters
    @in e: binary eccentricity
    @in M: binary's total mass in kg

    @out : GW frequency in Hz
    '''

    return (1+e)**(1.1954)/(1-e**2)**(3/2)*np.sqrt(G_Newt*M)/np.pi/a**(3/2)

def aEj(m1,m2,m3,Hrate,v_esc):
    '''
    Binary-single ejection semimajor axis

    @in m1   : primary   binary component in kg
    @in m2   : secondary binary component in kg
    @in m3   : third single body in kg
    @in Hrate: binary-single hardening rate
    @in v_esc: escape velocity from cluster environment in m/s

    @out     : critical semimajor axis for ejection in meters
    '''

    # total binary mass:
    m12  = m1+m2
    
    # total mass of the binary-single system:
    m123 = m12+m3
    
    # reduced mass of the binary-single system:
    mu   = m12*m3/m123
    
    return Hrate/2/np.pi*m1*m2*m3/m12**3*G_Newt*mu/v_esc**2

def veloDisp(m,xi,mAvg,Mcl,rh):
    '''
    Velocity dispersion of mass component m.
    
    @in m   : massive component mass in kg
    @in xi  : temperature ratio
    @in mAvg: average mass ik kg
    @in Mcl : cluster mass in kg
    @in rh  : cluster half-mass radius in m
    
    @out    : velocity dispersion in m/s
    '''

    # rms velocity of stars:
    vStar = np.sqrt(0.4*G_Newt*Mcl/rh)
    
    return np.sqrt(mAvg/m*xi)*vStar

def veloDispRel(m1,m2,xi,mAvg,Mcl,rh):
    '''
    Relative velocity dispersion between mass components m1 and m2.

    @in m1  : first  mass component in kg
    @in m2  : second mass component in kg
    @in xi  : equipartition parameter
    @in mAvg: average mass ik kg
    @in Mcl : cluster mass in kg
    @in rh  : cluster half-mass radius in m
    
    @out    : relative velocity dispersion in m/s
    '''

    return np.sqrt(veloDisp(m1,xi,mAvg,Mcl,rh)**2+veloDisp(m2,xi,mAvg,Mcl,rh)**2)

def Rseg(m,MBH,xi,mAvg,Mcl,rh):
    '''
    Segregation radius.

    @in m   : massive component mass in kg
    @in MBH : total mass in black holes
    @in xi  : temperature ratio
    @in mAvg: average mass ik kg
    @in Mcl : cluster mass in kg
    @in rh  : cluster half-mass radius in meters
    
    @out    : segregation radius in meters
    '''

    return 1/xi*MBH/Mcl*m/mAvg*rh

def Vseg(m,MBH,xi,mAvg,Mcl,rh):
    '''
    Segregation volume.

    @in m   : massive component mass in kg
    @in MBH : total mass in black holes
    @in xi  : temperature ratio
    @in mAvg: average mass ik kg
    @in Mcl : cluster mass in kg
    @in rh  : cluster half-mass radius in m
    
    @out    : segregation volume in m^3
    '''

    return 4*np.pi/3*Rseg(m,MBH,xi,mAvg,Mcl,rh)**3

def mergerRemnant(m1,m2,chi1,chi2,theta1,theta2,dPhi):
    '''
    Final mass, final spin parameter and GW kick velocity of a merger remnant
    calculated with the PRECESSION package.
    Ref: D.Gerosa & M.Kesden, PRD 93 (2016), 124066.

    Parameter description:

    @in  m1      : first  BH mass in <units>
    @in  m2      : second BH mass in <units>
    @in  chi1    : first  BH spin parameter in [0,1]
    @in  chi2    : second BH spin parameter in [0,1]
    @in  theta1  : angle of spin 1 with angular momentum in radians
    @in  theta2  : angle of spin 2 with angular momentum in radians
    @in  dPhi    : angle between the orbital plane projections of the spins in radians
    
    @out mRem    : merger remnant BH mass in <units>
    @out chiRem  : merger remnant BH spin parameter in [0,1]
    @out vGWkick : merger remnant GW kick velocity in km/s
    '''

    # BBH mass ratio:
    q = np.min([m1,m2])/np.max([m1,m2])
    
    chi_p   = chi1
    chi_s   = chi2
    
    theta_p = theta1
    theta_s = theta2
    
    # order `1`-> primary, `2`-> secondary:
    if m2>m1:
        
        chi_p   = chi2
        chi_s   = chi1
        
        theta_p = theta2
        theta_s = theta1
    
    M,m_1,m_2,S_1,S_2 = pre.get_fixed(q,chi_p,chi_s) # units M=1
    
    # Final mass of merger remnant in solar masses:
    mRem = pre.finalmass(theta_p,theta_s,dPhi,q,S_1,S_2)*(m1+m2)

    # Final spin of merger remnant:
    chiRem = pre.finalspin(theta_p,theta_s,dPhi,q,S_1,S_2)

    # Final GW kick (converted to m/s):
    vGWkick = pre.finalkick(theta_p,theta_s,dPhi,q,S_1,S_2,maxkick=False,kms=True,more=False)*1e3

    return mRem,chiRem,vGWkick

# Binary-single hardening Rate: H = A*(1+a/a0)**gamma; best fit params:

q_sma     = np.array([ 1,1/3,1/9,1/27,1/81,1/243 ]) # mass ratio values

A_sma     = np.array([ 14.55,15.82,17.17,18.15,18.81,19.16 ])

a0_sma    = np.array([ 3.48,4.18,3.59,3.32,3.87,4.17 ])

gamma_sma = np.array([ -0.25,-0.90,-0.79,-0.77,-0.82,-0.86 ])

A_sma_interpol     = interpolate.interp1d(q_sma,A_sma    ,kind='linear')
a0_sma_interpol    = interpolate.interp1d(q_sma,a0_sma   ,kind='linear')
gamma_sma_interpol = interpolate.interp1d(q_sma,gamma_sma,kind='linear')

def H_rate(a,q):
    '''
    Hardening rate.
    
    @in a: binary's semimajor axis nomralized to its hardening sma
    @in q: binary's mass ratio

    @out : hardening rate
    '''

    try:
        Hrate = A_sma_interpol(q)*(1+a/a0_sma_interpol(q))**gamma_sma_interpol(q)
    except:
        Hrate = 15.0
      
    return Hrate

# Binary-single eccentricity growth rate: K = A*(1+a/a0)**gamma + B; best fit params:

# mass ratio values:
q_ecc     = np.array([ 1,1/3,1/9,1/27 ])

# eccentricity values:
e_ecc     = np.array([ 0, 0.15,0.3,0.45,0.6,0.75,0.9, 1 ])

A_ecc     = np.array([ [0,0,0,0] , [0.037,0.082,0.051,0.064] , [0.075,0.095,0.111,0.143] , [0.105,0.129,0.172,0.212] , \
    [0.121,0.166,0.181,0.216] , [0.134,0.159,0.179,0.173] , [0.082,0.095,0.117,0.129] , [0,0,0,0] ])

a0_ecc    = np.array([ [0.4,0.05,0.4,0.3] , [0.339,0.042,0.385,0.284] , [0.151,0.213,0.307,1.033] , [0.088,0.137,0.526,0.722] , \
    [0.090,0.081,0.251,0.430] , [0.064,0.079,0.195,0.771] , [0.085,0.122,0.400,0.329] , [0.09,0.2,0.5,0.2] ])

gamma_ecc = np.array([ [-4,-0.2,-0.9,-2] , [-3.335,-0.168,-0.891,-1.206] , [-1.548,-1.152,-1.167,-1.537] , \
    [-0.893,-0.655,-1.174,-1.257] , [-0.895,-0.546,-1.169,-1.163] , [-0.544,-0.497,-0.846,-1.934] , \
        [-0.663,-0.716,-1.170,-1.125] , [-0.7,-0.8,-1.2,-1.2] ])

B_ecc     = np.array([ [0,0,0,0] , [-0.012,-0.048,-0.011,0.021] , [-0.008,-0.012,-0.007,-0.021] , [-0.005,-0.006,-0.016,-0.022] ,\
    [-0.008,-0.006,-0.007,-0.014] , [-0.006,-0.010,-0.004,-0.014] , [-0.004,-0.008,-0.001,-0.020] , [0,0,0,0] ])

A_ecc_interpol     = interpolate.interp2d(q_ecc,e_ecc,A_ecc    ,kind='linear')
a0_ecc_interpol    = interpolate.interp2d(q_ecc,e_ecc,a0_ecc   ,kind='linear')
gamma_ecc_interpol = interpolate.interp2d(q_ecc,e_ecc,gamma_ecc,kind='linear')
B_ecc_interpol     = interpolate.interp2d(q_ecc,e_ecc,B_ecc    ,kind='linear')

def K_rate(a,e,q):
    '''
    Eccentricity growth rate.

    @in a: binary's semimajor axis normilized to its hardening sma
    @in e: binary's eccentricity
    @in q: binary's mass ratio
    
    @out : eccentricity growth rate
    '''
  
    try:
        Krate = A_ecc_interpol(q,e)*(1+a/a0_ecc_interpol(q,e))**gamma_ecc_interpol(q,e)+B_ecc_interpol(q,e)
    except:
        Krate = 0.0
    
    return Krate[0]

# 2-dim interpolation of the ZAMS Mass Remnant Mass from SEVN package (for faster performance)
# [Spera, M., Mapelli, M., Bressan, A., 2015, MNRAS, 451, 2086]
# Delayed core-collapse engine is assumed
# ------------------------------------------------------------------------------------------------------------------------------

# Reading files exported from SEVN code and stored according to metallicity for various ZAMS masses:
MzamsMrem1  = np.load('./MzamsMrem/MzamsMrem1.npz' ); Mrem_delayed_1  = MzamsMrem1 ['Mrem1' ]
MzamsMrem2  = np.load('./MzamsMrem/MzamsMrem2.npz' ); Mrem_delayed_2  = MzamsMrem2 ['Mrem2' ]
MzamsMrem3  = np.load('./MzamsMrem/MzamsMrem3.npz' ); Mrem_delayed_3  = MzamsMrem3 ['Mrem3' ]
MzamsMrem4  = np.load('./MzamsMrem/MzamsMrem4.npz' ); Mrem_delayed_4  = MzamsMrem4 ['Mrem4' ]
MzamsMrem5  = np.load('./MzamsMrem/MzamsMrem5.npz' ); Mrem_delayed_5  = MzamsMrem5 ['Mrem5' ]
MzamsMrem6  = np.load('./MzamsMrem/MzamsMrem6.npz' ); Mrem_delayed_6  = MzamsMrem6 ['Mrem6' ]
MzamsMrem7  = np.load('./MzamsMrem/MzamsMrem7.npz' ); Mrem_delayed_7  = MzamsMrem7 ['Mrem7' ]
MzamsMrem8  = np.load('./MzamsMrem/MzamsMrem8.npz' ); Mrem_delayed_8  = MzamsMrem8 ['Mrem8' ]
MzamsMrem9  = np.load('./MzamsMrem/MzamsMrem9.npz' ); Mrem_delayed_9  = MzamsMrem9 ['Mrem9' ]
MzamsMrem10 = np.load('./MzamsMrem/MzamsMrem10.npz'); Mrem_delayed_10 = MzamsMrem10['Mrem10']
MzamsMrem11 = np.load('./MzamsMrem/MzamsMrem11.npz'); Mrem_delayed_11 = MzamsMrem11['Mrem11']
MzamsMrem12 = np.load('./MzamsMrem/MzamsMrem12.npz'); Mrem_delayed_12 = MzamsMrem12['Mrem12']

# collect remnant masses with various metallicity values in a single array:
Mrem_delayed = np.array([Mrem_delayed_1,Mrem_delayed_2,Mrem_delayed_3,Mrem_delayed_4,Mrem_delayed_5,Mrem_delayed_6,Mrem_delayed_7,\
    Mrem_delayed_8,Mrem_delayed_9,Mrem_delayed_10,Mrem_delayed_11,Mrem_delayed_12])

# Metallicity should not be out of this range: [1e-4,1.7e-2]:
Zvalues = np.array([1.0e-4,2.0e-4,5.0e-4,1.0e-3,2.0e-3,4.0e-3,6.0e-3,8.0e-3,1.0e-2,1.4e-2,1.7e-2,2.0e-2])

# Mass should not be out of this range: [20,340] solar masses

Npoints = 100
Mzams = np.linspace(20,340,Npoints)

# interpolate:
MremInterpol = interpolate.interp2d(Mzams,Zvalues,Mrem_delayed,kind='linear',bounds_error=True)

def Mrem(M,Z):
    '''
    Remnant mass as a function of progenitor metallicity and ZAMS mass.

    @in M: ZAMS mass in solar masses ; in range [20,340]
    @in Z: absolute metallicity      ; in range [1e-4,1.7e-2]

    @out : remnant mass in solar masses (scalar or array depending on M)
    '''

    M_lowerEdge = 60  # absolute lower edge of the upper mass gap (in solar masses)
    M_upperEdge = 120 # absolute upper edge of the upper mass gap (in solar masses)
    
    # check if mass input is an array or not:
    if isinstance(M,np.ndarray): # M is array
        
        out = MremInterpol(M,Z)*(np.heaviside(M_lowerEdge*np.ones(M.size)-MremInterpol(M,Z),0)\
            + np.heaviside(MremInterpol(M,Z)-M_upperEdge*np.ones(M.size),0))
    else: # M is not array
        
        out = MremInterpol(M,Z)*(np.heaviside(M_lowerEdge-MremInterpol(M,Z),0)\
            + np.heaviside(MremInterpol(M,Z)-M_upperEdge,0))
    
    return out

#def tRelax(N,rh,m):
def tRelax(Mcl, N, rh, m):
    '''
    Half-mass relaxation timescale.

    @in Mcl: cluster total mass (stars+gas)
    @in N  : Number of objects
    @in rh : half-mass radius
    @in m  : average mass

    @out  : half-mass relaxation timescale in seconds
    '''
    
    #return 0.138*N**(1/2)*rh**(3/2)/m**(1/2)/G_Newt**(1/2)/np.log(0.4*N)
    #return 128 * Myr * (Mcl/1e5/Msun)**(1/2) * (rh/pc)**(3/2) * 0.6*Msun/m * 8/ np.log(0.02*N)
    return (rh**3/G_Newt/Mcl)**(1/2)*N/8/np.log(0.1*N)

def RateCapture(m1,m2,n1,n2,xi,mAvg,Mcl,rh):
    '''
    Dynamical two-body capture volumetric rate density.

    @in m1  : first  mass in kg
    @in m2  : second mass in kg
    @in n1  : local volumetric number density of m1 objects
    @in n2  : local volumetric number density of m2 objects
    @in xi  : equipartition parameter
    @in mAvg: average mass in cluster in kg, scalar
    @in Mcl : cluster mass in kg, scalar
    @in rh  : half-mass radius of cluster in meters, scalar

    @out    : capture rate density in SI
    '''

    # relative velocity dispersion in m/s:
    vRel = veloDispRel(m1,m2,xi,mAvg,Mcl,rh)
    
    # capture cross section:
    crossSectionVrel = 3**(11/14)*2**(3/14)*math.gamma(5/7)/np.sqrt(np.pi)\
        *2*np.pi*(85*np.pi/6/np.sqrt(2))**(2/7)*G_Newt**2/c_light**(10/7)\
            *(m1+m2)**(10/7)*(m1*m2)**(2/7)/vRel**(11/7)

    return n1*n2*crossSectionVrel

def RateInter(m1,m2,n1,n2,rp,xi,mAvg,Mcl,rh):
    '''
    Dynamical two-body close encounter volumetric rate density.

    @in m1  : first  mass in kg
    @in m2  : second mass in kg
    @in n1  : local volumetric number density of m1 objects
    @in n2  : local volumetric number density of m2 objects
    @in rp  : pericenter of interaction
    @in xi  : equipartition parameter
    @in mAvg: average mass in cluster in kg, scalar
    @in Mcl : cluster mass in kg, scalar
    @in rh  : half-mass radius of cluster in meters, scalar

    @out    : interaction rate density in SI
    '''

    # relative velocity dispersion in m/s:
    vRel = veloDispRel(m1,m2,xi,mAvg,Mcl,rh)
    
    # two-body interaction cross section = geomtrical + focusing terms:
    crossSection = np.pi*rp**2*(1+2*G_Newt*(m1+m2)/rp/vRel**2)
    
    return n1*n2*vRel*crossSection

def RateExchange(m1,m2,m3,n12,n3,a12,xi,mAvg,Mcl,rh):
    '''
    Dynamical binary-single exchange volumetric rate density.
    Exchange scheme: 1-2 -> 3-2, object 1 is substituted.

    @in m1  : to-be-swapped binary member mass in kg
    @in m2  : to-be-retained in the binary mass in kg
    @in m3  : substitutor third mass in kg
    @in n12 : number density of binaries 1-2 in SI
    @in n3  : number densiry of third bodies 3 in SI
    @in a12 : semimajor axis of binary 1-2
    @in xi  : equipartition parameter
    @in mAvg: average mass in cluster in kg, scalar
    @in Mcl : cluster mass in kg, scalar
    @in rh  : half-mass radius of cluster in meters, scalar
    
    @out      : exchange rate density in SI
    '''

    m12 = m1+m2 # old binary mass
    m23 = m2+m3 # new binary mass
    m13 = m1+m3
    
    m123 = m1+m2+m3 # total mass
    
    xx = m1/m12
    yy = m3/m123
    
    # relative velocity dispersion is m/s:
    vRel = veloDispRel(m1+m2,m3,xi,mAvg,Mcl,rh)
    
    # cross section [D.G.Heggie, P.Hut, S.L.W.McMillan, Astrophys.J. 467 (1996), 359-369.]:
    crossSection = 1.39*AU**2*(a12/0.1/AU)*(10e3/vRel)**2*(m123/Msun)*(m23/m123)**(1/6)\
            *(m3/m13)**(7/2)*(m123/m12)**(1/3)*(m13/m123)*np.exp(3.7+7.49*xx-1.89*yy-15.49*xx**2-2.93*xx*yy\
                -2.92*yy**2+3.07*xx**3+13.15*xx**2*yy-5.23*xx*yy**2+3.12*yy**3)
    
    return n12*n3*vRel*crossSection

def Rate3body(m,n,v,etaMin):
    '''
    Three-body binary formation volumetric rate density.
    Ref. [M.Morscher et al, Astrophys. J. 800, 9 (2015)]

    @in m     : mass scale of  in kg
    @in n     : number density of masses in 1/m^3
    @in v     : velocity dispersion of masses in m/s
    @in etaMin: minimum hardness ratio

    @out      : three-body rate density in SI
    '''

    return 3**(9/2)*np.pi**(13/2)/2**(25/2)*etaMin**(-11/2)*(1+2*etaMin)*(1+3*etaMin)*n**3*(G_Newt*m)**5/v**9

def hardness_sampler(r, m1, m2, m3, mmean, x, etamin=5):
    
    # function to sample the hardness of a binary formed by interaction with a third object
    # the sampling is based on inverting CDF (which reduces to solving a poly equation in this case)
    
    # INPUT: r -- a sample from uniform distribution on [0,1]
    # PARS: m1,m2 -- masses of a forming binary; m3 -- mass of an intermediary;\
    #       mmean -- mean mass; x -- equipartition parameter; etamin -- minimum hardness
    # OUTPUT: eta0[0] -- hardness of the formed binary
    
    rmass = m1*m2/(m1+m2)
    M1,M2,M3,RMass = m1/mmean,m2/mmean,m3/mmean,rmass/mmean
    
    A = (1/M2 + 1/M2)/(M1**(-2*x) + M2**(-2*x))
    B = M3**(-x)*np.sqrt(2*RMass)
    
    Anorm = etamin**-5.5 * (1 + 2*etamin*A) * \
                (1 + B/np.sqrt(etamin))
    
    coefs = [Anorm*(-1 + r), 0, 0, 0, 0, 0, 0, 0, 0, 2*A, 2*A*B, 1, B]
    
    roots = np.roots(coefs)
    eta0 = np.real(roots[np.logical_and(np.imag(roots)==0.,np.real(roots)>0)])**2
    assert len(eta0) == 1, "hardness sampler: Found more than one sample (>1 real positive root to CDF == r)"
    
    return eta0[0]

def vEscape(Mcl,rh):
    '''
    Escape velocity.

    @in Mcl: total cluster mass in kg
    @in rh : half-mass relaxation timescale

    @out   : escape velocity in m/s
    '''

    return 2*np.sqrt(0.4*G_Newt*Mcl/rh)

def f_fb(Mzams):
    '''
    Fraction of ejected supernova mass that falls back onto the newly-borned proto-comapct object.
    
    @in Mzams: ZAMS star mass in solar masses
    
    @out     : fall-back fraction
    '''
    
    # Chandrasekhar limit in solar masses:
    Mch = 1.4
    
    b_36 = 4.36e-4
    b_37 = 5.22
    b_38 = 6.84e-2
    
    # core mass at the Base of the Asymptotic Giant Branch:
    McBAGB = (b_36*Mzams**b_37 + b_38)**(1/4)
    
    # Carbon/Oxygen core mass:
    M_CO = np.max([Mch,0.773*McBAGB-0.35])
    
    # Determine proto-compact object mass:
    if   M_CO<=3.5:
        
        M_proto = 1.2
    
    elif M_CO>=3.5 and M_CO<6.0:
        
        M_proto = 1.3
    
    elif M_CO>=6.0 and M_CO<11.:
        
        M_proto = 1.4
    
    else:
        
        M_proto = 1.6
    
    # Determine fall-back fraction:
    a2 = 0.133 - 0.093/(Mzams - M_proto)
    b2 = -11*a2 + 1
    
    if   M_CO<2.5:
        
        Mfb = 0.2
        ffb = Mfb/(Mzams-M_proto)
    
    elif M_CO>=2.5 and M_CO<3.5:
        
        Mfb = 0.5*M_CO - 1.05
        ffb = Mfb/(Mzams-M_proto)
    
    elif M_CO>=3.5 and M_CO<11:
        
        ffb = a2*M_CO + b2
    
    else:
        
        ffb = 1
    
    return ffb

# From Hurley et al. (2000):
   
alpha_a1 = 1.593890e3
beta_a1  = 2.053038e3
gamma_a1 = 1.231226e3
eta_a1   = 2.327785e2
mu_a1    = 0.0e0

alpha_a2 = 2.706708e3
beta_a2  = 1.483131e3
gamma_a2 = 5.772723e2
eta_a2   = 7.411230e1
mu_a2    = 0.0e0

alpha_a3 = 1.466143e2
beta_a3  = -1.048442e2
gamma_a3 = -6.795374e1
eta_a3   = -1.391127e1
mu_a3    = 0.0e0

alpha_a4 = 4.141960e-2
beta_a4  = 4.564888e-2
gamma_a4 = 2.958542e-2
eta_a4   = 5.571483e-3
mu_a4    = 0.0e0

alpha_a5 = 3.426349e-1
beta_a5  = 0.0e0
gamma_a5 = 0.0e0
eta_a5   = 0.0e0
mu_a5    = 0.0e0

alpha_a6 = 1.949814e1
beta_a6  = 1.758178e0
gamma_a6 = -6.008212e0
eta_a6   = -4.470533e0
mu_a6    = 0.0e0

alpha_a7 = 4.903830e0
beta_a7  = 0.0e0
gamma_a7 = 0.0e0
eta_a7   = 0.0e0
mu_a7    = 0.0e0

alpha_a8 = 5.212154e-2
beta_a8  = 3.166411e-2
gamma_a8 = -2.750074e-3
eta_a8   = -2.271549e-3
mu_a8    = 0.0e0

alpha_a9 = 1.312179e0
beta_a9  = -3.294936e-1
gamma_a9 = 9.231860e-2
eta_a9   = 2.610989e-2
mu_a9    = 0.0e0

alpha_a10 = 8.073972e-1
beta_a10  = 0.0e0
gamma_a10 = 0.0e0
eta_a10   = 0.0e0
mu_a10    = 0.0e0

alpha_b36_ = 1.445216e-1
beta_b36_  = -6.180219e-2
gamma_b36_ = 3.093878e-2
eta_b36_   = 1.567090e-2
mu_b36_    = 0.0e0

alpha_b37_ = 1.304129e0
beta_b37_  = 1.395919e-1
gamma_b37_ = 4.142455e-3
eta_b37_   = -9.732503e-3
mu_b37_    = 0.0e0

alpha_b38_ = 5.114149e-1
beta_b38_  = -1.160850e-2
gamma_b38_ = 0.0e0
eta_b38_   = 0.0e0
mu_b38_    = 0.0e0   
   
def tBGB(M, Z):
    '''
    @in M: ZAMS mass in kg
    @in Z: metallicity
    '''
  
    M = M/Msun
    
    zeta = np.log10(Z/Zsun)
    a1 = alpha_a1 + beta_a1*zeta + gamma_a1*zeta**2 + eta_a1*zeta**3 + mu_a1*zeta**4
    a2 = alpha_a2 + beta_a2*zeta + gamma_a2*zeta**2 + eta_a2*zeta**3 + mu_a2*zeta**4
    a3 = alpha_a3 + beta_a3*zeta + gamma_a3*zeta**2 + eta_a3*zeta**3 + mu_a3*zeta**4
    a4 = alpha_a4 + beta_a4*zeta + gamma_a4*zeta**2 + eta_a4*zeta**3 + mu_a4*zeta**4
    a5 = alpha_a5 + beta_a5*zeta + gamma_a5*zeta**2 + eta_a5*zeta**3 + mu_a5*zeta**4
    
    return (a1 + a2*M**4 + a3*M**(5.5) + M**7) / (a4*M**2 + a5*M**7) * Myr   
   
def tMS(M, Z):
    '''
    Main-sequence lifetime in sec.
    
    @in M: ZAMS mass in kg
    @in Z: metallicity
    '''
    
    M_ = M
    M = M/Msun
    
    zeta = np.log10(Z/Zsun)
    a6  = alpha_a6  + beta_a6 *zeta + gamma_a6 *zeta**2 + eta_a6 *zeta**3 + mu_a6 *zeta**4
    a7  = alpha_a7  + beta_a7 *zeta + gamma_a7 *zeta**2 + eta_a7 *zeta**3 + mu_a7 *zeta**4
    a8  = alpha_a8  + beta_a8 *zeta + gamma_a8 *zeta**2 + eta_a8 *zeta**3 + mu_a8 *zeta**4
    a9  = alpha_a9  + beta_a9 *zeta + gamma_a9 *zeta**2 + eta_a9 *zeta**3 + mu_a9 *zeta**4
    a10 = alpha_a10 + beta_a10*zeta + gamma_a10*zeta**2 + eta_a10*zeta**3 + mu_a10*zeta**4
    
    mu = np.max([0.5, 1.0 - 0.01*np.max([a6/M**a7, a8+a9/M**a10])])
    x = np.max([0.95, np.min([0.95 - 0.03*(zeta + 0.30103), 0.99])])
    tHook = mu*tBGB(M_, Z)/Myr
    
    return np.max([tHook, x*tBGB(M_, Z)/Myr]) * Myr   
   
N_grid = 500
M_grid = np.linspace(20, 340, N_grid)
Z_grid = np.logspace(np.log10(1e-4), np.log10(2e-2), N_grid)

Mremnants_F12d = np.loadtxt('MzamsMrem_F12d.txt', unpack=True)

Mremnants_F12d = np.transpose(Mremnants_F12d)

MremInterpol_F12d = interpolate.interp2d(M_grid, Z_grid, Mremnants_F12d, kind='linear', bounds_error=True)

def Mrem_F12d(M, Z):
    return MremInterpol_F12d(M, Z)
   
def Mrem_F12d_simul(M, Z):

    flag_command = "echo "+str(M)+" "+str(Z)+" 4 | ./sse_new.exe | tail -2 | head -1 | awk '{print $1 $2}'"
    Mrem_command = "echo "+str(M)+" "+str(Z)+" 4 | ./sse_new.exe | tail -2 | head -1 | awk '{print $NF}'"
    
    flag = str((Popen(flag_command, shell=True, stdout=PIPE).stdout).read())
    
    Mrem = float((Popen(Mrem_command, shell=True, stdout=PIPE).stdout).read())
    
    if flag=="b'BlackHole\\n'":
        return Mrem+0
    else:
        return 0.0
  
# end of file
