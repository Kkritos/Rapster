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

def remnant_kick(m1, m2, chi1, chi2, theta1, theta2, dPhi):
    """
    Merger remnant kick velocity in km/s.
    Assumes isotropic spin directions.
    From Gerosa & Kesdsen (2016) and references therein.
    
    Inputs:
    @in m1: primary mass
    @in m2: secondary mass
    @in chi1: dimensionless spin parameter of primary
    @in chi2: dimensionless spin parameter of secondary
    @in theta1: angle between orbital ang. mom. and primary spin vector
    @in theta2: angle between orbital ang. mom. and secondary spin vector
    @in dPhi: angle between spin projections in orbital plane
    """
    
    q = m2 / m1 if m2 < m1 else m1 / m2
    
    eta = q / (1 + q)**2
    
    A, B, H, V11, VA, VB, VC, C2, C3, zeeta = \
    1.2e4, -0.93, 6.9e3, 3677.76, 2481.21, 1792.45, 1506.52, 1140, 2481, 145*np.pi/180
    
    cost1, cost2 = np.cos(theta1), np.cos(theta2)
    sint1, sint2 = np.sqrt(1 - cost1**2), np.sqrt(1 - cost2**2)
    
    chi1_para, chi2_para = chi1 * cost1, chi2 * cost2
    
    chi1_perp, chi2_perp = chi1 * sint1, chi2 * sint2
    
    vm = A * eta**2 * (1 - q) / (1 + q) * (1 + B * eta)
    
    Delta_para = (chi1_para - q * chi2_para               ) / (1 + q)
    Delta_perp = (chi1_perp - q * chi2_perp * np.cos(dPhi)) / (1 + q)
    
    chi_para = (chi1_para + q**2 * chi2_para) / (1 + q)**2
    chi_perp = (chi1_perp + q**2 * chi2_perp * np.cos(dPhi)) / (1 + q)**2
    
    vs_perp = H * eta**2 * Delta_para
    
    vs_para = 16 * eta**2 * (Delta_perp * (V11 + 2 * VA * chi_para + 4 * VB * chi_para**2 + 8 * VC * chi_para**3) \
                             + 2 * chi_perp * Delta_para * (C2 + 2 * C3 * chi_para)) \
    * np.cos(np.random.rand() * 2 * np.pi)
    
    vk = np.sqrt(vm**2 + 2 * vm * vs_perp * np.cos(zeeta) + vs_perp**2 + vs_para**2)
    
    return vk

def remnant_spin(m1, m2, chi1, chi2, theta1, theta2, dPhi):
    """
    Dimensionles merger remnant spin parameter in [0, 1].
    From Gerosa & Kesdsen (2016) and references therein.
    
    Inputs:
    @in m1: primary mass
    @in m2: secondary mass
    @in chi1: dimensionless spin parameter of primary
    @in chi2: dimensionless spin parameter of secondary
    @in theta1: angle between orbital ang. mom. and primary spin vector
    @in theta2: angle between orbital ang. mom. and secondary spin vector
    @in dPhi: angle between spin projections in orbital plane
    """
    
    q = m2 / m1 if m2 < m1 else m1 / m2
    
    eta = q / (1 + q)**2
    
    cost1, cost2 = np.cos(theta1), np.cos(theta2)
    sint1, sint2 = np.sqrt(1 - cost1**2), np.sqrt(1 - cost2**2)
    
    chi1_para, chi2_para = chi1 * cost1, chi2 * cost2
    
    chi1_perp, chi2_perp = chi1 * sint1, chi2 * sint2
    
    Delta_para = (chi1_para - q * chi2_para               ) / (1 + q)
    Delta_perp = (chi1_perp - q * chi2_perp * np.cos(dPhi)) / (1 + q)
    
    chi_para = (chi1_para + q**2 * chi2_para) / (1 + q)**2
    chi_perp = (chi1_perp + q**2 * chi2_perp * np.cos(dPhi)) / (1 + q)**2
    
    chi1_chi2 = chi1_para * chi2_para + chi1_perp * chi2_perp * np.cos(dPhi)
    
    chi_squared = (chi1**2 + q**4 * chi2**2 + 2 * q**2 * chi1_chi2) / (1 + q)**4
    
    t0, t2, t3, s4, s5 = -2.8904, -3.51712, 2.5763, -0.1229, 0.4537
    
    el = 2 * np.sqrt(3) + t2 * eta + t3 * eta**2 + s4 * (1 + q)**4 / (1 + q**2) * chi_squared \
    + (s5 * eta + t0 + 2) * (1 + q)**2 / (1 + q**2) * chi_para
    
    chi_f = np.sqrt(q**2 * el**2 / (1 + q)**4 + chi_squared + 2 * q * el * chi_para / (1 + q)**2)
    
    chi_f = chi_f if chi_f < 1 else 1
    
    return chi_f

def remnant_mass(m1, m2, chi1, chi2, theta1, theta2, dPhi):
    """
    Merger remnant mass.
    From Gerosa & Kesdsen (2016) and references therein.
    
    Inputs:
    @in m1: primary mass
    @in m2: secondary mass
    @in chi1: dimensionless spin parameter of primary
    @in chi2: dimensionless spin parameter of secondary
    @in theta1: angle between orbital ang. mom. and primary spin vector
    @in theta2: angle between orbital ang. mom. and secondary spin vector
    @in dPhi: angle between spin projections in orbital plane
    """
    
    q = m2 / m1 if m2 < m1 else m1 / m2
    
    eta = q / (1 + q)**2
    
    cost1, cost2 = np.cos(theta1), np.cos(theta2)
    sint1, sint2 = np.sqrt(1 - cost1**2), np.sqrt(1 - cost2**2)
    
    chi1_para, chi2_para = chi1 * cost1, chi2 * cost2
    
    chi_para = (chi1_para + q**2 * chi2_para) / (1 + q)**2
    
    Z1 = 1 + (1 - chi_para**2)**(1 / 3) * ((1 + chi_para)**(1 / 3) + (1 - chi_para)**(1 / 3))
    Z2 = np.sqrt(3 * chi_para**2 + Z1**2)
    
    r_isco = 3 + Z2 - np.sign(chi_para) * np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))
    
    E_isco = np.sqrt(1 - 2 / 3 / r_isco)
    
    p0, p1 = 0.04827, 0.01707
    
    m_f = (m1 + m2) * (1 - eta * (1 - 4 * eta) * (1 - E_isco) - 16 * eta**2 * (p0 + 4 * p1 * chi_para * (chi_para + 1)))
    
    return m_f

def merger_remnant(m1, m2, chi1, chi2, theta1, theta2, dPhi):
    """
    Final mass, final spin parameter and GW kick velocity of a merger remnant.
    Ref: D.Gerosa & M.Kesden, PRD 93 (2016), 124066.
    
    @in  m1: first  BH mass [<units>]
    @in  m2: second BH mass [<units>]
    @in  chi1: first  BH spin parameter in [0,1]
    @in  chi2: second BH spin parameter in [0,1]
    @in  theta1: angle of spin 1 with angular momentum [rad]
    @in  theta2: angle of spin 2 with angular momentum [rad]
    @in  dPhi: angle between the orbital plane projections of the spins [rad]
    
    @out mRem    : merger remnant BH mass [<units>]
    @out chiRem  : merger remnant BH spin parameter in [0,1]
    @out vGWkick : merger remnant GW kick velocity [km/s]
    """
    
    m_rem = remnant_mass(m1, m2, chi1, chi2, theta1, theta2, dPhi)
    chi_rem = remnant_spin(m1, m2, chi1, chi2, theta1, theta2, dPhi)
    vGW_kick = remnant_kick(m1, m2, chi1, chi2, theta1, theta2, dPhi)
    
    return m_rem, chi_rem, vGW_kick

# end of file.
