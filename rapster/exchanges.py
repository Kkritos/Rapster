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
from .tidal_disruptions import *
from .fastrandom import choice_weighted, first_index_of
from .bhpool import FAM


def _choose(values, p):
    """np.random.choice(values, p=p) and the first index holding the drawn value
    (as np.where(values == value)[0][0]); same draw, same random numbers."""
    value, k = choice_weighted(values, p)
    if k is None:
        k = int(np.atleast_1d(np.squeeze(np.where(values==value))+0)[0])
    else:
        k = first_index_of(values, k)
    return value, k

def StarStar_to_BHstar(seed, t, z, k_ex1, N_ex1, m_avg, BH, ab, pairs, N_BHstar, N_tdeBHstar, v_star, vBH, tdes, binaries, m_min, m_max, with_tdes, f_accreted, EoS):
    """
    Creates BH-star binaries from star-star pairs, or a TDE occurs during the binary-single interaction.

    @in seed: simulation seed
    @in t: current simulation time
    @in z: current redshift
    @in k_ex1: current number of star-star -> BH-star exchanges
    @in N_ex1: total number of star-star -> BH-star exchanges
    @in m_avg: average mass
    @in BH: pool of single BHs (bhpool.BHPool: masses, spins, generations, TDE counts)
    @in ab: array of star-star semimajor axes
    @in pairs: array of BH-star pairs
    @in N_BHstar: number of BH-star pairs
    @in N_tdeBHstar: number of BH-star TDEs
    @in v_star: stellar velocity dispersion (3D)
    @in vBH: BH velocity dispersion (3D)
    @in tdes: tdes array 
    @in binaries: binaries array 
    @in m_min: min ZAMS star mass 
    @in m_max: max ZAMS star mass
    @in with_tdes: if =1 allow TDEs, else do not allow them
    @in f_accreted: fraction of the disrupted star accreted by the compact object
    @in EoS: NS equation of state; a name registered in compact_accretion.EOS_TABLES

    @out: all inputs
    """
 
    if k_ex1>0: # perform star-star -> BH-star exchange(s)
        
        for i in range(k_ex1):

            if ab.size == 0:  # no star-star binaries left to exchange, exit loop
                break

            if BH.size == 0:  # no single BHs left, exit loop
                break

            # sample the BH that substitutes one of the stars (probability ~ its mass):
            k = BH.draw(FAM['mass'])
            m, s, g, h = BH.get(k)
                
            a, kss = _choose(ab, ab/np.sum(ab))

            # delete star-star:
            ab = delete_at(ab, kss)

            # probability for TDE:
            p_TDE = min(2*R_sun*(m/m_avg)**(1/3)/a, 1.0) if with_tdes==1 else 0.0

            if np.random.rand() < p_TDE:
                # then TDE occurs and binary disrupts:
                
                tde_type = 2 # TDE during star-star+BH encounter

                # Sample star mass from evolving mass function:
                m_star, R_star = get_star(t, tBH_form, m_min, m_max)

                N_tdeBHstar, m, s, g, h, tdes = single_BH_TidalDisruption(seed, t, z, N_tdeBHstar, tde_type, m_avg, m_star, R_star, m, s, g, h, v_star, vBH, tdes, f_accreted, EoS)

                # update the BH, which stays in the single population:
                BH.replace(k, m, s, g, h)

            else:
                # BH substitutes a star and no TDE occurs:

                # update semimajor axis from energy conservation:
                a = a * m / m_avg
            
                # append pair:
                pairs = np.append(pairs, [[a, m, s, g, h]], axis=0)
            
                N_ex1+=1 # update number of star-star -> BH-star exchanges
                N_BHstar+=1
                
                BH.remove(k)

    return seed, t, z, k_ex1, N_ex1, m_avg, BH, ab, pairs, N_BHstar, N_tdeBHstar, v_star, vBH, tdes, binaries, m_min, m_max, with_tdes, f_accreted, EoS

def BHstar_to_BBH(seed, t, z, k_ex2, N_ex2, m_avg, BH, pairs, binaries, N_BBH, N_BHstar, N_tdeBHstar, v_star, vBH, tdes, m_min, m_max, with_tdes, f_accreted, EoS):
    """
    Creates BH-BH binaries from BH-star pairs, or a TDE occurs during the binary-single interaction.

    @out: all inputs
    """

    if k_ex2>0: # perform BH-star -> BH-BH exchange(s)
        
        for i in range(k_ex2):

            if BH.size == 0: # no single BHs left for exchange, exit loop
                break

            # draw a single BH that will substitute the star in the BH-star pair, with
            # probability ~ (mean pair BH mass + m): a mixture of uniform and mass-weighted draws:
            m_pairs = np.mean(pairs[:, 1])
            if np.random.rand() * (m_pairs * BH.size + BH.total(FAM['mass'])) < m_pairs * BH.size:
                k2 = BH.draw_uniform()
            else:
                k2 = BH.draw(FAM['mass'])
            m2, s2, g2, h2 = BH.get(k2)  # mass, spin, generation and number of TDEs of the second BH
            
            # draw a BH-star pair:
            if np.sum(pairs[:, 0]) == 0:
                continue
            ap, kp = _choose(pairs[:, 0], pairs[:, 0] / np.sum(pairs[:, 0]))
                
            m1 = pairs[kp][1]
            s1 = pairs[kp][2]
            g1 = pairs[kp][3]
            h1 = pairs[kp][4]

            # delete pair:
            pairs = delete_row(pairs, kp)
            
            # probability for TDE:
            p_TDE = min(R_sun*(m2/m_avg)**(1/3)/ap, 1.0) if with_tdes==1 else 0.0

            if np.random.rand() < p_TDE:
                # then TDE occurs and BH-star binary disrupts:

                tde_type = 3 # TDE during BH-star+BH encounter

                # Sample star mass from evolving mass function:
                m_star, R_star = get_star(t, tBH_form, m_min, m_max)

                N_tdeBHstar, m2, s2, g2, h2, tdes = single_BH_TidalDisruption(seed, t, z, N_tdeBHstar, tde_type, m_avg, m_star, R_star, m2, s2, g2, h2, v_star, vBH, tdes, f_accreted, EoS)

                # update secondary BH:
                BH.replace(k2, m2, s2, g2, h2)

                # release BH from BH-star pair into the single population:
                BH.add(m1, s1, g1, h1)

            else:
                # BH substitutes the star and no TDE occurs:

                # semimajor axis from energy conservation:
                sma = ap * m2 / m_avg
            
                # eccentricity:
                eccen = np.sqrt(np.random.rand())
            
                # append binary:
                binaries = np.append(binaries, [[np.random.randint(0, 999999999), 1, sma, eccen, m1, m2, s1, s2, g1, g2, t, z, 0, h1, h2]], axis=0)
            
                BH.remove(k2)
            
                N_BHstar = N_BHstar - 1
            
                N_ex2+=1 # update number of BH-star -> BH-BH exchanges
                N_BBH+=1
            
    return seed, t, z, k_ex2, N_ex2, m_avg, BH, pairs, binaries, N_BBH, N_BHstar, N_tdeBHstar, v_star, vBH, tdes, m_min, m_max, with_tdes, f_accreted, EoS

# end of file
