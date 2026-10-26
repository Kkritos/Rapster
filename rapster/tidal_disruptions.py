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
from .compact_accretion import *
from .bhpool import FAM

def draw_TDE_BH(BH, m_star):
    """
    Slot of the single BH that disrupts a star of mass m_star, drawn with weights
    p = (m_star + m) m^(1/3) (as np.random.choice(mBH, p=p/np.sum(p)) did): an exact
    mixture of the pool weight families m^(1/3) and m^(4/3).
    """
    T13 = m_star * BH.total(FAM['c13'])
    T43 = BH.total(FAM['c43'])
    return BH.draw(FAM['c13']) if np.random.rand() * (T13 + T43) < T13 else BH.draw(FAM['c43'])

def disrupt_star(seed, t, z, tde_type, m_avg, m_star, R_star, m, s, g, h, m_mean, vSTAR, vBH, tdes, f_accreted, EoS):
    """
    Tidal disruption of a star by the compact object (m, s, g, h): samples the orbit,
    evolves the spin while the accreted debris is added, and records the TDE.

    @in m_mean: mean single-BH mass (Msun), for the relative velocity
    (other inputs as in BH_TidalDisruptions)

    @out m_new, s_new, h_new: mass, spin and TDE count after the disruption
    """

    # orbital inclination:
    iota = np.arccos(np.random.uniform(-1, 1)) # radians

    # marginally bound orbit R_mb(m, s, iota) (functions.py): not needed during the run,
    # so it is computed for all TDE rows at once when the output is written (fill_R_mb;
    # the row stores the placeholder R_MB_PENDING until then):
    r_mb = R_MB_PENDING

    # tidal radius:
    r_t = R_star * (m/m_star)**(1/3)

    # pericenter radius:
    r_p = np.random.rand()*r_t

    # penetration parameter:
    beta = r_t/r_p

    # mass increment:
    dm = f_accreted*m_star

    # direction of accretion:
    prograde = True if np.cos(iota)>0 else False

    # maximum NS mass (TOV) for the chosen EoS:
    M_NS_max = get_eos(EoS).M_TOV

    # nature of compact object:
    NS = True if m<M_NS_max else False

    # final spin after disk accretion (spin_after_accretion: evolve_spin_during_accretion's
    # final spin, integrated with RK4 steps for BHs that grow by <= 2%):
    s_new = spin_after_accretion(Mi=m, Mf=m+dm, NS=NS, chi=s, dM=dm/100, eos=EoS, prograde=prograde)

    # relative velocity:
    v_rel = np.sqrt(vSTAR**2*m_avg/m_star + m_mean/m*vBH**2)

    # append tde:
    tdes.append([seed, t, z, tde_type, m_star, R_star, m, s, g, r_t, r_p, beta, iota, r_mb, dm, s_new, v_rel, h])

    return m + dm, s_new, h + 1

def BH_TidalDisruptions(seed, t, z, k_tde, N_tde, tde_type, m_avg, m_star, R_star, BH, vSTAR, vBH, tdes, binaries, pairs, f_accreted, EoS):
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
    @in BH: pool of single BHs (bhpool.BHPool: masses, spins, generations, TDE counts)
    @in vSTAR: velocity dispersion of stars (km/s)
    @in vBH: velocity dispersion of BHs (km/s)
    @in tdes: tdes array [seed, t, z, type, mstar, Rstar, m_BH, s_BH, g_BH, r_t, r_p, beta, iota, r_mb, dm, s_new, v_rel, h_BH]
    @in binaries: [ind, channel, a, e, m1, m2, s1, s2, g1, g2, t_form, z_form, Nex, h1, h2]
    @in pairs: [a, m, s, g, h]
    @in f_accreted: fraction of the disrupted star accreted
    @in EoS: NS equation of state; a name registered in compact_accretion.EOS_TABLES

    @out: all inputs
    """

    if BH.size==0:
        k_tde = 0

    if k_tde>0: # perform BH TDE(s)
        
        for i in range(k_tde):
            
            if BH.size==0: # BH pool exhausted mid-step
                break

            N_tde+=1 # update number of BH TDEs
            
            # sample the disrupting single BH, with weights p = (m_star + m)*m**(1/3):
            k = draw_TDE_BH(BH, m_star)
            m, s, g, h = BH.get(k)

            m_new, s_new, h_new = disrupt_star(seed, t, z, tde_type, m_avg, m_star, R_star, m, s, g, h,
                                               BH.mean_mass(), vSTAR, vBH, tdes, f_accreted, EoS)

            # update BH mass, spin and TDE count:
            BH.set_mass(k, m_new)
            BH.set_spin(k, s_new)
            BH.set_tde_count(k, h_new)

    return seed, t, z, k_tde, N_tde, tde_type, m_avg, m_star, R_star, BH, vSTAR, vBH, tdes, binaries, pairs, f_accreted, EoS

def single_BH_TidalDisruption(seed, t, z, N_tde, tde_type, m_avg, m_star, R_star, m, s, g, h, vSTAR, vBH, tdes, f_accreted, EoS):
    """
    Tidal disruption of a star by one given compact object (m, s, g, h), e.g. a BH
    in a BH-star pair (as BH_TidalDisruptions with a single BH).

    @out: N_tde, m, s, g, h (updated), tdes
    """
    m, s, h = disrupt_star(seed, t, z, tde_type, m_avg, m_star, R_star, m, s, g, h, m, vSTAR, vBH, tdes, f_accreted, EoS)
    return N_tde + 1, float(m), float(s), float(g), float(h), tdes

def try_BBH_star_disruption(rp, a, m1, m2, s1, s2, g1, g2, h1, h2, m_star, R_star, f_accreted, EoS, seed, t, z, tdes, m_avg, vBH, v_star):
    """
    Checks whether a single passage of a field star past a BBH (at pericenter rp)
    results in tidal disruption, and if so, performs the disruption (Channel A or B).
    Designed to be called once per passage, so it can be reused for both a single
    (non-resonant) encounter and for each individual passage of a resonant encounter.

    @in rp: pericenter of this passage (pc)
    @in a: BBH semimajor axis (pc)
    @in m1, m2: BBH component masses (Msun)
    @in s1, s2: BBH component dimensionless spins
    @in g1, g2: BBH component generations
    @in h1, h2: BBH component tdes counts
    @in m_star: mass of the field star (Msun)
    @in R_star: radius of the field star (pc)
    @in f_accreted: fraction of disrupted star mass accreted
    @in EoS: NS equation of state; a name registered in compact_accretion.EOS_TABLES
    @in seed, t, z: simulation seed, time, redshift
    @in tdes: tdes array [seed, t, z, type, mstar, Rstar, m_BH, s_BH, g_BH, r_t, r_p, beta, iota, r_mb, dm, s_new, v_rel, h_BH]
    @in m_avg: average stellar mass (Msun)
    @in vBH: 3D BH velocity dispersion (km/s)
    @in v_star: 3D star velocity dispersion (km/s)

    @out: disrupted (bool), m1, m2, s1, s2, h1, h2, tdes (updated if disrupted=True, unchanged otherwise)
    """

    # whole-binary (monopole) tidal radius: treats the BBH as a single point mass
    # of total mass (m1+m2), valid when the star's pericenter stays well outside
    # the binary's own separation a:
    r_t_bin = R_star * ((m1+m2)/m_star)**(1/3)

    if rp > r_t_bin:
        # pericenter is outside the binary's tidal radius: no disruption this passage.
        return False, m1, m2, s1, s2, h1, h2, tdes

    # maximum NS mass (TOV) for the chosen EoS:
    M_NS_max = get_eos(EoS).M_TOV

    if rp > a:
        # CHANNEL A: "whole-binary" disruption (monopole regime, rp between a and r_t_bin).
        # tde_type = 22.

        # debris mass fraction accreted by each BH, weighted by its own mass:
        f1 = m1 / (m1+m2)
        f2 = m2 / (m1+m2)

        # mass increment accreted by each BH:
        dm1 = f1 * f_accreted * m_star
        dm2 = f2 * f_accreted * m_star

        NS1 = True if m1<M_NS_max else False
        NS2 = True if m2<M_NS_max else False

        # evolve the spin of each BH/NS as it accretes its share of the debris:
        s1_new = spin_after_accretion(Mi=m1, Mf=m1+dm1, NS=NS1, chi=s1, dM=dm1/100, eos=EoS, prograde=True)
        s2_new = spin_after_accretion(Mi=m2, Mf=m2+dm2, NS=NS2, chi=s2, dM=dm2/100, eos=EoS, prograde=True)

        # relative velocity between the star and the binary's center of mass:
        v_rel_tde = np.sqrt(v_star**2*m_avg/m_star + (m1+m2)/m_star*vBH**2)

        # record one TDE entry per accreting BH; beta = r_t_bin/rp is the penetration parameter:
        tdes.append([seed, t, z, 22, m_star, R_star, m1, s1, g1, r_t_bin, rp, r_t_bin/rp, 0.0, 0.0, dm1, s1_new, v_rel_tde, h1])
        tdes.append([seed, t, z, 22, m_star, R_star, m2, s2, g2, r_t_bin, rp, r_t_bin/rp, 0.0, 0.0, dm2, s2_new, v_rel_tde, h2])

        m1 = m1 + dm1; s1 = s1_new; h1 = h1 + 1
        m2 = m2 + dm2; s2 = s2_new; h2 = h2 + 1

        return True, m1, m2, s1, s2, h1, h2, tdes

    else:
        # CHANNEL B: "single-component" disruption (rp <= a). tde_type = 21.

        # select which BH is the disruptor, weighted by mass:
        p1 = (m_star + m1) * m1**(1/3)
        p2 = (m_star + m2) * m2**(1/3)
        disrupt_1 = np.random.rand() < p1/(p1+p2)

        m_d = m1 if disrupt_1 else m2
        s_d = s1 if disrupt_1 else s2
        g_d = g1 if disrupt_1 else g2
        h_d = h1 if disrupt_1 else h2

        # standard single-BH tidal radius for the chosen disruptor:
        r_t = R_star * (m_d/m_star)**(1/3)

        dm = f_accreted * m_star

        NS = True if m_d<M_NS_max else False

        s_new = spin_after_accretion(Mi=m_d, Mf=m_d+dm, NS=NS, chi=s_d, dM=dm/100, eos=EoS, prograde=True)

        v_rel_tde = np.sqrt(v_star**2*m_avg/m_star + np.mean([m1, m2])/m_star*vBH**2)

        # record the TDE entry; beta = r_t/rp is the penetration parameter:
        tdes.append([seed, t, z, 21, m_star, R_star, m_d, s_d, g_d, r_t, rp, r_t/rp, 0.0, 0.0, dm, s_new, v_rel_tde, h_d])

        if disrupt_1:
            m1 = m_d + dm; s1 = s_new; h1 = h_d + 1
        else:
            m2 = m_d + dm; s2 = s_new; h2 = h_d + 1

        return True, m1, m2, s1, s2, h1, h2, tdes

# End of file.
