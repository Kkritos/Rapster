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

def fast_sample_3bodyBinary(mBH, n_samples=1):
    """
    Fast sampling for using rejection sampling with complexity O(N).

    @in mBH: array of BH masses
    @in n_samples: number of samples (m1, m2, m3)
    """
    
    # Pre-calculate proposal weights
    # m1 and m2 are dominated by m^4
    weights_binary = mBH**4
    p_binary = weights_binary / weights_binary.sum()
    
    # m3 is dominated by m^2.5
    weights_third = mBH**2.5
    p_third = weights_third / weights_third.sum()
    
    # Define the 'Interaction' term for the rejection criteria
    # Since (m1+m2)^-0.5 and (m1+m2+m3)^-0.5 are DECREASING functions,
    # the maximum value occurs at the smallest possible masses.
    min_m = np.min(mBH)
    max_interaction_val = (min_m + min_m)**(-0.5) * (min_m + min_m + min_m)**(-0.5)
    
    results = []
    while len(results) < n_samples:
        # Sample candidates
        # We need 3 unique indices
        indices = np.random.choice(len(mBH), size=3, replace=False)
        m1, m2 = mBH[indices[0]], mBH[indices[1]]
        m3 = mBH[indices[2]]
        
        # Calculate the interaction/correction factor
        # P_actual / P_proposal_separable
        interaction_val = (m1 + m2)**(-0.5) * (m1 + m2 + m3)**(-0.5)
        
        # Accept/Reject
        if np.random.rand() < (interaction_val / max_interaction_val):
            results.append((m1, m2, m3))
            
    return results[0] if n_samples == 1 else results

def three_body_binary(t, z, k_3bb, mBH_avg, binaries, mBH, sBH, gBH, hBH, vBH, N_3bb, N_BBH, random_pairing=False):
    """
    @in t: simulation time
    @in z: simulation redshift
    @in k_3bb: number of 3bbs in current step
    @in mBH_avg: average BH mass
    @in binaries: array of BBHs
    @in mBH: array of single BH masses
    @in sBH: array of single BH spins
    @in gBH: array of single BH generations
    @in hBH: array of BH tdes count
    @in vBH: 3D BH velocity dispersion
    @in N_3bb: number of 3bbs
    @in N_BBH: number of BBHs
    @in random_pairing: if True, use uniform random pairing instead of mass-weighted (m^5)

    @out: all inputs
    """
    
    if k_3bb>0:
        
        for i in range(k_3bb):

            N_3bb+=1

            # sample masses that form the 3bb:
            if random_pairing:
                m1, m2 = np.random.choice(mBH, size=2, replace=False)
                m3 = mBH_avg
            else:
                m1, m2, m3 = fast_sample_3bodyBinary(mBH, n_samples=1)
            
            # find index locations of the sampled BHs:
            k1 = np.squeeze(np.where(mBH==m1))+0
            k2 = np.squeeze(np.where(mBH==m2))+0
            
            if isinstance(k1, np.ndarray):
                k1 = k1[0]
            if isinstance(k2, np.ndarray):
                k2 = k2[0]
                
            # initial hardness of newly formed 3bb:
            eta = sample_hardness()
            
            # semimajor axis:
            sma = G_Newton * m1 * m2 / eta / m3 / vBH**2
            
            # eccentricity (thermal):
            eccen = np.sqrt(np.random.rand())
            
            # append binary:
            binaries = np.append(binaries, [[np.random.randint(0, 999999999), 3, sma, eccen, m1, m2, sBH[k1], sBH[k2], gBH[k1], gBH[k2], t, z, 0, hBH[k1], hBH[k2]]], axis=0)
            
            # update number of in-cluster BBHs:
            N_BBH+=1
            
            # delete single BHs:
            mBH = np.delete(mBH, [k1, k2])
            sBH = np.delete(sBH, [k1, k2])
            gBH = np.delete(gBH, [k1, k2])
            hBH = np.delete(hBH, [k1, k2])
            
    return t, z, k_3bb, mBH_avg, binaries, mBH, sBH, gBH, hBH, vBH, N_3bb, N_BBH

# End of file.
