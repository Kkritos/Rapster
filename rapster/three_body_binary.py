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

def p_3bodyBinary(m1, m2, m3):
    """
    Joint probability density function for masses m1, m2, and m3 to participate in a three-body interaction.
    Masses m1 and m2 form the binary, while m3 carries away the required binding energy in the form of heat.
    The result is not normalized.
    """
    
    return (m1*m2)**(4)*m3**(5/2)/(m1+m2)**(1/2)/(m1+m2+m3)**(1/2)
p_3bodyBinary = np.vectorize(p_3bodyBinary)

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
            else:
                p1 = np.array([p_3bodyBinary(m_1, mBH[mBH!=m_1], np.mean(mBH)).sum() for m_1 in mBH]) # marginalized probability
                p1 /= p1.sum() # normalize margninalized probability
                m1 = np.random.choice(mBH, size=1, replace=False, p=p1)[0] # sample first mass
                p2 = p_3bodyBinary(m1, mBH[mBH!=m1], np.mean(mBH)) # conditional probability
                p2 /= p2.sum() # normalize conditional probability
                m2 = np.random.choice(mBH[mBH!=m1], size=1, replace=False, p=p2)[0] # sample second mass
            
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
            sma = G_Newton * m1 * m2 / eta / mBH_avg / vBH**2
            
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
