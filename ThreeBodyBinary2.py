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

#@njit
def three_body_binary(t, z, k_3bb, mBH_avg, binaries, mBH, sBH, gBH, vBH, N_3bb, N_BBH):
    """
    @in t: simulation time
    @in z: simulation redshift
    @in k_3bb: number of 3bbs in current step
    @in mBH_avg: average BH mass
    @in binaries: array of BBHs
    @in mBH: array of single BH masses
    @in sBH: array of single BH spins
    @in gBH: array of single BH generations
    @in vBH: 3D BH velocity dispersion
    @in N_3bb: number of 3bbs
    @in N_BBH: number of BBHs
    
    @out: all inputs
    """
    
    if k_3bb>0:
        
        for i in range(k_3bb):

            N_3bb+=1

            m1, m2 = np.random.choice(mBH, size=2, replace=False, p=mBH**(5) / np.sum(mBH**(5)))
            
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
            binaries = np.append(binaries, [[np.random.randint(0, 999999999), 3, sma, eccen, m1, m2, sBH[k1], sBH[k2], gBH[k1], gBH[k2], t, z, 0]], axis=0)
            
            # update number of in-cluster BBHs:
            N_BBH+=1
            
            # delete single BHs:
            mBH = np.delete(mBH, [k1, k2])
            sBH = np.delete(sBH, [k1, k2])
            gBH = np.delete(gBH, [k1, k2])
            
    return t, z, k_3bb, mBH_avg, binaries, mBH, sBH, gBH, vBH, N_3bb, N_BBH

# end of file
