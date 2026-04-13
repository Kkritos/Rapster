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

def StarStar_to_BHstar(k_ex1, N_ex1, m_avg, mBH, sBH, gBH, hBH, ab, pairs, N_BHstar):
    """
    @in k_ex1: current number of star-star -> BH-star exchanges
    @in N_ex1: total number of star-star -> BH-star exchanges
    @in m_avg: average mass
    @in mBH: array of single BH masses
    @in sBH: array of single BH spins
    @in gBH: array of single BH generations
    @in hBH: array of BH tdes count
    @in ab: array of star-star semimajor axes
    @in pairs: array of BH-star pairs
    @in N_BHstar: number of BH-star pairs
    
    @out: all inputs
    """
    
    if k_ex1>0: # perform star-star -> BH-star exchange(s)
        
        for i in range(k_ex1):
            
            N_ex1+=1
            
            # sample mass of the BH that substitutes one of the stars:
            m = np.random.choice(mBH, p=mBH/np.sum(mBH))
            
            k = np.squeeze(np.where(mBH==m))+0
            
            if isinstance(k, np.ndarray):
                k=k[0]
                
            s = sBH[k]
            g = gBH[k]
            h = hBH[k]
                
            a = np.random.choice(ab, p=ab/np.sum(ab))
                
            kss = np.squeeze(np.where(ab==a))+0
            
            if isinstance(kss, np.ndarray):
                kss=kss[0]
                
            a = a * m / m_avg
            
            # append pair:
            pairs = np.append(pairs, [[a, m, s, g, h]], axis=0)
            
            # delete star-star:
            ab = np.delete(ab, kss)
            
            N_BHstar+=1
                
            mBH = np.delete(mBH, k)
            sBH = np.delete(sBH, k)
            gBH = np.delete(gBH, k)
            hBH = np.delete(hBH, k)

    return k_ex1, N_ex1, m_avg, mBH, sBH, gBH, hBH, ab, pairs, N_BHstar

def BHstar_to_BBH(t, z, k_ex2, N_ex2, m_avg, mBH, sBH, gBH, hBH, pairs, binaries, N_BBH, N_BHstar):
    """
    Creates BH-BH binaries from BH-star pairs.    

    @out: all inputs
    """
    
    if k_ex2>0: # perform BH-star -> BH-BH exchange(s)
        
        for i in range(k_ex2):
            
            N_ex2+=1
            
            # draw a single BH that will substitute the star in the BH-star pair:
            m2 = np.random.choice(mBH, p=(np.mean(pairs[:, 1]) + mBH)/np.sum(np.mean(pairs[:, 1]) + mBH))
            
            # location of the sampled BH:
            k2 = np.squeeze(np.where(mBH==m2))+0
            
            if isinstance(k2, np.ndarray):
                k2=k2[0]
                
            s2 = sBH[k2] # spin of the second BH
            g2 = gBH[k2] # generation of the second BH
            h2 = hBH[k2] # number of second BH's tdes
            
            # draw a BH-star pair:
            ap = np.random.choice(pairs[:, 0], p=pairs[:, 0] / np.sum(pairs[:, 0]))
            
            kp = np.squeeze(np.where(pairs[:, 0]==ap))+0
            
            if isinstance(kp, np.ndarray):
                kp=kp[0]
                
            m1 = pairs[kp][1]
            s1 = pairs[kp][2]
            g1 = pairs[kp][3]
            h1 = pairs[kp][4]
            
            # semimajor axis:
            sma = ap * m2 / m_avg
            
            # eccentricity:
            eccen = np.sqrt(np.random.rand())
            
            # append binary:
            binaries = np.append(binaries, [[np.random.randint(0, 999999999), 1, sma, eccen, m1, m2, s1, s2, g1, g2, t, z, 0, h1, h2]], axis=0)
            
            # delete pair:
            pairs = np.delete(pairs, kp, axis=0)
            
            mBH = np.delete(mBH, k2)
            sBH = np.delete(sBH, k2)
            gBH = np.delete(gBH, k2)
            hBH = np.delete(hBH, k2)
            
            N_BHstar = N_BHstar - 1
            
            N_BBH+=1
            
    return t, z, k_ex2, N_ex2, m_avg, mBH, sBH, gBH, hBH, pairs, binaries, N_BBH, N_BHstar

# end of file
