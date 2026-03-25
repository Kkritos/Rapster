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

import numpy as np
import pickle
from scipy.optimize import fsolve
from scipy.interpolate import CloughTocher2DInterpolator
import os

def x_mb(a_BH, iota):
    """
    Marginally bound orbit normalized to gravitational radius.
    """
    def fun(x):
        # Ensure x is positive and physically sensible
        if x <= 0: return 1e10 
        
        inner = x**5 - a_BH**2 * x**3 * np.cos(iota)**2
        # If the sqrt interior is negative, return a high penalty instead of NaN
        if inner < 0: return 1e10
        
        return x**4 - 4*x**3 - a_BH**2 * (1 - 3*np.cos(iota)**2)*x**2 + a_BH**4 * np.cos(iota)**2 + 4 * a_BH * np.sin(iota) * np.sqrt(x**5 - a_BH**2 * x**3 * np.cos(iota)**2)
        
    x = fsolve(fun, 4.0)
    return x
x_mb = np.vectorize(x_mb)

# bulk of domain:
N_points_bulk = 10**3

a_BH_vec = np.random.uniform(-1, 1, N_points_bulk)
iota_vec = np.random.uniform(0, np.pi/2, N_points_bulk)

# boundary of parameter space:
N_points_boundary = 10**3

a_BH_vec = np.concatenate([a_BH_vec, np.ones(N_points_boundary)])
iota_vec = np.concatenate([iota_vec, np.random.uniform(0, np.pi/2, N_points_boundary)])

a_BH_vec = np.concatenate([a_BH_vec, -np.ones(N_points_boundary)])
iota_vec = np.concatenate([iota_vec, np.random.uniform(0, np.pi/2, N_points_boundary)])

a_BH_vec = np.concatenate([a_BH_vec, np.random.uniform(-1, 1, N_points_boundary)])
iota_vec = np.concatenate([iota_vec, np.pi/2*np.ones(N_points_boundary)])

a_BH_vec = np.concatenate([a_BH_vec, np.random.uniform(-1, 1, N_points_boundary)])
iota_vec = np.concatenate([iota_vec, np.zeros(N_points_boundary)])

# corners of domain:
a_BH_vec = np.concatenate([a_BH_vec, np.array([-1, -1, 1, 1])])
iota_vec = np.concatenate([iota_vec, np.array([0, np.pi/2, 0, np.pi/2])])

# compute values in the scatter grid:
x_mb_vec = x_mb(a_BH_vec, iota_vec)

if __name__=="__main__":

    # Define the path to your Data folder:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    pickle_path = os.path.join(BASE_DIR, 'x_mb_vec.pkl')

    # Bundle variables into a dictionary:
    data_to_save = {
        'x_mb_vec': x_mb_vec,
        'a_BH_vec': a_BH_vec,
        'iota_vec': iota_vec
    }

    # Save (Dump) the vector:
    with open(pickle_path, 'wb') as f:
        pickle.dump(data_to_save, f)

    print(f"Vector successfully saved to {pickle_path}")

# End of file.
