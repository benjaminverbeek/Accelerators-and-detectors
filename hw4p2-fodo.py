# Homework 4 problem 2: FODO

import numpy as np

# function returning the map for drift space of length l
def drift_space(l):
    # return 2x2 matrix
    return np.array([[1, l], [0, 1]])

# function returning the map for quadropoles of length l, constant k
def quad_map(k, l, focusing=True, thin=False):
    if not focusing:    # defocusing
        l = -l  # negative focus point
    if thin or not thin: # not thin non functioning... need imaginary nums?
        return np.array([[1, 0], [-(k*l), 1]])
     # return 2x2 matrix
    '''
    if k < 0:
        sqK = np.sqrt(abs(k))   #  NOTE: correct?
        return np.array([[np.cos(l*sqK), 1/sqK * np.sin(l*sqK) ], [-sqK*np.sin(l*sqK), np.cos(l*sqK)]])
    else:
        sqK = np.sqrt(k)
        return np.array([[np.cos(l*sqK), 1/sqK * np.sin(l*sqK) ], [-sqK*np.sin(l*sqK), np.cos(l*sqK)]])
        '''

k1, k2 = 3.5, 3.5
M1 = drift_space(2.5)
M2 = quad_map(k1, 0.15, thin=True)
M3 = drift_space(0.002)
M4 = quad_map(k2, 0.30, focusing=False, thin=True)
M5 = drift_space(0.002) 
M6 = quad_map(k1, 0.15, thin=True)
M7 = drift_space(2.0)

sx = 0.60 # mm
sxp= 0.50 # mrad
sx12 = 0
sy = 0.25 # mm
syp= 1.30 # mrad
sy12 = 0.22*10**-6

# sigma matrix for x
sigma_x = np.array([[sx**2, sx12], [sx12, sxp**2]])
# sigma matrix for y
sigma_y = np.array([[sy**2, sy12], [sy12, syp**2]])

system_matrix = M7 @ M6 @ M5 @ M4 @ M3 @ M2 @ M1
print(system_matrix, "system matrix")
sigma_x_prop = system_matrix @ sigma_x @ system_matrix.T
print(sigma_x_prop, "sigma_x_prop")

k1, k2 = 3.5, 3.5
M1 = drift_space(2.5)
M2 = quad_map(k1, 0.15, thin=True, focusing=False)
M3 = drift_space(0.002)
M4 = quad_map(k2, 0.30, focusing=True, thin=True)
M5 = drift_space(0.002) 
M6 = quad_map(k1, 0.15, thin=True, focusing=False)
M7 = drift_space(2.0)
sigma_y_prop = system_matrix @ sigma_y @ system_matrix.T
print(sigma_y_prop, "sigma_y_prop")