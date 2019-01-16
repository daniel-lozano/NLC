#!/usr/bin/env python

# ==============================================================================
# WRITTEN BY: Jeremy Goh (for OFYP 18/19)
# SUPERVISOR: Prof. Michel Gingras
#
# Stores the constants pertaining to the pyrochlore lattice, i.e., variables
#   independent of the NLC order being considered, or parameters that should
#   *not* be changed across different runs/trials.
#
# Imported into:
#   - neutron_scattering_main.py
#   - thermodynamic_main.py
# ==============================================================================

import numpy as np


def rotation(axis, angle):
    """
        Generates matrix for rotating coordinates in 3D about the given axis.
        
        :param axis: Rotation axis: np.ndarray
        :param angle: Rotation angle in radians: float
        :return: The 3-D rotation matrix R  -> np.ndarray
        """
    
    # axis must be normalized first
    axis /= np.sqrt(axis.dot(axis))
    
    R = np.array([
                  [np.cos(angle) + axis[0] ** 2 * (1 - np.cos(angle)),
                   axis[0] * axis[1] * (1 - np.cos(angle)) - axis[2] * np.sin(angle),
                   axis[0] * axis[2] * (1 - np.cos(angle)) + axis[1] * np.sin(angle)],
                  
                  [axis[1] * axis[0] * (1 - np.cos(angle)) + axis[2] * np.sin(angle),
                   np.cos(angle) + axis[1] ** 2 * (1 - np.cos(angle)),
                   axis[1] * axis[2] * (1 - np.cos(angle)) - axis[0] * np.sin(angle)],
                  
                  [axis[2] * axis[0] * (1 - np.cos(angle)) - axis[1] * np.sin(angle),
                   axis[2] * axis[1] * (1 - np.cos(angle)) + axis[0] * np.sin(angle),
                   np.cos(angle) + axis[2] ** 2 * (1 - np.cos(angle))]
                  ])
                  
    return R

# physical constants
KB = 0.086173  # Boltzmann constant [meV/K]
U_B = 5.788E-2  # Bohr magneton [meV/T]

# pyrochlore coordinates in global coord
# XYZ = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # global cubic axes (not basis vectors)
# SITE_R = np.array([[0, 0, 0], [1, 1, 0], [1, 0, 1], [0, 1, 1],  # pyrochlore 0 (sites 0-3)
#                    [2, 2, 0], [3, 3, 0], [3, 2, 1], [2, 3, 1],  # pyrochlore 1 (sites 4-7)
#                    [2, 0, 2], [3, 1, 2], [3, 0, 3], [2, 1, 3],  # pyrochlore 2 (sites 8-11)
#                    [0, 2, 2], [1, 3, 2], [1, 2, 3], [0, 3, 3]  # pyrochlore 3 (sites 12-15)
#                    ]) / 4  # displacement vectors for all 16 sites in the conventional unit cell

# local z axes for the 4 site types of pyrochlore lattice
Z_DIR = np.array([[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]]) / np.sqrt(3)

# cluster-specific look-up dictionaries
# number of sites for each cluster type
N_DICT = {'0': 1, '1': 4, '2': 7, '3': 10, '4Y': 13, '4I': 13, '4L': 13}

# site types for each cluster type
ST_DICT = {'0': [0],
    '1': [0, 1, 2, 3],
    '2': [0, 1, 2, 3, 3, 0, 1],
    '3': [0, 1, 2, 3, 3, 0, 1, 1, 2, 3],
    '4Y': [0, 1, 2, 3, 3, 0, 1, 1, 2, 3, 2, 3, 0],
    '4I': [0, 1, 2, 3, 3, 0, 1, 1, 2, 3, 3, 0, 1],
    '4L': [0, 1, 2, 3, 3, 0, 1, 1, 2, 3, 0, 1, 2]
            }
# nearest neighbour pairings for each cluster type
NN_PAIRS_DICT = {'0': [],
    '1': [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
    '2': [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
          (2, 4), (2, 5), (2, 6), (4, 5), (4, 6), (5, 6)],
    '3': [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
        (2, 4), (2, 5), (2, 6), (4, 5), (4, 6), (5, 6),
        (5, 7), (5, 8), (5, 9), (7, 8), (7, 9), (8, 9)],
    '4Y': [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
        (2, 4), (2, 5), (2, 6), (4, 5), (4, 6), (5, 6),
        (5, 7), (5, 8), (5, 9), (7, 8), (7, 9), (8, 9),
        (6, 10), (6, 11), (6, 12), (10, 11), (10, 12), (11, 12)],
    '4I': [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
          (2, 4), (2, 5), (2, 6), (4, 5), (4, 6), (5, 6),
          (5, 7), (5, 8), (5, 9), (7, 8), (7, 9), (8, 9),
          (8, 10), (8, 11), (8, 12), (10, 11), (10, 12), (11, 12)],
    '4L': [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
          (2, 4), (2, 5), (2, 6), (4, 5), (4, 6), (5, 6),
          (5, 7), (5, 8), (5, 9), (7, 8), (7, 9), (8, 9),
          (9, 10), (9, 11), (9, 12), (10, 11), (10, 12), (11, 12)]
          }
# position vectors to sites in NLCE clusters
R_DICT = {'0': np.array([[0, 0, 0]]),
    '1': np.array([[0, 0, 0],[1, 1, 0], [1, 0, 1], [0, 1, 1]]) / 4,
    '2': np.array([[0, 0, 0],[1, 1, 0], [1, 0, 1], [0, 1, 1],
                [2, -1, 1], [2, 0, 2], [1, -1, 2]]) / 4,
    '3': np.array([[0, 0, 0],[1, 1, 0], [1, 0, 1], [0, 1, 1],
                   [2, -1, 1], [2, 0, 2], [1, -1, 2],
                    [3, 1, 2], [3, 0, 3], [2, 1, 3]]) / 4,
    '4Y': np.array([[0, 0, 0],[1, 1, 0], [1, 0, 1], [0, 1, 1],
                    [2, -1, 1], [2, 0, 2], [1, -1, 2],
                    [3, 1, 2], [3, 0, 3], [2, 1, 3],
                    [1, -2, 3], [0, -1, 3], [0, -2, 2]]) / 4,
    '4I': np.array([[0, 0, 0],[1, 1, 0], [1, 0, 1], [0, 1, 1],
                    [2, -1, 1], [2, 0, 2], [1, -1, 2],
                    [3, 1, 2], [3, 0, 3], [2, 1, 3],
                    [4, -1, 3], [4, 0, 4], [3, -1, 4]]) / 4,
    '4L': np.array([[0, 0, 0],[1, 1, 0], [1, 0, 1], [0, 1, 1],
                    [2, -1, 1], [2, 0, 2], [1, -1, 2],
                    [3, 1, 2], [3, 0, 3], [2, 1, 3],
                    [2, 2, 4], [1, 1, 4], [1, 2, 3]]) / 4
        }

# bond-dependent phase factors
PHASE = 0.5 + 1j * np.sqrt(3) / 2  # exp(i*pi/3)
PHASEC = np.conj(PHASE)  # exp(-i*pi/3)
GAMMA = np.array([[0, 1, -PHASEC, -PHASE],
                  [1, 0, -PHASE, -PHASEC],
                  [-PHASEC, -PHASE, 0, 1],
                  [-PHASE, -PHASEC, 1, 0]])

# All symmetry allowed rotation matrices for O_h point group symmetry for pyrochlores
SYM = np.zeros((48, 3, 3))

SYM[0] = rotation(np.array([0, 0, 1.0]), 0.0)
SYM[1] = rotation(np.array([+1.0, +1.0, +1.0]), 2 * np.pi / 3)
SYM[2] = rotation(np.array([+1.0, -1.0, -1.0]), 2 * np.pi / 3)
SYM[3] = rotation(np.array([-1.0, +1.0, -1.0]), 2 * np.pi / 3)
SYM[4] = rotation(np.array([-1.0, -1.0, +1.0]), 2 * np.pi / 3)
SYM[5] = rotation(np.array([+1.0, +1.0, +1.0]), 4 * np.pi / 3)
SYM[6] = rotation(np.array([+1.0, -1.0, -1.0]), 4 * np.pi / 3)
SYM[7] = rotation(np.array([-1.0, +1.0, -1.0]), 4 * np.pi / 3)
SYM[8] = rotation(np.array([-1.0, -1.0, +1.0]), 4 * np.pi / 3)
SYM[9] = rotation(np.array([1.0, 0, 0]), np.pi / 2)
SYM[10] = rotation(np.array([0, 1.0, 0]), np.pi / 2)
SYM[11] = rotation(np.array([0, 0, 1.0]), np.pi / 2)
SYM[12] = rotation(np.array([1.0, 0, 0]), 2 * np.pi / 2)
SYM[13] = rotation(np.array([0, 1.0, 0]), 2 * np.pi / 2)
SYM[14] = rotation(np.array([0, 0, 1.0]), 2 * np.pi / 2)
SYM[15] = rotation(np.array([1.0, 0, 0]), 3 * np.pi / 2)
SYM[16] = rotation(np.array([0, 1.0, 0]), 3 * np.pi / 2)
SYM[17] = rotation(np.array([0, 0, 1.0]), 3 * np.pi / 2)
SYM[18] = rotation(np.array([0, +1.0, -1.0]), np.pi)
SYM[19] = rotation(np.array([+1.0, 0, -1.0]), np.pi)
SYM[20] = rotation(np.array([+1.0, -1.0, 0]), np.pi)
SYM[21] = rotation(np.array([0, +1.0, +1.0]), np.pi)
SYM[22] = rotation(np.array([+1.0, 0, +1.0]), np.pi)
SYM[23] = rotation(np.array([+1.0, +1.0, 0]), np.pi)

SYM[24] = -rotation(np.array([0, 0, 1.0]), 0.0)
SYM[25] = -rotation(np.array([+1.0, +1.0, +1.0]), 2 * np.pi / 3)
SYM[26] = -rotation(np.array([+1.0, -1.0, -1.0]), 2 * np.pi / 3)
SYM[27] = -rotation(np.array([-1.0, +1.0, -1.0]), 2 * np.pi / 3)
SYM[28] = -rotation(np.array([-1.0, -1.0, +1.0]), 2 * np.pi / 3)
SYM[29] = -rotation(np.array([+1.0, +1.0, +1.0]), 4 * np.pi / 3)
SYM[30] = -rotation(np.array([+1.0, -1.0, -1.0]), 4 * np.pi / 3)
SYM[31] = -rotation(np.array([-1.0, +1.0, -1.0]), 4 * np.pi / 3)
SYM[32] = -rotation(np.array([-1.0, -1.0, +1.0]), 4 * np.pi / 3)
SYM[33] = -rotation(np.array([1.0, 0, 0]), np.pi / 2)
SYM[34] = -rotation(np.array([0, 1.0, 0]), np.pi / 2)
SYM[35] = -rotation(np.array([0, 0, 1.0]), np.pi / 2)
SYM[36] = -rotation(np.array([1.0, 0, 0]), 2 * np.pi / 2)
SYM[37] = -rotation(np.array([0, 1.0, 0]), 2 * np.pi / 2)
SYM[38] = -rotation(np.array([0, 0, 1.0]), 2 * np.pi / 2)
SYM[39] = -rotation(np.array([1.0, 0, 0]), 3 * np.pi / 2)
SYM[40] = -rotation(np.array([0, 1.0, 0]), 3 * np.pi / 2)
SYM[41] = -rotation(np.array([0, 0, 1.0]), 3 * np.pi / 2)
SYM[42] = -rotation(np.array([0, +1.0, -1.0]), np.pi)
SYM[43] = -rotation(np.array([+1.0, 0, -1.0]), np.pi)
SYM[44] = -rotation(np.array([+1.0, -1.0, 0]), np.pi)
SYM[45] = -rotation(np.array([0, +1.0, +1.0]), np.pi)
SYM[46] = -rotation(np.array([+1.0, 0, +1.0]), np.pi)
SYM[47] = -rotation(np.array([+1.0, +1.0, 0]), np.pi)


