

import itertools as it
import logging as log
import os
from datetime import timedelta
from time import time

from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian, quantum_LinearOperator
from scipy.stats import cauchy

from constants import *

# -----------------------------------------------------
#  <model parameters>  !!! (CHANGE PARAMETERS HERE) !!!
# -----------------------------------------------------
# Parameters are stored as a list of dictionaries.
# Append new dictionaries for new parameters you want to test
# Each dictionary should contain the following test parameters
#   as key-value pairs.
# :name      for debugging + name of file to save results in
# :clusters  list of clusters to consider for NLCE: ['0','1','2','3','4Y','4I','4L']
#            the order of NLC-4 clusters must NOT be modified
# :gz        Lande g-factor
# :Jzz       Ising coupling [meV]
# :J+-       symmetric planar exchange [meV]
# :J+-+-     anisotropic planar exchange [meV]
# :temp      temperature [K]
# :B_ext     external magnetic field in global coord [T]
# :hhl       wavevectors (hh0) and (00l)
#            check that 0 is not included --> otherwise might have divide by 0 error
# :trials    number of times to repeat each cluster calculation
#            set to *no* disorder averaging (Recommended for >NLC-2)
#            any Int>1 will activate disorder fields in the Hamiltonian
# :Gamma     width for disorder TransVerse fields Lorentzian [meV]

parameters_list = []
parameters_list.append({
                       'name': 'ns_final0',
                       'clusters': ['0', '1', '2', '3', '4Y', '4I', '4L'],
                       'gz': 4.32,
                       'Jzz': 1,
                       'J+-': 0.,
                       'J+-+-': 0.,
                       'temp': 0.1,
                       'B_ext': [0.0, 0.0, 0],
                       'hhl': np.arange(-2.501, 2.501, 0.025),
                       'trials': 1,
                       'Gamma': 0.19,
                       })
print "----%----------%----------%------STARTING CALCULUS----%----------%----------%----------%------"

def get_thermal_average(eigenvals, eigenvecs,uu_op, temp):t
    """
        Calculates thermal average of the moment-moment correlators,
        for both (SF, NSF) neutron scattering.
        
        :param eigenvals: array of eigenvalues for a given Hamiltonian (ndarray)
        :param eigenvecs: array of eigenvectors for a given Hamiltonian (ndarray)
        :param uu_op: OPerator to calculate thermal average of (site-dependent)
        :param temp: temperature to calculate thermal average at
        :return: averaged correlation
        """
    ss_avg=0
    
    Z=0
    
    # shift eigenvalue spectrum to avoid negative E and large Z values;
    # leaves thermal avg unchanged
    eigenvals = eigenvals-eigenvals.min()
    
    for E_n, n in zip(eigenvals, eigenvecs.T):
        ss_avg += np.dot(n.conj(), uu_op.matvec(n)) * np.exp(-E_n / (KB * temp))
        Z += np.exp(-E_n / (KB * temp))
    
    return ss_avg.real / Z


def construct_hamiltonian(nl, basis,J_zz, J_pm, J_pmpm, gz, B_ext, width= 0.19, tv_field= False):
    """
        Constructs the symmetry-allowed, nearest-neighbour, effective spin-1/2 Hamiltonian for
        non-Kramers;
        Allows for random transverse field model and external magnetic fields to be applied;
        
        :param nlc:      Choice of which NLC-? calculation to perform;
        Corresponds to number and arrangement of tetrahedra (1, 2, 3, 4I, 4L, 4Y);
        :param basis:    Basis corresponding to nlc cluster
        :param J_zz:     Ising spin-ice coupling [meV]
        :param J_pm:     XY-like exchange [meV]
        :param J_pmpm:   In-plane anisotropic exchange [meV]
        :param gz:       Lande g-factor
        :param B_ext:    External magnetic field in global coordinates [T]
        :param width:    width for disorder TransVerse fields Lorentzian [meV]
        :param tv_field: Boolean to activate random TransVerse fields at each site (default False)
        :return:         Hamiltonian object (sparse matrix)
        """
    # generate indices of Nearest Neighbour pairs + Site Types (for Z_DIR) of the N sites
    # cf. pyrochlore_NLC_numbering.pdf
    N = N_DICT[nlc]  # number of sites
    ST = ST_DICT[nlc]  # Site Types
    NN_pairs = NN_PAIRS_DICT[nlc]  # Nearest Neighbour pairings
    
    # define coefficients of operators using site-coupling lists
    Jzz = [[J_zz, i, j] for (i, j) in NN_pairs]
    Jpm = [[-J_pm, i, j] for (i, j) in NN_pairs]
    Jpp = [[J_pmpm * GAMMA[ST[i], ST[j]], i, j] for (i, j) in NN_pairs]
    Jmm = [[J_pmpm * np.conj(GAMMA[ST[i], ST[j]]), i, j] for (i, j) in NN_pairs]
    m_z = [[-(U_B * gz) * Z_DIR[ST[i]].dot(B_ext), i] for i in range(N)]
    
    if tv_field:
        # generate Lorentz distribution for disorder fields
        tv_fields = cauchy(scale=width)
        h_x = [[-np.abs(tv_fields.rvs()), i] for i in range(N)]
    else:
        h_x = []

    # static and dynamic lists
    static = [["zz", Jzz], ["+-", Jpm], ["-+", Jpm], ["++", Jpp], ["--", Jmm],
              ["x", h_x], ["z", m_z]]
    dynamic = []

    # compute the Hamiltonian
    H = hamiltonian(static, dynamic, basis=basis,
                check_herm=False, check_symm=False)
    # print(f"Hamiltonian:\n{H.toarray()}")
    
    return H
















