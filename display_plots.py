
import itertools as it
import logging as log
import os
import matplotlib.pyplot as plt
from datetime import timedelta
from time import time

#from quspin.basis import spin_basis_1d
#from quspin.operators import hamiltonian, quantum_LinearOperator
#from scipy.stats import cauchy

from constants import *
from sums import *



print("----%----------%----------%------STARTING CALCULUS----%----------%----------%----------%------")


Jzz=1. #First neighbors interaction constant
Jpm=0.0
B_field=[0.0,0.0,0.0] # field
T=0.1/KB # Kelvin
gz=4.32 #Lande factor

q=2*np.pi*np.arange(-2.501, 2.501, 0.1)#0.025

cluster=['0','1','2','3','4Y','4I','4L']

#Arrays to store the different scattering results
c_SF_intensity = np.zeros((len(cluster), q.size, q.size))
c_NSF_intensity = np.zeros((len(cluster), q.size, q.size))
print c_NSF_intensity[0][0,0]


print("First Neighbor interaction constant=",Jzz)
print("First Neighbor exchange=",Jpm)
print("Magnetic field=",B_field)

SF=np.load("data.npz")['SF']
NSF=np.load("data.npz")['NSF']
print SF.shape[0]

for i in range(SF.shape[0]):
    
    c_SF_intensity=SF[i]
    c_NSF_intensity=NSF[i]
    
    print "variation",abs(np.amax(c_NSF_intensity)-np.amin(c_NSF_intensity))
    plt.figure(figsize=(12,5))
    plt.subplot(121)
    im1=plt.imshow(c_NSF_intensity,cmap="gist_heat")
    plt.colorbar(im1,orientation="vertical")
    plt.title("NSF,Cluster type="+str(i))

    plt.subplot(122)
    im2=plt.imshow(c_SF_intensity,cmap="gist_heat")
    plt.colorbar(im2,orientation="vertical")
    plt.title("SF,Cluster type="+str(i))

    plt.show()


Tot_contributions_NSF=cluster_contributions(NSF,4)
Tot_contributions_SF=cluster_contributions(NSF,4)
print SF.shape
print Tot_contributions_SF.shape
print NSF.shape
print Tot_contributions_NSF.shape



