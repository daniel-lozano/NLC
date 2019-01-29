
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

hh,l=np.meshgrid(q,q)

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
plt.close()


Tot_contributions_NSF=cluster_contributions(NSF,4)
Tot_contributions_SF=cluster_contributions(SF,4)

Bare_sums_SF=bare_sum(Tot_contributions_SF)
Bare_sums_NSF=bare_sum(Tot_contributions_NSF)

Euler_sums_SF=euler_sum(Tot_contributions_SF,3,4)
Euler_sums_NSF=euler_sum(Tot_contributions_NSF,3,4)

#print SF.shape
#print Tot_contributions_SF.shape
#print Bare_sums_SF.shape
#print NSF.shape
#print Tot_contributions_NSF.shape
#print Bare_sums_NSF.shape
#
#levels=np.linspace(np.amin(Euler_sums_SF[1]),np.amax(Euler_sums_SF[1]),50)
#im=plt.contour(hh,l,Euler_sums_SF[1],cmap="gist_heat",levels=levels)
#plt.colorbar(im,orientation="vertical")
#plt.title("Bare sum NSF order=" +str(1))
#plt.show()
#


for i in range(1,Bare_sums_SF.shape[0]):
    
    plt.figure(figsize=(12,5))
    plt.subplot(121)
    im2=plt.imshow(Bare_sums_NSF[i],cmap="gist_heat")
    plt.colorbar(im2,orientation="vertical")
    plt.title("Bare sum NSF order= " +str(i))
    
    levels=np.linspace(np.amin(Bare_sums_SF[i]),np.amax(Bare_sums_SF[i]),50)
    
    plt.subplot(122)
    im1=plt.contour(hh,l,Bare_sums_SF[i],cmap="gist_heat",levels=levels)
    plt.colorbar(im1,orientation="vertical")
    plt.title("Bare sum SF order=" +str(i))

    plt.show()

for i in range(Euler_sums_SF.shape[0]):
    
    plt.figure(figsize=(7,7))
    plt.subplot(221)
    im1=plt.imshow(Bare_sums_NSF[i-2],cmap="gist_heat")
    plt.colorbar(im1,orientation="vertical")
    plt.title("Bare sum NSF order=" +str(i+3))
    
    plt.subplot(222)
    levels=np.linspace(np.amin(Bare_sums_SF[i-2]),np.amax(Bare_sums_SF[i-2]),50)
    im2=plt.contour(hh,l,Bare_sums_SF[i-2],cmap="gist_heat",levels=levels)
    plt.colorbar(im2,orientation="vertical")
    plt.title("Bare sum SF order= " +str(i+3))

    plt.subplot(223)
    im1=plt.imshow(Euler_sums_NSF[i],cmap="gist_heat")
    plt.colorbar(im1,orientation="vertical")
    plt.title("Euler sum NSF order=" +str(i+3))
    
    plt.subplot(224)
    levels=np.linspace(np.amin(Euler_sums_SF[i]),np.amax(Euler_sums_SF[i]),50)
    im2=plt.contour(hh,l,Euler_sums_SF[i],cmap="gist_heat",levels=levels)
    plt.colorbar(im2,orientation="vertical")
    plt.title("Euler sum SF order= " +str(i+3))

    plt.savefig("Euler_vs_Bare_sums_order"+str(i+3)+".png")
    plt.show()




