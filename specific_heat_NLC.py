import itertools as it
import logging as log
import os
import matplotlib.pyplot as plt
from datetime import timedelta
from time import time

from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian, quantum_LinearOperator
from scipy.stats import cauchy

from constants import *
from functions import *
from sums import *


print("----%----------%----------%------STARTING CALCULUS----%----------%----------%----------%------")


Jzz=1 #First neighbors interaction constant
Jpm=0
Jppmm=0.
Jzpm=0
h=0. # Field strength, should be in Teslas
B_field=np.array([1.0,1.0,0.0])*h # Field direction
gz=4.32 #Lande factor
gzz=1.80
gxy=4.32
Eta=[[0,-1,np.exp(np.pi*1j/3),np.exp(-np.pi*1j/3)],
     [-1,0,np.exp(-np.pi*1j/3),np.exp(np.pi*1j/3)],
     [np.exp(np.pi*1j/3),np.exp(-np.pi*1j/3),0,-1],
     [np.exp(-np.pi*1j/3),np.exp(np.pi*1j/3),-1,0]]

Gamma=-np.conj(Eta)

cluster=['0','1','2','3','4Y','4I','4L']

#Arrays to store the different scattering results

T_array=np.logspace(-2,3,50)
### Array for the magnetization ###
#Magnetization=np.zeros((len(cluster),len(T_array)))
### Array for the specific heat ###
C_h=np.zeros((len(cluster),len(T_array)))
C_E=np.zeros((len(cluster),len(T_array)))


print("First Neighbor interaction constant(Jzz)=",Jzz)
print("First Neighbor exchange(Jpm)=",Jpm)
print("Double excitation(Jppmm)=",Jppmm)
print("Mixed(Jzpm)=",Jzpm)
print("Magnetic field=",B_field)

Z_scatt=np.array([1, -1, 0]) / np.sqrt(2) # Scattering Polarization direction
Z_magne=np.array([1, 1, 0]) / np.sqrt(2)

for i in range(len(T_array)):
    T=T_array[i]
#    print("Temperature=%f" %T)

    c_magnetization=[]
    c_heat=[]
    c_Energy=[]
    
    for c in range(len(cluster)-4):
        
        magnetization=0 ### define variable that will contain the magnetization of the cluster
        
        N=N_DICT[cluster[c]] #Size of the system
        
        ST=ST_DICT[cluster[c]] #Site types
        
        NN=NN_PAIRS_DICT[cluster[c]] #Pairs in the system
        
        Positions=R_DICT[cluster[c]] #Positions of the sites in the lattice
        
        basis=spin_basis_1d(L=N,S='1/2',pauli=True)#Basis of the system

        H=Hamiltonian(basis,N,Jzz,Jpm,Jppmm,Jzpm,B_field,NN,ST)
        #Obtaining the eigen values and eigen vectors of the Hamiltonian
        eigenvals,eigenvect=H.eigh()
        
        c_heat.append(specific_heat(eigenvals,eigenvect,T)/N)
        c_Energy.append(Energy(eigenvals,eigenvect,T))

    ### Getting the contribution of each cluster for the magnetization ###

#    contributions=cluster_contributions_array(np.array(c_magnetization),min(c,4))
    contributions_heat=cluster_contributions_array(np.array(c_heat),min(c,4))
    contributions_energy=cluster_contributions_array(np.array(c_Energy),min(c,4))

    for k in range(len(contributions_heat)):
#        Magnetization[k,i]=(sum(contributions[:(k+1)]))/h   ### usual order k ###
        C_h[k,i]=(sum(contributions_heat[:(k+1)])) ### divided by h
        C_E[k,i]=(sum(contributions_energy[:(k+1)]))

marker=["--",".-","^-","*-"]
for i in range(1,len(contributions_heat)):
    plt.subplot(121)
    plt.semilogx(T_array,C_E[i,:],marker[i-1],label=str(i))
    if(i==1):
        plt.xlabel("$ T/J $")
        plt.title(" $ E/N $")
    plt.legend()
    plt.subplot(122)
    plt.semilogx(T_array,C_h[i,:],marker[i-1],label=str(i))


plt.xlabel("$ T/J $")
plt.legend()
plt.title(" $ C_v $")
plt.show()


for i in range(1,len(contributions_heat)):
    plt.subplot(121)
    plt.semilogx(T_array,abs(C_E[i,:]),marker[i-1],label=str(i))
    plt.legend()
    if(i==1):
        plt.xlabel("$ T/J $")
        plt.title(" $ E/N $")
    
    plt.subplot(122)
    plt.semilogx(T_array,abs(C_h[i,:]),marker[i-1],label=str(i))

plt.xlabel("$ T/J $")
plt.legend()
plt.title(" $ C_v $")
#plt.ylim(0.,100)
#plt.ylim(0.2,2.6)
plt.show()


