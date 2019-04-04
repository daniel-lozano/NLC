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


Jzz=0.17 #First neighbors interaction constant
Jpm=0.05
Jppmm=0.05
Jzpm=-0.14
h=0.2 # Field strength, should be in Teslas
B_field=np.array([1.0,0.0,0.0])*h # Field direction
gz=4.32 #Lande factor
gzz=1.80
gxy=4.32
Eta=[[0,-1,np.exp(np.pi*1j/3),np.exp(-np.pi*1j/3)],
     [-1,0,np.exp(-np.pi*1j/3),np.exp(np.pi*1j/3)],
     [np.exp(np.pi*1j/3),np.exp(-np.pi*1j/3),0,-1],
     [np.exp(-np.pi*1j/3),np.exp(-np.pi*1j/3),-1,0]]

Gamma=-np.conj(Eta)

cluster=['0','1','2','3','4Y','4I','4L']

#Arrays to store the different scattering results

T_array=np.linspace(1,10,20)
Magnetization=np.zeros((len(cluster),len(T_array)))

print("First Neighbor interaction constant(Jzz)=",Jzz)
print("First Neighbor exchange(Jpm)=",Jpm)
print("Double excitation(Jppmm)=",Jppmm)
print("Mixed(Jzpm)=",Jzpm)
print("Magnetic field=",B_field)

Z_scatt=np.array([1, -1, 0]) / np.sqrt(2) # Scattering Polarization direction
Z_magne=np.array([1, 1, 0]) / np.sqrt(2)

for i in range(len(T_array)):
    T=T_array[i]
    print("Temperature=%f" %T)
    
    c_magnetization=[]
    
    for c in range(len(cluster)-5):
        
        magnetization=0 ### define variable that will contain the magnetization of the cluster
        
        N=N_DICT[cluster[c]] #Size of the system
        
        ST=ST_DICT[cluster[c]] #Site types
        
        NN=NN_PAIRS_DICT[cluster[c]] #Pairs in the system
        
        Positions=R_DICT[cluster[c]] #Positions of the sites in the lattice
        
        basis=spin_basis_1d(L=N,S='1/2',pauli=True)#Basis of the system

        H=Hamiltonian(basis,N,Jzz,Jpm,Jppmm,Jzpm,B_field,NN,ST)
        #Obtaining the eigen values and eigen vectors of the Hamiltonian
        eigenvals,eigenvect=H.eigh()
        
        for s1 in range(N):
            
            ### Generating the coefficient of the quantum operator ###
            coefficient=[[gzz,s1]]#U_B*gz
            
            ### Creating the quantum operator to be averaged ###
            linear_op=quantum_LinearOperator([['z',coefficient]], basis=basis, check_herm=False, check_symm=False)
            
            #Average value of the operator
            Op_T_average=thermal_average(eigenvals,eigenvect,linear_op,T)# Defining the operator to be calculated for a specific couple of sites
            
            ### Adding all the possible configurations due to symmetry transformations ###
            
            # Defining the type site orientation
            Tot_proyection=0
            proyection=abs(np.dot(Z_DIR[ST[s1]],Z_magne)) #Proyection factor to measure magnetization

            magnetization+=(Op_T_average.real*proyection/N)
#            if(c!=0):
#                for sym in range(SYM.shape[0]):
#
#                    z_1=np.dot(SYM[sym],Z_DIR[ST[s1]]) #Defining orientation for the spin configurations
#                    proyection=np.dot(z_1,Z_magne) #Proyection factor to measure magnetization
#                    Tot_proyection+=  proyection
#                    magnetization+=Op_T_average.real*proyection/48. #Divided in all the possible symmetry related configurations
#            print("Total proyection is " , Tot_proyection)
        print("Magnetization %f in cluster %d" %(magnetization,c))
        c_magnetization.append(magnetization/N)# Normalizing in the number of sites of the given cluster

    ### Getting the contribution of each cluster for the magnetization ###
#    print("magnetization per cluster type")
#    print(c_magnetization)
    contributions=cluster_contributions_array(np.array(c_magnetization),min(c,4))
#    print("Contributions for temperature %f" %T)
    print("Contributions")
    print(contributions)
    print("Bare sums")
    print(bare_sum_array(contributions))
    for k in range(len(contributions)):
        Magnetization[k,i]=(sum(contributions[:(k+1)]))/h# usual order k

for i in range(len(contributions)):
    plt.semilogx(T_array,Magnetization[i,:],"o-",label=str(i))

plt.xlabel("$ T\ (K) $")
plt.legend()
plt.xlim(1,10)
plt.title("h=" +str(h))
#plt.ylim(0.2,2.6)
plt.ylim(0.2,2.6)
plt.show()

