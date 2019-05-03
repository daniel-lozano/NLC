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


h=float(input("Enter a value for the magnetic field="))# Field strength, should be in Teslas
B_direct=np.array([1.0,0.0,0.0])
B_field=B_direct*h # Field direction

cluster=['0','1','2','3','4Y','4I','4L']

#Arrays to store the different scattering results

T_array=np.sort(np.random.lognormal(0.5, 1,20))#np.linspace(1,10,20)
Magnetization=np.zeros((len(cluster),len(T_array)))

print("First Neighbor interaction constant(Jzz)=",Jzz)
print("First Neighbor exchange(Jpm)=",Jpm)
print("Double excitation(Jppmm)=",Jppmm)
print("Mixed(Jzpm)=",Jzpm)
print("Magnetic field=",B_field)


Z_magne=np.array([1, 1, 1]) / np.sqrt(3)

for i in range(len(T_array)):
    T=T_array[i]
    if(i%4==0):
        print("Temperature=%f" %T)
    
    c_magnetization=[]
    
    for c in range(len(cluster)-5):
        
        magnetization=0 ### define variable that will contain the magnetization of the cluster
        
        N=N_DICT[cluster[c]] ### Size of the system
        
        ST=ST_DICT[cluster[c]] ### Site types
        
        NN=NN_PAIRS_DICT[cluster[c]] ### Pairs in the system
        
        Positions=R_DICT[cluster[c]] ### Positions of the sites in the lattice
        
        basis=spin_basis_1d(L=N,S='1/2',pauli=True) ### Basis of the system

        H=Hamiltonian(basis,N,Jzz,Jpm,Jppmm,Jzpm,B_field,NN,ST)
        ### Obtaining the eigen values and eigen vectors of the Hamiltonian
        eigenvals,eigenvect=H.eigh()
        
        ### Building coefficients for the matrix ###
        coefficient=[]
        num_of_sym=0
        
        for sym in (SYM.shape[0]):
            B_rot=np.dot(SYM[sym],B_direct)
            dot_product=abs(np.dot(B_rot,B_direct))
            if(dot_product==1):
                num_of_sym+=1
        
        
        print("There are %d symmetries that fullfil the condiction" %num_of_sym)
        for s1 in range(N):
            
            Z_dir=Z_DIR[ST[s1]]
            ### Generating the coefficient of the quantum operator ###
            proyection=abs(np.dot(Z_DIR[ST[s1]],Z_magne)) ###proyecting onto the easy z axis define by the magnetic field
            coefficient.append([gzz*abs(np.dot(Z_DIR[ST[s1]],Z_magne)),s1])#[[gzz,s1]]
            
        ### Creating the quantum operator to be averaged ###
        linear_op=quantum_LinearOperator([['z',coefficient]], basis=basis, check_herm=False, check_symm=False)
        
        #Average value of the operator
        Op_T_average=thermal_average(eigenvals,eigenvect,linear_op,T)# Defining the operator to be calculated for a specific couple of sites
        
        #Proyection factor to measure magnetization

        magnetization+=(Op_T_average.real)# Normalizing in the number of sites of the given cluster
    
        c_magnetization.append(magnetization)

    ### Getting the contribution of each cluster for the magnetization ###
#    print("Magnetization per cluster type")
#    print(c_magnetization)
    contributions=cluster_contributions_array(np.array(c_magnetization),min(c,4))

#    print("Contributions")
#    print(contributions)

    for k in range(len(contributions)):
        Magnetization[k,i]=(sum(contributions[:(k+1)]))/h# usual order k

for i in range(len(contributions)):
    plt.semilogx(T_array,abs(Magnetization[i,:]),"-",label=str(i))

plt.xlabel("$ T\ (K) $")
plt.legend()
plt.xlim(1,10)
plt.title("h=" +str(h))
#plt.ylim(0.2,2.6)
#plt.ylim(0.2,2.6)
plt.show()

