
import itertools as it
import logging as log
import os
from datetime import timedelta
from time import time

from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian, quantum_LinearOperator
from scipy.stats import cauchy

from constants import *


print("----%----------%----------%------STARTING CALCULUS----%----------%----------%----------%------")

Jzz=1. #First neighbors interaction constant
Jpm=0.0
h=0.0 # field
T=10 # Kelvin
gz=4.32 #Lande factor

print("First Neighbor interaction constant=",Jzz)
print("First Neighbor exchange=",Jpm)
print("Magnetic field=",h)

N=2 #Size of the system

basis=spin_basis_1d(L=N,S='1/2',pauli=True)#Basis of the system

def gen_hamiltonian(basis,N,Jzz,Jpm,h):

    print("Number of sites=",N)
    print("Printing basis of the system")
    print(basis)

    J_zz=[[Jzz,i,(i+1)%N] for i in range(N)]
    J_pm=[[-Jpm,i,(i+1)%N] for i in range(N)]
    z_field=[[(U_B * gz)*h,i] for i in range(N)] #Carefull with the directions!!! it should be Z\dotB

    #Time independe parameters of the Hamiltonian
    static=[["zz",J_zz],["z",z_field],["+-",J_pm],["-+",J_pm]]
    
    #Time dependent parameters of the Hamiltonian
    dynamic=[]

    print(static)

    H=hamiltonian(static,dynamic,dtype=np.float64,basis=basis,check_herm=False, check_symm=False)
    return H

def get_thermal_average(eigenvals,eigenvect,linear_op,Temp):

    average=0 #Average of the linear operator
    Z=0 #Partition function
    
    for val,vect in zip(eigenvals,eigenvect.T):
        average+=np.dot(vect.conj(),linear_op.matvec(vect))*np.exp(-val/(KB*Temp))
        Z+=np.exp(-val/(KB*Temp))
    
    return average/Z





H=gen_hamiltonian(basis,N,Jzz,Jpm,h)
#Obtaining the eigen values and eigen vectors of the Hamiltonian
eigenvals,eigenvect=H.eigh()

#Generating the coefficient of the quantum operator

coefficient=[[(U_B*KB)**2,0,0]]#   [[value,site1,site2 ]]

linear_op=quantum_LinearOperator([['zz',coefficient]], basis=basis, check_herm=False, check_symm=False)
Op_T_average=get_thermal_average(eigenvals,eigenvect,linear_op,T)

print(Op_T_average)



print("Hamiltonian")
print(H.toarray())

print("Eigen vectors")
size=len(H.eigh()[0])

for i in range(size):
    print H.eigh()[0][i],H.eigh()[1].T[i]




