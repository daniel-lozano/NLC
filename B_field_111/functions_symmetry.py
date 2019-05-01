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
from sums import *

Jzz=0.17 #First neighbors interaction constant
Jpm=0.05
Jppmm=0.05
Jzpm=-0.14

gz=4.32 #Lande factor
gzz=1.80
gxy=4.32

KB = 0.086173  # Boltzmann constant [meV/K]
U_B = 5.788E-2  # Bohr magneton [meV/T]


Eta=[[0,-1.,np.exp(np.pi*1j/3.),np.exp(-np.pi*1j/3.)],
     [-1.,0,np.exp(-np.pi*1j/3.),np.exp(np.pi*1j/3.)],
     [np.exp(np.pi*1j/3.),np.exp(-np.pi*1j/3.),0,-1.],
     [np.exp(-np.pi*1j/3),np.exp(np.pi*1j/3.),-1,0]]

Gamma=-np.conj(Eta)

def Hamiltonian(basis,N,Jzz,Jpm,Jppmm,Jzpm,B_field,h_field,NN,ST,coef):
    
    
    J_zz=[[Jzz,i,j] for (i,j) in NN]
    J_pm=[[-Jpm,i,j] for (i,j) in NN]
    
    #    Double exitation terms
    J_pp=[[Jppmm*Gamma[ST[i],ST[j]],i,j] for (i,j) in NN]
    J_mm=[[Jppmm*np.conj(Gamma[ST[i],ST[j]]),i,j] for (i,j) in NN]
    
    #    Exitation plus magnetization term
    
    J_zpj=[[Jzpm*Eta[ST[i]][ST[j]],i,j] for (i,j) in NN]
    J_zmj=[[Jzpm*np.conj(Eta[ST[i]][ST[j]]),i,j] for (i,j) in NN]
    
    J_zpi=[[Jzpm*Eta[ST[j]][ST[i]],j,i] for (i,j) in NN]
    J_zmi=[[Jzpm*np.conj(Eta[ST[j]][ST[i]]),j,i] for (i,j) in NN]
    
    ### multiplying the constants too the coefficients ###
    
    COEF=coef.copy()
    for i in range(len(COEF)):
        COEF[i][0]*=-U_B*gz*h_field

    z_field=COEF#[[-(U_B * gz)*np.dot(B_field,Z_DIR[ST[i]]),i] for i in range(N)] #Carefull with the directions!!! it should be Z\dotB

    #Time independe parameters of the Hamiltonian
    static=[["zz",J_zz],
            ["z",z_field],
            ["+-",J_pm],["-+",J_pm], #Exchange term
            ["++",J_pp],["--",J_mm], #Double flip term
            ["z+",J_zpj],["z+",J_zpi],["z-",J_zmj],["z-",J_zmi]] #Non krammers term

    #Time dependent parameters of the Hamiltonian
    dynamic=[]
    
    H=hamiltonian(static,dynamic,basis=basis,check_herm=False, check_symm=False)#dtype=np.float64,
    return H

def thermal_average(eigenvals,eigenvect,linear_op,Temp):
    
    average=0 #Average of the linear operator
    Z=0 #Partition function
    eigenvals-=np.ones(len(eigenvals))*min(eigenvals)#Set minimum eigen value at 0
    for val,vect in zip(eigenvals,eigenvect.T):
        average+=np.dot(vect.conj(),linear_op.matvec(vect))*np.exp(-val/(KB*Temp))
        Z+=np.exp(-val/(KB*Temp))
    return average/Z

def specific_heat(eigenvals,eigenvect,Temp):
    
    average=0
    average_2=0 #Average of the linear operator
    Z=0 #Partition function
    eigenvals-=np.ones(len(eigenvals))*min(eigenvals)#Set minimum eigen value at 0
    for val,vect in zip(eigenvals,eigenvect.T):
        average+=(val)*np.exp(-val/(KB*Temp))
        average_2+=(val**2)*np.exp(-val/(KB*Temp))
        Z+=np.exp(-val/(KB*Temp))
    
    c_h=(average_2-average**2)/(Z*KB*Temp**2)*6.022E23 * 1.602E-22
    
    
    return c_h



