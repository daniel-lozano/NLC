
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


print("----%----------%----------%------STARTING CALCULUS----%----------%----------%----------%------")


Jzz=1#0.17 #First neighbors interaction constant
Jpm=0#0.05
Jppmm=0#0.05
Jzpm=0#-0.14
B_field=[0.0,0.0,0.0] # field
T=10.0/KB # Kelvin
gz=4.32 #Lande factor

Eta=[[0,-1,np.exp(np.pi*1j/3),np.exp(-np.pi*1j/3)],
     [-1,0,np.exp(-np.pi*1j/3),np.exp(np.pi*1j/3)],
     [np.exp(np.pi*1j/3),np.exp(-np.pi*1j/3),0,-1],
     [np.exp(-np.pi*1j/3),np.exp(-np.pi*1j/3),-1,0]]

Gamma=-np.conj(Eta)


size=int(input("Enter number of points for the scattering structure factor (int): "))
q=2*np.pi*np.arange(-2.501, 2.501, 0.1)#0.025
ql=2*np.pi*np.linspace(0.01, 5.01, size)
qh=2*np.pi*np.linspace(0.01, 2.501, size)
cluster=['0','1','2','3','4Y','4I','4L']

#Arrays to store the different scattering results
c_SF_intensity = np.zeros((len(cluster), ql.size, qh.size))
c_NSF_intensity = np.zeros((len(cluster), ql.size, qh.size))
print(c_NSF_intensity[0][0,0])


print("First Neighbor interaction constant=",Jzz)
print("First Neighbor exchange=",Jpm)
print("Magnetic field=",B_field)



def gen_hamiltonian(basis,N,Jzz,Jpm,B_field,NN):

    print("Number of sites=",N)


    J_zz=[[Jzz,i,j] for (i,j) in NN]
    J_pm=[[-Jpm,i,j] for (i,j) in NN]
    
    #    Double exitation terms
    J_pp=[[Jppmm*Gamma[ST[i],ST[j]],i,j] for (i,j) in NN]
    J_mm=[[Jppmm*np.conj(Gamma[ST[i],ST[j]]),i,j] for (i,j) in NN]
    
    #    Exitation plus magnetization term
    
    J_zpj=[[Jzpm*Eta[ST[i]][ST[j]],i,j] for (i,j) in NN]
    J_zpi=[[Jzpm*Eta[ST[j]][ST[i]],j,i] for (i,j) in NN]
    
    J_zmj=[[Jzpm*np.conj(Eta[ST[i]][ST[j]]),i,j] for (i,j) in NN]
    J_zmi=[[Jzpm*np.conj(Eta[ST[j]][ST[i]]),j,i] for (i,j) in NN]

    
    
    z_field=[[-(U_B * gz)*np.matmul(B_field,Z_DIR[ST[i]]),i] for i in range(N)] #Carefull with the directions!!! it should be Z\dotB

    #Time independe parameters of the Hamiltonian
    static=[["zz",J_zz],["z",z_field],["+-",J_pm],["-+",J_pm],["++",J_pp],["--",J_mm],["z+",J_zpj],["z+",J_zpi],["z-",J_zmj],["z-",J_zmi]]
    
    #Time dependent parameters of the Hamiltonian
    dynamic=[]


    H=hamiltonian(static,dynamic,basis=basis,check_herm=False, check_symm=False)#dtype=np.float64,
    return H

def get_thermal_average(eigenvals,eigenvect,linear_op,Temp):

    average=0 #Average of the linear operator
    Z=0 #Partition function
    eigenvals-=np.ones(len(eigenvals))*min(eigenvals)
    
    for val,vect in zip(eigenvals,eigenvect.T):
        average+=np.dot(vect.conj(),linear_op.matvec(vect))*np.exp(-val/(KB*Temp))
        Z+=np.exp(-val/(KB*Temp))
    
    return average/Z



#Defining cluster to be use

for c in range(len(cluster)):

    q_time = time()
    
    print("Initialicing calculus for cluster type ", cluster[c])
    
    N=N_DICT[cluster[c]] #Size of the system
    
    NN=NN_PAIRS_DICT[cluster[c]] #Pairs in the system
    
    ST=ST_DICT[cluster[c]] #Site types
    
    Positions=R_DICT[cluster[c]] #Positions of the sites in the lattice
    
    print("N=",N)
#    print "Nearest neighbors", NN
#    print "Site types",ST
#    print "Positions",Positions

    Z_scatt=np.array([1, -1, 0]) / np.sqrt(2) # Scattering Polarization direction

    basis=spin_basis_1d(L=N,S='1/2',pauli=True)#Basis of the system
#    print("Basis of the system")
#    print(basis)
    H=gen_hamiltonian(basis,N,Jzz,Jpm,B_field,NN)
    
    #Obtaining the eigen values and eigen vectors of the Hamiltonian
    eigenvals,eigenvect=H.eigh()
    if(Jpm==0):
        for i in eigenvect:
            if(list(i).count(1)+list(i).count(-1)!=1):
                print("There's a non classical contribution!!!")
        #Finding all the correlation coefficients
    
    for s1,s2 in it.product(range(N),range(N)):
        
       
        #Generating the coefficient of the quantum operator
        coefficient=[[(U_B*gz)**2,s1,s2]] #   [[value,site1,site2 ]]
        linear_op=quantum_LinearOperator([['zz',coefficient]], basis=basis, check_herm=False, check_symm=False)
        
        #Average value of the operator
        Op_T_average=get_thermal_average(eigenvals,eigenvect,linear_op,T)
        
       
        

    #At this point the prefactor due to the direction and the exponential factor must be multiply by the average factor previously found to have the total scattering function, the different type of cluster must be included aswell

        for sym in range(SYM.shape[0]):
            
            #Relative position vector and Z directions for both spin in the correlation factor
        
            r_ij=np.dot(SYM[sym],Positions[s1]-Positions[s2])
            
            z_1=np.dot(SYM[sym],Z_DIR[ST[s1]])
            z_2=np.dot(SYM[sym],Z_DIR[ST[s2]])
            
            dir_s1=np.cross(z_1, Z_scatt)
            dir_s2=np.cross(z_2, Z_scatt)

            Projection_factor_NSF=np.dot(z_1,Z_scatt)*np.dot(z_2,Z_scatt)
           
            for h,l in it.product(range(qh.size),range(ql.size)):#
    
                q_vector=np.array([qh[h],qh[h],ql[l]])
                
                X_scatt=np.array([qh[h],qh[h],ql[l]])/np.sqrt(np.dot(q_vector,q_vector)) #np.sqrt(2.*q[h]**2+q[l]**2)
                
                Projection_factor_SF=np.dot(dir_s1,X_scatt)*np.dot(dir_s2,X_scatt)

                c_NSF_intensity[c][l,h]+=Op_T_average.real*np.cos(q_vector.dot(r_ij))*Projection_factor_NSF/48.

                c_SF_intensity[c][l,h]+=Op_T_average.real*np.cos(q_vector.dot(r_ij))*Projection_factor_SF/48.

    print("time used for cluster"+ str(c)+ "=", timedelta(seconds=time()-q_time))
    print("")



np.savez_compressed("data_T"+str(round(T,2))+"_J"+str(Jzz)+"_Jpm"+str(Jpm),SF=c_SF_intensity, NSF=c_NSF_intensity, QH=qh, QL=ql)



