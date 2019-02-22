
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


Jzz=1. #First neighbors interaction constant
Jpm=0.0
B_field=[0.0,0.0,0.0] # field
T=0.1/KB # Kelvin
gz=4.32 #Lande factor

q=2*np.pi*np.arange(-4.01, 4.01, 0.025)#np.arange(-2.501, 2.501, 0.1)

cluster=['0','1','2','3','4Y','4I','4L']

#Arrays to store the different scattering results
c_SF_intensity = np.zeros((len(cluster), q.size, q.size))
c_NSF_intensity = np.zeros((len(cluster), q.size, q.size))
print(c_NSF_intensity[0][0,0])


print("First Neighbor interaction constant=",Jzz)
print("First Neighbor exchange=",Jpm)
print("Magnetic field=",B_field)



def gen_hamiltonian(basis,N,Jzz,Jpm,B_field,NN):

    print("Number of sites=",N)


    J_zz=[[Jzz,i,j] for (i,j) in NN]
    J_pm=[[-Jpm,i,j] for (i,j) in NN]
    z_field=[[(U_B * gz)*np.matmul(B_field,Z_DIR[ST[i]]),i] for i in range(N)] #Carefull with the directions!!! it should be Z\dotB

    #Time independe parameters of the Hamiltonian
    static=[["zz",J_zz],["z",z_field],["+-",J_pm],["-+",J_pm]]
    
    #Time dependent parameters of the Hamiltonian
    dynamic=[]


    H=hamiltonian(static,dynamic,dtype=np.float64,basis=basis,check_herm=False, check_symm=False)
    return H

def get_thermal_average(eigenvals,eigenvect,linear_op,Temp,aiao_OP,s1,s2):

    average=0 #Average of the linear operator
    Z=0 #Partition function
    eigenvals-=np.ones(len(eigenvals))*min(eigenvals)
    contributing=0
    sites=[]
    
    for i in range(N):
        sites.append(i)
    sites=np.array(sites)
    for val,vect in zip(eigenvals,eigenvect.T):
        use=1
        
        
        for operator in aiao_OP: #Get states in the all in all out state
                
            coefficient=[[operator[i],i] for i in range(N)]

            linear_op=quantum_LinearOperator([['z',coefficient]], basis=basis, check_herm=False, check_symm=False)
                
            val=np.dot(vect.conj(),linear_op.matvec(vect))
                
            if(abs(val)==4):
                use=0
                break

        
        if(use==1):
            contributing+=1
            average+=np.dot(vect.conj(),linear_op.matvec(vect))*np.exp(-val/(KB*Temp))
            Z+=np.exp(-val/(KB*Temp))
    if(s1==0 and s2==0):
        print(str(contributing)+" out of " + str(len(eigenvals)) )
    return average/Z



def find_aiao(c,N,basis): #Find the configurations that have an all-in or an all-out state
    
    OPERATORS=[] #Stores the sites to get the all in all out state
    L=2**N
    print("L=",L)
#    print "Nearest neighbors", NN
#    print "Site types",ST
#    print "Positions",Positions
#    
    if(c>0):
        min_num=min(c,4)
        
        print("finding positions")
        for t in range(min_num):
            sites_of_z=np.zeros(N)
            
            for index in range(4):
                sites_of_z[index+3*t]=1
#            if(t>0):
#                sites_of_z[3*t]=-1
            OPERATORS.append(sites_of_z)
    
        indices=[]
        print(OPERATORS)
       
        return np.array(OPERATORS)
    return []






#Defining cluster to be use

for c in range(len(cluster)):
    time_c=time()
    print("")
    print("Initialicing calculus for cluster type ", cluster[c])
    
    N=N_DICT[cluster[c]] #Size of the system
    
    NN=NN_PAIRS_DICT[cluster[c]] #Pairs in the system
    
    ST=ST_DICT[cluster[c]] #Site types
    
    Positions=R_DICT[cluster[c]] #Positions of the sites in the lattice
    
    print("N=",N)
    
    OPERATORS=[] #Stores the sites to get the all in all out state
    L=2**N
    print("L=",L)
#    print "Nearest neighbors", NN
#    print "Site types",ST
#    print "Positions",Positions

    Z_scatt=np.array([1, -1, 0]) / np.sqrt(2) # Scattering Polarization direction

    basis=spin_basis_1d(L=N,S='1/2',pauli=True)#Basis of the system
    
    aiao_OP=find_aiao(c,N,basis)
    

    H=gen_hamiltonian(basis,N,Jzz,Jpm,B_field,NN)
    
    #Obtaining the eigen values and eigen vectors of the Hamiltonian
    eigenvals,eigenvect=H.eigh()
#    print(eigenvals)

    

    for i in eigenvect:
        if(list(i).count(1)+list(i).count(-1)!=1):
            print("%%%%%%%%%%%%%%%%%%%%%%There's a non classical contribution!!!%%%%%%%%%%%%%%%%%%%%%%")
    #Finding all the correlation coefficients
    
    for s1,s2 in it.product(range(N),range(N)):

        #Generating the coefficient of the quantum operator
        coefficient=[[(U_B*gz)**2,s1,s2]]#   [[value,site1,site2 ]]
        linear_op=quantum_LinearOperator([['zz',coefficient]], basis=basis, check_herm=False, check_symm=False)
        
        #Average value of the operator
        #Op_T_average=get_thermal_average(eigenvals,eigenvect,linear_op,T) Old line
        Op_T_average=get_thermal_average(eigenvals,eigenvect,linear_op,T,aiao_OP,s1,s2)
        

    #At this point the prefactor due to the direction and the exponential factor must be multiply by the average factor previously found to have the total scattering function, the different type of cluster must be included aswell

        for sym in range(SYM.shape[0]):
            
            #Relative position vector and Z directions for both spin in the correlation factor
        
            r_ij=np.matmul(SYM[sym],Positions[s1]-Positions[s2])
            
            z_1=np.matmul(SYM[sym],Z_DIR[ST[s1]])
            z_2=np.matmul(SYM[sym],Z_DIR[ST[s2]])
            
            dir_s1=np.cross(z_1, Z_scatt)
            dir_s2=np.cross(z_2, Z_scatt)

            Projection_factor_NSF=np.dot(z_1,Z_scatt)*np.dot(z_2,Z_scatt)
           


            for h,l in it.product(range(q.size),range(q.size)):#
    
                q_vector=np.array([q[h],q[h],q[l]])
                
                X_scatt=np.array([q[h],q[h],q[l]])/np.sqrt(2.*q[h]**2+q[l]**2)
                
                Projection_factor_SF=np.dot(dir_s1,X_scatt)*np.dot(dir_s2,X_scatt)

                c_NSF_intensity[c][l,h]+=Op_T_average.real*np.cos(q_vector.dot(r_ij))*Projection_factor_NSF/48.

                c_SF_intensity[c][l,h]+=Op_T_average.real*np.cos(q_vector.dot(r_ij))*Projection_factor_SF/48.


    print("Time used for cluster "+ str(c) +"="+str(time()-time_c))


np.savez_compressed("NLC1_no_aiao_data_T"+str(round(T,2))+"_J"+str(Jzz)+"_Jpm"+str(Jpm),SF=c_SF_intensity, NSF=c_NSF_intensity)



