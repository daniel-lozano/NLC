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
from functions_symmetry import *
from sums import *
from visualice import find_symm


print("----%----------%----------%------STARTING CALCULUS----%----------%----------%----------%------")


h_max=float(input("Enter max value for the magnetic field="))# Field strength, should be in Teslas
H_field=np.linspace(0.1,h_max,100)
dh=abs(H_field[1]-H_field[0])

ans=input("Plot cluster (y) or (n): ")
B_direct=np.array([1.0,1.0,1.0])/np.sqrt(3)

cluster=['0','1','2','3']#,'4Y','4I','4L']

#Arrays to store the different scattering results

T_array=np.linspace(1, 2,1)#np.linspace(1,10,20)
Magnetization=np.zeros((len(cluster),len(T_array)))

print("\nFirst Neighbor interaction constant(Jzz)=",Jzz)
print("\nFirst Neighbor exchange(Jpm)=",Jpm)
print("\nDouble excitation(Jppmm)=",Jppmm)
print("\nMixed(Jzpm)=",Jzpm)
print("\nMagnetic field direction=",B_direct)
#####################################################################
########### Cluster type with symmetry axis (111) ###################
#####################################################################
W_values={'0' : [1,0,0,0,0,0],
          '1' : [-4,1,0,0,0,0],
          '2p': [1,-2,0,1,0,0],
          '2' : [1,-2,1,0,0,0],
          '3p': [0,1,0,-2,0,1],
          '3' : [0,1,-1,-1,1,0]
          }

###   L_multiplicity=[ L(0),L(1),L(2),L(2p),L(3),L(3p)] ###
L_multiplicity=np.array([1.,1./2,1./7,6./7,6./5,9./5])

cluster_field=['0','1','2','2p','3','3p']

####################################################################
######## Symmetry axis set by the magnetic field applied ###########
####################################################################

Z_magne=np.array([1, 1, 1]) / np.sqrt(3)

### Perpendicular axis to the symmetry axis ###
X_magne=np.array([1, -1, 0]) / np.sqrt(2)
Y_magne=np.array([-1, 1, 0]) / np.sqrt(2)

########################################################################
###### transformations that keep the symmetry and that don't ###########
########################################################################
index_sym,index_no_sym=find_symm(B_direct,SYM)

Energy_h=np.zeros(len(H_field))
for i in range(len(H_field)):
    
    h=H_field[i]
    B_field=B_direct*h # Field direction
    T=T_array[0]
    if(i%4==0):
        print("H_field=%f @ T=%f" %(h,T))
        print("B_field=",B_field)
    
    
    types_of_c=[]

    Energy_cluster=[]

    for c in range(len(cluster)):
        
        N=N_DICT[cluster[c]] ### Size of the system
        ST=ST_DICT[cluster[c]] ### Site types
        NN=NN_PAIRS_DICT[cluster[c]] ### Pairs in the system
        Positions=R_DICT[cluster[c]] ### Positions of the sites in the lattice
        
     
        ### Symmetry used to get other type of cluster ###
        sym_to_use=0
        if(c==3):
            sym_to_use=1

        ### Arrays for the projections of the easy axis directions in the magnetization directions ###
        projection_para_z=[]
        projection_para_x=[]
        projection_para_y=[]

        ### Coefficients for the directions of the magnetization ###
        coef_para_z=[]
        coef_para_x=[]
        coef_para_y=[]

        ### ndArray for the direction of spins ###
        direc_para=np.zeros((3,N),dtype=float)

        if(c>1):
            projection_perp_z=[]
            projection_perp_x=[]
            projection_perp_y=[]
            
            coef_perp_z=[]
            coef_perp_x=[]
            coef_perp_y=[]
            
            direc_perp=np.zeros((3,N))

        ### Getting proyection factor for cluster c for every site ###

        for sp in range(N):
            
            ### Direction of the easy spin axis for the site sp ###
            z_par=Z_DIR[ST[sp]]
            ### Projections of the axis in the different directions ###
            projection_para_z.append(np.dot(Z_magne,z_par))
            projection_para_x.append(np.dot(X_magne,z_par))
            projection_para_y.append(np.dot(Y_magne,z_par))
            
            ### Coefficient for the magnetization in each direction
            coef_para_z.append([projection_para_z[-1],sp])
            coef_para_x.append([projection_para_x[-1],sp])
            coef_para_y.append([projection_para_y[-1],sp])
            
            direc_para[0][sp]=z_par[0]
            direc_para[1][sp]=z_par[1]
            direc_para[2][sp]=z_par[2]
            
            
            if(c>1):
                z_per=np.dot(SYM[index_no_sym[sym_to_use]],z_par) ### Another cluster type accesed
                
               
                projection_perp_z.append(np.dot(Z_magne,z_per))
                projection_perp_x.append(np.dot(X_magne,z_per))
                projection_perp_y.append(np.dot(Y_magne,z_per))
                
                coef_perp_z.append([abs(projection_perp_z[-1]),sp])
                coef_perp_x.append([abs(projection_perp_x[-1]),sp])
                coef_perp_y.append([abs(projection_perp_y[-1]),sp])
                
                direc_perp[0][sp]=z_per[0]
                direc_perp[1][sp]=z_per[1]
                direc_perp[2][sp]=z_per[2]

        ### Basis of the system and Hamiltonian of the system ###
        basis=spin_basis_1d(L=N,S='1/2',pauli=True)
        ### NOTE THAT THE COEFFICIENT FOR THE B FIELD ARE DIFFERENTE NOW! ###
        H_para=Hamiltonian(basis,N,Jzz,0,0,0,B_field,h,NN,ST,coef_para_z)### Constant imported from file
#        H_para=Hamiltonian(basis,N,Jzz,Jpm,Jppmm,Jzpm,B_field,h,NN,ST,coef_para_z)### Constant imported from file
    
        eigenvals,eigenvect=H_para.eigh()

        E_para=average_energy(eigenvals,eigenvect,T)
        Energy_cluster.append(E_para)

        types_of_c.append(str(c))


        if(c>1):
            
            H_perp=Hamiltonian(basis,N,Jzz,0,0,0,B_field,h,NN,ST,coef_perp_z)### Constant imported from file
            #            H=Hamiltonian(basis,N,Jzz,Jzz,Jpm,Jppmm,B_field,h,NN,ST,coef_perp_z)
            eigenvals_perp,eigenvect_perp=H_perp.eigh()

            E_perp=average_energy(eigenvals_perp,eigenvect_perp,T)
            Energy_cluster.append(E_perp)


        ### Plotting clusters only one when i=0 ###
        if(ans=="y" and i==0):
            
            Positions_transform=np.dot(SYM[index_no_sym[sym_to_use]],Positions.T).T

            fig=plt.figure()

            ax=fig.gca(projection='3d')
            ax.plot(Positions[:,0],Positions[:,1],Positions[:,2],"o-",label="cluster "+ clusters[c])
            ax.quiver(Positions[:,0],Positions[:,1],Positions[:,2],direc_para[0,:],direc_para[1,:],direc_para[2,:],length=0.1,color="r")
            ax.plot(Positions_transform[:,0],Positions_transform[:,1],Positions_transform[:,2],"o-",label="cluster perp "+ clusters[c])
            if(c>1): ax.quiver(Positions_transform[:,0],Positions_transform[:,1],Positions_transform[:,2],direc_perp[0,:],direc_perp[1,:],direc_perp[2,:],length=0.1,color="g")



            x=np.linspace(-1,1,10)
            y=x.copy()
            z=x.copy()
            ax.plot(x,y,z,"k--",linewidth=1)
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            plt.legend()
            plt.show()
        if(i==0):
            print(types_of_c)

    Energy_cluster=np.array(Energy_cluster)


    Weights=np.zeros(len(cluster_field))


    for w in range(len(cluster_field)):
        t=cluster_field[w]
        factors=np.array(W_values[t])
        
        Weights[w]=np.dot(factors,Energy_cluster)

    if(i==0):
        print(types_of_c)
        print("Energy of clusters")
        print(Energy_cluster)
        print("Weights vector")
        print(Weights)
        print("Factors")
        print(W_values)

    Energy_h[i]=np.dot(Weights,L_multiplicity)
#    print("--L(c)-- multiplicity")
#    print(L_multiplicity)
#    print("E_K=3")
#    print(Energy_h)
#    input("PARAR PARA REVISAR!")

plt.plot(H_field,Energy_h,label="Energy")

plt.plot(H_field[1:], (Energy_h[1:]-Energy_h[:len(H_field)-1])/dh, label="Derivative")
plt.ylabel("Energy")
plt.xlabel("H_field")
plt.legend()
plt.show()




