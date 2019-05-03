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
from visualice import find_symm


print("----%----------%----------%------STARTING CALCULUS----%----------%----------%----------%------")


h=float(input("Enter a value for the magnetic field="))# Field strength, should be in Teslas
B_direct=np.array([1.0,1.0,1.0])/np.sqrt(3)
B_field=B_direct*h # Field direction

cluster=['0','1','2','3']#,'4Y','4I','4L']

#Arrays to store the different scattering results

T_array=np.linspace(0.5, 1,1)#np.linspace(1,10,20)
Magnetization=np.zeros((len(cluster),len(T_array)))

print("First Neighbor interaction constant(Jzz)=",Jzz)
print("First Neighbor exchange(Jpm)=",Jpm)
print("Double excitation(Jppmm)=",Jppmm)
print("Mixed(Jzpm)=",Jzpm)
print("Magnetic field=",B_field)

### Cluster type with symmetry axis (111) ###

W_values={'0' : [1,0,0,0,0,0],
          '1' : [-4,1,0,0,0,0],
          '2p': [1,-2,0,1,0,0],
          '2' : [1,-2,1,0,0,0],
          '3p': [0,1,0,-2,0,1],
          '3' : [0,1,-1,-1,1,0]
          }

###   L_multiplicity=[ L(0),L(1),L(2),L(2p),L(3),L(3p)] ###
L_multiplicity=np.array([1.,1./2,1./7,6./7,6./5,9./5])

for i in ['0','1','2','2p','3','3p']:
    print(i,W_values[i])

### Symmetry axis set by the magnetic field applied ###
Z_magne=np.array([1, 1, 1]) / np.sqrt(3)

### Perpendicular axis to the symmetry axis ###
X_magne=np.array([1, -1, 0]) / np.sqrt(2)
Y_magne=np.array([-1, 1, 0]) / np.sqrt(2)

### transformations that keep the symmetry and that don't ###
index_sym,index_no_sym=find_symm(B_field,SYM)

for i in range(len(T_array)):
    
    T=T_array[i]
    
    if(i%4==0):
        print("Temperature=%f" %T)
    
    types_of_c=[]

    magnetization_cluster=[]

    for c in range(len(cluster)):
        
        N=N_DICT[cluster[c]] ### Size of the system
        ST=ST_DICT[cluster[c]] ### Site types
        NN=NN_PAIRS_DICT[cluster[c]] ### Pairs in the system
        Positions=R_DICT[cluster[c]] ### Positions of the sites in the lattice
        
        ### Basis of the system and Hamiltonian of the system ###
        basis=spin_basis_1d(L=N,S='1/2',pauli=True)
        H=Hamiltonian(basis,N,Jzz,Jpm,Jppmm,Jzpm,B_field,NN,ST)### Constant imported from file

        eigenvals,eigenvect=H.eigh()
        
        ### Getting the projection factors and the coefficients for the quantum operator ###
        
        ### Symmetry used to get other type of cluster ###
        sym_to_use=0
        if(c==3):
            sym_to_use=1
        
        projection_para_z=[]
        coef_para_z=[]
        direc_para=np.zeros((3,N),dtype=float)
        if(c>1):
            projection_perp_z=[]
            coef_perp_z=[]
            direc_perp=np.zeros((3,N))

        print("\nGetting proyection factor for cluster %d " %c)


        for sp in range(N):
    
            z_par=Z_DIR[ST[sp]]
            projection_para_z.append(np.dot(Z_magne,z_par))
            coef_para_z.append([projection_para_z[-1],sp])
            
            direc_para[0][sp]=z_par[0]
            direc_para[1][sp]=z_par[1]
            direc_para[2][sp]=z_par[2]
            
            
            if(c>1):
                z_per=np.dot(SYM[index_no_sym[sym_to_use]],z_par)### Another configuration accesed
### Checking is the vector direction is properly transformed
#                for vec in Z_DIR:
#                    if(1-np.dot(z_per,vec)<1E-2):
#                        print("Condition fulfilled")
#                        print(vec)
#                        print(z_per)
                projection_perp_z.append(np.dot(Z_magne,z_per))
                coef_perp_z.append([abs(projection_perp_z[-1]),sp])
                
                direc_perp[0][sp]=z_per[0]
                direc_perp[1][sp]=z_per[1]
                direc_perp[2][sp]=z_per[2]



#        if(c>1):
#            print(coef_perp)
#            print("Perp Sum is %f" %sum(projection_perp))
#            print(coef_para)
#            print("Para Sum is %f" %sum(projection_para))
        print("directions")
        print(direc_para)
        print(direc_para[0,:])

        ### Getting the operator for the magnetization ###
        Sz_para=quantum_LinearOperator([['z',coef_para_z]], basis=basis, check_herm=False, check_symm=False)
        types_of_c.append(str(c))
        Sz_average_para=thermal_average(eigenvals,eigenvect,Sz_para,T)#

        if(Sz_average_para.imag<1E-15):
            magnetization_cluster.append(Sz_average_para.real)
        else:
            print("Error due to imaginary values!")
            break


        if(c>1):
            Sz_perp=quantum_LinearOperator([['z',coef_perp_z]], basis=basis, check_herm=False, check_symm=False)
            types_of_c.append(str(c)+"p")
            Sz_average_perp=thermal_average(eigenvals,eigenvect,Sz_perp,T)
            
            if(Sz_average_perp.imag<1E-15):
                magnetization_cluster.append(Sz_average_perp.real)
            else:
                print("Error due to imaginary values!")
                break


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

        print(types_of_c)
        print(magnetization_cluster)

        


