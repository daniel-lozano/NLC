import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

### Import all the info for the clusters
from constants import *

### Line of symmetry ###
B=np.array([1,1,1])/np.sqrt(3)

def plot_sym_line():
    x=np.linspace(-1,1,10)
    y=x.copy()
    z=x.copy()

    ax.plot(x,y,z,"k--")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.legend()

def find_symm(B,SYM):
    index_sym=[]
    index_no_sym=[]
    
    ### if there is no symmetry axis then all the symmetries are returned ###
    if(np.dot(B,B)==0):
        return np.arange(np.shape(SYM)[0])

    for i in range(np.shape(SYM)[0]):
        transform=np.dot(SYM[i],B)
        norm_t=np.dot(transform,B)
        if(round((norm_t),6)==1):
            index_sym.append(i)
        else:
            index_no_sym.append(i)
    return index_sym,index_no_sym


### Clusters to be called ###
answer=input("Check symmetries with \hat{B}=(111)?: (y) or (n) ")
if(answer=="y"):

    clusters=['0','1','2','3']
    c=2
    r_sites=R_DICT[clusters[c]]

    fig=plt.figure()
    ax=fig.gca(projection='3d')
    plot_sym_line()

    ax.plot(r_sites[:,0],r_sites[:,1],r_sites[:,2],"o-")
    plt.show()

    ### Indexes where the symmetry is conserved and where is not ###
    index_sym,index_no_sym=find_symm(B,SYM)

    print("Number of symmetries that are conserved is %d out of %d" %(len(index_sym),np.shape(SYM)[0]) )
    print(index_sym)

    for i in index_sym:
        r_transform=np.dot(SYM[i],r_sites.T)
        ### Plotting ###
        fig=plt.figure()
        ax=fig.gca(projection='3d')
        ax.plot(r_transform.T[:,0],r_transform.T[:,1],r_transform.T[:,2],"o-",label="cluster "+ clusters[c])
        plot_sym_line()

        ax.set_title("SYM[%d]" %i)# +str(i)+"]")
        plt.show()
    #    print(r_transform.T)

    ### plotting new cluster ###

    fig=plt.figure()
    ax=fig.gca(projection='3d')
    plot_sym_line()
    r_transform=np.dot(SYM[index_no_sym[0]],r_sites.T)
    ax.plot(r_transform.T[:,0],r_transform.T[:,1],r_transform.T[:,2],"o-",label="cluster "+ clusters[c])
    plt.show()







