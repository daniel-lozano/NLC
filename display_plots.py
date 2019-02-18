
import matplotlib.pyplot as plt
from datetime import timedelta
from time import time
from sys import argv

#from quspin.basis import spin_basis_1d
#from quspin.operators import hamiltonian, quantum_LinearOperator
#from scipy.stats import cauchy

from constants import *
from sums import *



print("----%----------%----------%------STARTING CALCULUS----%----------%----------%----------%------")


Jzz=1. #First neighbors interaction constant
Jpm=0.0
B_field=[0.0,0.0,0.0] # field
T=0.1/KB # Kelvin
gz=4.32 #Lande factor

q=2*np.pi*np.arange(-4.01, 4.01, 0.1)#np.arange(-2.501, 2.501, 0.1)#

hh,l=np.meshgrid(q,q)

cluster=['0','1','2','3','4Y','4I','4L']

#Arrays to store the different scattering results
c_SF_intensity = np.zeros((len(cluster), q.size, q.size))
c_NSF_intensity = np.zeros((len(cluster), q.size, q.size))
print c_NSF_intensity[0][0,0]


print("First Neighbor interaction constant=",Jzz)
print("First Neighbor exchange=",Jpm)
print("Magnetic field=",B_field)
name=argv[1]
SF=np.load(name)['SF']
NSF=np.load(name)['NSF']
print SF.shape[0]

ans1=raw_input("Plot cluster constributions? (yes) or (no):")

keyword="viridis"

for i in range(SF.shape[0]):
    
    c_SF_intensity=SF[i]
    c_NSF_intensity=NSF[i]
    Tot=c_SF_intensity+c_NSF_intensity
    
    print "variation",abs(np.amax(c_NSF_intensity)-np.amin(c_NSF_intensity))
    
    if(ans1=="yes"):
    
        plt.figure(figsize=(12,5))
        plt.subplot(131)
        im1=plt.imshow(c_NSF_intensity,cmap=keyword)
        cbar=plt.colorbar(im1,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("NSF,Cluster type="+str(i))

        plt.subplot(132)
        im2=plt.imshow(c_SF_intensity,cmap=keyword)
        cbar=plt.colorbar(im2,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("SF,Cluster type="+str(i))
        
        plt.subplot(133)
        im2=plt.imshow(Tot,cmap=keyword,extent=(min(q)/(2*np.pi), max(q)/(2*np.pi), min(q)/(2*np.pi), max(q)/(2*np.pi)))
        cbar=plt.colorbar(im2,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Total,Cluster type="+str(i))
        plt.show()


Tot_contributions_NSF=cluster_contributions(NSF,4)
Tot_contributions_SF=cluster_contributions(SF,4)

Bare_sums_SF=bare_sum(Tot_contributions_SF)
Bare_sums_NSF=bare_sum(Tot_contributions_NSF)

Euler_sums_SF=euler_sum(Tot_contributions_SF,3,4)
Euler_sums_NSF=euler_sum(Tot_contributions_NSF,3,4)

#print SF.shape
#print Tot_contributions_SF.shape
#print Bare_sums_SF.shape
#print NSF.shape
#print Tot_contributions_NSF.shape
#print Bare_sums_NSF.shape
#
#levels=np.linspace(np.amin(Euler_sums_SF[1]),np.amax(Euler_sums_SF[1]),50)
#im=plt.contour(hh,l,Euler_sums_SF[1],cmap="gist_heat",levels=levels)
#plt.colorbar(im,orientation="vertical")
#plt.title("Bare sum NSF order=" +str(1))
#plt.show()
#

ans2=raw_input("Plot Bare sums? (yes) or (no):")
ans3=raw_input("Plot Euler vs Bare sums? (yes) or (no):")
if (ans2=="yes"):
    
    for i in range(1,Bare_sums_SF.shape[0]):
    
        plt.figure(figsize=(12,5))
        plt.subplot(131)
        im2=plt.imshow(Bare_sums_NSF[i],cmap=keyword,extent=(min(q)/(2*np.pi), max(q)/(2*np.pi), min(q)/(2*np.pi), max(q)/(2*np.pi)))
        cbar=plt.colorbar(im2,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Bare sum NSF order= " +str(i))
    
        levels=np.linspace(np.amin(Bare_sums_SF[i]),np.amax(Bare_sums_SF[i]),200)

        plt.subplot(132)
        im1=plt.imshow(Bare_sums_SF[i],cmap=keyword,extent=(min(q)/(2*np.pi), max(q)/(2*np.pi), min(q)/(2*np.pi), max(q)/(2*np.pi)))
        cbar=plt.colorbar(im1,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Bare sum SF order=" +str(i))
        
        levels=np.linspace(np.amin(Bare_sums_SF[i]+Bare_sums_NSF[i]),np.amax(Bare_sums_SF[i]+Bare_sums_NSF[i]),200)
        
        plt.subplot(133)
        im1=plt.imshow(Bare_sums_SF[i]+Bare_sums_NSF[i],cmap=keyword,extent=(min(q)/(2*np.pi), max(q)/(2*np.pi), min(q)/(2*np.pi), max(q)/(2*np.pi)))
        cbar=plt.colorbar(im1,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Bare sum Total order=" +str(i))

        plt.show()

if (ans3=="yes"):
    for i in range(Euler_sums_SF.shape[0]):
    
        plt.figure(figsize=(7,7))
        plt.subplot(231)
        im1=plt.imshow(Bare_sums_NSF[i-2],cmap=keyword)
        cbar=plt.colorbar(im1,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Bare sum NSF order=" +str(i+3))
    
        plt.subplot(232)
#        levels=np.linspace(np.amin(Bare_sums_SF[i-2]),np.amax(Bare_sums_SF[i-2]),200)
#        im2=plt.contour(hh,l,Bare_sums_SF[i-2],cmap=keyword,levels=levels)
        im2=plt.imshow(Bare_sums_SF[i-2],cmap=keyword,extent=(min(q)/(2*np.pi), max(q)/(2*np.pi), min(q)/(2*np.pi), max(q)/(2*np.pi)))
        cbar=plt.colorbar(im2,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Bare sum SF order= " +str(i+3))
        
        plt.subplot(233)
#        levels=np.linspace(np.amin(Bare_sums_SF[i-2]+Bare_sums_NSF[i-2]),np.amax(Bare_sums_SF[i-2]+Bare_sums_NSF[i-2]),200)
#        im2=plt.contour(hh,l,Bare_sums_SF[i-2]+Bare_sums_NSF[i-2],cmap=keyword,levels=levels)
        im3=plt.imshow(Bare_sums_SF[i-2]+Bare_sums_NSF[i-2],cmap=keyword,extent=(min(q)/(2*np.pi), max(q)/(2*np.pi), min(q)/(2*np.pi), max(q)/(2*np.pi)))
        cbar=plt.colorbar(im3,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Bare sum total order= " +str(i+3))

        plt.subplot(234)
        im1=plt.imshow(Euler_sums_NSF[i],cmap=keyword)
        cbar=plt.colorbar(im1,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Euler sum NSF order=" +str(i+3))
    
        plt.subplot(235)
#        levels=np.linspace(np.amin(Euler_sums_SF[i]),np.amax(Euler_sums_SF[i]),200)
#        im2=plt.contour(hh,l,Euler_sums_SF[i],cmap=keyword,levels=levels)
        im2=plt.imshow(Euler_sums_SF[i],cmap=keyword,extent=(min(q)/(2*np.pi), max(q)/(2*np.pi), min(q)/(2*np.pi), max(q)/(2*np.pi)) )
        cbar=plt.colorbar(im2,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Euler sum SF order= " +str(i+3))
        
        plt.subplot(236)
#        levels=np.linspace(np.amin(Euler_sums_SF[i]+Euler_sums_NSF[i]),np.amax(Euler_sums_SF[i]+Euler_sums_NSF[i]),200)
#        im3=plt.contour(hh,l,Euler_sums_SF[i]+Euler_sums_NSF[i],cmap=keyword,levels=levels)
        im3=plt.imshow(Euler_sums_SF[i]+Euler_sums_NSF[i],cmap=keyword,extent=(min(q)/(2*np.pi), max(q)/(2*np.pi), min(q)/(2*np.pi), max(q)/(2*np.pi)) )
        cbar=plt.colorbar(im3,orientation="horizontal",pad=0.1)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        plt.title("Euler sum Tot order= " +str(i+3))

        plt.savefig("Euler_vs_Bare_sums_order"+str(i+3)+".png")
        plt.show()



from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

ans4=raw_input("Plot 3D plot of SF? (yes) or (no):")

if(ans4=="yes"):
    fig1=plt.figure()
    ax = fig1.gca(projection='3d')


    surf=ax.plot_surface(hh,l,Euler_sums_SF[-1],cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    ax.set_zlim(np.amin(Euler_sums_SF[-1]), np.amax(Euler_sums_SF)*2)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_xlabel("$[h,h,0]$")
    ax.set_ylabel("$[0,0,l]$")

    fig1.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()


ans5=raw_input("Plot 3D plot of NSF? (yes) or (no):")
if(ans5=="yes"):
    fig2=plt.figure()
    ax = fig2.gca(projection='3d')


    surf=ax.plot_surface(hh,l,Euler_sums_NSF[-1],cmap=cm.coolwarm,
                         linewidth=0, antialiased=False)
    ax.set_zlim(np.amin(Euler_sums_NSF[-1]), np.amax(Euler_sums_NSF))
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_xlabel("$[h,h,0]$")
    ax.set_ylabel("$[0,0,l]$")

    fig2.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
#
#N=len(q)
#number=5
#
#for i in range(number):
#    pos=str(number)+"1"+str(i+1)
#    plt.subplot(pos)
#    plt.plot(q,Euler_sums_NSF[-1][:,(5+i*10)%N],"k-",label="hh="+str(q[(5+i*10)%N]))
#
#    plt.legend()
#    plt.ylim(np.amin(Euler_sums_NSF[-1]),np.amax(Euler_sums_NSF[-1]))
#
#
#plt.show()
#
##index=np.where()












