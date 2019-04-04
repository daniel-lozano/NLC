import numpy as np
import math as mt

'''
DESCRIPTION: This function returns the contribution of each singular term defined as  multiplicity times the weigth of the cluster type.
PARAMETERS: The parameter of the function is an array with the values of the mean values of an operator.
OUTPUT: The function returns an array with the expression for each contribution up to some order k defined by the size of the input array.
'''

def cluster_contributions(weights ,order):
#weights=[O0,O1,O2,O3,O4Y,O4I,O4L,....] ndarray
    print("order=", order)
    if(order+1>5):
        raise ValueError("The order entered is to high! max order=4")
    contribution=np.zeros((order+1,weights.shape[1],weights.shape[1]))

    for i in range(order+1):
        ### multiplicity factor is the prefactor number to every combination of weigths ###
        if i==0: contribution[i]=1.*(weights[0])
        
        if i==1: contribution[i]=0.5*(weights[1]-4*weights[0])
        
        if i==2: contribution[i]=1.0*(weights[2]-2*weights[1]+weights[0])
        
        if i==3: contribution[i]=3.0*(weights[3]-2*weights[2]+weights[1])
        
        if i==4: contribution[i]=2.0*(weights[4]-3*weights[3]+3*weights[2]-1*weights[0])+\
                                3.0*(weights[5]-2*weights[3]+1*weights[2])+\
                                6.0*(weights[6]-2*weights[3]+1*weights[2])
    return contribution

def cluster_contributions_array(weights ,order):
    #weights=[O0,O1,O2,O3,O4Y,O4I,O4L,....] ndarray
    if(order+1>5):
        raise ValueError("The order entered is to high! max order=4")
    contribution=np.zeros(order+1)

    for i in range(order+1):
        
        if i==0: contribution[i]=1*(weights[0])
        
        if i==1: contribution[i]=0.5*(weights[1]-4*weights[0])
        
        if i==2: contribution[i]=1*(weights[2]-2*weights[1]+weights[0])
        
        if i==3: contribution[i]=3*(weights[3]-2*weights[2]+weights[1])
        
        if i==4: contribution[i]=2*(weights[4]-3*weights[3]+3*weights[2]-1*weights[0])+\
            3*(weights[5]-2*weights[3]+1*weights[2])+\
            6*(weights[6]-2*weights[3]+1*weights[2])
    return contribution



'''
DESCRIPTION: This function returns the bare sums of a series
PARAMETERS: The parameter of the function is an array with the terms of the series
OUTPUT: The function returns an array with the bare sums of the series, the bare sum order corresponds to the index of the list.
'''


def bare_sum(serie_terms):
    
    sum_terms=np.zeros(serie_terms.shape)
    Contribution=np.zeros((serie_terms.shape[1],serie_terms.shape[2]))
    
    for i in range(len(serie_terms)):
        Contribution+=serie_terms[i]
        sum_terms[i]+=Contribution
    
    return sum_terms

def bare_sum_array(serie_terms):
    
    sum_terms=np.zeros(len(serie_terms))
    Contribution=0
    
    for i in range(len(serie_terms)):
        Contribution+=serie_terms[i]
        sum_terms[i]+=Contribution
    
    return sum_terms
#

'''
    DESCRIPTION: This function returns an array with the pascal triangle numbers
    PARAMETERS: The parameter of the function are the order to which the numbers are calculated and the number of terms 
    OUTPUT: The function returns an array pascal triangle numbers
'''


def pascal_triangle(order,number_of_terms):#number>order+1
    
    array=np.zeros(number_of_terms)
    
    if(number_of_terms<=order):
        raise ValueError("the number of terms is smaller than the order required")
    
    for i in range(order+1):
        array[i]=mt.factorial(order)/((mt.factorial(order-i))*(mt.factorial(i)))
    
    return array


'''
    DESCRIPTION:
    PARAMETERS:
    OUTPUT:
'''
def product(array_like,ndarray_like):
    product=np.zeros((ndarray_like.shape[1],ndarray_like.shape[2]))
    for i in range(len(array_like)):
        product+=array_like[i]*ndarray_like[i]
    
    return product

def euler_sum(series_terms,start,stop):

    #Adding first terms before starting the Euler transformation

    first_terms=np.zeros((series_terms.shape[1],series_terms.shape[2]))
    ### Adding first terms S_0+S_1+S_2
    for i in range(start):
        first_terms+=series_terms[i]
    ### Terms S_3 and S_4, [S_3,S_4]
    terms_to_add=series_terms[start:(stop+1)]

    sum_tot=np.zeros((stop-start+1,series_terms.shape[1],series_terms.shape[2]))

    for i in range(stop-start+1):#

        print("Pascal triangle", pascal_triangle(i,len(terms_to_add)),2.**(i+1))
        
        if(i==0):### Adding Sum_3=S_0+S_1+S_2+S_3/2
            sum_tot[i]=first_terms+product(pascal_triangle(i,len(terms_to_add)),terms_to_add)/(2.**(i+1))
        
        else:    ### Adding Sum_4=S_0+S_1+S_2+S_3/2+(S_3+S_4)/4=Sum_3+(S_3+S_4)/4
            sum_tot[i]=sum_tot[i-1]+product(pascal_triangle(i,len(terms_to_add)),terms_to_add)/(2.**(i+1))

    return sum_tot





