import numpy as np
import math as mt

'''
DESCRIPTION: This function returns the contribution of each singular term defined as  multiplicity times the weigth of the cluster type.
PARAMETERS: The parameter of the function is an array with the values of the mean values of an operator.
OUTPUT: The function returns an array with the expression for each contribution up to some order k defined by the size of the input array.
'''

def cluster_contributions(weights ,order):
#weights=[W0,W1,W2,W3,W4Y,W4I,W4L,....]
    print("order=", order)
    if(order+1>5):
        raise ValueError("The order entered is to high! max order=4")
    contribution=np.zeros((order+1,weights.shape[1],weights.shape[1]))

    for i in range(order+1):
        
        if i==0: contribution[i]=1*(weights[0])
        
        if i==1: contribution[i]=0.5*(weights[1]-4*weights[0])
        
        if i==2: contribution[i]=1*(weights[2]-2*weights[1]+weights[0])
        
        if i==3: contribution[i]=3*(weights[3]-2*weights[2]+weights[1])
        
        if i==4: contribution[i]=2*(weights[4]-3*weights[3]+3*weights[2]-1*weights[0])+\
                                3*(weights[5]-2*weights[3]+1*weights[2])+\
                                6*(weights[6]-2*weights[3]+1*weights[2])
    return contribution

#W=np.ones(2)
#print W
#print cluster_contributions(W)

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

#for i in range(4):
#    print i
#    print pascal_triangle(i,4)

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
    print(type(series_terms))
    
    first_terms=np.zeros((series_terms.shape[1],series_terms.shape[2]))
    
    for i in range(start):
        first_terms+=series_terms[i]
        
    terms_to_add=series_terms[start:(stop+1)]

    sum_tot=np.zeros((stop-start+1,series_terms.shape[1],series_terms.shape[2]))

    for i in range(stop-start+1):
        
        print(pascal_triangle(i,len(terms_to_add)))
        
        sum_tot[i]=first_terms+product(pascal_triangle(i,len(terms_to_add)),terms_to_add)/2.**(i+1)

    return sum_tot




