#!/usr/bin/env python

# ==============================================================================
# WRITTEN BY: Jeremy Goh (for OFYP 18/19)
# SUPERVISOR: Prof. Michel Gingras
#
# This module contains the following functions and algorithms:
# - gen_contributions: Generates contributions at each NLC order (\Lambda)
# - gen_baresums: Computes (bare) partial sums from contributions
# - euler: Euler's transformation (for resummation)
# ==============================================================================

import numpy as np

def gen_contributions(terms: np.ndarray, max_order: int)->np.ndarray:
    """
        Generates contributions at each NLC order in a numpy array, given the individual cluster averages.
        Thus, this incorporates the multiplicities L(c)  and weights W(c) for the pyrochlore lattice.
        Primarily for Euler summation.
        
        :param terms: array containing cluster averages (raw output from main scripts)
        :param max_order: max NLC order should correspond to the size of the 'terms' array
        :return: array of contributions up to max NLC order
        """
    
    contributions_arr = np.zeros((max_order+1, *terms[0].shape))  # init array of contributions
    for order in range(max_order+1):
        if order == 0:
            contributions_arr[order] = terms[0]
        elif order == 1:
            contributions_arr[order] = 0.5*terms[1] - 2*terms[0]
        elif order == 2:
            contributions_arr[order] = terms[2] - 2*terms[1] + terms[0]
        elif order == 3:
            contributions_arr[order] = 3*terms[3] - 6*terms[2] + 3*terms[1]
        elif order == 4:
            # ensure index match cluster type -> [4]: 4Y, [5]: 4I, [6]: 4L
            contributions_arr[order] = 2*terms[4] + 3*terms[5] + 6*terms[6] - \
                24*terms[3] + 15*terms[2] - 2*terms[1]
        else:
            raise ValueError("Invalid max_order provided!", max_order)

    return contributions_arr


def gen_baresums(terms: np.ndarray, return_last: bool=True):
    """
        Computes the partial sums naively (bare sums; without resummation techniques).
        Primarily for Wynn algorithm.
        
        :param terms: array of contributions at each NLC order (cf. gen_contributions function)
        :param return_last: whether to return the array of partial sums or just the last partial sum
        :return: array of (bare) partial sums (if return_last is False),
        the final partial sum (otherwise)
        """
    sum_tmp = 0
    baresums_arr = np.zeros_like(terms)
    for order in range(terms.shape[0]):
        sum_tmp += terms[order]
        baresums_arr[order] += sum_tmp
    if return_last: return baresums_arr[-1]
    else: return baresums_arr


# Euler's transformation algorithm
def euler(cont: np.ndarray, start: int, stop: int):
    """
        Evaluates analytical expression for Euler transformation (generic, up to NLC-4)
        
        :param cont: array of contributions at each NLC order (cf. gen_contributions function)
        :param start: order to start Euler transformation with
        :param stop: order to stop Euler transformation with (here, should always be <= 4)
        :return: Euler transformed property
        """
    p_sum = np.zeros_like(cont[0])
    for i in range(start):
        p_sum += cont[i]
    for i in range(stop-start+1):
        # manually insert forward-differencing operator
        if i == 0: fdo = cont[start]
        elif i == 1: fdo = cont[start+1] + cont[start]
        elif i == 2: fdo = cont[start+2] + 2*cont[start+1] + cont[start]
        elif i == 3: fdo = cont[start+3] + 3*cont[start+2] + 3*cont[start+1] + cont[start]
        else: raise ValueError(f'Invalid combination of start ({start}) and stop ({stop}) values.')
        
        p_sum += 1/(2**(i+1)) * fdo
    return p_sum
