"""
a set of different functions for specifying planning rates 
"""
import numpy as np

def exp(capacity=3, param=[1/3, 0.01]):
    """ exponential decay function for planning rates

    args:
        capacity (int): capacity of planning
        param (list/np.array): [slope, scale]
    """
    x = np.arange(capacity)
    theta = param[1]*np.exp(-x*param[0])

    return theta

def inv(capacity=3, param=[1/3, 0.01]):
    """ inverse decay function for planning rates

    args:
        capacity (int): capacity of planning
        param (list/np.array): [slope, scale]
    """
    x = np.arange(capacity)
    theta = param[0]/(x+1)

    return theta

def arb():
    """ arbitrary function for planning rates
    """
    pass