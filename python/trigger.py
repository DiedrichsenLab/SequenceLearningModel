"""
a set of different functions for specifying triggering rule 
"""
import numpy as np

def indep(X,ss,w,param=[0.72]):
    """ 
    Args:
        X (np.array): 5*win_trig, hidden state
        w (int): triggering window
        ss (list/np.array): steady state of hidden state
        param (list): triggering parameter
    """
    horizoncheck = False
    nHorizon = X.shape[1]
    element_done = np.full(w,0)

    for i in range(nHorizon):
        if any(X[:,i] > (ss[i]*param[0])):
            element_done[i] = 1

    if np.sum(element_done)==nHorizon:
        horizoncheck = True

    return horizoncheck, element_done



def sum(X,ss,w,param=[0.45]):
    """ 
    Args:
        X (np.array): 5*win_trig, hidden state
        w (int): triggering window
        ss (list/np.array): steady state of hidden state
        param (list): triggering parameter
    """
    horizoncheck = False
    nHorizon = X.shape[1]
    element_done = np.full(w,0)

    if nHorizon==0:
        horizoncheck = True
        return horizoncheck, element_done

    X = np.sum(X,axis=1)
    ss = np.sum(ss)

    if any(X>(ss*param[0])):
        horizoncheck = True
        element_done[:]=1

    return horizoncheck, element_done