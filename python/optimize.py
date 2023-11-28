from experiment import Exp
from sim_trial import do_sim
import numpy as np
from scipy.stats import logistic
from model import Model
from scipy.optimize import minimize,fmin
from plot import plot_IPI
import matplotlib.pyplot as plt


def loss(param,M,fit,T,ipi_ideal):
    
    # check if we have Aintegrate and Ainhibit
    for i,p in enumerate(fit):
        if p=='Aintegrate' or p=='Ainhibit' or p=='coeff':
            param[i] = logistic.cdf(param[i])

    M.reinit(fit,param)

    T,SIM = do_sim(M,T)

    ipi = np.diff(T.pressTime)

    if len(T.pressTime)<T.numPress:
        return np.inf

    loss = np.linalg.norm(ipi-ipi_ideal)
    #print(loss)
    return loss

def optimize(M,fit,T,ipi_ideal):

    param0 = M.get_param(fit)
    for i,p in enumerate(fit):
        if p=='Aintegrate' or p=='Ainhibit' or p=='coeff': # these parameters are bounded between 0 and 1
            param0[i] = logistic.ppf(param0[i])

    if param0:
        out = fmin(loss,param0,args=(M,fit,T,ipi_ideal),maxiter=10000,disp=False)
        #print(out.success)
    
    return M
    

if __name__ == "__main__":

    # model
    M = Model(SigEps=0,wFuture=2,coeff=0.7,Bound=1.4675,capacity=3.5,Aintegrate=0.9956,power=1)
    
    # experiment
    T = Exp()
    T.seqShow(RT=6000)
    
    # ideal IPI
    ipi_ideal = [100,130,130,130,130,130,130,130,100]
    ipi_ideal = [150, 180, 185, 193, 196, 199, 202, 200, 150]
    ipi_ideal = [92,132,208,208,208,208,208,208,132]

    # what parameters to fit
    fit = ['coeff','capacity']

    # optimize
    M = optimize(M,fit,T,ipi_ideal)
    # test
    T,SIM = do_sim(M,T)
    ipi = np.diff(T.pressTime)
    IPI=[]
    IPI.append(ipi)
    IPI.append(ipi_ideal)
    label=['model','ideal']
    fig,ax = plot_IPI(IPI,label)
    plt.show()  

