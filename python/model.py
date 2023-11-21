import numpy as np

class Model:
    def __init__(self, numOptions=5, Aintegrate=0.995, Ainhibit=0, theta=0.01, dT_motor=90, dT_visual=70, SigEps=0.01, Bound=1, capacity=3, coeff=0.5, wFuture=0, power=1):

        # Model parameters
        self.numOptions = numOptions
        self.Aintegrate = Aintegrate
        self.Ainhibit = Ainhibit
        self.theta = theta
        self.dT_motor = dT_motor
        self.dT_visual = dT_visual
        self.SigEps = SigEps
        self.Bound = Bound
        self.capacity = capacity

        # newly added parameters
        self.coeff = coeff
        self.wFuture = wFuture

    def get_param(self,fit):
        param = []
        for i, p in enumerate(fit):
            param.append(self.__dict__[p])
        return param

    def reinit(self,fit,param):
        for i, p in enumerate(fit):
            self.__dict__[p] = param[i]

    def plan_exp(self,dec,nDecision,power=1):
        mult = np.exp(-np.power(dec - nDecision,power) / self.capacity)  # Distribution of weight onto future decisions
        mult[dec < nDecision] = 0  # Made decisions will just decay
        rate = self.theta * mult
        return rate
    
    def plan_inv(self,dec,nDecision):
        denom = dec-nDecision+1
        mask = denom>0
        rate = np.zeros_like(denom,dtype=float)
        rate[mask] = self.theta*(1/denom[mask])
        return rate


    def trig_indep(self,X_subset,rate_subset):

        horizoncheck = False
        nHorizon = X_subset.shape[1]
        element_done = np.full(self.wFuture,0)


        for i in range(nHorizon):
            if any(X_subset[:,i] > (rate_subset[i]/(1-self.Aintegrate)) * self.coeff):
                element_done[i] = 1
        if np.sum(element_done)==nHorizon:
            horizoncheck = True

        return horizoncheck, element_done
    
    def trig_sum(self,X_subset,rate_subset):

        horizoncheck = False
        nHorizon = X_subset.shape[1]
        element_done = np.full(self.wFuture,0)

        if nHorizon==0:
            horizoncheck = True
            return horizoncheck, element_done

        X = np.sum(X_subset,axis=1)
        rate = np.sum(rate_subset)

        if any(X>(rate/(1-self.Aintegrate))*self.coeff):
            horizoncheck = True
            element_done[:]=1

        return horizoncheck, element_done