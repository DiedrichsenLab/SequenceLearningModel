import numpy as np
import plan
import trigger

class Model:
    def __init__(self, numOptions=5, Aintegrate=0.996, Ainhibit=0, dT_motor=90, dT_visual=70, SigEps=0, 
                 plan_dict = {'func':'exp', 'capacity':10, 'param':[0.01,0.34]},
                 trigger_dict = {'func':'indep', 'win_trig':1, 'param':[0.6]}):
        """
        args:
            param (list/np.array): [scale,....]
            win_trig (int): window size for triggering (1 means naive triggering)
        """

        # Model parameters
        self.numOptions = numOptions
        self.Aintegrate = Aintegrate
        self.Ainhibit = Ainhibit
        self.dT_motor = dT_motor
        self.dT_visual = dT_visual
        self.SigEps = SigEps

        # relate to planning and triggering
        self.plan_dict = plan_dict
        self.trigger_dict = trigger_dict


    def get_theta(self):
        self.theta = getattr(plan,self.plan_dict['func'])(capacity = self.plan_dict['capacity'],param = self.plan_dict['param'])
        return self.theta
    
    def get_trigger(self,X):
        w = self.trigger_dict['win_trig']
        theta = self.theta
        A = self.Aintegrate

        # calculate the steady state of hidden state
        ss = theta[:w]/(1-A)

        trigger_check, element_done = getattr(trigger,self.trigger_dict['func'])(X,ss,w,param=self.trigger_dict['param'])
        return trigger_check, element_done
    


    #def get_bound(self,perc=0.587):
    #    self.Bound = perc*(self.plan_dict['param'][0])/(1-self.Aintegrate)

    # def get_param(self,fit):
    #     param = []
    #     for i, p in enumerate(fit):
    #         param.append(self.__dict__[p])
    #     return param

    # def reinit(self,fit,param):
    #     for i, p in enumerate(fit):
    #         self.__dict__[p] = param[i]