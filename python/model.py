import numpy as np
import plan
import trigger

class Model:
    def __init__(self, numOptions=5, Aintegrate=0.995, Ainhibit=0, dT_motor=90, dT_visual=70, SigEps=0, Bound=1, 
                 plan_dict = {'func':'exp', 'capacity':10, 'param':[0.34, 0.01]},
                 trigger_dict = {'func':'indep', 'win_trig':0, 'param':[0.72]}):

        # Model parameters
        self.numOptions = numOptions
        self.Aintegrate = Aintegrate
        self.Ainhibit = Ainhibit
        self.dT_motor = dT_motor
        self.dT_visual = dT_visual
        self.SigEps = SigEps
        self.Bound = Bound

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
        ss = theta[1:(w+1)]/(1-A)

        horizoncheck, element_done = getattr(trigger,self.trigger_dict['func'])(X,ss,w,param=self.trigger_dict['param'])
        return horizoncheck, element_done

    # def get_param(self,fit):
    #     param = []
    #     for i, p in enumerate(fit):
    #         param.append(self.__dict__[p])
    #     return param

    # def reinit(self,fit,param):
    #     for i, p in enumerate(fit):
    #         self.__dict__[p] = param[i]