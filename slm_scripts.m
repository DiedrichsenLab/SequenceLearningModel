% simulations with the deaults of the function 
%%things to work on:
% does the nmber of opions raise the reaction time
% what is the best paratmet combo for MT /RT

[IPIs, Exam] =    slm_diagModel();

% fliplr([0.2:.2:3])
%% simulation with the 'boxcar' function
close all

 [IPIs, Exam] =    slm_diagModel('DecayFunc' , 'boxcar' , 'numSimulations' , 50,...
     'Aintegrate' , [0.96:0.002:0.98] , 'theta_stim' , [.01:0.002:0.03] , 'DecayParam' , [4], 'SigEps' , 0.0135,'SeqLength' , 14 , 'Horizons' , [1:6 , 14]);
 sum(Exam.R.isError)/length(Exam.R.isError)
 slm_diagModelViz(IPIs , 'MT_RT_vs_Horizon')
%% 
% load('slm_IPIsExp.mat')
load('slm_IPIsBoxcar.mat')
% separate the full horizons only
% A = getrow( IPIs , IPIs.Horizon == max(IPIs.Horizon));
% A = getrow( IPIs , IPIs.Horizon == 14);
slm_diagModelViz(IPIs , 'MT_RT_vs_theta' , 'Aintegrate_select' , 0.975)
slm_diagModelViz(IPIs , 'MT_RT_vs_Aintegrate' , 'theta_select' , 0.01)
%%
slm_diagModelViz(IPIs , 'MT_RT_vs_Horizon_constantAiTs', 'Aintegrate_select' , 0.9950)
%%
A = getrow( IPIs , IPIs.Aintegrate == 0.98 & IPIs.theta_stim == 0.01);
slm_diagModelViz(IPIs , 'MT_RT_vs_Horizon','Aintegrate_select' , 0.975)
slm_diagModelViz(IPIs , 'MT_RT_vs_Horizon')
%%
% A = getrow( IPIs , IPIs.Horizon == 10);
slm_diagModelViz(IPIs , 'IPIs_vs_theta')

slm_diagModelViz(A, 'IPIs_vs_Aintegrate')

