% simulations with the deaults of the function 
%%things to work on:
% does the nmber of opions raise the reaction time
% what is the best paratmet combo for MT /RT

[IPIs, Exam] =    slm_diagModel();

% fliplr([0.2:.2:3])
parameters:

%% simulation with the 'boxcar' function
 [IPIs, Exam] =    slm_diagModel('DecayFunc' , 'boxcar' , 'numSimulations' , 30,...
     'Aintegrate' , 0.975 , 'theta_stim' , 0.01 , 'DecayParam' , [2], 'SigEps' , 0.005,'SeqLength' , 14);
%% 
% load('slm_IPIsExp.mat')
load('slm_IPIsBoxcar.mat')
% separate the full horizons only
% A = getrow( IPIs , IPIs.Horizon == max(IPIs.Horizon));
A = getrow( IPIs , IPIs.Horizon == 14);
slm_diagModelViz(A , 'MT_RT_vs_theta' , 'Aintegrate_select' , 0.975)
slm_diagModelViz(IPIs , 'MT_RT_vs_Aintegrate')
%%
slm_diagModelViz(IPIs , 'MT_RT_vs_Horizon_constantAiTs', 'Aintegrate_select' , 0.9950)
%%
A = getrow( IPIs , IPIs.Aintegrate == 0.98 & IPIs.theta_stim == 0.01);
slm_diagModelViz(IPIs , 'MT_RT_vs_Horizon','Aintegrate_select' , 0.975)
slm_diagModelViz(IPIs , 'MT_RT_vs_Horizon')
%%
A = getrow( IPIs , IPIs.Horizon == 10);
slm_diagModelViz(A , 'IPIs_vs_theta')

slm_diagModelViz(A, 'IPIs_vs_Aintegrate')

