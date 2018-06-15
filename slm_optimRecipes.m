%% 1- optimize the decision boundary and decayParam with Giacomos parameters Aintegrate ana Theta_stim with noise set to 0
% clear param parName par day R
parName = { 'bAll' , 'Box'};
loBound = [];
hiBound = [];
h=0;
load('param0_0_4  5.mat')
initParam = [.5 , 2];
day = [4 5];
h = 0;
[Param Fval] = slm_optimize(Dall ,  initParam , 'parName' , parName,'runNum' ,['1_box_',num2str(h),'_',num2str(day)] , 'cycNum' , 7 ,'samNum'  , [5] ,...
    'ItrNum' , 1000 , 'loBound' , loBound , 'hiBound' , hiBound , 'Day' , day , 'Horizon' , [1:5] , 'poolHorizons' , [5:13],...
    'noise' , 0 ,  'subjNum' , [1:15] , 'desiredField' , {'MT'} , 'noisefreeRep' , [],  ...
    'MsetField' , {'PlanningCurve' , 'logistic'  ,'theta_stim' ,0.0084,'Aintegrate' , 0.985});
close all

%% 2- optimize the initial decision boundary to get the RTs for every window size
% clear param parName par day R
parName = { 'binit' , 'B_coef'};
loBound = [];
hiBound = [];
h=0;
load('param1_1_0_4  5.mat.mat')
initParam = [.3 , 1];
day = [4 5];
h = 0;
[Param Fval] = slm_optimize(Dall ,  initParam , 'parName' , parName,'runNum' ,['1_1_',num2str(h),'_',num2str(day)] , 'cycNum' , 7 ,'samNum'  , [5] ,...
    'ItrNum' , 1000 , 'loBound' , loBound , 'hiBound' , hiBound , 'Day' , day , 'Horizon' , [1:5] , 'poolHorizons' , [5:13],...
    'noise' , 0 ,  'subjNum' , [1:15] , 'desiredField' , {'MT'} , 'noisefreeRep' , [],  ...
    'MsetField' , {'PlanningCurve' , 'logistic'  ,'theta_stim' ,0.0084,'Aintegrate' , 0.985});
close all
%% 2 - create the noise-free simulation
day = [4 5];
load('param1_1_0_4  5.mat')
parName = param.parName(end,:); % [0.492,2.32]
par = param.par(end , :);
[R] = slm_optimSimulate(Dall , par  , 'parName' , parName,'samNum'  , 100 ,...
        'Day' , day, 'Horizon' , [1:5] , 'poolHorizons' , [5:13] , 'noise' ,0, 'subjNum' , [1:15],...
        'MsetField' , ...
        {'PlanningCurve' , 'ramp' ,  'SigEps' , 0.0 ,'theta_stim' ,0.0084,'Aintegrate' , 0.985});% ,'NumPresses',1 , 'stimulus' , [3]);


%% 3- optimize the noise level using the parametrs of the noise free optimization as presets
% in the noisey round, we fit not to the actual data, but to the niseless
% simulation

% clear param parName par day R
parName = { 'SigEps'};%,'theta_stim' ,'Aintegrate'};
loBound = [];
hiBound = [];
h=0;
load('param1_0_4  5.mat')
initParam = .0089;%param.par(end , 1);
day = [4 5];
h = 0;
[Param Fval] = slm_optimize(Dall ,  initParam , 'parName' , parName,'runNum' ,['2_',num2str(h),'_',num2str(day)] , 'cycNum' , 7 ,'samNum'  , [5] ,...
    'ItrNum' , 1000 , 'loBound' , loBound , 'hiBound' , hiBound , 'Day' , day , 'Horizon' , [1:5] , 'poolHorizons' , [5:13],...
    'noise' , 1 ,  'subjNum' , [1:15] , 'desiredField' , {'MT'} , 'noisefreeRep' , R,  ...
    'MsetField' , {'PlanningCurve' , 'exp'  ,'theta_stim' ,0.0084 ,'Aintegrate' , 0.985 , 'bAll' , param.par(end , 1) , 'DecayParam' ,  param.par(end , 2)});
close all



%% 4-  noisey simulation  

day = [4 5];

load('param1_0_4  5.mat')
b = param.par(end , 1);
dp = param.par(end , 2);
load('param2_0_4  5.mat')
parName = param.parName(end,:); 
par = .022;%param.par(end , :);% [0.00856624999917112]
[RN] = slm_optimSimulate(Dall , par  , 'parName' , parName,'samNum'  , 100 ,...
        'Day' , day, 'Horizon' , [1:5] , 'poolHorizons' , [5:13] , 'noise' ,1, 'subjNum' , [1:15],...
        'MsetField' ,...
        {'PlanningCurve' , 'exp'  ,'theta_stim' ,0.0084 ,'Aintegrate' , 0.985 , 'bAll' , b, 'DecayParam',dp});%,'NumPresses',1 , 'stimulus' , [3]);
    
%% single finger visual    
% R_Nsing = RN;
% R_sing = R;

R_Nsing = getrow(R_Nsing , ~R_Nsing.isError);
figure('color' , 'white')
for tn = 1:length(R_Nsing.MT)
    plot([1:2:length(R_Nsing.X{tn})*2] , R_Nsing.X{tn}(3,:) , 'color' , [.6 .6 .6]);
    hold on
    plot([R_Nsing.decisionTime(tn) R_Nsing.decisionTime(tn)] , [0 nanmean(R_sing.B{1})] , 'color' , 'b')
%     plot([1:2:length(R_Nsing.X{tn})*2] , R_Nsing.B{tn} , 'color' , 'r')
end
plot([nanmedian(R_Nsing.decisionTime) nanmedian(R_Nsing.decisionTime)] , [0 nanmean(R_sing.B{1})] , ':' , 'color' , 'm','LineWidth' , 3)
line([0 length(R_sing.X{1})*2] , [mean(R_sing.B{1}) mean(R_sing.B{1})] , 'color' , 'r' , 'LineWidth' , 3)
hold on
plot([1:2:length(R_sing.X{1})*2] , R_sing.X{1}(3,:) , 'color' , 'k' , 'LineWidth' , 3);
line([R_sing.decisionTime R_sing.decisionTime] , [0 mean(R_sing.B{1})] , 'color' , 'b' , 'LineWidth' , 3)
line([0 length(R_sing.X{1})*2] , [mean(R_sing.B{1}) mean(R_sing.B{1})] , 'color' , 'r' , 'LineWidth' , 3)
    
%% sequence visual   
% R_seq = R;
% RN_seq = RN;
Dall.RT = Dall.AllPressTimes(:,1)-1500;
A = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day , [4 5]) &...
ismember(Dall.Horizon , [1:5]));
figure('color' , 'white')
subplot(211)
hold on
plot(R_seq.MT , 'o-', 'color' , [0 0 1] )  
lineplot(A.Horizon  , A.MT ,  'plotfcn','nanmedian' , 'linecolor' ,  [0 1 0 ],...
                    'errorcolor' , [0 1 0]) 
                
lineplot(RN_seq.singleH , RN_seq.MT , 'subset' , ~RN_seq.isError, 'plotfcn','nanmedian' , 'linecolor' , [1 0 0 ],...
                    'errorcolor' , [1 0 0])


                
subplot(212)
hold on
plot(R_seq.RT , 'o-', 'color' , [0 0 1] )  
lineplot(RN_seq.singleH , RN_seq.RT , 'subset' , ~RN_seq.isError, 'plotfcn','nanmedian' , 'linecolor' , [1 0 0 ],...
                    'errorcolor' , [1 0 0])

lineplot(A.Horizon  , A.RT ,  'plotfcn','nanmedian' , 'linecolor' ,  [0 1 0 ],...
                    'errorcolor' , [0 1 0]) 
              

    
    %% 3 - Simulate with levels of noise
se = [0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09];
h = 0;
clear param parName par day R
day = [4 5];
filename = 'paramexp.mat';
% filename = 'paramlog.mat';
% filename = ['param'  , '6_', num2str(h),'_',num2str(day), '.mat'];
load(['/Users/nedakordjazi/Documents/GitHub/SequenceLearningModel/' , filename]);
parName = param.parName(end,:);
allFit = [];
parName = [param.parName(end,:) , 'SigEps'];
for i = 1:length(se)
    day = [4 5];
    close all
    clear par R
    par = [param.par(end , :) se(i)];
    [R] = slm_optimSimulate(Dall , par  , 'parName' , parName,'samNum'  , 100 ,...
        'Day' , day, 'Horizon' , [1:13] , 'poolHorizons' , [5:13] , 'noise' ,1 , 'subjNum' , [1:15],...
        'MsetField' ,...
        {'PlanningCurve' , 'exp', 'SigEps' , 0.034 ,'theta_stim' ,0.0084 ,'Aintegrate' , 0.985});
    R.SigEps = se(i)*ones(size(R.MT));
    allFit = addstruct(allFit , R);
end

A = getrow(allFit , ~allFit.isError);
figure('color' , 'white')
lineplot(A.singleH , A.MT , 'plotfcn' , 'nanmedian',...
                'split', A.SigEps  , 'leg' , 'auto');