%% 1- optimize 
% clear param parName par day R
parName = { 'bAll' ,'theta_stim' ,'Aintegrate','SigEps'};
loBound = [];
hiBound = [];
h=0;
initParam = param.par(end , :);
day = [4 5];
h = 0;
[Param Fval] = slm_optimize(Dall ,  initParam , 'parName' , parName,'runNum' ,['3_',num2str(h),'_',num2str(day)] , 'cycNum' , 7 ,'samNum'  , [20] ,...
    'ItrNum' , 1000 , 'loBound' , loBound , 'hiBound' , hiBound , 'Day' , day , 'Horizon' , [1:13] , 'poolHorizons' , [5:13],...
    'noise' , 1 ,  'subjNum' , [1:15] , 'desiredField' , {'MT'} , ...
    'MsetField' , {'PlanningCurve' , 'exp'});
close all




%% 2- Simulate

clear param parName par day R
day = [4 5];
filename = 'param3_0_4  5.mat';
% filename = 'paramlog.mat';
% filename = ['param'  , '6_', num2str(h),'_',num2str(day), '.mat'];
load(['/Users/nkordjazi/Documents/GitHub/SequenceLearningModel/' , filename]);
parName = param.parName(end,:);
par = param.par(end , :);
[RN] = slm_optimSimulate(Dall , par  , 'parName' , parName,'samNum'  , 100 ,...
        'Day' , day, 'Horizon' , [1:13] , 'poolHorizons' , [5:13] , 'noise' ,1, 'subjNum' , [1:15],...
        'MsetField' , {'PlanningCurve' , 'exp'} , 'NumPresses' , 5);%  , 'stimulus' , [3]);
    
figure('color' , 'white')
for tn = 1:length(R_N.MT)
    plot([1:2:length(R_N.X{tn})*2] , R_N.X{tn}(3,:) , 'color' , [.6 .6 .6]);
    hold on
    plot([R_N.pressTime(tn) R_N.pressTime(tn)] , [0 nanmean(R.B{1})] , 'color' , 'b')
%     plot([1:2:length(R_N.X{tn})*2] , R_N.B{tn} , 'color' , 'r')
end
line([0 length(R.X{1})*2] , [mean(R.B{1}) mean(R.B{1})] , 'color' , 'r' , 'LineWidth' , 3)
hold on
plot([1:2:length(R.X{1})*2] , R.X{1}(3,:) , 'color' , 'k' , 'LineWidth' , 3);
line([R.pressTime R.pressTime] , [0 mean(R.B{1})] , 'color' , 'b' , 'LineWidth' , 3)
    
    
    
    
    
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
        'MsetField' , {'PlanningCurve' , 'exp'});
    R.SigEps = se(i)*ones(size(R.MT));
    allFit = addstruct(allFit , R);
end

A = getrow(allFit , ~allFit.isError);
figure('color' , 'white')
lineplot(A.singleH , A.MT , 'plotfcn' , 'nanmedian',...
                'split', A.SigEps  , 'leg' , 'auto');