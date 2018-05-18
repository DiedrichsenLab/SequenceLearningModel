% simulations with the deaults of the function
%%things to work on:
% does the nmber of opions raise the reaction time
% what is the best paratmet combo for MT /RT

[IPIs, Exam] =    slm_diagModel();

% fliplr([0.2:.2:3])
%% Recepie for model diagnosis
close all

[IPIs, Exam] = slm_diagModel( 'numSimulations' , 5,...
    'SigEps' , 0.01 ,'DecayParam' , 2 ,'Aintegrate' , 0.9785 , 'theta_stim' , .0084 , 'Capacity' , 4 ,...
    'SeqLength' , 14,'Horizons' , [1:2:14],'Ainhibit' , [0.0]);


%%


%%
[IPIs, Exam] =    slm_diagModel( 'numSimulations' , 50,...
    'SigEps' , 0.01 ,'DecayParam' , 2 ,'Aintegrate' , 0.976 , 'theta_stim' , .01 , 'Capacity' , 3 ,...
    'SeqLength' , 14,'Horizons' , [1:2:14 , 14],'Ainhibit' , [0.0],'DecayParam' , 2);


[R,SIM,Trials,Models]=slm_testModel('simpleSeq','SigEps' , 0.01 ,'DecayParam' , 2 , 'Aintegrate' , 0.976 , 'theta_stim' , .0084 ,...
    'Horizons' , [1:14] , 'numSimulations' , 100 , 'Capacity' , 1,'VizProgress',0,'SeqLength' , 14);
%% Capacity
A = getrow(IPIs , IPIs.singleH==14);
A = IPIs;
figure('color' , 'white')
subplot(211)
lineplot(A.Capacity , A.MT , 'style_thickline')
title('MT')
grid on
set(gca , 'FontSize' , 16)
xlabel('Buffer Size')
subplot(212)
lineplot(A.Capacity , A.RT,'style_thickline')
title('RT')
xlabel('Buffer Size')
grid on
set(gca , 'FontSize' , 16)
%% Horizon
A = IPIs;
figure('color' , 'white')
subplot(211)
lineplot(A.singleH , A.MT , 'style_thickline')
title('MT')
grid on
set(gca , 'FontSize' , 16)
xlabel('Horizon')
subplot(212)
lineplot(A.singleH , A.RT,'style_thickline')
title('RT')
xlabel('Horizon')
grid on
set(gca , 'FontSize' , 16)


%% IPI

figure('color' , 'white')
IPIs.ipiNum = repmat([1:size(IPIs.pressTime , 2)-1] , size(IPIs.ipi , 1) , 1);
H = repmat(IPIs.singleH , 1 , size(IPIs.pressTime , 2)-1);
A = IPIs;%getrow(IPIs,IPIs.singleH==H(h));
index = fliplr([reshape(A.ipiNum , numel(A.ipiNum) , 1) reshape(H, numel(H) , 1)]);
data  = reshape(A.ipi , numel(A.ipi) , 1);
lineplot(index , data , 'style_shade');
title('IPIs')
xlabel('IPIs number')
grid on
set(gca , 'FontSize' , 16)
%%
A = IPIs;
IPIs = getrow(A , A.a == 1);

ipiMat = zeros(5 , 5);
for tn = 1:length(IPIs.MT)
    stim = Exam.T.stimulus(tn,:);
    for pr = 1:size(stim , 2)-1
        ipiMat(stim(pr) , stim(pr+1)) = ipiMat(stim(pr) , stim(pr+1)) + IPIs.ipi(pr);
    end
end
ipiMat= ipiMat/tn;
imagesc(ipiMat)
axis square

IPIs.a = ones(length(IPIs.MT) , 1);
IPIs1.a = zeros(length(IPIs1.MT) , 1);

IPIs = addstruct(IPIs , IPIs1);







A = IPIs1;
figure('color' , 'white')
subplot(211)
lineplot(A.Ainhibit , A.MT)
title('MT')
xlabel('Ainhibit')
subplot(212)
lineplot(A.Ainhibit , A.RT)
title('RT')
xlabel('Ainhibit')






A = getrow(IPIs , IPIs.singleH==14);
figure('color' , 'white')
subplot(211)
lineplot(A.Capacity , A.MT)
title('MT')
xlabel('Capacity')
subplot(212)
lineplot(A.Capacity , A.RT)
title('RT')
xlabel('Capacity')

A14 = getrow(IPIs , IPIs.singleH == 14 & IPIs.RT>=700 &IPIs.RT<=800 & IPIs.MT>=4500 & IPIs.MT<=5000);
plot(A14.MT , A14.RT)


A1 = getrow(IPIs , IPIs.singleH == 1 & IPIs.RT>=630 &IPIs.RT<=670 & IPIs.MT>=6300 & IPIs.MT<=6800);
plot(A1.MT , A1.RT)

% visually inspect
A = A14;
figure
for i = 1:length(A.MT)
    a = A.Aintegrate(i);
    t = A.theta_stim(i);
    s = A.SigEps(i);
    c = A.Capacity(i);
    d = A.DecayParam(i);
    T = getrow(A1 , A1.Aintegrate == a & A1.theta_stim == t & A1.Capacity == c &...
        A1.SigEps == s & A1.DecayParam == d);
    if ~isempty(T.Aintegrate)
        Target = getrow(IPIs , IPIs.Aintegrate == a & IPIs.theta_stim == t & ...
            IPIs.Capacity == c & IPIs.SigEps == s & IPIs.DecayParam == d);
    end
    subplot(211)
    lineplot(Target.singleH , Target.MT)
    subplot(212)
    lineplot(Target.singleH , Target.RT)
    title(num2str(i))
    pause()
    drawnow
end

% best i = 37 ~ 41
i = 41;
a = A14.Aintegrate(i);
t = A14.theta_stim(i);
s = A14.SigEps(i);
c = A14.Capacity(i);
d = A14.DecayParam(i);
Target = getrow(IPIs , IPIs.Aintegrate == a & IPIs.theta_stim == t & ...
    IPIs.Capacity == c & IPIs.SigEps == s & IPIs.DecayParam == d);

A = Target;
figure('color' , 'white')
subplot(211)
lineplot(A.singleH , A.MT)
title('MT')
xlabel('Horizon')
subplot(212)
lineplot(A.singleH , A.RT)
title('RT')
xlabel('Horizon')

%% recipie for learning

Sequences = [3 4 2 2 4 4 2 4 5 2 4 1 2 4;
    5 1 2 4 5 1 2 4 1 2 3 2 4 1;
    1 3 1 2 4 2 3 1 2 4 4 5 1 2];

[R,SIM,Trials,Models]=slm_testModel('seqLearn','Sequences' , Sequences , 'SigEps' , 0.01 ,'DecayParam' , 2 , 'Aintegrate' , 0.976 , 'theta_stim' , .0084 ,...
    'Horizons' , [1:14] , 'Capacity' , 4,'plotSim',0,'NumTrials' , 100);
[R,SIM,Trials,Models]=slm_testModel('simpleSeq','SigEps' , 0.01 ,'DecayParam' , 2 , 'Aintegrate' , 0.976 , 'theta_stim' , .0084 ,...
    'Horizons' , [1:14] , 'numSimulations' , 10 , 'Capacity' , 1,'plotSim',0,'SeqLength' , 14);

[R,SIM,Trials,Models]=slm_testModel('simpleSeq' , 'SigEps' , 0.01 ,'DecayParam' , 2 ,...
    'Aintegrate' , 0.98 , 'theta_stim' , .0084 , 'Capacity' , 3,'Ainhibit', 0.001,'SeqLength',14);


[IPIs, Exam] =    slm_diagModel( 'Ainhibit' , [0.0],'DecayParam' , 2);



% [IPIs, Exam] =    slm_diagModel( 'numSimulations' , 50,...
%     'SigEps' , 0.01 ,'DecayParam' , 2 ,'Aintegrate' , 0.976 , 'theta_stim' , .01 , 'Capacity' , 3 ,...
%     'SeqLength' , 14,'Horizons' , [1:2:14 , 14],'Ainhibit' , [0.0],'DecayParam' , 2);

%% Optmization recipes
% optimization step number 1
% we want to keep the number of parameters under 6 to maintin good performnace in fminsearch. 
% in the first round of optimization, set the boundry value to groups of
% presses and mainly optimize for 'theta_stim'  'Aintegrate' 'SigEps' on
% day 1  - window size = 1

parName = {'theta_stim'  'Aintegrate' 'SigEps' 'Bound(1)' 'Bound(2:3)' 'Bound(4:12)' , 'Bound(13:14)'};
initParam = [.0084 0.976 0.01 .45 .45 .45 .45];
loBound = [0.007 , 0.75 , 0.01 .1 .1 .1 .1];
hiBound = [0.02 , 0.988 0.05 .6 .6 .6 .6];
 
[Param Fval] = slm_optimize(Dall , initParam , 'parName' , parName,'runNum' ,1 , 'cycNum' , 1 ,'samNum'  , [] ,...
    'ItrNum' , 300 , 'loBound' , loBound , 'hiBound' , hiBound , 'Day' , [5] , 'Horizon' , [1] , 'poolHorizons' , [7:13]);

%% 
% optimization step number 2
% now keep 'theta_stim'  'Aintegrate' 'SigEps' consant within the
% slm_optimSimTrial as the last iteration of the previous optimization and
% optimize for more fine grained boundry values
% day 1  - window size = 1


parName = {'Bound(1)' 'Bound(2)' 'Bound(3)' 'Bound(4:10)' 'Bound(11)' 'Bound(12)' 'Bound(13)' 'Bound(14)' 'dtGrowth'};  


loBound = [.1 .1 .1 .1 .1 .1 .1 .1 1];
hiBound = [.6 .6 .6 .6 .6 .6 .6 .6 5] ;
H = {[1] [2] [3] [4] [5] [6] [7:13]};
for h = 1:length(H)
    if h == 1
        initParam = [0.45 0.45 0.45 0.45 0.45  0.45 0.45 0.45 1]; 
    else
        initParam = [0.45 0.45 0.45 0.45 0.45  0.45 0.45 0.45 3]; 
    end
    [Param Fval] = slm_optimize(Dall , initParam , 'parName' , parName,'runNum' ,3+.1*h , 'cycNum' , 5 ,'samNum'  , [] ,...
        'ItrNum' , 20 , 'loBound' , loBound , 'hiBound' , hiBound , 'Day' , [1] , 'Horizon' , H{h} , 'poolHorizons' , [7:13],...
        'customizeInitParam' , 1);
    close all
end
%%
M = [];
par = [];
parName = {'Bound(1)' 'Bound(2)' 'Bound(3)' 'Bound(4:10)' 'Bound(11)' 'Bound(12)' 'Bound(13)' 'Bound(14)'};  
H = {[1] [2] [3] [4] [5] [6] [7:13]};
for h = 1:length(H)
    filename = ['param2.'  , num2str(h) , '.mat'];
    load(['/Users/nedakordjazi/Documents/GitHub/SequenceLearningModel/optim2/' , filename]);
    par(h,:) = param.par(end , :);
    R = slm_optimSimulate(Dall , par(h,:)  , 'parName' , parName,'samNum'  , [] , 'Day' , 1, 'Horizon' , H{h} , 'poolHorizons' , [7:13]);
    T.RT     = R(:,1);
    T.IPI    = R(:,2:14);
    T.MT     = R(:,15);
    T.window = h*ones(size(T.RT));
    M = addstruct(M , T);
    clear T
    close all
end

%% MT RT



All  = M;

figure('color' , 'white')
subplot(211)
colorz = {[0.840000000000000,0.360000000000000,0.501176470588235],[0.360000000000000,0.456470588235294,0.760000000000000]};
lineplot(All.window , All.MT , 'plotfcn' , 'nanmean',...
    'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'} );

title('MT')
grid on
set(gca , 'FontSize' , 16)
xlabel('Horizon')
subplot(212)
lineplot(All.window , All.RT , 'plotfcn' , 'nanmean',...
    'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'});
title('RT')
xlabel('Horizon')
grid on
set(gca , 'FontSize' , 16)


%% IPI
clear Fit Act
Fit = M;
Fit.ipiNum = repmat([1:13] , length(Fit.window) , 1);
Fit.IPI = reshape(Fit.IPI , numel(Fit.IPI) , 1);

Fit.window  = repmat(Fit.window , 1 , 13);
Fit.window  = reshape(Fit.window , numel(Fit.IPI) , 1);
Fit.ipiNum = reshape(Fit.ipiNum , numel(Fit.ipiNum) , 1);

All  = Fit ;

figure('color' , 'white')

lineplot(All.ipiNum , All.IPI , 'plotfcn' , 'nanmean',...
    'linecolor' , colorz,'split'  , All.window , ...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'});

title('IPIs')
xlabel('IPIs number')
grid on
set(gca , 'FontSize' , 16)


