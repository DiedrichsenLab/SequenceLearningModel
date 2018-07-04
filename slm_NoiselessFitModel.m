function slm_NoiselessFitModel(what , Dall , varargin)
planFunc = 'exp';
initParam = [.55 , 7];
parName = { 'bAll' , 'DecayParam'};
NumPresses = size(Dall.AllPress , 2);
stimulus = [];
ItrNum = 1000;
cycNum = 1;
day = [4 5]; % just fitting days 4, 5
Horizon = [1:5];
initPlan = 'uniform';
c = 1;
NameExt = [];
input_initalParam = [];
input_parName = [];
MsetField = {};
loBound = [];
hiBound = [];
MSF = [];
optimizeIPINumber = [1:3];
includeMT = 0;
includeRT = 1;
diffMT = 0;
if length(varargin)==1
    varargin = varargin{1}; % coming from a highr loop, so nneds to be unpacked
end
while(c<=length(varargin))
    switch(varargin{c})
        case {'runNum'}
            % the number of run
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'ItrNum'}
            % the number of times to sample the data and go through the optmization process
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'cycNum'}
            % the number of times to sample the data and go through the optmization process
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'NumPresses'}
            % the number of presses to include
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'stimulus'}
            % for when you want constant stimulus
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'planFunc'}
            % 'exp', 'logistic', 'box' , 'ramp' 'logistic+ramp'
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'NameExt'}
            % extension for fimename for fitting more than once
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'MsetField'}
            % the names and values of the fields we want to set in M
            % has to be cell of value names, followed by their values
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Horizon'}
            % what Horizon to include
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'initPlan'}
            % for the arbitrary planning function this defines the starting point
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'input_initalParam'}
            % for the arbitrary initial parameters to not use the hardcoded ones
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'input_parName'}
            % for the arbitrary parameter names to not use the hardcoded ones
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'loBound'}
            % lower bound for parameters
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'hiBound'}
            % Higher bound for parameters
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'includeRT'}
            % include RT in fitting or not
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'optimizeIPINumber'}  % for the 'FitIPIRT' case
            % IPI numbers to include in the optimization
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'includeMT'}  % for the 'FitIPIRT' case
            % include MT in theoptimization or not
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'diffMT'}  % for the 'FitIPIRT' case
            % use the difference between MTs instead of MTs
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end


switch planFunc
    case 'logistic'
        parName = { 'bAll' , 'B_coef1' 'B_coef2'};
        initParam = [.45 3 1];
    case 'ramp'
        parName = { 'bAll' , 'rampDecay1' , 'rampDecay2'};
        initParam = [.50 , 7 0];
    case 'box'
        parName = { 'bAll' , 'Box'};
        initParam = [.55 , 5];
    case 'box_logistic'
        parName = { 'bAll' , 'B_coef1' 'B_coef2'};
        initParam = [.5 , 3 3];
    case 'box_exp'
        parName = { 'bAll' , 'DecayParam'};
        initParam = [.5 , 7];
    case 'box_ramp'
        parName = { 'bAll' , 'rampDecay'};
        initParam = [.50 , 7 ];
    case 'arbitrary'
        parName = { 'bAll' , 'planFunc'};
        switch initPlan % initialize with one of the following options
            case 'ramp'
                initParam = [.4 linspace(1,0,NumPresses)];             % start from RAMP
            case 'exp'
                initParam = [.45 exp(-([1:NumPresses]-1)./4)];     % start from EXP
            case 'logistic'
                Xdomain = [-floor(NumPresses/2):30];        % start from logistic
                PF = 1./(1+1*exp(Xdomain));
                PF = PF(1:NumPresses);%.*[ones(1,3) zeros(1,NumPresses-3)];
                initParam = [.51 , PF];
            case 'box_exp'
                initParam = [.51 exp(-([1:NumPresses]-1)./4).*([ones(1,4) zeros(1,NumPresses-4)])];     % start from Box-EXP
            case 'uniform'
                initParam = [.5 ones(1,NumPresses)];     % start from all equal uniform
            case 'random'
                initParam = [.1 rand(1,NumPresses)];     % start from all random
        end
end
if ~isempty(input_initalParam)
    initParam = input_initalParam;
end
if ~isempty(input_parName)
    parName = input_parName;
end
% baseDir = '/Users/nedakordjazi/Documents/GitHub/SequenceLearningModel/';
baseDir = '/Users/nkordjazi/Documents/GitHub/SequenceLearningModel/';
switch what
    case 'stepwiseWindowPlan'
        windo = {1 2 3 4 5};
        MsetField = {'planFunc(1)'    [1]};
        for h = 1:length(windo)
            NameExt{h} = ['MTRT_H' , num2str(windo{h}(1))];
        end
        %% get bound and theta stim from w  = 1;
        input_parName = {'bAll'    'theta_stim' };%    'planFunc(3)'    'planFunc(4)'    'planFunc(5)'};
        input_initalParam = [0.602	0.00857];%	0.45	0.20	0.0888];
%         input_parName = {'bAll'   'planFunc(2)'   'planFunc(3)'    'planFunc(4)'    'planFunc(5)'};
%         input_initalParam = [0.7000  0.7000    0.5000    0.3000    0.1000];
        slm_NoiselessFitModel('FitMTRT' , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [1],'input_initalParam' , input_initalParam,...
                'input_parName' , input_parName,'NameExt' , 'MTall_H1','loBound' , loBound , 'hiBound',hiBound,'includeRT',1,...
                'MsetField' ,MsetField,'optimizeIPINumber' , [1],'includeMT' , 1 , 'diffMT' , 0);
        slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary','NameExt' , 'MTall_H1' , 'Horizon' , [1:5],'MsetField' ,MsetField)
        %% fit all Windows for  MT using H=1:5 and initial conditions
        load([baseDir ,'arbitraryMTall_H1_4  5/param_arbitraryMTall_H1_4  5.mat'])
        MsetField = {'planFunc(1)'    [1]};
        for pp = 1:2
            MsetField = [MsetField , param.parName{end , pp} ,  {param.par(end,pp)}];
        end 
        
        input_parName = {'planFunc(2)'};%    'planFunc(3)'    'planFunc(4)'    'planFunc(5)'};
        input_initalParam = [.7];%	0.45	0.20	0.0888];
%         input_parName = {'bAll'   'planFunc(2)'   'planFunc(3)'    'planFunc(4)'    'planFunc(5)'};
%         input_initalParam = [0.7000  0.7000    0.5000    0.3000    0.1000];
        slm_NoiselessFitModel('FitMTRT' , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [1:2],'input_initalParam' , input_initalParam,...
                'input_parName' , input_parName,'NameExt' , 'MTall_H1-2','loBound' , loBound , 'hiBound',hiBound,'includeRT',1,...
                'MsetField' ,MsetField,'optimizeIPINumber' , [1],'includeMT' , 1 , 'diffMT' , 1);
        slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary','NameExt' , 'MTall_H1-2' , 'Horizon' , [1:5],'MsetField' ,MsetField)
        %%
        load([baseDir ,'arbitraryMTall_Hall1-2_4  5/param_arbitraryMTall_Hall1-2_4  5.mat'])
        MsetField = [MsetField , param.parName{end , end} ,  {param.par(end,end)}];
         MsetField = {'planFunc(1)'    [1]};
        
        input_parName = {'bAll'   'theta_stim'   'planFunc(2)' 'planFunc(3)' 'planFunc(4)' 'planFunc(5)'};
        input_initalParam = [0.55	0.0086	0.7	0.47	0.32	0.05];%	0.20	0.0888];
%         input_parName = {'bAll'   'planFunc(2)'   'planFunc(3)'    'planFunc(4)'    'planFunc(5)'};
%         input_initalParam = [0.7000  0.7000    0.5000    0.3000    0.1000];
        slm_NoiselessFitModel('FitIPIRT' , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [1:5],'input_initalParam' , input_initalParam,...
                'input_parName' , input_parName,'NameExt' , 'MTHall_H','loBound' , loBound , 'hiBound',hiBound,'includeRT',1,...
                'MsetField' ,MsetField,'optimizeIPINumber' , [1],'includeMT' , 1 , 'diffMT' , 1);
        slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary','NameExt' , 'MTall_H2-3' , 'Horizon' , [1:5],'MsetField' ,MsetField)
        %%
        load([baseDir ,'arbitraryMTall_H2-3_4  5/param_arbitraryMTall_H2-3_4  5.mat'])
        
        input_parName = {'planFunc(2)' 'planFunc(3)' 'planFunc(4)' 'planFunc(5)'};
        input_initalParam = param.par(end,:);%	0.20	0.0888];
%         input_parName = {'bAll'   'planFunc(2)'   'planFunc(3)'    'planFunc(4)'    'planFunc(5)'};
%         input_initalParam = [0.7000  0.7000    0.5000    0.3000    0.1000];
        slm_NoiselessFitModel('FitMTRT' , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [2 5],'input_initalParam' , input_initalParam,...
                'input_parName' , input_parName,'NameExt' , 'MTdiff_H2-3','loBound' , loBound , 'hiBound',hiBound,'includeRT',1,...
                'MsetField' ,MsetField,'optimizeIPINumber' , [1],'includeMT' , 1 , 'diffMT' , 1);
         slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary','NameExt' , 'MTdiff_H2-3' , 'Horizon' , [1:5],'MsetField' ,MsetField)
        %%
        %% fit W = 1 ,2  for IPIs 1 and MT
        h = [1 2];
        W  = 'FitIPIRT';
        load([baseDir ,'arbitraryMTall_Hall1_4  5/param_arbitraryMTall_Hall1_4  5.mat'])
%         input_parName = {'bAll' , 'theta_stim' ,'Aintegrate','planFunc(2)'};
%         input_initalParam = [0.52000    0.0083    0.9852    0.7000];
        input_parName = {'bAll' ,'planFunc(2)'};
        input_initalParam = [0.52000    0.7000];
        slm_NoiselessFitModel(W , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [1 2],'input_initalParam' , input_initalParam,...
            'input_parName' , input_parName,'NameExt' , 'IPI_H1-2','loBound' , loBound , 'hiBound',hiBound,'includeRT',1,...
            'optimizeIPINumber' , [1:2],'includeMT' , 1 , 'MsetField' ,MsetField);
        slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary','NameExt' , 'IPI_H1-2' , 'Horizon' , [1:5],'MsetField' ,MsetField)
        %% fit W = 3 : 5  for IPIs 1:3 using H=1,2 and initial conditions
        load([baseDir ,'arbitraryIPI_H1-2_4  5/param_arbitraryIPI_H1-2_4  5.mat'])
        MsetField = {'planFunc(1)'    [1]};
        for pp = 1:3
            MsetField = [MsetField , param.parName{end , pp} ,  {param.par(end,pp)}];
        end   
        
        input_initalParam = [];
        for h   = 3:5
            input_parName = {};
            input_initalParam = [];
            for hi = 2:h
                input_parName =  [input_parName {['planFunc(' , num2str(hi) , ')']}];
            end
            if h >3
                load([baseDir , 'arbitrary',NameExt{h-1},'_4  5/param_arbitrary',NameExt{h-1},'_4  5.mat'])
                input_initalParam = param.par(end,:);
                input_initalParam = [input_initalParam-.01 max(input_initalParam(end)-.1 , 0)];
            else
                load([baseDir ,'arbitraryIPI_H1-2_4  5/param_arbitraryIPI_H1-2_4  5.mat'])
                input_initalParam = param.par(end,end); %planFunc(2)
                input_initalParam = [input_initalParam .5];
            end
            slm_NoiselessFitModel(W , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [2:5],'input_initalParam' , input_initalParam,...
                'input_parName' , input_parName,'NameExt' , NameExt{h},'loBound' , loBound , 'hiBound',hiBound,'includeRT',1,...
                'optimizeIPINumber' , [1 2],'includeMT' , 1 ,'MsetField' ,MsetField);
            slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary','NameExt' , NameExt{h} , 'Horizon' , [1:windo{h}],'MsetField' ,MsetField)
        end
        %% fit all Windows for IPIs 1:3 5:9 using H=1:5 and initial conditions
        MsetField = {'planFunc(1)'    [1]};
        
        load([baseDir , 'arbitrary',NameExt{h},'_4  5/param_arbitrary',NameExt{h},'_4  5.mat'])
        input_initalParam = param.par(end,:)-.01;
        input_parName     = param.parName(end , :);
        load([baseDir ,'arbitraryIPI_H1-2_4  5/param_arbitraryIPI_H1-2_4  5.mat'])
        input_initalParam = [input_initalParam, param.par(end,1:3)];
        input_parName     = [input_parName    , param.parName(end,1:3)];
        
        slm_NoiselessFitModel(W , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [1:5],'input_initalParam' , input_initalParam,...
                'input_parName' , input_parName,'NameExt' , 'IPI_Hall','loBound' , loBound , 'hiBound',hiBound,'includeRT',1,...
                'optimizeIPINumber' ,[1:3 5:9],'MsetField' ,MsetField);
        slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary','NameExt' , 'IPI_Hall' , 'Horizon' , [1:5],'MsetField' ,MsetField)
        %% fit all Windows for IPIs 1 and MT using H=1:5 and initial conditions
        load([baseDir ,'arbitraryIPI_Hall_4  5/param_arbitraryIPI_Hall_4  5.mat']) 
        input_parName = {'bAll'    'theta_stim'    'Aintegrate'    'planFunc(2)'    'planFunc(3)'    'planFunc(4)'    'planFunc(5)'};
        input_initalParam = [0.6000    0.0082    0.9862    0.5000    0.5000    0.5000    0.5000];
%         for pp = 1:size(param.par , 2)
%             input_parName = [input_parName , param.parName{end , pp} ];
%             input_initalParam = [input_initalParam ,  param.par(end,pp)];
%         end
        slm_NoiselessFitModel(W , Dall , 'planFunc' , 'arbitrary', 'Horizon' , [1:5],'input_initalParam' , input_initalParam,...
                'input_parName' , input_parName,'NameExt' , 'IPIMTall_Hall','loBound' , loBound , 'hiBound',hiBound,'includeRT',1,...
                'MsetField' ,MsetField,'optimizeIPINumber' , [1],'includeMT' , 1);
        slm_NoiselessFitModel('Simulate' , Dall , 'planFunc' , 'arbitrary','NameExt' , 'IPIMTall_Hall' , 'Horizon' , [1:5],'MsetField' ,MsetField)
        
    case 'FitIPIRT'
        %% STEP 1 - fit the ball and the planning function parametrs to get MT
        MSF = {'PlanningCurve' , planFunc  ,'theta_stim' ,0.0084,'Aintegrate' , 1};
        MSF = [MSF , MsetField];
        
        if includeMT
            if diffMT
                optim = {'IPI' , 'diffMT'};
            else
                optim = {'diffIPI' , 'MT'};
            end
        else
            optim= {'IPI'};
        end
        [Param Fval] = slm_optimize(Dall ,  initParam , 'parName' , parName,'runNum' ,['_',planFunc,NameExt,'_',num2str(day)],...
            'Horizon' , Horizon , 'noise' , 0 ,  'subjNum' , [1:15] , 'desiredField' , optim ,'MsetField' , MSF ,...
            'NumPresses' , NumPresses,'loBound' , loBound , 'hiBound',hiBound,'optimizeIPINumber',optimizeIPINumber , 'diffMT' , diffMT);
        if includeRT
            slm_NoiselessFitModel('FitRT' , Dall , varargin);
        end
    case 'FitMTRT'
        if diffMT 
            optim = {'diffMT'};
        else
            optim = {'MT'};
        end
        %% STEP 1 - fit the ball and the planning function parametrs to get MT
        MSF = {'PlanningCurve' , planFunc  ,'theta_stim' ,0.0084,'Aintegrate' ,1};
        MSF = [MSF , MsetField];
        [Param Fval] = slm_optimize(Dall ,  initParam , 'parName' , parName,'runNum' ,['_',planFunc,NameExt,'_',num2str(day)],...
            'Horizon' , Horizon , 'noise' , 0 ,  'subjNum' , [1:15] , 'desiredField' , optim ,'MsetField' ,...
            MSF , 'NumPresses' , NumPresses,'loBound' , loBound , 'hiBound',hiBound, 'diffMT' , diffMT);
        if includeRT
            slm_NoiselessFitModel('FitRT' , Dall , varargin);
        end
    case 'FitRT'
        %% STEP 2 - with parameters of STEP 1 fit the initial decision boundary to get the RTs for every window size
        Horizon = [1:5];
        MSF = {'PlanningCurve' , planFunc  ,'theta_stim' ,0.0084,'Aintegrate' ,1};
        MSF = [MSF , MsetField];
        saveDir = [planFunc,NameExt,'_',num2str(day)];
        cd([baseDir , saveDir]);
        parName = { 'bInit'};
        load(['param_',planFunc,NameExt,'_',num2str(day),'.mat'])
        parcount = 1;
        for p = 1:size(param.parName,2)
            if strcmp(param.parName{end,p} , 'planFunc')
                MSF = [MSF , param.parName{end,p},param.par(end,parcount:parcount+(NumPresses-1))];
                parcount = parcount  +NumPresses;
            else
                MSF = [MSF , param.parName{end,p},param.par(end,parcount)];
                parcount = parcount  +1;
            end
        end
        
        for  h = 1:length(Horizon)
            [Param Fval] = slm_optimize(Dall ,  .49 , 'parName' , parName,'runNum' ,['_',planFunc,NameExt,'Binit_',num2str(h),'_',num2str(day)],...
                'samNum'  , [5] ,'Horizon' , [Horizon(h)] ,'noise' , 0 ,  'subjNum' , [1:15] , 'desiredField' , {'RT'} ,...
                'MsetField' , MSF ,'saveDir' , saveDir, 'NumPresses' , NumPresses,'loBound' , loBound , 'hiBound',hiBound);
        end
    case 'Simulate'
        %% STEP 3 - create the noise-free simulation
        saveDir = [planFunc,NameExt,'_',num2str(day)];
        cd([baseDir , saveDir]);
        MSF = {'PlanningCurve' , planFunc  ,'theta_stim' ,0.0084,'Aintegrate' , 1};
        MSF = [MSF , MsetField];
        load(['param_',planFunc,NameExt,'_',num2str(day),'.mat'])
        parcount = 1;
        for p = 1:size(param.parName,2)
            if strcmp(param.parName{end,p} , 'planFunc')
                MSF = [MSF , param.parName{end,p},param.par(end,parcount:parcount+(NumPresses-1))];
                parcount = parcount  +NumPresses;
            else
                MSF = [MSF , param.parName{end,p},param.par(end,parcount)];
                parcount = parcount  +1;
            end
        end
        AllR = [];
        for  h = 1:length(Horizon)
            clear  param
            fname = ['param_',planFunc,NameExt,'Binit_',num2str(h),'_',num2str(day) , '.mat'];
            load(fname)
            parName = param.parName(end,:);
            par = param.par(end , :);
            [R] = slm_optimSimulate(Dall , par  , 'parName' , parName,'samNum'  , 100 ,...
                'Day' , day, 'Horizon' , [Horizon(h)] , 'poolHorizons' , [5:13] , 'noise' ,0, 'subjNum' , [1:15],'MsetField' , MSF, 'NumPresses' , NumPresses);
            AllR = addstruct(AllR , R);
        end
        %% STEP 4 - visualize
        view = 1; % 1 if you just want to view so everything is plotted in one figure - 0 for publications, separate figures
        c1 = [255, 153, 179]/255; % Random Red Tones
        ce = [153, 0, 51]/255;
        for rgb = 1:3
            tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 6);
        end
        for i = 1:length(tempcol)
            colz{i,1} = tempcol(: , i)';
        end
        
        clear tempcol
        c1 = [153, 194, 255]/255; % structured blue tones
        ce = [0, 0, 102]/255;
        for rgb = 1:3
            tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 6);
        end
        for i = 1:length(tempcol)
            colz{i,2} = tempcol(: , i)';
            avgCol{i} = mean([colz{i,2} ; colz{i,1}],1);
        end
        
        R_seq = AllR;
        Dall.RT = Dall.AllPressTimes(:,1)-1500;
        A = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day , [4 5]) &...
            ismember(Dall.Horizon , [1:5]));
        figure('color' , 'white')
        if view
            subplot(231)
        else
            subplot(121)
        end
        % MT
        hold on
        plot(R_seq.MT , 'o-', 'color' , [0 0 1] )
        lineplot(A.Horizon  , A.MT ,  'plotfcn','nanmedian' ,'style_thickline',...
            'linecolor' ,  'm','errorcolor' , 'm')
        title('MT')
        set(gca,'FontSize' , 18,'GridAlpha' , .2 , 'Box' , 'off','YLim' , [3000 7000] , 'YTick' , [3000:1000: 6000],...
            'YTickLabel' , [3:6] , 'XTick', [1:5] ,'XTickLabel' ,{'1' , '2' , '3' , '4' , '5 - 13'} );
        ylabel('Execution time [s]','FontSize' , 20)
        xlabel('Viewing window size (W)' ,'FontSize' , 20)
        
        % RT
        if view
            subplot(232)
        else
            subplot(122)
        end
        hold on
        plot(R_seq.RT , 'o-')
        lineplot(A.Horizon  , A.RT ,  'plotfcn','nanmedian' , ...
            'style_thickline','linecolor' ,  'm','errorcolor' , 'm')
        title('RT')
        set(gca,'FontSize' , 18,'GridAlpha' , .2 , 'Box' , 'off',...
            'YLim' , [400 800] , 'YTick' , [400:100: 800],...
            'YTickLabel' , [0.4:.1:0.8] , 'XTick', [1:5] ,'XTickLabel' ,{'1' , '2' , '3' , '4' , '5 - 13'} );
        ylabel('Initial reaction time [s]','FontSize' , 20)
        xlabel('Viewing window size (W)' ,'FontSize' , 20)
        
        % IPI
        
        Fit.IPI = AllR.IPI;
        Fit.IPI = reshape(Fit.IPI , numel(Fit.IPI) , 1);
        Act.IPI = A.IPI;
        Act.IPI = reshape(Act.IPI , numel(Act.IPI) , 1);
        
        Fit.singleH  = repmat(nanmean(AllR.Horizon, 2), 1 , size(A.AllPress,2)-1);
        Fit.singleH  = reshape(Fit.singleH , numel(Fit.IPI) , 1);
        Act.singleH  = repmat(nanmean(A.Horizon, 2) , 1 , size(A.AllPress,2)-1);
        Act.singleH  = reshape(Act.singleH , numel(Act.IPI) , 1);
        
        
        Fit.ipiNum = repmat(1:size(AllR.stimulus,2)-1 , size(AllR.stimulus,1) , 1);
        Fit.ipiNum = reshape(Fit.ipiNum , numel(Fit.ipiNum) , 1);
        Act.ipiNum = repmat(1:size(AllR.stimulus,2)-1 , length(A.AllPress) , 1);
        Act.ipiNum = reshape(Act.ipiNum , numel(Act.ipiNum) , 1);
        
        Act.fitoract = ones(size(Act.ipiNum));
        Fit.fitoract = zeros(size(Fit.ipiNum));
        
        All  = addstruct(Fit , Act);
        colorz = colz(:,1);
        if view
            subplot(234)
        else
            figure('color' , 'white')
            subplot(121)
        end
        
        H = unique(All.singleH);
        for h= 1:length(H)
            A = getrow(Fit , Fit.singleH==H(h));
            plot(A.ipiNum  ,A.IPI, '-o' , 'Color' , colorz{h} , ...
                'MarkerEdgeColor' , colorz{h} , 'MarkerFaceColor' , colorz{h},...
                'LineWidth' , 1.5 , 'MarkerSize' , 5)
            hold on
        end
        legend({'W = 1' , 'W = 2' , 'W = 3' , 'W = 4' , 'W = 5-13'} , 'Box' , 'off')
        xlabel('IPI number')
        ylabel('Inter-press interval time [s]')
        set(gca , 'FontSize' , 16 , 'Box' , 'off' , 'YLim' , [150 650],'YTick' , [200:100: 600],...
            'YTickLabel' , [0.2:.1:0.6] )
        
        if view
            subplot(235)
        else
            subplot(122)
        end
        colorz = colz(:,2);
        lineplot(All.ipiNum , All.IPI , 'plotfcn' , 'nanmedian',...
            'split', All.singleH  , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
            'linewidth' , 1.5 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 5, 'markercolor' , colorz , 'leg' , {'W = 1' , 'W = 2' , 'W = 3' , 'W = 4' , 'W = 5-13'} , ...
            'subset' , All.fitoract == 1);
        
        ylabel('Inter-press interval time [s]')
        xlabel('IPIs number')
        set(gca , 'FontSize' , 16 , 'Box' , 'off' , 'YLim' , [150 650],'YTick' , [200:100: 600],...
            'YTickLabel' , [0.2:.1:0.6] )
        
        if view
            subplot(2,3,[3 6])
        else
            figure('color' , 'white')
        end
        plot(R_seq.planFunc(1,:), '-o' ,'LineWidth' , 2 , 'MarkerSize' , 5)
        set(gca , 'XLim' , [1 size(R_seq.stimulus,2)] , 'Box' , 'off')
        title([planFunc , ' planning function'])
end
