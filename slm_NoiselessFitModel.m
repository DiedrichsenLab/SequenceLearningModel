function slm_NoiselessFitModel(what , Dall , varargin)
planFunc = 'exp';
initParam = [.55 , 7];
parName = { 'bAll' , 'DecayParam'};
NumPresses = [];
stimulus = [];
ItrNum = 1000;
cycNum = 1;
day = [4 5]; % just fitting days 4, 5
c = 1;
NameExt = [];
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
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end

switch planFunc
    case 'logistic'
        parName = { 'bAll' , 'B_coef1' 'B_coef2'};
        initParam = [.52 3 1];
    case 'ramp'
        parName = { 'bAll' , 'rampDecay' };
        initParam = [.50 , 5 4];
    case 'box'
        parName = { 'bAll' , 'Box'};
        initParam = [.55 , 5];
    case 'box_logistic'
        parName = { 'bAll' , 'B_coef1' 'B_coef2'  'Box'};
        initParam = [.4 , 2 3 4];
    case 'box_exp'
        parName = { 'bAll' , 'DecayParam'  'Box'};
        initParam = [.5 , 7 5];
    case 'box_ramp'
        parName = { 'bAll' , 'rampDecay' ,'Box'};
        initParam = [.50 , 7  3];
        
end
baseDir = '/Users/nedakordjazi/Documents/GitHub/SequenceLearningModel/';
% baseDir = '/Users/nkordjazi/Documents/GitHub/SequenceLearningModel/';
switch what
    case 'Fit'
        %% STEP 1 - fit the ball and the planning function parametrs to get MT
        MSF = {'PlanningCurve' , planFunc  ,'theta_stim' ,0.0084,'Aintegrate' , 0.985};
        [Param Fval] = slm_optimize(Dall ,  initParam , 'parName' , parName,'runNum' ,['_',planFunc,NameExt,'_',num2str(day)],...
            'Horizon' , [1:5] , 'noise' , 0 ,  'subjNum' , [1:15] , 'desiredField' , {'MT'} ,'MsetField' , MSF);
        
        %% STEP 2 - with parameters of STEP 1 fit the initial decision boundary to get the RTs for every window size
        saveDir = [planFunc,NameExt,'_',num2str(day)];
        cd([baseDir , saveDir]);
        parName = { 'bInit'};
        load(['param_',planFunc,NameExt,'_',num2str(day),'.mat'])
        for p = 1:size(param.par,2)
            MSF = [MSF , param.parName{end,p},param.par(end,p)];
        end
        for  h = 1:5
            [Param Fval] = slm_optimize(Dall ,  initParam , 'parName' , parName,'runNum' ,['_',planFunc,NameExt,'Binit_',num2str(h),'_',num2str(day)],...
                'samNum'  , [5] ,'Horizon' , [h] ,'noise' , 0 ,  'subjNum' , [1:15] , 'desiredField' , {'RT'} ,  'MsetField' , MSF ,...
                'saveDir' , saveDir);
        end
    case 'Simulate'
        %% STEP 3 - create the noise-free simulation
        saveDir = [planFunc,NameExt,'_',num2str(day)];
        cd([baseDir , saveDir]);
        MSF = {'PlanningCurve' , planFunc  ,'theta_stim' ,0.0084,'Aintegrate' , 0.985};
        load(['param_',planFunc,NameExt,'_',num2str(day),'.mat'])
        for p = 1:size(param.par,2)
            MSF = [MSF , param.parName{end,p},param.par(end,p)];
        end
        AllR = [];
        for  h = 1:5
            clear  param
            fname = ['param_',planFunc,NameExt,'Binit_',num2str(h),'_',num2str(day) , '.mat'];
            load(fname)
            parName = param.parName(end,:); 
            par = param.par(end , :);
            [R] = slm_optimSimulate(Dall , par  , 'parName' , parName,'samNum'  , 100 ,...
                'Day' , day, 'Horizon' , [h] , 'poolHorizons' , [5:13] , 'noise' ,0, 'subjNum' , [1:15],'MsetField' , MSF);
            AllR = addstruct(AllR , R);
        end
        %% STEP 4 - visualize
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
        subplot(231)
        % MT
        hold on
        plot(R_seq.MT , 'o-', 'color' , [0 0 1] )
        lineplot(A.Horizon  , A.MT ,  'plotfcn','nanmedian' , 'linecolor' ,  [0 1 0 ],'errorcolor' , [0 1 0])
        title('MT')
        
        % RT
        subplot(234)
        hold on
        plot(R_seq.RT , 'o-', 'color' , [0 0 1] )
        lineplot(A.Horizon  , A.RT ,  'plotfcn','nanmedian' , 'linecolor' ,  [0 1 0 ],'errorcolor' , [0 1 0])
        title('RT')
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
        subplot(232)
        
        H = unique(All.singleH);
        for h= 1:length(H)
            A = getrow(Fit , Fit.singleH==H(h));
            plot(A.ipiNum  ,A.IPI, '-o' , 'Color' , colorz{h} , ...
                'MarkerEdgeColor' , colorz{h} , 'MarkerFaceColor' , colorz{h},...
                'LineWidth' , 1.5 , 'MarkerSize' , 5)
            hold on
        end
        legend({'H = 1' , 'H = 2' , 'H = 3' , 'H = 4' , 'H = 5-13'} , 'Box' , 'off')
        title('IPIs - fitted')
        xlabel('IPIs number')
        grid on
        set(gca , 'FontSize' , 16 , 'Box' , 'off' , 'YLim' , [150 700])
        
        subplot(235)
        colorz = colz(:,2);
        lineplot(All.ipiNum , All.IPI , 'plotfcn' , 'nanmedian',...
            'split', All.singleH  , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
            'linewidth' , 1.5 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 5, 'markercolor' , colorz , 'leg' , {'H = 1' , 'H = 2' , 'H = 3' , 'H = 4' , 'H = 5-13'} , ...
            'subset' , All.fitoract == 1);
        
        title('IPIs - Actual')
        xlabel('IPIs number')
        grid on
        set(gca , 'FontSize' , 16 , 'YLim' , [150 700])
        
        subplot(2,3,[3 6])
        plot(R_seq.planFunc(1,:), '-o' ,'LineWidth' , 1.5 , 'MarkerSize' , 5)
        set(gca , 'XLim' , [1 size(R_seq.stimulus,2)] , 'Box' , 'off')
        title([planFunc , ' planning function'])
end
