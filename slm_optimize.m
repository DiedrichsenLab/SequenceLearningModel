function [Param Fval] = slm_optimize(Dall , initParam , varargin)
cycNum = 1;
samNum = 20;
tol = [0.5 0.03];
ItrNum = 1000;
NumPresses = size(Dall.AllPress , 2);
loBound = [];
hiBound = [];
Day = [4 5];
poolHorizons = [5:13];
noisefreeRep = [];
optimizeIPINumber = [1:3];
c = 1;
diffMT = 0;



while(c<=length(varargin))
    switch(varargin{c})
        
        case {'parName'}
            % names of paramters
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'runNum'}
            % the number of run
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'cycNum'}
            % the number of times to sample the data and go through the optmization process
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'samNum'}
            % number of samples to pick from the data - leave [] to get all
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'tol'}
            % tolerance margines for std and mean respectively
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'ItrNum'}
            % number of iterations within a cycle
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
        case {'Day'}
            % what days to include
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Horizon'}
            % what Horizon to include
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'poolHorizons'}
            % [] or the horizons you want to pool together
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'noise'}
            % 1 or 0 -
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'subjNum'}
            % subjects to include in modeling
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'desiredField'}
            % the name of the field in the data to minimize the error on e.g. MT TR IPI...
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'MsetField'}
            % the names and values of the fields we want to set in M
            % has to be cell of value names, followed by their values
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
        case {'noisefreeRep'}
            % the noise free response for situations where we want the
            % noisey response to be around the median of the noise -free
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'saveDir'}
            % Directory to save the mat files to
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'optimizeIPINumber'}
            % IPI numbers to include in the optimization
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'diffMT'}  % for the 'FitIPIRT' case
            % use the difference between MTs instead of MTs
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
         case {'mORi'}
            % m for mcbook, i for iMac
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end
switch mORi
    case 'm'
        mainDir = '/Users/nedakordjazi/Documents/GitHub/SequenceLearningModel/';
    case 'i'
        mainDir = '/Users/nkordjazi/Documents/GitHub/SequenceLearningModel/';
end
if ~exist('saveDir')
    saveDir = [mainDir , runNum(2:end)];
else
    saveDir = [mainDir , saveDir];
end

%% optimization
if ~isempty(poolHorizons)
    Dall.Horizon(ismember(Dall.Horizon , poolHorizons)) = poolHorizons(1);
end

Dall = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day , Day) &...
    ismember(Dall.Horizon , Horizon) & ismember(Dall.SN , subjNum));



Horizon = unique(Dall.Horizon);
for i = 1:cycNum
    %% in cycle numbers bigger than 1, retrieve the last param set from the last cycle and use as initials
    disp(['Initializing optimization cycle number ' , num2str(i) , '/', num2str(cycNum) , ' with ' , num2str(ItrNum) , ' iterations...'])
    if i==1
        mkdir(saveDir);
    else
        load([saveDir , '/param' , runNum , '.mat'])
        initParam = param.par(end , :);
    end
    Dall.RT = Dall.AllPressTimes(:,1)-1500;
    %% set up the desired output
    G = tapply(Dall , {'Horizon'} , {'MT' , 'nanmean'},{'RT' , 'nanmean'},{'IPI' , 'nanmean'});
    G.MT = G.MT;%/NumPresses;
    if diffMT
        G.diffMT = [G.MT ; -diff(10*G.MT,1,1)];
    end
    x_desired = [];
    x_d = [];
    if exist('noisefreeRep') & ~isempty(noisefreeRep)
        for xd = 1:length(desiredField)
            eval(['x_desired = [x_desired ; noisefreeRep.' , desiredField{xd} , '];'] )
        end
    else
        for xd = 1:length(desiredField)
            if sum(strcmp(desiredField{xd} , 'IPI')) | sum(strcmp(desiredField{xd} , 'diffIPI'))
                for hh = 1:length(Horizon)
                    F = getrow(G , G.Horizon==Horizon(hh));
                    eval(['x_d = F.' , desiredField{xd} ,';'] )
                    x_desired = [x_desired  [x_d(optimizeIPINumber)]];% x_desired(12:13)];
                end
            else
                eval(['x_desired = [x_desired  G.' , desiredField{xd} , '''];'] )
            end
        end
        
    end

    x_desired = reshape(x_desired , 1,numel(x_desired));

    %% subsample the data
    A = [];
    if ~noise
        samNum = 1;
    end
    A = [];
    for h = 1:length(Horizon)
        N = getrow(Dall , Dall.Horizon == Horizon(h)); % when the noise is off all the trials will turn out identical
        if ~isempty(samNum)
            N =  getrow(N , randperm(length(N.TN) , samNum));
        end
        A = addstruct(A , N);
    end
    ANA = A;
    %% set up the otimization handle object
    model = @(param,T, M, opts) slm_optimSimTrial(param , T ,M,opts); % Model Function
    
    %% set up the inputs to the model funtion T , M
    if ~isempty(NumPresses)
        SeqLength = NumPresses;
    else
        SeqLength = size(Dall.AllPress , 2);
    end
    T.TN = ANA.TN;
    T.Horizon = ANA.Horizon;
    T.numPress = SeqLength.*ones(length(ANA.TN) , 1);
    T.stimTime = zeros(length(ANA.TN) , SeqLength);
    if exist('stimulus') & ~isempty(stimulus)
        T.stimulus = repmat(stimulus ,length(ANA.TN) , 1); % constant stimulus
    else
        T.stimulus = ANA.AllPress(:,1:SeqLength);
    end
    T.forcedPressTime = nan(length(ANA.TN) , SeqLength);
    T.Horizon =repmat(T.Horizon , 1, SeqLength) .*(ones(length(T.TN),SeqLength));
    for tn = 1:length(T.TN)
        T.Horizon(tn , 1:T.Horizon(tn)) = NaN;
    end
    
    %% set the default M values
    M.numOptions    = 5;
    M.dT_visual     = 100;
    M.Ainhibit      = 0;
    M.Capacity      = 1;
    M.dT_motor      = 120;
    M.dtGrowth      = 1;
    M.TSDecayParam  = 3;
    M.Aintegrate    = 0.98;
    M.bAll          = 0.5;     % press boundary for 5 window sizes
    M.PlanningCurve = 'exp';   % other options: 'logistic', 'box' , 'ramp'
    M.DecayParam    = 7;       % the decay constant for the 'exp' option of PlanningCurve
    M.B_coef1       = 1;       % for the 'logistic' option of PlanningCurve
    M.B_coef2       = 0;       % for the 'logistic' option of PlanningCurve
    M.Box           = 1;       % box size for the 'boxcar' option of PlanningCurve
    M.rampDecay1     = size(T.stimulus , 2);   % for the 'ramp' option of PlanningCurve
    M.rampDecay2     = 0;   % for the 'ramp' option of PlanningCurve
    M.theta_stim    = 0.01;
    M.parName       = parName;
    M.planFunc = zeros(1,NumPresses); % for the arbitrary option of the PlanningCurve
    if~noise
        M.SigEps    = 0;
    else
        if sum(ismember(Day , [4 5]))
            M.SigEps    = [0.0035 0.0045 0.005]; % 0.0035 for window 1  --- 0.0045 for the rest;
        else
            M.SigEps    = [0.0029 0.0036 0.0036];
        end
        
    end
    
    c = 1;
    %% re-set the fields that have been defined in input
    while(c<=length(MsetField))
        eval(['M.',MsetField{c} '= MsetField{c+1};']);
        c=c+2;
    end
    M.bInit         = M.bAll;  % initial bound for 5 window sizes
    %% set up the cost function
    opts.runNum       = runNum;
    opts.cycNum       = i;
    opts.mode         = 'optim';
    opts.desiredField = desiredField;
    opts.saveDir      = saveDir;
    opts.optimizeIPINumber = optimizeIPINumber;
    opts.diffMT = diffMT;
    OLS = @(param) nansum(nansum((model(param,T, M ,opts) - x_desired).^2));
    %% optimization
    opts = optimset('MaxIter', ItrNum ,'TolFun',1,'Display','iter' , 'TolX' , 1e-6);
    if isempty(loBound) | isempty(hiBound)
        [Param Fval] = fminsearch(OLS,initParam,opts);
    else
        [Param Fval] = fminsearchbnd(OLS,initParam,loBound,hiBound, opts);
    end
end


