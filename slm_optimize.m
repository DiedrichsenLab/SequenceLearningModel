function [Param Fval] = slm_optimize(Dall , initParam , varargin)
cycNum = 50;
samNum = 20;
tol = [0.5 0.03];
ItrNum = 20;
NumPresses = size(Dall.AllPress , 2);
c = 1;
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
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end
mainDir = '/Users/nkordjazi/Documents/GitHub/';
% mainDir = '/Users/nedakordjazi/Documents/GitHub/';
%% optimization
Dall = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day , Day) &...
    ismember(Dall.Horizon , Horizon) & ismember(Dall.SN , subjNum));
if ~isempty(poolHorizons)
    Dall.Horizon(ismember(Dall.Horizon , poolHorizons)) = poolHorizons(1);
end


Horizon = unique(Dall.Horizon);
for i = 1:cycNum
    %% in cycle numbers bigger than 1, retrieve the last param set from the last cycle and use as initials
    disp(['Initializing optimization cycle number ' , num2str(i) , '/', num2str(cycNum) , ' with ' , num2str(ItrNum) , ' iterations...'])
    if i>1
        load([mainDir , 'SequenceLearningModel/param' , runNum , '.mat'])
        initParam = param.par(end , :);
    end
    ANA = getrow(Dall , ismember(Dall.Horizon , Horizon));
    ANA.RT = ANA.AllPressTimes(:,1)-1500;
    %% set up the desired output
    G = tapply(ANA , {'Horizon'} , {'MT' , 'nanmedian'},{'RT' , 'nanmedian'},{'IPI' , 'nanmedian'});
    x_desired = [];
    for xd = 1:length(desiredField)
        eval(['x_desired = [x_desired ; G.' , desiredField{xd} , '];'] )
    end
    x_desired = x_desired';
    %% subsample the data
    A = [];
    if ~noise
        samNum = 1;
    end
    A = [];
    for h = 1:length(Horizon)
        N = getrow(ANA , ANA.Horizon == Horizon(h)); % when the noise is off all the trials will turn out identical
        if ~isempty(samNum)
            N =  getrow(N , randperm(length(N.TN) , samNum));
        end
        A = addstruct(A , N);
    end
    ANA = A;
    %% set up the otimization handle object
    model = @(param,T, M, opts) slm_optimSimTrial(param , T ,M,opts); % Model Function
    
    %% set up the inputs to the model funtion T , M
    SeqLength = NumPresses;
    T.TN = ANA.TN;
    T.Horizon = ANA.Horizon;
    T.numPress = SeqLength.*ones(length(ANA.TN) , 1);
    T.stimTime = zeros(length(ANA.TN) , SeqLength);
    if exist('stimulus')
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
    M.numOptions = 5;
    M.dT_visual  = 90;
    M.Ainhibit   = 0;
    M.Capacity   =1;
    M.dT_motor   = 150;
    M.dtGrowth = 1;
    M.TSDecayParam  = 7.75;
    M.Aintegrate  = 0.98;
    M.bAll = 0.6;
    M.Bound = M.bAll.*ones(1,size(T.stimulus , 2)); % boundry is a vector of length maxPresses
    M.PlanningCurve = 'exp'; % other options: 'logistic', 'box' , 'ramp'
    M.DecayParam   = 7; % the decay constant for the 'exp' option of PlanningCurve
    M.B_coef = 1;       % for the 'logistic' option of PlanningCurve
    M.Box = 1;          % box size for the 'boxcar' option of PlanningCurve
    M.rampDecay = size(T.stimulus , 2);   % number of steps between 1 and 0 for the 'ramp' option of PlanningCurve
    M.theta_stim = 0.01;
    if~noise
        M.SigEps      = 0;
    else
        M.SigEps      = 0.02;
    end
    M.parName = parName;
    c = 1;
    %% re-set the fields that have been defined in input
    while(c<=length(MsetField))
        eval(['M.',MsetField{c} '= MsetField{c+1};']);
        c=c+2;
    end
    
    %% set up the cost function
    opts.runNum       = runNum;
    opts.cycNum       = i;
    opts.mode         = 'optim';
    opts.desiredField = desiredField;
    OLS = @(param) nansum(nansum((model(param,T, M ,opts) - x_desired).^2));
    %% optimization
    opts = optimset('MaxIter', ItrNum ,'TolFun',1e+03,'Display','iter' , 'TolX' , 1e-4);
    if isempty(loBound) | isempty(hiBound)
        [Param Fval] = fminsearch(OLS,initParam,opts);
    else
        [Param Fval] = fminsearchbnd(OLS,initParam,loBound,hiBound, opts);
    end
end


