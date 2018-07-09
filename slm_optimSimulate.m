function [R] = slm_optimSimulate(Dall, par  , varargin)


samNum = 20;
tol = [0.5 0.03];
c = 1;
NumPresses = size(Dall.AllPress , 2);
while(c<=length(varargin))
    switch(varargin{c})
        case {'parName'}
            % names of paramters
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

if ~isempty(poolHorizons)
    Dall.Horizon(ismember(Dall.Horizon , poolHorizons)) = poolHorizons(1);
end

Dall = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day , Day) &...
    ismember(Dall.Horizon , Horizon) & ismember(Dall.SN , subjNum));


ANA = getrow(Dall , ismember(Dall.Horizon , Horizon));
ANA.RT = ANA.AllPressTimes(:,1)-1500;
Horizon = unique(Dall.Horizon);
%% subsample the data
A = [];
if ~noise
    samNum = 1;
end
A = [];

for h = 1:length(Horizon)
    N = getrow(ANA , ANA.Horizon == Horizon(h)); % when the noise is off all the trials will turn out identical
    idx = randi([1 length(N.TN)],1,samNum);
    N =  getrow(N , idx);
    A = addstruct(A , N);
end
ANA = A;
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
M.numOptions    = 5;
M.dT_visual     = 100;
M.Ainhibit      = 0;
M.Capacity      = 1;
M.dT_motor      = 120;
M.dtGrowth      = 1;
M.TSDecayParam  = 3;
M.Aintegrate    = 0.98;
M.bAll          = 0.5;     % press boundary for 5 window sizes
M.PlanningCurve = 'exp'; % other options: 'logistic', 'box' , 'ramp'
M.DecayParam    = 7; % the decay constant for the 'exp' option of PlanningCurve
M.B_coef1       = 1;       % for the 'logistic' option of PlanningCurve
M.B_coef2       = 0;       % for the 'logistic' option of PlanningCurve
M.Box           = 1;          % box size for the 'boxcar' option of PlanningCurve
M.rampDecay1     = size(T.stimulus , 2);   % for the 'ramp' option of PlanningCurve
M.rampDecay2     = 0;   % for the 'ramp' option of PlanningCurve
M.theta_stim    = 0.0001;
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

%% re-set the fields that have been defined in input
c = 1;
for c = 1:length(parName)
    eval(['M.',parName{c} '= par(c);']);
end
c = 1;
while(c<=length(MsetField))
    eval(['M.',MsetField{c} '= MsetField{c+1};']);
    c=c+2;
end
M.bInit         = M.bAll;  % initial bound for 5 window sizes
%% simulate
opts.runNum       = [];
opts.cycNum       = [];
opts.mode         = 'sim';
opts.desiredField = [];

R = slm_optimSimTrial(par , T , M , opts);

end
