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


Dall = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day , Day) & ...
    ismember(Dall.Horizon , Horizon) & ismember(Dall.SN , subjNum));
if ~isempty(poolHorizons)
    Dall.Horizon(ismember(Dall.Horizon , poolHorizons)) = poolHorizons(1);
end


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
    if ~isempty(samNum)
        N =  getrow(N , randperm(length(N.TN) , samNum));
    end
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
M.dT_visual     = 90;
M.Ainhibit      = 0;
M.Capacity      = 1;
M.dT_motor      = 150;
M.dtGrowth      = 1;
M.TSDecayParam  = 3;
M.Aintegrate    = 0.98;
M.bAll          = 0.5;     % press boundary for 5 window sizes
M.bInit         = M.bAll; % initial bound for 5 window sizes
M.PlanningCurve = 'exp'; % other options: 'logistic', 'box' , 'ramp'
M.DecayParam    = 7; % the decay constant for the 'exp' option of PlanningCurve
M.B_coef1       = 1;       % for the 'logistic' option of PlanningCurve
M.B_coef2       = 0;       % for the 'logistic' option of PlanningCurve
M.Box           = 1;          % box size for the 'boxcar' option of PlanningCurve
M.rampDecay1     = size(T.stimulus , 2);   % for the 'ramp' option of PlanningCurve
M.rampDecay2     = 0;   % for the 'ramp' option of PlanningCurve
M.theta_stim    = 0.01;
M.parName       = parName;
M.planFunc = zeros(1,NumPresses); % for the arbitrary option of the PlanningCurve
if~noise
    M.SigEps    = 0;
else
    M.SigEps    = 0.02;
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

%% simulate
opts.runNum       = [];
opts.cycNum       = [];
opts.mode         = 'sim';
opts.desiredField = [];

R = slm_optimSimTrial(par , T , M , opts);

%% Horizon
plt = 0;
if plt
    ANA = getrow(Dall , ismember(Dall.Horizon , Horizon));
    ANA.RT = ANA.AllPressTimes(:,1)-1500;
    C = R;
    R = getrow(R , ~R.isError);
    All = [];
    Act = [];
    Fit = R;
    Act.singleH = ANA.Horizon;
    Fit.singleH  = nanmean(R.Horizon, 2);
    Act.MT = ANA.MT;
    Act.RT = ANA.AllPressTimes(:,1) - 1500;
    
    Act.fitoract = ones(size(Act.MT));
    Fit.fitoract = zeros(size(Fit.MT));
    
    All  = addstruct(Fit , Act);
    
    figure('color' , 'white')
    subplot(211)
    colorz = colz(3,:);{[0.840000000000000,0.360000000000000,0.501176470588235],[0.360000000000000,0.456470588235294,0.760000000000000]};
    lineplot(All.singleH , All.MT , 'plotfcn' , 'nanmedian',...
        'split', All.fitoract  , 'linecolor' , colorz,...
        'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
        'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'});
    title('MT')
    grid on
    set(gca , 'FontSize' , 16, 'YLim' , [2000 9000])
    xlabel('Horizon')
    
    
    subplot(212)
    lineplot(All.singleH , All.RT , 'plotfcn' , 'nanmedian',...
        'split', All.fitoract  , 'linecolor' , colorz,...
        'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
        'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'});
    title('RT')
    xlabel('Horizon')
    grid on
    set(gca , 'FontSize' , 16 , 'YLim' , [300 750])
    
    MTRT = All;
    %% IPI
    All = [];
    Act = [];
    Fit = [];
    clear Fit Act
    Fit.IPI = R.IPI;
    Fit.IPI = reshape(Fit.IPI , numel(Fit.IPI) , 1);
    Act.IPI = ANA.IPI;
    Act.IPI = reshape(Act.IPI , numel(Act.IPI) , 1);
    
    Fit.singleH  = repmat(nanmean(R.Horizon, 2), 1 , size(ANA.AllPress,2)-1);
    Fit.singleH  = reshape(Fit.singleH , numel(Fit.IPI) , 1);
    Act.singleH  = repmat(nanmean(ANA.Horizon, 2) , 1 , size(ANA.AllPress,2)-1);
    Act.singleH  = reshape(Act.singleH , numel(Act.IPI) , 1);
    
    
    Fit.ipiNum = repmat(1:size(T.stimulus,2)-1 , size(T.stimulus,1) , 1);
    Fit.ipiNum = reshape(Fit.ipiNum , numel(Fit.ipiNum) , 1);
    Act.ipiNum = repmat(1:size(T.stimulus,2)-1 , length(ANA.AllPress) , 1);
    Act.ipiNum = reshape(Act.ipiNum , numel(Act.ipiNum) , 1);
    
    Act.fitoract = ones(size(Act.ipiNum));
    Fit.fitoract = zeros(size(Fit.ipiNum));
    
    
    All  = addstruct(Fit , Act);
    colorz = colz(:,1);
    figure('color' , 'white')
    subplot(211)
    if noise
        lineplot(All.ipiNum , All.IPI , 'plotfcn' , 'nanmedian',...
            'split', All.singleH  , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
            'linewidth' , 1.5 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 5, 'markercolor' , colorz , 'leg' , {'H = 1' , 'H = 2' , 'H = 3' , 'H = 4' , 'H = 5-13'} , ...
            'subset' , All.fitoract == 0);
    else
        H = unique(All.singleH);
        for h= 1:length(H)
            A = getrow(Fit , Fit.singleH==H(h));
            plot(A.ipiNum  ,A.IPI, '-o' , 'Color' , colorz{h} , ...
                'MarkerEdgeColor' , colorz{h} , 'MarkerFaceColor' , colorz{h},...
                'LineWidth' , 1.5 , 'MarkerSize' , 5)
            hold on
        end
        legend({'H = 1' , 'H = 2' , 'H = 3' , 'H = 4' , 'H = 5-13'} , 'Box' , 'off')
    end
    title('IPIs - fitted')
    xlabel('IPIs number')
    grid on
    set(gca , 'FontSize' , 16 , 'Box' , 'off' , 'YLim' , [150 700])
    
    subplot(212)
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
    IPI = All;
end
