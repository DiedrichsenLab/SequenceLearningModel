function [T,SIM]=slm_simTrial(M,T,varargin)
%% function [T,SIM]=slm_simTrial(M,T,varargin);
%
% Incoporates horizon size (T.Horizon) as well as buffer size (M.capacity), and preparation time (M.prepTime)
% Multiple planning possible as long as M.capacity = 1, this funcion should work exacly same as [T,SIM]=slm_simTrial(M,T);
% Simulates a trial using a parallel evidence-accumulation model for sequential response tasks

M.capacity = 1; % always set capacity to 1 since it's soft capacity
dT = 2;     % delta-t in ms
maxTime = 4000; % Maximal time for trial simulation
c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'DecayFunc'}
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
            % define the parameters for the decay function
        case {'DecayParam'}
            % for 'exp' this would be the time constant (defaul = 1)
            % for 'linear' this would be a negative slope (default = -1/seqlength)
            % for 'boxcar' this would be the number of 1s in a row (default = 5)
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end

% Determine length of the trial
maxPresses = max(T.numPress);

% number of decision steps is defined by the M.capacity and T.Horizon, whichever results in more decision steps
dec=1:maxPresses;  % Number of decision
maxPlan = M.capacity; % Number of digits being planned in every decision step
if isfield(T,'Horizon')
    % set the stimTime for presses that are not shown at time 0 to NaN.
    % These will be filled with press times
    T.stimTime(~isnan(T.Horizon)) = NaN;
else
    T.Horizon=nan(size(T.stimulus,1),size(T.stimulus,2)); % default to full horizon
end
T.pressTime = nan(size(T.stimulus,1),size(T.stimulus,2)); % to be filled with press times

% initialize variables
X = zeros(M.numOptions,maxTime/dT,maxPresses); % Hidden state
S = zeros(M.numOptions,maxTime/dT,maxPresses); % Stimulus present
t = (1:maxTime/dT)*dT-dT; % Time in ms

% implement forced-RT collapsing decision boundary (logistic decay)
if ~isnan(T.forcedPressTime(1,1))
    a = 0.005;%0.05; % the decay constant: the bigger the faster the decay --> reaches zero faster
    b = T.forcedPressTime + 100; % in ms, the inflexion point
    B = (M.Bound ./ (1 + exp ( a * (t - b) ) ) ) - M.Bound / 2; % Decision Bound - logistic decay
    B(B<=0)=0;
else
    B = ones(1,maxTime/dT)*M.Bound; % Decision Bound - constant value
end

i = 1;                   % Index of simlation
nDecision = 1;           % Current decision to make
%numPresses = 0;         % Number of presses produced
isPressing = 0;          % Is the motor system currently occupied?
remPress = maxPresses;   % Remaining presses. This variable will be useful if/whenever multiple digits are being planned in every decsion step
%collapsedBound = 0;

% Set up parameters
A  = eye(M.numOptions)*(M.Aintegrate)+...
    (~eye(M.numOptions)*M.Ainhibit); % A defined the matrix of autoregressive coefficients

% Use logistic growth
% T.stimTime=T.stimTime+250;
a = .1; %0.09; % the growth constant: the bigger the faster the growth --> reaches Bound faster
b = 400; %200; %400; % in ms, how long it takes for the function to reach max
G = (M.Bound/2) ./ (1 + exp ( -a * (t - (T.stimTime(1)+b/2) ) ) ); % logistic growth

%% Start time-by-time simulation
while remPress && i<maxTime/dT
    mult = zeros(1,length(dec));
    % Update the stimulus: Fixed stimulus time
    indx = find(t(i)>(T.stimTime+M.dT_visual)); % Index of which stimuli are present T.
    % stimuli of greater than Horizon size will be unavailable
    for j=indx
        if T.stimulus(j)>0 && T.stimulus(j)<=5 % in single resp task, some stimuli are distractors
            S(T.stimulus(j),i,j)=1;
        end;
    end
    
    % Update the evidence state
    eps = randn([M.numOptions 1 maxPresses]) * M.SigEps;
    switch DecayFunc
        case 'exp'
            if ~exist('DecayParam','var')
                DecayParam = maxPlan;
            end
            mult=exp(-(dec-nDecision)./DecayParam);      % How much stimulus exponentia decay
        case 'linear'
            if ~exist('DecayParam','var')
                DecayParam = 1/(1-max(dec));
            end
            linDecay=DecayParam*(dec-max(dec));  % How much stimulus linear decay
            mult(nDecision:end) = linDecay(1:max(dec)-nDecision+1);  % How much stimulus linear decay
        case 'boxcar'
            if ~exist('DecayParam','var')
                DecayParam = 5;
            end
            bc = min(max(dec) , nDecision+DecayParam);
            mult(nDecision:bc) = 1;  % How much stimulus linear decay
    end
    mult(dec<nDecision)=0;                    % Made decisions will just decay
    for j =1:maxPresses
        %X(:,i+1,j) = (A*X(:,i,j)) + (M.theta_stim.*mult(j).*S(:,i,j)) + dT*eps(:,1,j);
        X(:,i+1,j) = (A*X(:,i,j)) + (M.theta_stim.*mult(j).*S(:,i,j).*G(i)) + dT*eps(:,1,j);
    end
    
    %     % implement forced-RT collapsing decision boundary (exponential decay)
    %     if ~isnan(T.forcedPressTime(1,1)) && collapsedBound==0
    %         if t(i+1)>=0%800%T.forcedPressTime-T.respWindow*2
    %             %cti=i+1:(i+1+T.respWindow*3)-1; % collapsing time interval
    %             cti=i+1:(i+1+(2500/dT)-1); % collapsing time interval
    %             b=150; %0.01; % decay constant
    %             B(cti) = M.Bound .* -expm1( ( -(max(cti) - cti) ./ b ) ); % Boundary update: exponential decay
    %             collapsedBound=1; B(cti(end):end)=0;
    %         end
    %     end
    
    % find the press indecies that have to be planed in this decision cycle
    % doing it this way will be useful if/whenever multiple digits are being planned in every decsion step
    indx1 = nDecision * maxPlan - (maxPlan-1):min(nDecision * maxPlan , maxPresses);
    
    % determine if we issue a decision
    % motor command is not issued unless all the presses that have to planned in one step hit the boundary
    if ~isPressing && any(squeeze(X(:,i+1,indx1))>B(i+1)) && t(i+1)>=1750%(T.forcedPressTime-T.respWindow*3)
        count = 1;
        for prs = indx1
            [~,T.response(1,prs)]=max(X(:,i+1,prs));
            T.decisionTime(1,prs) = t(i+1);                            % Decision made at this time
            T.pressTime(prs) = T.decisionTime(prs)+count*M.dT_motor; % Press time delayed by motor delay
            % if there are any stumuli that have not appeared yet, set their stimTime to press time of Horizon presses before
            if sum(isnan(T.stimTime))
                idx2  = find(isnan(T.stimTime));
                for k = 1:length(idx2)
                    if ~isnan(T.pressTime(idx2(k) - T.Horizon(idx2(k))))
                        T.stimTime(idx2(k)) = T.pressTime(prs);
                    end
                end
            end
            isPressing = 1;                % Motor system engaged
            count = count+1;
        end
    end;
    
    % Update the motor system: Checking if current movement is done
    if (isPressing)
        if (t(i+1))>=T.pressTime(prs)
            isPressing = 0;
            % update the remaining presses
            remPress = max(0 , remPress - maxPlan);
            nDecision = nDecision+1;       % Waiting for the next decision
        end;
    end
    i=i+1;
end;

%% because the muber of presses to be planned is high, sometimes the trial times out and the decisionis not reached, so we need to account for that
if ~isfield(T , 'response') || (length(T.response) < maxPresses && i >= maxTime/dT)
    T.decisionTime(1 : maxPresses) = NaN;
    T.response(1 : maxPresses) = NaN;
    T.MT = NaN;
    T.timingError=1;
    T.isError=1;
    if (nargout>1)
        SIM.X = NaN;     % Hidden state
        SIM.S = NaN;     % Stimulus present
        SIM.B = NaN;     % Bound
        SIM.t = NaN;     % Time
        SIM.bufferSize = M.capacity;
    end;
else
    T.isError = zeros(size(T.TN));
    for i = 1:size(T.response,1)
        T.isError(i) = ~isequal(T.stimulus(i,1:T.numPress) , T.response(i,1:T.numPress));
    end
    % add timing errors for single resp exp
    if T.numPress==1 % single resp exp
        T.MT = NaN;
        if T.pressTime(1)<(T.forcedPressTime-T.respWindow) || T.pressTime(1)>(T.forcedPressTime+T.respWindow)
            T.timingError=1;
            T.isError=1;
        else % simple seq exp
            T.timingError=0;
        end
    else
        T.MT = max(T.pressTime,[],2);
    end
    tmax = T.pressTime(maxPresses);
    i = find(t == tmax);
    if (nargout>1)
        SIM.X = X(:,1:i-1,:); % Hidden state
        SIM.S = S(:,1:i-1,:); % Stimulus present
        SIM.B = B(1,1:i-1);     % Bound
        SIM.t = t(1,1:i-1);    % Time
        SIM.bufferSize = M.capacity;
    end;
end