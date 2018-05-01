function [T,SIM]=slm_simTrial(M,T,varargin)
% function [T,SIM]=slm_simTrial(M,T,varargin);
% incoporates horizon size (T.Horizon) as well as buffer size (M.capacity)
% multiple planning possible
% Simulates a trial using a parallel evidence-accumulation model for sequential response tasks 


% the model is basically an ARX (auroregressive with exogenious input).
% X(i) = A*X(i-1) + Theta*Stimulus(i) + Noise      X defined the tensor of current state of decision making horseraces
%%
dT = 2;            %delta-t in ms
maxTime = 50000;            % Maximal time for trial simulation
M.capacity = min(M.capacity , max(T.Horizon)); % this controls for situations where horizon size is smalled thatn capacity
c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'DecayParam'}  
            % for 'exp' this would be the time constant (defaul = 1)
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

maxPresses = max(T.numPress);    % Determine length of the trial

% number of decision steps is defined by the M.capacity and T.Horizon, whichever results in more decision steps
% calculate the number of decision steps as total number of presses - capacity
% this is because the first "capacity" presses would be planned in one decision step
dec=1: max(maxPresses - M.capacity+1,M.capacity);  % Number of decision steps


% for the first decision step, "capacity" digits are planned and for the rest, the shift is 1 by 1.
maxPlan = ones(1 , length(dec)); % Number of digits being planned in every decision step
if length(dec)>2
    maxPlan(1) = M.capacity;
end

% set the stimTime for presses that are not shown at time 0 to NaN (depending on the horizon size)
% These will be filled with press times (depending on the horizon size)
T.stimTime(~isnan(T.Horizon)) = NaN;

T.pressTime = NaN*ones(size(T.Horizon)); % to be filled with press times

% initialize variables 
X = zeros(M.numOptions,maxTime/dT,maxPresses); % Hidden state 
S = zeros(M.numOptions,maxTime/dT,maxPresses); % Stimulus present 
% implement forced-RT collapsing decision boundary (logistic decay)
if ~isnan(T.forcedPressTime(1,1))
    a = 0.005;%0.05; % the decay constant: the bigger the faster the decay --> reaches zero faster
    b = T.forcedPressTime + 100; % in ms, the inflexion point
    B = (M.Bound ./ (1 + exp ( a * (t - b) ) ) ) - M.Bound / 2; % Decision Bound - logistic decay
    B(B<=0)=0;
else
    B = ones(1,maxTime/dT)*M.Bound; % Decision Bound - constant value
end
t = [1:maxTime/dT]*dT-dT;   % Time in ms  
i = 1;                   % Index of simlation 
nDecision = 1;           % Current decision to male 
numPresses = 0;          % Number of presses produced 
isPressing = 0;          % Is the motor system currently occupied? 
remPress = maxPresses;   % Remaining presses. This variable will be useful if/whenever multiple digits are being planned in every decsion step
% Set up parameters 

A  = eye(M.numOptions)*(M.Aintegrate)+(~eye(M.numOptions)).*M.Ainhibit;      % A defined the matrix of autoregressive coefficients
% A(3,1) = 0.05;
prs = 0; % indexes the pressesd digits

% find the press indecies that have to be planed in the first decision cycle
PlanIndx= prs+1 : prs+1+(maxPlan(nDecision)-1);


% multiplier funstion for the stimulus evidence intake 
cap_mult = ones(1,M.capacity-1);
mult = [cap_mult , zeros(1,length(dec))];
mult = [mult(logical(mult)) , exp(-[dec(1:end)-nDecision]./DecayParam)];      % How much stimulus exponentia decay
decPressCount = 1;
%% Forced Prep time_____________ Use logistic growth for stimulus horse race
a = .1; %0.09; % the growth constant: the bigger the faster the growth --> reaches Bound faster
b = 400; %200; %400; % in ms, how long it takes for the function to reach max
G = (M.Bound/2) ./ (1 + exp ( -a * (t - (T.stimTime(1)+b/2) ) ) ); % logistic growth
%% Linear growth for dt_motor to start faster and slow down to steady state
% implimenting the idea of making dT a function of the percentage of the M.capacity that you have planned ahead
dtgrowth = linspace(M.dT_motor ,M.dT_motor*.9, M.capacity);
plannedAhead = zeros(1,maxPresses); % the number of digits planned ahead on each press
%%
while nDecision<=length(dec) && i<maxTime/dT    
    % Update the stimulus: Fixed stimulus time
    indx = find(t(i)>(T.stimTime+M.dT_visual)); % Index of which stimuli are present T.
 
    % stimuli of greater than Horizon size will be unavailable
    if sum(indx)
        for j=indx
            if T.stimulus(j)>0 && T.stimulus(j)<=5 % in single resp task, some stimuli are distractors
                S(T.stimulus(j),i,j)=1;
            end
        end
    end
    
    
    
    % Update the evidence state
    eps = randn([M.numOptions 1 maxPresses]) * M.SigEps;
    for j = 1:maxPresses
        if ~isnan(T.forcedPressTime(1,1))
            X(:,i+1,j) = (A*X(:,i,j)) + (M.theta_stim.*mult(j).*S(:,i,j).*G(i)) + dT*eps(:,1,j);
        else
            X(:,i+1,j)= A * X(:,i,j) + M.theta_stim .* mult(j) .* S(:,i,j) + dT*eps(:,1,j);
        end
    end
    % if the system is in the first decision cycle, it wont start pressing
    % until all the digits within the buffer size have been planned
    if ~isPressing && sum(sum(squeeze(X(:,i+1,PlanIndx(decPressCount)))>B(i+1))) >= 1
        plannedAhead(PlanIndx) = length(PlanIndx):-1:1;
        if decPressCount <= length(PlanIndx)
            T.decisionTime(1,PlanIndx(decPressCount)) = t(i+1);                            % Decision made at this time
            [~,T.response(1,PlanIndx(decPressCount))]=max(X(:,i+1,PlanIndx(decPressCount)));
            decPressCount = decPressCount + 1;
        end
        % after filling the buffer, start pressing
        if decPressCount == length(PlanIndx)+1
            dtcount = 0;
            for prs = PlanIndx
                T.pressTime(prs) = T.decisionTime(PlanIndx(end))+sum(dtgrowth(plannedAhead(PlanIndx(1):PlanIndx(1+dtcount)))); % Press time delayed by motor delay
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
                dtcount = dtcount+1;
            end
        end
    end
    % Update the motor system: Checking if current movement is done
    if (isPressing)
        if (t(i+1))>=T.pressTime(prs)
            isPressing = 0;
            mult(prs+1 : end) = mult(prs:end-1);
            mult(1:prs) = 0;
            % update the remaining presses
            remPress = max(0 , remPress - maxPlan(nDecision));
            nDecision = nDecision+1;       % Waiting for the next decision
            decPressCount = 1;
            if nDecision<=length(dec)
                % update the press indecies that have to be planed in next decision cycle
                PlanIndx= prs+1 : prs+1+(maxPlan(nDecision)-1);
                decPressCount = 1;
            end
        end
    end
    i=i+1;
end;
% because the muber of presses to be planned is high, sometime the trial
% times out and the decisionis not reached, so we need to account for that
if ~isfield(T , 'pressTime') || ~isfield(T , 'response') || (length(T.pressTime) < maxPresses && i >= maxTime/dT)
    T.decisionTime(length(T.pressTime)+1 : maxPresses) = NaN;
    T.response(length(T.pressTime)+1 : maxPresses) = NaN;
    T.pressTime(length(T.pressTime)+1 : maxPresses) = NaN;
    tmax = NaN;
    if (nargout>1)
        SIM.X = NaN;     % Hidden state
        SIM.S = NaN;     % Stimulus present
        SIM.B = NaN;     % Bound
        SIM.t = NaN;     % Time
        SIM.bufferSize = M.capacity;
    end;
else
    T.MT = max(T.pressTime,[],2);
    T.isError = zeros(size(T.TN));
    for i = 1:length(T.isError)
        T.isError(i) = ~isequal(T.stimulus(i,:) , T.response(i,:));
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