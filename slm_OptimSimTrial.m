function R = slm_optimSimTrial(par , T , runNum ,cycNum , parName , mode)

%% parameters in the model (M) to be optimized
% M.capacity   = par(1);
% M.theta_stim = par(2);
% M.Aintegrate = par(3);
% M.Ainhibit   = par(4);
% M.dtGrowth   = par(5);
% M.SigEps     = par(6);
% DecayParam   = par(7);

%% outpu    [Raction time - 13 IPIS]
cd('/Users/nkordjazi/Documents/GitHub/SequenceLearningModel')
D = dir;
N = {};
for i = 1:length(D)
    N = [N {D(i).name}];
end
% save a new emty variable to ammend with optimization iterations
if~isempty(runNum) % is runNum is empty it means we are just simulating with a set of parametrs
    if ~sum(strcmp(N , ['param' , num2str(runNum) , '.mat']))
        param.cycNum = []; % number of resampling cycles within a run
        param.itrNum = []; % optimization iteration number within a cycle
        param.par = [];      % paramets
        param.parName = {};
        save(['/Users/nkordjazi/Documents/GitHub/SequenceLearningModel/param' , num2str(runNum) , '.mat'] ,'param' );
    end
    % ammend the ptimization matrix
    load(['/Users/nkordjazi/Documents/GitHub/SequenceLearningModel/param' , num2str(runNum) , '.mat'])
    
    P.cycNum = cycNum; % run number
    if ~isempty(param.itrNum)
        P.itrNum = param.itrNum(end)+1; % optimization cycle number within a run
    else
        P.itrNum = 1;
    end
    P.par = par;
    P.parName = parName;
    param = addstruct(param , P);
    save(['/Users/nkordjazi/Documents/GitHub/SequenceLearningModel/param' , num2str(runNum) , '.mat'] ,'param' );
end
%% we are going to hardcode tha parametrs in the model that we want to keep constant

% M.Bound      = 0.45;
M.numOptions = 5;
M.dT_visual  = 70;
M.Ainhibit   = 0;
M.capacity   = 1;
DecayParam   = 4;
M.dT_motor   = 120;
M.dtGrowth = 1;
M.theta_stim  = .0084;
M.Aintegrate  = 0.976;
M.SigEps      = 0.01;
M.Bound(1)    = .45;
M.Bound(2:3)  = .45;
M.Bound(4:12) = .45;
M.Bound(13:14)= .45;

for pn = 1:length(parName)
    eval(['M.' , parName{pn} , ' = par(pn);'] )
end

%%
% function [T,SIM]=slm_simTrial(M,T,varargin);
% incoporates horizon size (T.Horizon) as well as buffer size (M.capacity)
% multiple planning possible
% Simulates a trial using a parallel evidence-accumulation model for sequential response tasks


% the model is basically an ARX (auroregressive with exogenious input).
% X(i) = A*X(i-1) + Theta*Stimulus(i) + Noise      X defined the tensor of current state of decision making horseraces

%%
AllT = T;
AllR = [];
R = [];
for trls = 1:length(T.TN)
    T = getrow(AllT , trls);
    dT = 2;            %delta-t in ms
    maxTime = 50000;            % Maximal time for trial simulation
    M.capacity = min(M.capacity , nanmean(T.Horizon)); % this controls for situations where horizon size is smalled thatn capacity
    
    
    maxPresses = max(T.numPress);    % Determine length of the trial
    
    % number of decision steps is defined by the M.capacity and T.Horizon, whichever results in more decision steps
    % calculate the number of decision steps as total number of presses - capacity
    % this is because the first "capacity" presses would be planned in one decision step
    dec=1: max(maxPresses - M.capacity+1,M.capacity);  % Number of decision steps
%     dec=1: maxPresses;
    
    
    % for the first decision step, "capacity" digits are planned and for the rest, the shift is 1 by 1.
    maxPlan = ones(1 , length(dec)); % Number of digits being planned in every decision step
    if length(dec)>2 & length(dec)<maxPresses
        maxPlan(1) = M.capacity;
    end
    
    % set the stimTime for presses that are not shown at time 0 to NaN (depending on the horizon size)
    % These will be filled with press times (depending on the horizon size)
    T.stimTime(~isnan(T.Horizon)) = NaN;
    
    T.pressTime = NaN*ones(size(T.Horizon)); % to be filled with press times
    T.pressTime = NaN*ones(size(T.Horizon)); % to be filled with press times
    T.decisionTime = NaN*ones(size(T.Horizon));
    T.response = NaN*ones(size(T.Horizon));
    % initialize variables
    X = zeros(M.numOptions,maxTime/dT,maxPresses); % Hidden state
    S = zeros(M.numOptions,maxTime/dT,maxPresses); % Stimulus present
    % implement forced-RT collapsing decision boundary (logistic decay)

    B = ones(size(T.stimulus,2),maxTime/dT).*repmat(M.Bound'  ,1, maxTime/dT); % Decision Bound - constant value
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
    %% Linear growth for dt_motor to start faster and slow down to steady state
    % implimenting the idea of making dT a function of the percentage of the M.capacity that you have planned ahead
    dtgrowth = linspace(M.dT_motor ,M.dT_motor* M.dtGrowth, M.capacity);
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
        if ~isPressing && sum(sum(squeeze(X(:,i+1,PlanIndx(decPressCount)))>B(prs+1 , i+1))) >= 1
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
        if size(T.pressTime , 2)~=size(T.stimulus , 2) || size(T.decisionTime , 2)~=size(T.stimulus , 2)
            T.pressTime = nan(1 , size(T.stimulus , 2));
            T.decisionTime = nan(1 , size(T.stimulus , 2));
            T.response = nan(1 , size(T.stimulus , 2));
        end
    AllR = addstruct(AllR , T);
    
end
AllR.MT = AllR.pressTime(:,end) - AllR.pressTime(:,1);
AllR.RT =  AllR.pressTime(:,1);
AllR.singleH = nanmean(AllR.Horizon , 2);
AllR.IPI = diff(AllR.pressTime , [], 2);
switch mode
    case {'optim'}
        R =  [AllR.RT mean(AllR.IPI(:,1:2) , 2) , mean(AllR.IPI(:,4:9),2), mean(AllR.IPI(:,12:13),2)];
    case{'sim'}
        R =  [AllR.RT AllR.IPI AllR.MT];
end

