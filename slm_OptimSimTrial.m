function R = slm_optimSimTrial(par , T , runNum ,cycNum , parName , mode , noise)

%% parameters in the model (M) to be optimized
% M.capacity   = par(1);
% M.theta_stim = par(2);
% M.Aintegrate = par(3);
% M.Ainhibit   = par(4);
% M.dtGrowth   = par(5);
% M.SigEps     = par(6);
% DecayParam   = par(7);

%% outpu    [Raction time - 13 IPIS]

D = dir;
N = {};
for i = 1:length(D)
    N = [N {D(i).name}];
end
% mainDir = '/Users/nkordjazi/Documents/GitHub/';
mainDir = '/Users/nedakordjazi/Documents/GitHub/';
% save a new emty variable to ammend with optimization iterations
if~isempty(runNum) % is runNum is empty it means we are just simulating with a set of parametrs
    cd([mainDir , 'SequenceLearningModel'])
    if ~sum(strcmp(N , ['param' , runNum , '.mat']))
        param.cycNum = []; % number of resampling cycles within a run
        param.itrNum = []; % optimization iteration number within a cycle
        param.par = [];      % paramets
        param.parName = {};
        save([mainDir , 'SequenceLearningModel/param' , runNum , '.mat'] ,'param' );
    end
    % ammend the ptimization matrix
    load([mainDir , 'SequenceLearningModel/param' , runNum , '.mat'])
    
    P.cycNum = cycNum; % run number
    if ~isempty(param.itrNum)
        P.itrNum = param.itrNum(end)+1; % optimization cycle number within a run
    else
        P.itrNum = 1;
    end
    P.par = par;
    P.parName = parName;
    param = addstruct(param , P);
    save([mainDir , 'SequenceLearningModel/param' , runNum , '.mat'] ,'param' );
end
%% we are going to hardcode tha parametrs in the model that we want to keep constant

% M.Bound      = 0.45;
M.numOptions = 5;
M.dT_visual  = 90;
M.Ainhibit   = 0;
M.Capacity   =5;% 7;
M.DecayParam   = 7;
M.dT_motor   = 150;
M.dtGrowth = 1;
M.TSDecayParam  = 7.75;%5;
M.TSmin = 0.03775;
M.Aintegrate  = 0.94173;
if~noise
    M.SigEps      = 0;%0.01;
else
    M.SigEps      = 0.02;
end
bAll = 0.6;
M.Bound = bAll.*ones(M.Capacity,size(T.stimulus , 2));
M.Bound(:,1) = [.60042 ;.601295 ;.606129; .6037 ;.602804];
M.B0 = 4;%3.85478167568902;
M.B_coef = 0.351509146252228;
origCap = M.Capacity; % to preserve the original M.Capacity in designs that the capacity/horizon keeps changing
% for pn = 1:length(parName)
%     eval(['M.' , parName{pn} , ' = par(pn);'] )
% end
%% modulating theta_stim with capacity

% ========= creating an exponential decay
% M.theta_stim  = M.TSmin*exp(-([M.Capacity:-1:1]-1)./M.TSDecayParam);%linspace(0.01,0.04,M.Capacity);

% ========= modulating with actual MTs
load([mainDir ,'SequenceLearningModel/MTRTday5.mat'])
% M.theta_stim = M.TSmin + (K.MT_norm).* M.TSmin;
M.theta_stim = [0.034989 0.04042 0.047 0.054 0.0622];
%%
% function [T,SIM]=slm_simTrial(M,T,varargin);
% incoporates horizon size (T.Horizon) as well as buffer size (M.Capacity)
% multiple planning possible
% Simulates a trial using a parallel evidence-accumulation model for sequential response tasks


% the model is basically an ARX (auroregressive with exogenious input).
% X(i) = A*X(i-1) + Theta*Stimulus(i) + Noise      X defined the tensor of current state of decision making horseraces

%%
AllT = T;
AllR = [];
R = [];
clear Growth
H = [1:6];
Xdomain = [-8:8];
for h = 1:length(H) % the decay constant. the bigger the b the faster the decay --> reached 0 faster
    B1 = M.B0 - (max(H)-h+1)*M.B_coef;
    Decay(h,:) = 1./(1+1*exp(B1*Xdomain));
end
for trls = 1:length(T.TN)
    M.Capacity   = origCap;
    T = getrow(AllT , trls);
    dT = 2;            %delta-t in ms
    maxTime = 50000;            % Maximal time for trial simulation
    M.Capacity = floor(min(M.Capacity , nanmean(T.Horizon))); % this controls for situations where horizon size is smalled thatn capacity
    
    %     after setting the dtgrowth set the capacity to 1
    
    maxPresses = max(T.numPress);    % Determine length of the trial
    
    % number of decision steps is defined by the M.Capacity and T.Horizon, whichever results in more decision steps
    % calculate the number of decision steps as total number of presses - capacity
    % this is because the first "capacity" presses would be planned in one decision step
    dec=1: max(maxPresses - M.Capacity+1,M.Capacity);  % Number of decision steps
%       dec=1: maxPresses;
    
    
    % for the first decision step, "capacity" digits are planned and for the rest, the shift is 1 by 1.
    maxPlan = ones(1 , length(dec)); % Number of digits being planned in every decision step
    if length(dec)>2 & length(dec)<maxPresses
        maxPlan(1) = M.Capacity;
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
    
    B = ones(size(T.stimulus,2),maxTime/dT).*repmat(M.Bound(M.Capacity , :)'  ,1, maxTime/dT); % Decision Bound - constant value
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
    
    
    %% multiplier funstion for the stimulus evidence intake
    % ============== exponential decay   1
    %     cap_mult = ones(1,M.Capacity-1);
    %     mult = [cap_mult , zeros(1,length(dec))];
    %     mult = [mult(logical(mult)) , exp(-[dec-nDecision]./M.DecayParam)];      % How much stimulus exponentia decay
    %
    % ============== exponential decay   2
    mult = exp(-([1:maxPresses]-1)./M.DecayParam);
    mult(M.Capacity+1:end) = 0;
    T.mult = mult;
    
    % ============== logistic decay
    %     mult = Decay(M.Capacity , 1:maxPresses);
    
    % ============== all ones --> when modulating theta_stim
    %     mult = ones(1,maxPresses);
    
    %%
    
    dtgrowth = linspace(M.dT_motor ,M.dT_motor* M.dtGrowth, M.Capacity);
    decPressCount = 1;
    %% Linear growth for dt_motor to start faster and slow down to steady state
    % implimenting the idea of making dT a function of the percentage of the M.Capacity that you have planned ahead
    
    plannedAhead = zeros(1,maxPresses); % the number of digits planned ahead on each press
    %%
    eps = randn([M.numOptions maxTime/dT maxPresses]) * M.SigEps;
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
        
        for j = 1:maxPresses
            if ~isnan(T.forcedPressTime(1,1))
                X(:,i+1,j) = (A*X(:,i,j)) + (M.theta_stim(M.Capacity).*mult(j).*S(:,i,j).*G(i)) + eps(:,i+1,j);
            else
                X(:,i+1,j)= A * X(:,i,j) + M.theta_stim(M.Capacity).* mult(j) .* S(:,i,j) + eps(:,i+1,j);
            end
        end
        % if the system is in the first decision cycle, it wont start pressing
        % until all the digits within the buffer size have been planned
        if ~isPressing && sum(squeeze(X(:,i+1,PlanIndx(decPressCount)))>B(prs+1 , i+1)) >= 1
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
        if i == 250000
            keyboard
        end
    end;
    T.X{1} = X(:,1:i,:);
    
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
        if length(unique(AllR.singleH))==1
%             R =  [AllR.RT AllR.IPI(: , 1:13)];% mean(AllR.IPI(: ,4:10) , 2) AllR.IPI(: , 11:13)];
            R =  [AllR.RT]';
            R(isnan(R))=10e+10;
        else
            R =  [AllR.RT]';% mean(AllR.IPI(: ,4:10) , 2) AllR.IPI(: , 11:13)];
            R(isnan(R))=10e+10;
        end
    case{'sim'}
        R =  AllR;
end

