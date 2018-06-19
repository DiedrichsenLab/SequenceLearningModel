function R = slm_optimSimTrial(par , T , M ,opts)
%% Intro
% function R = slm_optimSimTrial(par , T , opts.runNum ,opts.cycNum , parName , mode , noise)
% incoporates horizon size (T.Horizon) as well as buffer size (M.Capacity)
% multiple planning possible
% Simulates a trial using a parallel evidence-accumulation model for sequential response tasks


% the model is basically an ARX (auroregressive with exogenious input).
% X(i) = A*X(i-1) + Theta*Stimulus(i) + Noise      X defined the tensor of current state of decision making horseraces

%% housekeeping

D = dir;
N = {};
for i = 1:length(D)
    N = [N {D(i).name}];
end
% mainDir = '/Users/nkordjazi/Documents/GitHub/SequenceLearningModel';
mainDir = '/Users/nedakordjazi/Documents/GitHub/SequenceLearningModel/';
% save a new emty variable to ammend with optimization iterations
if~isempty(opts.runNum) % is opts.runNum is empty it means we are just simulating with a set of parametrs
    cd(opts.saveDir )
    if ~sum(strcmp(N , ['param' , opts.runNum , '.mat']))
        param.cycNum = []; % number of resampling cycles within a run
        param.itrNum = []; % optimization iteration number within a cycle
        param.par = [];      % paramets
        param.parName = {};
        save([opts.saveDir ,'/param' , opts.runNum , '.mat'] ,'param' );
    end
    % ammend the ptimization matrix
    load([opts.saveDir ,'/param' , opts.runNum , '.mat'])
    
    P.cycNum = opts.cycNum; % run number
    if ~isempty(param.itrNum)
        P.itrNum = param.itrNum(end)+1; % optimization cycle number within a run
    else
        P.itrNum = 1;
    end
    P.par = par;
    P.parName = M.parName;
    param = addstruct(param , P);
    save([opts.saveDir ,'/param' , opts.runNum , '.mat'] ,'param' );
end

origCap = M.Capacity; % to preserve the original M.Capacity in designs that the capacity/horizon keeps changing
%% substitute the pre-set parameters with optimization parameters
for pn = 1:length(M.parName)
    eval(['M.' , M.parName{pn} , ' = par(pn);'] )
end
M.Bound = [M.bInit ones(1 ,size(T.stimulus,2)-1)*M.bAll]; % boundry is a vector of length maxPresses
%% 
AllT = T;
AllR = [];
R = [];
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
    
    B = ones(size(T.stimulus,2),maxTime/dT).*repmat(M.Bound'  ,1, maxTime/dT); % Decision Bound - constant value
    t = [1:maxTime/dT]*dT-dT;   % Time in ms
    i = 1;                   % Index of simlation
    nDecision = 1;           % Current decision to male
    numPresses = 0;          % Number of presses produced
    isPressing = 0;          % Is the motor system currently occupied?
    remPress = maxPresses;   % Remaining presses. This variable will be useful if/whenever multiple digits are being planned in every decsion step
    % Set up parameters
    
    A  = eye(M.numOptions)*(M.Aintegrate)+(~eye(M.numOptions)).*M.Ainhibit;      % A defined the matrix of autoregressive coefficients
    prs = 0; % indexes the pressesd digits
    
    % find the press indecies that have to be planned in the first decision cycle
    PlanIndx= prs+1 : prs+1+(maxPlan(nDecision)-1);
    
    
    %% planning curve funsction for the stimulus evidence intake
    switch M.PlanningCurve
        case 'exp'
            %% ============== exponential decay   1 --> ones for the
            % capacity and then exp decay
            cap_mult = ones(1,M.Capacity-1);
            planFunc = [cap_mult , zeros(1,length(dec))];
            planFunc = [planFunc(logical(planFunc)) , exp(-[dec-nDecision]./M.DecayParam)];      % How much stimulus exponentia decay
            %% ============== exponential decay   2 --> expdecay and then zeros from capacity onwards
%             planFunc = exp(-([1:maxPresses]-1)./M.DecayParam);
%             planFunc(M.Capacity+1:end) = 0;
%             T.planFunc = planFunc;
        case 'logistic'
            Xdomain = [-M.B_coef2:20];
            planFunc = 1./(1+1*exp(M.B_coef1*(Xdomain)));
            planFunc = planFunc(1:size(T.stimulus , 2));
        case 'box'
            planFunc = zeros(1,maxPresses);
            planFunc(1:floor(M.Box)) = 1;
        case 'ramp'
            planFunc = zeros(1,maxPresses);
            planFunc(1:floor(M.rampDecay)) = linspace(1,0,floor(M.rampDecay));
        case 'box_logistic'
            Xdomain = [-M.B_coef2:20];
            planFunc1 = 1./(1+1*exp(M.B_coef1*(Xdomain)));
            planFunc1 = planFunc1(1:size(T.stimulus , 2));
            
            planFunc2 = zeros(1,maxPresses);
            planFunc2(1:floor(M.rampDecay)) = linspace(1,0,floor(M.rampDecay));
            
            planFunc = planFunc1.*planFunc2;
        case 'box_exp'
            cap_mult = ones(1,M.Capacity-1);
            planFunc1 = [cap_mult , zeros(1,length(dec))];
            planFunc1 = [planFunc1(logical(planFunc1)) , exp(-[dec-nDecision]./M.DecayParam)];
            
            planFunc2 = zeros(1,maxPresses);
            planFunc2(1:floor(M.Box)) = 1;
            
            planFunc = planFunc1.*planFunc2;
        case 'box_ramp'
            planFunc1 = zeros(1,maxPresses);
            planFunc1(1:floor(M.Box)) = 1;
            
            planFunc2 = zeros(1,maxPresses);
            planFunc2(1:floor(M.rampDecay)) = linspace(1,0,floor(M.rampDecay));
            
            planFunc = planFunc1.*planFunc2;
    end
    T.planFunc = planFunc;
    pltmult = 0;
    if pltmult
        figure('color' , 'white')
        plot(planFunc , '-o' ,'LineWidth' , 1.5 , 'MarkerSize' , 5)
        set(gca , 'XLim' , [1 maxPresses] , 'Box' , 'off')
    end

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
                X(:,i+1,j) = (A*X(:,i,j)) + (M.theta_stim.*planFunc(j).*S(:,i,j).*G(i)) + eps(:,i+1,j);
            else
                X(:,i+1,j)= A * X(:,i,j) + M.theta_stim.* planFunc(j) .* S(:,i,j) + eps(:,i+1,j);
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
                planFunc(prs+1 : end) = planFunc(prs:end-1);
                planFunc(1:prs) = 0;
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
    T.X{1} = X(:,1:i,:);
    T.B{1} = M.Bound;
    if size(T.pressTime , 2)~=size(T.stimulus , 2) || size(T.decisionTime , 2)~=size(T.stimulus , 2)
        T.pressTime = nan(1 , size(T.stimulus , 2));
        T.decisionTime = nan(1 , size(T.stimulus , 2));
        T.response = nan(1 , size(T.stimulus , 2));
    end
    pltTrial = 0;
    if pltTrial
        SIM.t = t(1,1:i);
        SIM.B = repmat(T.B{1} , length(SIM.t) , 1);
        SIM.X = T.X{1};
        slm_plotTrial('TrialHorseRace' , SIM , T)
    end
    AllR = addstruct(AllR , T);
    
end
AllR.MT = AllR.pressTime(:,end) - AllR.pressTime(:,1);
AllR.RT =  AllR.pressTime(:,1);
AllR.singleH = nanmean(AllR.Horizon , 2);
AllR.IPI = diff(AllR.pressTime , [], 2);
for tn = 1:size(AllR.stimulus,1)
    AllR.isError(tn,1) = ~isequal(AllR.stimulus(tn, :) , AllR.response(tn  ,:));
end
R = [];
switch opts.mode
    case {'optim'}
        for xd = 1:length(opts.desiredField)
            H = unique(AllR.singleH);
            for hh = 1:length(H)
                Temp = getrow(AllR , ~AllR.isError & AllR.singleH == H(hh));
                if ~isempty(Temp.MT)
                    Temp = tapply(Temp ,{'singleH'}, {opts.desiredField{xd} , 'median'} );
                    eval(['R = [R ; Temp.' , opts.desiredField{xd} , '];'] )
                else
                    R = [R;NaN];
                end
            end
        end
        R =  R';
        R(isnan(R))=10e+10;
    case{'sim'}
        R =  AllR;
end

