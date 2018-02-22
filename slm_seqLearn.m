function [T_seqLearn,SIM_seqLearn]=slm_seqLearn(M,AllT,varargin)
% function [T_seqLearn,SIM_seqLearn]=slm_seqLearn(M,AllT,varargin)
% incoporates horizon size (T.Horizon) as well as buffer size (M.capacity)
% multiple planning possible
% Simulates sequence learning using a parallel evidence-accumulation model for sequential response tasks
% AllT contains all the trials to be trained on

% the model is basically an ARX (auroregressive with exogenious input).
% X(i) = A*X(i-1) + Theta*Stimulus(i) + Noise      X defined the tensor of current state of decision making horseraces
%%
dT = 2;     % delta-t in ms
maxTime = 50000; % Maximal time for trial simulation
c = 1;

while(c<=length(varargin))
    switch(varargin{c})
        case {'DecayFunc'}
            % defines the decay type
            % can be 'exp', 'linear' or 'boxcar'
            % default is exponential
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'DecayParam'}
            % for 'exp' this would be the time constant (defaul = 1)
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'SeqLearningType'}
            % 'prob' for associative learning
            % 'chunk' for chunking mode
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'MapLearning'}
            % 0 or 1
            % Represents the general learnig of finger mappings that happens even in random sequences. could be implimented by lowering
            % the noise level, and can co-occur with either 'prob' , 'chunk' in parallel
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'AssLearnConst'}
            % Associative Learning Constant by which first order probabilities affect performance
            %(the higher the number the bigger the effect) Default  = 2
            % can set to 0 to onserve just the effect of MapLearning
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'MapLearnConst'}
            % Mapping Learning Constant by which the repetitions affect performance
            %(the higher the number the bigger the effect) Default  = 2
            % can set to 0 to onserve just the effect of just associations
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'VizProgress'}
            % Visualization Default 0
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end

%VizProgress = 0;
if exist('VizProgress','var')
    figure('color' , 'white');
end

% initilize the press counter to establish mapping learning
OptionsPressed = zeros(1, M.numOptions);

% ititialize the logistic growth function that impliments mapping learning.
% this will produce a coefficient based on number of presses made per digit
% this coefficient will increase Theta_stim accordingly
% the 3000 assumes that after having pressed each digits 3000 time,
% (approx.3 days, with blocks of 36 trials, 10 blocksper day, selength of 14),
% mapping will reach full speed
%%  logistic growth
% b = .009; % the growth constant. the bigger the b the faster the growth --> reached 1 faster
% Growth = 1./(1+3000*exp(-b*[1:length(AllT.TN)]));
%% exponential growth
% b1 = 15000; % the growth constant. the bigger the b the slower the growth --> reached 1 slower
% ThetaGrowth  = 1-exp(-((1:length(AllT.TN))-1)./b1);
b2 =2000;
ProbGrowth  = 1-exp(-((1:length(AllT.TN))-1)./b2);
%%
T_seqLearn = [];
% SIM_ammended = [];
A  = eye(M.numOptions)*(M.Aintegrate)+(~eye(M.numOptions)).*M.Ainhibit;      % A defined the matrix of autoregressive coefficients
SIM_seqLearn.firstTrans = zeros(M.numOptions , M.numOptions); %initialize the matrix of first order conditional probabilities
for tn = 1:length(AllT.TN)
    T = getrow(AllT , tn);
    M.capacity = min(M.capacity , max(T.Horizon));
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
    t = (1:maxTime/dT)*dT-dT;   % Time in ms
    
    % implement forced-RT collapsing decision boundary (logistic decay)
    if ~isnan(T.forcedPressTime(1,1))
        a = 0.005;%0.05; % the decay constant: the bigger the faster the decay --> reaches zero faster
        b = T.forcedPressTime + 100; % in ms, the inflexion point
        B = (M.Bound ./ (1 + exp ( a * (t - b) ) ) ) - M.Bound / 2; % Decision Bound - logistic decay
        B(B<=0)=0;
    else
        B = ones(1,maxTime/dT)*M.Bound; % Decision Bound - constant value
    end
    
    i = 1;                      % Index of simlation
    nDecision = 1;              % Current decision to male
    % numPresses = 0;             % Number of presses produced
    isPressing = 0;             % Is the motor system currently occupied?
    remPress = maxPresses;      % Remaining presses. This variable will be useful if/whenever multiple digits are being planned in every decsion step
    
    % Set up parameters
    prs = 0; % indexes the pressesd digits
    
    % find the press indecies that have to be planed in the first decision cycle
    indx1= prs+1 : prs+1+(maxPlan(nDecision)-1);
    
    cap_mult = ones(1,M.capacity-1);
    mult = [cap_mult , zeros(1,length(dec))];
    mult = [mult(logical(mult)) , exp(-(dec(1:end)-nDecision)./DecayParam)];      % How much stimulus exponentia decay
    decPressCount = 1;
    
    % initialize learning related parametrs
    %TempTrans = zeros(M.numOptions,M.numOptions); %initialize the matrix of first order conditional probabilities
    
    % Use logistic growth for stimulus horse race
    a = .1; %0.09; % the growth constant: the bigger the faster the growth --> reaches Bound faster
    b = 400; %200; %400; % in ms, how long it takes for the function to reach max
    G = (M.Bound/2) ./ (1 + exp ( -a * (t - (T.stimTime(1)+b/2) ) ) ); % logistic growth
    
    % update A matrix after each trial to account for assiciative learning
    while nDecision<=length(dec) && i<maxTime/dT
        
        % Update the stimulus: Fixed stimulus time
        indx = find(t(i)>(T.stimTime+M.dT_visual)); % Index of which stimuli are present T.
        % stimuli of greater than Horizon size will be unavailable
        if sum(indx)
            for j=indx
                S(T.stimulus(j),i,j)=1;
            end;
        end
        
        %% restore the A matrix for every press
        if T.seqNum > 1
            A  = eye(M.numOptions)*(M.Aintegrate)+(~eye(M.numOptions)).*M.Ainhibit;      % A defined the matrix of autoregressive coefficients
            % update the A matrix for each particular transition
            % the off diagonals are affected by the transition probability from
            % the previous press to the current press
            currentTransP = zeros(1,M.numOptions);
            if prs<maxPresses
                for opt = 1: M.numOptions
                    nextPress = T.stimulus(prs+1);
                    currentTransP(opt) = SIM_seqLearn.firstTrans(nextPress,opt);
                    A(nextPress , opt) = A(nextPress , opt) + (ProbGrowth(tn)*.1*currentTransP(opt));
                end
            end
        end
        %%
        
        % Update the evidence state
        eps = randn([M.numOptions 1 maxPresses]) * M.SigEps;
        
        for j =1:maxPresses
            if ~isnan(T.forcedPressTime(1,1))
                X(:,i+1,j) = (A*X(:,i,j)) + (M.theta_stim.*mult(j).*S(:,i,j).*G(i)) + dT*eps(:,1,j);
            else
                X(:,i+1,j)= A * X(:,i,j) + M.theta_stim .* mult(j) .* S(:,i,j) + dT*eps(:,1,j);
            end
        end
        
        % if the system is in the first decision cycle, it wont start pressing
        % until all the digits within the buffer size have been planned
        if ~isPressing && sum(sum(squeeze(X(:,i+1,indx1(decPressCount)))>B(i+1))) >= 1
            if decPressCount <= length(indx1)
                T.decisionTime(1,indx1(decPressCount)) = t(i+1);                            % Decision made at this time
                [~,T.response(1,indx1(decPressCount))]=max(X(:,i+1,indx1(decPressCount)));
                decPressCount = decPressCount + 1;
            end
            % after filling the buffer, start pressing
            if decPressCount == length(indx1)+1
                dtcount = 1;  % motor deltT counter
                for prs = indx1
                    T.pressTime(prs) = T.decisionTime(indx1(end))+dtcount*M.dT_motor; % Press time delayed by motor delay
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
            end;
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
                if nDecision<=length(dec)
                    % update the press indecies that have to be planed in next decision cycle
                    indx1= prs+1 : prs+1+(maxPlan(nDecision)-1);
                    decPressCount = 1;
                end
            end
        end
        i=i+1;
        
    end
    disp([num2str(tn) , ' Trials complete'])
    if ~isfield(T , 'pressTime') || ~isfield(T , 'response') || (length(T.pressTime) < maxPresses && i >= maxTime/dT)
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
    
    T.IPI = diff(T.pressTime);
    TempTrans = zeros(M.numOptions , M.numOptions); %Re-initialize and update the matrix of first order probabilities
    for tn1 = 1:tn
        for prs = 1:length(AllT.stimulus(tn1,:))-1
            TempTrans(AllT.stimulus(tn1 , prs) ,AllT.stimulus(tn1 , prs+1)) = TempTrans(AllT.stimulus(tn1 , prs) ,AllT.stimulus(tn1 , prs+1))+1;
        end
    end
    % find the conditional probabilities
    for cond  = 1:M.numOptions
        TempTrans(cond , :) = TempTrans(cond , :)/sum(TempTrans(cond , :));
    end
    TempTrans(isnan(TempTrans)) = 0;
    % set the diagonal to zero. otherwise Aintegrate would become too big
    SIM_seqLearn.firstTrans = TempTrans;% - diag(diag(TempTrans));
    
    % maplearning: find the number of times each finger has been pressed to modify the noise std
    for nop = 1:M.numOptions
        OptionsPressed(nop) = OptionsPressed(nop) + sum(T.stimulus == nop);
    end
    T_seqLearn = addstruct(T_seqLearn , T);
    T_seqLearn.MT = T_seqLearn.pressTime(:,end) - T_seqLearn.pressTime(:,1);
    S = getrow(T_seqLearn , T_seqLearn.seqNum > 0);
    R = getrow(T_seqLearn , T_seqLearn.seqNum ==0);
    if isempty(R.MT)
        R.MT = 0;
    elseif isempty(S.MT)
        S.MT = 0;
    end
    if exist('VizProgress','var')
        subplot(131)
        plot(S.MT , 'color' , 'b');
        hold on
        plot(R.MT , 'color' , 'r')
        legend({'Blue Trained Sequences' , 'Red Random Sequences'})
        title(['MT after ' , num2str(tn) , ' trials. Average Random_MT  = ' num2str(nanmean(R.MT)) , '    Average Trained_MT  = ' num2str(nanmean(S.MT))])
        drawnow()
        
        ipiMat = zeros(M.numOptions , M.numOptions);
        for tn1 = 1:length(T_seqLearn.TN)
            stim = T_seqLearn.stimulus(tn1,:);
            for pr = 1:size(stim , 2)-1
                ipiMat(stim(pr) , stim(pr+1)) = ipiMat(stim(pr) , stim(pr+1)) + T_seqLearn.IPI(tn1 , pr);
            end
        end
        ipiMat= ipiMat/tn1;
        subplot(132)
        imagesc(ipiMat)
        axis square
        title(['Average IPIs ' , num2str(tn) , ' trials.'])
        colorbar
        drawnow()
        
        subplot(133)
        imagesc(SIM_seqLearn.firstTrans)
        axis square
        title(['Probability Matrix after ' , num2str(tn) , ' trials.'])
        colorbar
        drawnow()
    end
    
end