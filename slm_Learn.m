function [T_ammended,SIM_ammended]=slm_Learn(M,T,varargin);
% function [T,SIM]=slm_simTrialCap(M,T);
% incoporates horizon size (T.Horizon)
% T is a batch of trials
% multiple planning possible
% as long as M.capacity = 1,this funcion should work exacly same as [T,SIM]=slm_simTrial(M,T);
% Simulates learning using a parallel evidence-accumulation model
% for sequential response tasks
M.capacity = 1; % always set capacity to 1 since it's soft capacity
dT = 2;     % delta-t in ms
maxTime = 50000; % Maximal time for trial simulation
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
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end
if ~exist('DecayFunc')
    DecayFunc = 'exp';
end


AllT = T;

% Initialize the trials structure to be ammended with response times, decision times...
T_ammended = [];
SIM_ammended = [];
firstTrans = zeros(M.numOptions , M.numOptions); %initialize the matrix of first order conditional probabilities
BoundMat = M.Bound * ones(M.numOptions , M.numOptions);
for tn = 1:length(AllT.TN)
    T = getrow(AllT , tn);
    maxPresses = T.numPress;
    % number of decision steps is defined by the M.capacity and T.Horizon, whichever results in more decision steps
    dec=1: maxPresses;  % Number of decision
    maxPlan = M.capacity; % Number of digits being planned in every decision step
    if isfield(T , 'Horizon')
        % set the stimTime for presses that are not shown at time 0 to NaN.
        % These will be filled with press times
        T.stimTime(~isnan(T.Horizon)) = NaN;
    end
    T.pressTime = NaN*ones(size(T.Horizon)); % to be fiied with press times
    decSteps = max(dec);
    
    % initialize variables
    X = zeros(M.numOptions,maxTime/dT,maxPresses); % Hidden state
    S = zeros(M.numOptions,maxTime/dT,maxPresses); % Stimulus present
    B = ones(maxTime/dT , M.numOptions)*M.Bound;                % Decision Bound - constant value
    t = [1:maxTime/dT]*dT-dT;   % Time in ms
    i = 1;                   % Index of simlation
    nDecision = 1;           % Current decision to male
    numPresses = 0;          % Number of presses produced
    isPressing = 0;          % Is the motor system currently occupied?
    remPress = maxPresses;   % Remaining presses. This variable will be useful if/whenever multiple digits are being planned in every decsion step
    % Set up parameters
    
    A  = eye(M.numOptions)*(M.Aintegrate-M.Ainhibit)+...
        ones(M.numOptions)*M.Ainhibit;
    
    % Start time-by-time simulation
    
    while remPress && i<maxTime/dT
        % extract the probability of the current transition, and apply it to the decision boundry
        currentTransP = zeros(1,M.numOptions);
        if remPress<maxPresses
            for opt = 1: M.numOptions
                currentTransP(opt) = firstTrans(T.stimulus(maxPresses-remPress) , opt);
            end
        end
        mult = zeros(1,length(dec));
        % Update the stimulus: Fixed stimulus time
        indx = find(t(i)>(T.stimTime+M.dT_visual)); % Index of which stimuli are present T.
        % stimuli of greater than Horizon size will be unavailable
        for j=indx
            S(T.stimulus(j),i,j)=1;
        end;
        
        % Update the evidence state
        eps = randn([M.numOptions 1 maxPresses]) * M.SigEps;
        switch DecayFunc
            case 'exp'
                if ~exist('DecayParam')
                    DecayParam = maxPlan;
                end
                mult=exp(-[dec-nDecision]./DecayParam);      % How much stimulus exponentia decay
            case 'linear'
                if ~exist('DecayParam')
                    DecayParam = 1/(1-max(dec));
                end
                linDecay=DecayParam*(dec-max(dec));  % How much stimulus linear decay
                mult(nDecision:end) = linDecay(1:max(dec)-nDecision+1);  % How much stimulus linear decay
            case 'boxcar'
                if ~exist('DecayParam')
                    DecayParam = 5;
                end
                bc = min(max(dec) , nDecision+DecayParam);
                mult(nDecision:bc) = 1;  % How much stimulus linear decay
        end
        mult(dec<nDecision)=0;                    % Made decisions will just decay
        for j =1:maxPresses
            X(:,i+1,j)= A * X(:,i,j) + M.theta_stim .* mult(j) .* S(:,i,j) + dT*eps(:,1,j);
        end
        
        % find the press indecies that have to be planed in this decision cycle
        % doing it this way will be useful if/whenever multiple digits are being planned in every decsion step
        indx1= nDecision * maxPlan - (maxPlan-1):min(nDecision * maxPlan , maxPresses);
        % Determine if we issue a decision
        % motor command is not issued unless all the presses that have to planned in one step hit the boundry
        if ~isPressing && sum(sum(squeeze(X(:,i+1,indx1))>((1-2*currentTransP).*B(i+1 , :))')) == length(indx1)
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
        % update the rransition probabilities with every new trials
        
    end;
    T.isError = zeros(size(T.TN));
    for i = 1:length(T.isError)
        T.isError(i) = ~isequal(T.stimulus(i,:) , T.response(i,:));
    end
    
    
    
    firstTrans = zeros(M.numOptions , M.numOptions); %Re-initialize and update the matrix of first order probabilities
    for tn1 = 1:tn
        for prs = 1:AllT.numPress(tn1)-1
            firstTrans(AllT.stimulus(tn1 , prs) ,AllT.stimulus(tn1 , prs+1)) = firstTrans(AllT.stimulus(tn1 , prs) ,AllT.stimulus(tn1 , prs+1))+1;
        end
    end
    % deivde the values in the firstTrans matrix by the total number of double tansitions
    num2Trans = sum(AllT.numPress(1:tn)) - tn; % a sequence of length N has N-1 double transitions
    firstTrans = firstTrans/num2Trans;
    
    T_ammended = addstruct(T_ammended , T);
    T_ammended.MT = T_ammended.pressTime(:,end) - T_ammended.pressTime(:,1);
    subplot(131)
    plot(T_ammended.MT)
    title(['MT after ' , num2str(tn) , ' trials.'])
    drawnow()
    
    subplot(132)
    imagesc((1-5*firstTrans).*BoundMat)
    axis square
    title(['Boundry Matrix after ' , num2str(tn) , ' trials.'])
    colorbar
    drawnow()
    
    subplot(133)
    imagesc(firstTrans)
    axis square
    title(['Probability Matrix after ' , num2str(tn) , ' trials.'])
    colorbar
    drawnow()
    
    
%     SIM_ammended = addstruct(SIM_ammended , SIM);
end
out = 1;
% because the muber of presses to be planned is high, sometime the trial
% times out and the decisionis not reached, so we need to account for that



