function [T,SIM]=slm_simTrial(M,T,varargin);
% function [T,SIM]=slm_simTrial(M,T);
% Simulates a trial using a parallel evidence-accumulation model 
% for sequential response tasks 

dT = 2;     % delta-t in ms
maxTime = 10000; % Maximal time for trial simulation 
vararginoptions(varargin,{'dT','maxTime'}); 

% Determine length of the trial 
maxPresses = max(T.numPress); 

% initialize variables 
X = zeros(M.numOptions,maxTime/dT,maxPresses); % Hidden state 
S = zeros(M.numOptions,maxTime/dT,maxPresses); % Stimulus present 
B = ones(1,maxTime/dT)*M.Bound;                % Decision Bound - constant value 
t = [1:maxTime/dT]*dT-dT;   % Time in ms  
i = 1;                    % Index of simlation 
nDecision = 1;           % Current decision to male 
numPresses = 0;          % Number of presses produced 
isPressing = 0;          % Is the motor system currently occupied? 

% Set up parameters 
dec=[1:maxPresses]; % Number of decision 
A  = eye(M.numOptions)*(M.Aintegrate-M.Ainhibit)+...
     ones(M.numOptions)*M.Ainhibit; 

% Start time-by-time simulation 
while numPresses<maxPresses && i<maxTime/dT
   
    % Update the stimulus: Fixed stimulus time  
    indx = find(t(i)>(T.stimTime+M.dT_visual)); % Index of which stimuli are present T.
    for j=indx' 
        S(T.stimulus(j),i,j)=1;
    end; 
    
    % Update the evidence state 
    eps = randn([M.numOptions 1 maxPresses]) * M.SigEps; 
    mult=exp(-[dec-nDecision]./M.capacity);  % How much stimulus 
    mult(dec<nDecision)=0;                  % Made decisions will just decay 
    for j=1:maxPresses 
        X(:,i+1,j)= A * X(:,i,j) + M.theta_stim .* mult(j) .* S(:,i,j) + dT*eps(:,1,j); 
    end; 
    % Determine if we issue a decision 
    if nDecision<=maxPresses && ~isPressing && any(X(:,i+1,nDecision)>B(i+1)) 
        [~,T.response]=max(X(:,i+1,nDecision));
        T.decisionTime(nDecision) = t(i+1);                            % Decision made at this time    
        T.pressTime(nDecision) = T.decisionTime(nDecision)+M.dT_motor; % Press time delayed by motor delay  
        isPressing = 1;                % Motor system engaged
        nDecision = nDecision+1;       % Waiting for the next decision 
    end; 
    
    % Update the motor system: Checking if current movement is done 
    if (isPressing) 
        if (t(i+1))>=T.pressTime(nDecision-1)
            isPressing = 0; 
            numPresses = numPresses+1; 
        end; 
    end 
    i=i+1; 
end; 

if (nargout>1) 
SIM.X = X(:,1:i-1,:); % Hidden state 
SIM.S = S(:,1:i-1,:); % Stimulus present 
SIM.B = B(1,1:i-1);     % Bound 
SIM.t = t(1,1:i-1);    % Time 
end; 

