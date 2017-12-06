function T=slm_simTrial(M,T,varargin);
% function T=slm_simTrial(M,T);
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

% Start time-by-time simulation 
while numPresses<maxPresses 
    
    % Update the stimulus: Fixed stimulus time  
    indx = find(t(i)>(T.stimTime+M.dT_visual)); % Index of which stimuli are present 
    for j=indx' 
        S(T.stimulus(j),i,j)=1;
    end; 
    
    % Update the evidence state 
    eps = randn([M.numOptions 1 maxPresses]) * M.SigEps; 
    X(:,i+1,:)= X(:,i,:) + M.theta_stim * S(:,i,:) + dT*eps; 
    
    % Determine if we issue a decision 
    if (any(X(:,i+1,nDecision)>B(i+1)) 
        [~,T.response]=max(X(:,i+1,nDecision));
        T.decisionTime(nDecision) = t(i+1);                            % Decision made at this time    
        T.pressTime(nDecision) = T.decisionTime(nDecision)+M.dt_motor; % Press time delayed by motor delay  
        
end; 




