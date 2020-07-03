function varargout=slm_testModel(what,varargin)
% Wrapper function to test different aspects of the sml toolbox 
switch(what)
   
    case 'singleResp'
        % Make Model 
        M.Aintegrate = 1;    % Diagnonal of A  
        M.Ainhibit = 0;      % Inhibition of A 
        M.theta = 0.01;  % Rate constant for integration of sensory information 
        M.dT_motor = 90;     % Motor non-decision time 
        M.dT_visual = 70;    % Visual non-decision time 
        M.SigEps    = 0.01;   % Standard deviation of the gaussian noise 
        M.Bound     = 0.45;     % Boundary condition 
        M.numOptions = 5;    % Number of response options 
        M.capacity   = 3;   % Exponential decay on theta into future 
        % Make experiment 
        T.TN = 1; 
        T.numPress = 1; 
        T.window = 1; 
        T.forcedPressTime = [NaN NaN]; 
        T.stimulus = 1; 
        
        R=[]; 
        for i=1:1000
            [TR,SIM]=slm_simTrial(M,T); 
            % slm_plotTrial(SIM,TR); 
            R=addstruct(R,TR); 
        end; 
            % slm_plotTrial(SIM,TR); 
        subplot(1,2,1); 
        histplot(R.pressTime,'split',R.stimulus==R.response,'style_bar1'); 
        subplot(1,2,2); 
        
        keyboard; 
    case 'simpleSeq'
          % Make Model 
        M.Aintegrate = 0.98;    % Diagnonal of A  
        M.Ainhibit = 0;      % Inhibition of A 
        M.theta = 0.01;  % Rate constant for integration of sensory information 
        M.dT_motor = 90;     % Motor non-decision time 
        M.dT_visual = 70;    % Visual non-decision time 
        M.SigEps    = 0.02;   % Standard deviation of the gaussian noise 
        M.Bound     = 0.45;     % Boundary condition 
        M.numOptions = 5;    % Number of response options 
        M.capacity   = 1;   % Exponential decay on theta into future  
        
        % Make experiment 
        T.TN = 1; 
        T.numPress = 5; 
        T.window = 1; 
        T.forcedPressTime = nan(5,2); 
        T.stimulus = [1;2;5;4;3];  
        
        R=[]; 
        for i=1:1000
            [TR,SIM]=slm_simTrial(M,T); 
            slm_plotTrial(SIM,TR); 
            R=addstruct(R,TR); 
        end; 
            % slm_plotTrial(SIM,TR); 
        subplot(1,2,1); 
        histplot(R.pressTime,'split',R.stimulus==R.response,'style_bar1'); 
        subplot(1,2,2); 
        
        keyboard; 
    case 'horizon'
        
        % Make Model 
        M.Aintegrate = 0.995;    % Diagnonal of A  
        M.Ainhibit = 0;      % Inhibition of A 
        x=[0:4]; 
        M.theta = 0.01;  % Rate constant for integration of sensory information 
        M.dT_motor = 100;     % Motor non-decision time 
        M.dT_visual = 100;    % Visual non-decision time 
        M.SigEps    = 0;   % Standard deviation of the gaussian noise 
        M.Bound     = 1;     % Boundary condition 
        M.numOptions = 5;    % Number of response options 
        M.capacity   = 3;   % Exponential decay on theta into future 
        
        % Make experiment 
        T.TN = 1; 
        T.numPress = 13; 
        T.window = 3;  
        T.RT = 600; 
        T.stimulus = [1 2 3 4 5 1 2 3 4 5 1 2 3];  
        
        T=vararginoptionsStruct(varargin,{'window','RT'},T);

        [R,SIM]=slm_simTrial(M,T); 
        slm_plotTrial(SIM,R); 
        
        varargout= {R}; 
    case 'windowExp' 
        
        % Make Model 
        M.Aintegrate = 0.995;    % Diagnonal of A  
        M.Ainhibit = 0;      % Inhibition of A 
        x=[0:4]; 
        M.theta = 0.01;  % Rate constant for integration of sensory information 
        M.dT_motor = 100;     % Motor non-decision time 
        M.dT_visual = 100;    % Visual non-decision time 
        M.SigEps    = 0;   % Standard deviation of the gaussian noise 
        M.Bound     = 1;     % Boundary condition 
        M.numOptions = 5;    % Number of response options 
        M.capacity   = 3;   % Exponential decay on theta into future 
        
        % Make experiment 
        T.TN = 1; 
        T.numPress = 13; 
        T.RT = NaN; 
        T.stimulus = [1 2 3 4 5 1 2 3 4 5 1 2 3];  
        
        RR=[];
        for w=[1:4]
            T.window=w; 
            [R,SIM]=slm_simTrial(M,T); 
            S.ITI=diff(R.pressTime)';
            S.num=[1:T.numPress-1]';
            S.window = ones(T.numPress-1,1)*w; 
            RR=addstruct(RR,S); 
        end;
        lineplot(RR.num,RR.ITI,'split',RR.window,'style_thickline')
        varargout= {RR}; 
        
end