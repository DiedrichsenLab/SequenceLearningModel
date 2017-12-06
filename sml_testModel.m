function vararout=sml_testModel(what,varargin)
% Wrapper function to test different aspects of the sml toolbox 
switch(what)
   
    case 'singleResp'
        % Make Model 
        M.Aintegrate = 1;    % Diagnonal of A  
        M.Ainhibit = 0;      % Inhibition of A 
        M.theta_stim = 0.01;  % Rate constant for integration of sensory information 
        M.dT_motor = 90;     % Motor non-decision time 
        M.dT_visual = 70;    % Visual non-decision time 
        M.SigEps    = 0.1;   % Standard deviation of the gaussian noise 
        M.Bound     = 1;     % Boundary condition 
        M.numOptions = 5;    % Number of response options 
        
        % Make experiment 
        T.TN = 1; 
        T.numPress = 1; 
        T.stimTime = 0; 
        T.forcedPressTime = [NaN NaN]; 
        
        R=[]; 
        for i=1:1000
            TR=sml_simTrial(M,T); 
            R=addstruct(R,TR); 
        end; 
        
    case 'simpleSeq'
        
end