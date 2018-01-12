function vararout=slm_testModel(what,varargin)
% Wrapper function to test different aspects of the sml toolbox
switch(what)
    
    case 'singleResp'
        % Make Model
        M.Aintegrate = 1;    % Diagnonal of A
        M.Ainhibit = 0;      % Inhibition of A
        M.theta_stim = 0.01;  % Rate constant for integration of sensory information
        M.dT_motor = 90;     % Motor non-decision time
        M.dT_visual = 70;    % Visual non-decision time
        M.SigEps    = 0.01;   % Standard deviation of the gaussian noise
        M.Bound     = 0.45;     % Boundary condition
        M.numOptions = 5;    % Number of response options
        M.capacity   = 3;   % Capacity for preplanning (buffer size)
        % Make experiment
        T.TN = 1;
        T.numPress = 1;
        T.stimTime = 0;
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
        %
        keyboard;
    case 'simpleSeq'
        seqLen = 10;
        R=[];
        tn = 1;
        % Make Models with defferent horizons
        cap= 1;
        for hrzn = 1:seqLen - 1
            M.Aintegrate = 0.98;    % Diagnonal of A
            M.Ainhibit = 0;      % Inhibition of A
            M.theta_stim = 0.01;  % Rate constant for integration of sensory information
            M.dT_motor = 90;     % Motor non-decision time
            M.dT_visual = 70;    % Visual non-decision time
            M.SigEps    = 0.01;   % Standard deviation of the gaussian noise
            M.Bound     = 0.45;     % Boundary condition
            M.numOptions = 5;    % Number of response options
            M.capacity   = cap;   % Capacity for preplanning (buffer size)
            
            % Make experiment
            T.TN = tn;
            T.bufferSize = cap;  % useful is we decide to maipulate buffersize (inherited from M)
            T.numPress = seqLen;
            T.stimTime = zeros(1 , seqLen);
            T.forcedPressTime = nan(1 , seqLen);
            % Horizon feature added. stimTime will be the actual time that the stimulus came on.
            T.Horizon = hrzn*(ones(1,seqLen));
            T.Horizon(1:hrzn) = NaN;
            for i=1:20
                % generate random stimuli every rep
                T.stimulus = randi(5 , 1, seqLen );
                [TR,SIM]=slm_simTrial(M,T);
                R=addstruct(R,TR);
            end;
            tn = tn +1;
        end
        slm_plotTrial('BlockMT' , SIM , R )
        slm_plotTrial('IPIHorizon' , SIM , R )
        keyboard;
        
end