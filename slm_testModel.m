function vararout=slm_testModel(what,varargin)
% Wrapper function to test different aspects of the sml toolbox
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
            % defines the parameters for the decay function
            % for 'exp' this would be the time constant (defaul = 1)
            % for 'linear' this would be a negative slope (default = -1/seqlength)
            % for 'boxcar' this would be the number of 1s in a row (default = 5)
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
            
        case {'SeqLength'}
            % defines the length of the sequence
            % default is 10
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'RandPortion'}
            % assign a proportion to the random sequences in the training block
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'numTrials'}
            % Number of trials in the training block
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Sequences'}
            % the sequences you want to train the model on
            % this is a matrix of size N*SeqLength, where N is the number of sequences
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

if ~exist('DecayFunc')
    DecayFunc = 'exp';
end
if ~exist('SeqLength')
    SeqLength = 10;
end
if ~exist('RandPortion')
    RandPortion = .5;
end
if ~exist('numTrials')
    numTrials = 1000;
end


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
        
        R=[];
        tn = 1;
        
        
        M.Aintegrate = 0.98;    % Diagnonal of A
        M.Ainhibit = 0;      % Inhibition of A
        M.theta_stim = 0.01;  % Rate constant for integration of sensory information
        M.dT_motor = 90;     % Motor non-decision time
        M.dT_visual = 70;    % Visual non-decision time
        M.SigEps    = 0.01;   % Standard deviation of the gaussian noise
        M.Bound     = 0.45;     % Boundary condition
        M.numOptions = 5;    % Number of response options
        M.capacity   = cap;   % Capacity for preplanning (buffer size)
        switch DecayFunc
            case 'exp'
                if ~exist('DecayParam')
                    DecayParam = cap;
                end
            case 'linear'
                if ~exist('DecayParam')
                    DecayParam = 1/(1-SeqLength);
                end
            case 'boxcar'
                if ~exist('DecayParam')
                    DecayParam = 5;
                end
        end
        % Make Models with defferent horizons
        for hrzn = 1:SeqLength - 1
            % Make experiment
            T.TN = tn;
            T.bufferSize = cap;  % useful is we decide to maipulate buffersize (inherited from M)
            T.numPress = SeqLength;
            T.stimTime = zeros(1 , SeqLength);
            T.forcedPressTime = nan(1 , SeqLength);
            % Horizon feature added. stimTime will be the actual time that the stimulus came on.
            T.Horizon = hrzn*(ones(1,SeqLength));
            T.Horizon(1:hrzn) = NaN;
            for i=1:20
                % generate random stimuli every rep
                T.stimulus = randi(5 , 1, SeqLength );
                [TR,SIM]=slm_simTrial(M,T,'DecayFunc' , DecayFunc,'DecayParam' , DecayParam);
                R=addstruct(R,TR);
            end;
            tn = tn +1;
        end
        slm_plotTrial('BlockMT' , SIM , R )
        slm_plotTrial('IPIHorizon' , SIM , R )
        keyboard;
        
    case 'SeqLearn'
        RndPor = ceil(numTrials*RandPortion);         % number of random sequences in the block
        SeqPor = numTrials - RndPor;                  % number of trained sequences in the block
        trPerSeq = ceil(SeqPor/size(Sequences, 1));   % number of trails per trained sequence
        
        
        R=[];
        cap= 1;
        M.Aintegrate = 0.98;    % Diagnonal of A
        M.Ainhibit = 0;      % Inhibition of A
        M.theta_stim = 0.01;  % Rate constant for integration of sensory information
        M.dT_motor = 90;     % Motor non-decision time
        M.dT_visual = 70;    % Visual non-decision time
        M.SigEps    = 0.01;   % Standard deviation of the gaussian noise
        M.Bound     = 0.45;     % Boundary condition
        M.numOptions = 5;    % Number of response options
        M.capacity   = cap;   % Capacity for preplanning (buffer size)
        switch DecayFunc
            case 'exp'
                if ~exist('DecayParam')
                    DecayParam = cap;
                end
            case 'linear'
                if ~exist('DecayParam')
                    DecayParam = 1/(1-SeqLength);
                end
            case 'boxcar'
                if ~exist('DecayParam')
                    DecayParam = 5;
                end
        end
        % Make Models with defferent horizons
        AllT = [];
        for hrzn = SeqLength
            for tn=1:RndPor
                % Random Sequences
                % Make experiment
                T.bufferSize = cap;  % useful is we decide to maipulate buffersize (inherited from M)
                T.numPress = SeqLength;
                T.seqType = 0;       % denoting it's random
                T.stimTime = zeros(1 , SeqLength);
                T.forcedPressTime = nan(1 , SeqLength);
                % Horizon feature added. stimTime will be the actual time that the stimulus came on.
                T.Horizon = hrzn*(ones(1,SeqLength));
                T.Horizon(1:hrzn) = NaN;
                T.stimulus = randi(5 , 1, SeqLength );
                AllT  = addstruct(AllT , T);
            end
            for ns = 1:size(Sequences, 1)
                for tn=1:trPerSeq
                    % Random Sequences
                    % Make experiment
                    T.bufferSize = cap;  % useful is we decide to maipulate buffersize (inherited from M)
                    T.numPress = SeqLength;
                    T.stimTime = zeros(1 , SeqLength);
                    T.forcedPressTime = nan(1 , SeqLength);
                    T.seqType = ns;       % denoting sequence number
                    % Horizon feature added. stimTime will be the actual time that the stimulus came on.
                    T.Horizon = hrzn*(ones(1,SeqLength));
                    T.Horizon(1:hrzn) = NaN;
                    T.stimulus = Sequences(ns , :);
                    AllT  = addstruct(AllT , T);
                end
            end
        end
        
        % Shuffle the trails to make them intermixed
        mixedTN = randperm(length(1:length(AllT.numPress)));
        AllT = getrow(AllT , mixedTN);
        AllT.TN = [1:length(AllT.numPress)]';
        
        [T,SIM]=slm_Learn(M,AllT);
        
        slm_plotTrial('BlockMT' , SIM , R )
        slm_plotTrial('IPIHorizon' , SIM , R )
        keyboard;
        
end