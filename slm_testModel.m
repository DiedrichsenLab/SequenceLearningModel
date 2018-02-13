function [R,SIM,Trials,Models]=slm_testModel(what,varargin)
% Wrapper function to test different aspects of the sml toolbox
c = 1;
%% Set default parameters
DecayFunc = 'exp';
SeqLength = 14;
numSimulations = 20;
Aintegrate = 0.98;
Ainhibit = 0.0;
theta_stim = 0.01;
dT_motor = 90;
SigEps = 0.01;
Bound = 0.45;
Horizons = [SeqLength];
Capacity = 1;
DecayParam = 1;
NumTrials = 1000;
RandPortion = .5;

%% manage varargin
while(c<=length(varargin))
    switch(varargin{c})
        
      
        case {'DecayParam'}
            % defines the parameters for the decay function
            % for 'exp' this would be the time constant (defaul = 1)
          
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'SeqLength'}
            % defines the length of the sequence
            % default is 10
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'numSimulations'}
            % defines the number of sequences to be simulated Default = 200
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Aintegrate'}
            % Diagonal of A: determines the efect of the previous value of each horserace on it's next value
            % (the one-back autoregressive coefficinnt) Default :  0.98
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Ainhibit'}
            % off-Diagonal of A: determines the effect of the previous
            % values of other horseraces on eachother's next values aka Lateral inhibition
            % Lateral inhibition Default 0
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'theta_stim'}
            % constant rate of information integration Default = 0.01
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'dt_motor'}
            % execution delta T Default = 90 ms
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'SigEps'}
            % Standard deviation of the gaussian noise Default = 0.01
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Bound'}
            % Decision Bound Default = 0.
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Horizons'}
            % The horizon sizes you want to include in the simulation
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Capacity'}
            % The plannign buffer size
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'NumTrials'}
            % the number of trials for the training block default  = 100
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'RandPortion'}
            % assign a proportion to the random sequences in the training block
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

%% the cases

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
        Trials = [];
        Models = [];
        tn = 1;
        % Make Models with defferent horizons
        for hrzn = Horizons
            M.Aintegrate = Aintegrate;  % Diagnonal of A
            M.Ainhibit = Ainhibit;      % Inhibition of A
            M.theta_stim = theta_stim;        % Rate constant for integration of sensory information
            M.dT_motor = dT_motor;            % Motor non-decision time
            M.dT_visual = 70;           % Visual non-decision time
            M.SigEps    = SigEps;         % Standard deviation of the gaussian noise
            M.Bound     = Bound;         % Boundary condition
            M.numOptions = 5;           % Number of response options
            M.capacity   = Capacity;         % Capacity for preplanning (buffer size)
            
            % Make experiment
            T.TN = tn;
            T.bufferSize = Capacity;  % useful is we decide to maipulate buffersize (inherited from M)
            T.numPress = SeqLength;
            T.stimTime = zeros(1 , SeqLength);
            T.forcedPressTime = nan(1 , SeqLength);
            % Horizon feature added. stimTime will be the actual time that the stimulus came on.
            T.Horizon = hrzn*(ones(1,SeqLength));
            T.Horizon(1:hrzn) = NaN;
            for i=1:numSimulations
                % generate random stimuli every rep
                T.stimulus = randi(5 , 1, SeqLength );
                [TR,SIM]=slm_simTrial(M,T,'DecayParam' , DecayParam);
                if sum(isnan(TR.pressTime))==0
                    R=addstruct(R,TR);
                    Trials = addstruct(Trials,T);
                    Models = addstruct(Models , M);
                end
                %                 slm_plotTrial('TrialHorseRace' , SIM , TR )
            end
            tn = tn +1;
        end
        %         slm_plotTrial('BlockMT' , SIM , R )
        %         slm_plotTrial('IPIHorizon' , SIM , R )
        
        
    case 'SeqLearn'
        if ~exist('Sequences')
            error('You have to provide training Sequences')
        end
       
        
        
        M.Aintegrate = Aintegrate;  % Diagnonal of A
        M.Ainhibit = Ainhibit;      % Inhibition of A
        M.theta_stim = theta_stim;        % Rate constant for integration of sensory information
        M.dT_motor = dT_motor;            % Motor non-decision time
        M.dT_visual = 70;           % Visual non-decision time
        M.SigEps    = SigEps;         % Standard deviation of the gaussian noise
        M.Bound     = Bound;         % Boundary condition
        M.numOptions = 5;           % Number of response options
        M.capacity   = Capacity;         % Capacity for preplanning (buffer size)
        AllT = [];
        
        
        TperH = ceil(NumTrials/length(Horizons));
        
        RndPor = ceil(TperH*RandPortion);         % number of random sequences in the block
        SeqPor = TperH - RndPor;                  % number of trained sequences in the block
        trPerSeq = ceil(SeqPor/size(Sequences, 1));   % number of trails per trained sequence
        
        for hrzn = Horizons
            for tn=1:RndPor
                % Random Sequences
                % Make experiment
                T.bufferSize = Capacity;  % useful is we decide to maipulate buffersize (inherited from M)
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
                    T.bufferSize = Capacity;  % useful is we decide to maipulate buffersize (inherited from M)
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
        [T,SIM]=slm_seqLearn(M,AllT,'DecayParam' , 2);
        
        
        
       
        %         slm_plotTrial('BlockMT' , SIM , R )
        %         slm_plotTrial('IPIHorizon' , SIM , R )
        
        
        
end