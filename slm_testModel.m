function [R,S,B,SIM,T,TR,M]=slm_testModel(what,varargin)
%% example call:
% [R,S,B,SIM,T,TR,M]=slm_testModel('singleResp','DecayFunc','exp','Aintegrate',0.975,'theta_stim',0.01,'DecayParam',2,'SigEps',0.01);

% Wrapper function to test different aspects of the sml toolbox
c = 1;
%% Set default parameters
subj = 1:20;
block = [1,2];
trial = 6:55;
plotSim = 0; % whether to plot each single trial simulation, or not
DecayFunc = 'exp';
DecayParam = 2;
SeqLength = 10;
%numSimulations = 200;
Aintegrate = 0.98;
Ainhibit = 0.0;
theta_stim = 0.01;
dT_motor = 90;
dT_visual = 70;
SigEps = 0.01;
Bound = 0.50; %0.90; %0.45;
numOptions = 5;
cap = 1;

%% manage varargin
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
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end

%% the cases

switch(what)
    
    case 'singleResp'
        % Make Model
        M.Aintegrate = Aintegrate;  % Diagnonal of A
        M.Ainhibit = Ainhibit;      % Inhibition of A
        M.theta_stim = theta_stim;        % Rate constant for integration of sensory information
        M.dT_motor = dT_motor;            % Motor non-decision time
        M.dT_visual = dT_visual;           % Visual non-decision time
        M.SigEps    = SigEps;         % Standard deviation of the gaussian noise
        M.Bound     = Bound;         % Boundary condition
        M.numOptions = numOptions;           % Number of response options
        M.capacity   = cap;         % Capacity for preplanning (buffer size)
        
        % Make experiment
        R=[]; S=[]; B=[];
        for s=subj
            for b=block
                [T]=slm_genTrial('sr2_rt',trial,s,b);
                for t=1:numel(trial)
                    [TR,SIM]=slm_simTrial(M,getrow(T,t),'DecayFunc',DecayFunc,'DecayParam',DecayParam);
                    if plotSim==1; slm_plotTrial('TrialHorseRace',SIM,TR); end
                    B=addstruct(B,TR);
                end
                B.BN=ones(numel(B.TN),1)*b;
                S=addstruct(S,B);
            end
            S.SN=ones(numel(S.TN),1)*s;
            R=addstruct(R,S);
        end
        R1=tapply(R,{'SN','prepTime'},{R.isError,'nanmean','name','ER','subset',R.timingError==0});
        plt.line(R1.prepTime,(1-R1.ER)*100,'errorbars','shade');
        xlabel('prep time'); ylabel('accuracy %');
        %keyboard;
        
    case 'simpleSeq'
        
        R=[];
        tn = 1;
        % Make Models with defferent horizons
        cap= 1;
        for hrzn = 1:2:SeqLength
            M.Aintegrate = Aintegrate;  % Diagnonal of A
            M.Ainhibit = Ainhibit;      % Inhibition of A
            M.theta_stim = theta_stim;        % Rate constant for integration of sensory information
            M.dT_motor = dT_motor;            % Motor non-decision time
            M.dT_visual = 70;           % Visual non-decision time
            M.SigEps    = SigEps;         % Standard deviation of the gaussian noise
            M.Bound     = Bound;         % Boundary condition
            M.numOptions = 5;           % Number of response options
            M.capacity   = cap;         % Capacity for preplanning (buffer size)
            switch DecayFunc
                case 'exp'
                    if ~exist('DecayParam','var')
                        DecayParam = cap;
                    end
                case 'linear'
                    if ~exist('DecayParam','var')
                        DecayParam = 1/(1-SeqLength);
                    end
                case 'boxcar'
                    if ~exist('DecayParam','var')
                        DecayParam = 5;
                    end
            end
            % Make experiment
            T.TN = tn;
            T.bufferSize = cap;  % useful is we decide to maipulate buffersize (inherited from M)
            T.numPress = SeqLength;
            T.stimTime = zeros(1 , SeqLength);
            T.forcedPressTime = nan(1 , SeqLength);
            % Horizon feature added. stimTime will be the actual time that the stimulus came on.
            T.Horizon = hrzn*(ones(1,SeqLength));
            T.Horizon(1:hrzn) = NaN;
            for i=1:numSimulations
                % generate random stimuli every rep
                T.stimulus = randi(5 , 1, SeqLength );
                [TR,SIM]=slm_simTrial(M,T,'DecayFunc' , DecayFunc,'DecayParam' , DecayParam);
                slm_plotTrial('TrialHorseRace' , SIM , TR )
                if i>1 && isequal(fields(TR) , fields(R))
                    R=addstruct(R,TR);
                elseif i == 1
                    R=addstruct(R,TR);
                end
            end
            tn = tn +1;
        end
%         slm_plotTrial('BlockMT' , SIM , R )
%         slm_plotTrial('IPIHorizon' , SIM , R )


        
end