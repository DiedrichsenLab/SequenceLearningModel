function  [IPIs, Exam] =    slm_diagModel(varargin)
% model examination with different parameters
% testing out different dacaye time constatnts for the exponential
% clear all
%% Set default parameters
DecayFunc = 'exp';
SeqLength = 14;
numSimulations = 200;
Aintegrate =.97:0.005:0.995;
Ainhibit = 0.0;
theta_stim = .01:.05:.5;
SigEps = 0.01;
Bound = 0.45;
Horizons = [1:SeqLength];

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
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end



switch DecayFunc
    case 'exp'
        if ~exist('DecayParam')
            DecayParam = .2:.2:2;
        end
    case 'linear'
        if ~exist('DecayParam')
            DecayParam = 1/(1-SeqLength);
        end
    case 'boxcar'
        if ~exist('DecayParam')
            DecayParam = 1:SeqLength;
        end
end
%% calculate the simulated data structure
c = 1;
for dp = 1:length(DecayParam)
    for aint = 1:length(Aintegrate)
        for ainh = 1:length(Ainhibit)
            for ts = 1:length(theta_stim)
                dp
                [Exam(c).R,Exam(c).SIM,Exam(c).T,Exam(c).M] = slm_testModel('simpleSeq','DecayParam',DecayParam(dp),...
                    'numSimulations' , numSimulations, 'Aintegrate' , Aintegrate(aint),'Ainhibit' , ...
                    Ainhibit(ainh),'theta_stim' , theta_stim(ts) , 'SeqLength' , SeqLength ,...
                    'numSimulations' , numSimulations,'DecayFunc' , DecayFunc , 'SigEps' , SigEps,'Bound' , Bound , 'Horizons' , Horizons);
                Exam(c).R.singleH = nanmean(Exam(c).R.Horizon,2);
                Exam(c).R.singleH(isnan(Exam(c).R.singleH)) = size(Exam(c).R.pressTime , 2);
                c = c+1;
            end
        end
    end
end

%% still unpacking....
IPIs = [];
c = 1;
fig = figure;
for dp = 1:length(DecayParam)
    for aint = 1:length(Aintegrate)
        for ainh = 1:length(Ainhibit)
            for ts = 1:length(theta_stim)
                A = Exam(c).R;
                c = c+1;
                H = unique(A.singleH);
                H(isnan(H)) = size(A.pressTime , 2);
                H = unique(H);
                ipis = [];
                for h = 1:length(H)
                    B = getrow(A , A.singleH == H(h));
                    ipis.Horizon = H(h);
                    ipis.DecayParam = DecayParam(dp);
                    ipis.Aintegrate = Aintegrate(aint);
                    ipis.Ainhibit   = Ainhibit(ainh);
                    ipis.theta_stim = theta_stim(ts);
                    Bi = B;
                    ipis.ipi{1} = reshape(B.pressTime , numel(B.pressTime),1);
                    ipis.label{1} = ones(size(ipis.ipi{1}));
                    ipis.label{1}(1:2*numSimulations) = 0; % mark the first two IPIS as beggining IPIs
                    ipis.MT{1} = Bi.MT;
                    ipis.RT{1} = Bi.pressTime(:, 1);
                    
                    [~ , P , E] = lineplot(ipis.label{1} , ipis.ipi{1});
                     
                    ipis.P_start = P(1);
                    ipis.P_end = P(2);
                    ipis.E_start = E(1);
                    ipis.E_end = E(2);
                    IPIs = addstruct(IPIs , ipis);
                end
            end
        end
    end
end
close(fig);
