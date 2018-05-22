function [Param Fval] = slm_optimize(Dall , rawinitParam , varargin)
cycNum = 50;
samNum = 20;
tol = [0.5 0.03];
ItrNum = 20;
c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'parName'}
            % names of paramters
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'runNum'}
            % the number of run
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'cycNum'}
            % the number of times to sample the data and go through the optmization process
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'samNum'}
            % number of samples to pick from the data - leave [] to get all
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'tol'}
            % tolerance margines for std and mean respectively
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'ItrNum'}
            % number of iterations within a cycle
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'loBound'}
            % lower bound for parameters
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'hiBound'}
            % Higher bound for parameters
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Day'}
            % what days to include
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Horizon'}
            % what Horizon to include
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'poolHorizons'}
            % [] or the horizons you want to pool together
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'customizeInitParam'}
            % 1 or 0 ---> for boundries where the initial parametrs
            % are best to be set as the proportions of the IPIs
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'noise'}
            % 1 or 0 -
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end
%% optimization
Dall = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day , Day) & ismember(Dall.Horizon , Horizon));
if ~isempty(poolHorizons)
    Dall.Horizon(ismember(Dall.Horizon , poolHorizons)) = poolHorizons(1);
end
Horizon = unique(Dall.Horizon);
initParam = rawinitParam;
for i = 1:cycNum
    disp(['Initializing optimization cycle number ' , num2str(i) , '/', num2str(cycNum) , ' with ' , num2str(ItrNum) , ' iterations...'])
    if i>1
        customizeInitParam = 0; % dont customize after the first cycle
        load(['/Users/nedakordjazi/Documents/GitHub/SequenceLearningModel/param' , runNum , '.mat'])
        initParam = param.par(end , :);
    end
    % Set up the T structure
    ANA0 = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day , Day) & ismember(Dall.Horizon , Horizon));
    A = [];
    for h = Horizon
        ANA = getrow(ANA0 , ANA0.Horizon == h);
        serr = std(ANA.MT);%/sqrt(length(ANA.MT));
        stdbound = [serr-tol(1)*serr serr+tol(1)*serr];
        meanbound = [mean(ANA.MT)-tol(2)*mean(ANA.MT) mean(ANA.MT)+tol(2)*mean(ANA.MT)];
        MT_std = 0;
        MT_mean = 0;
        while ~(MT_std & MT_mean)
            if ~isempty(samNum)
                temp =  getrow(ANA , randperm(length(ANA.TN) , samNum));
            else
                temp  = ANA;
            end
            serr = std(temp.MT);%/sqrt(length(temp.MT));
            MT_std = serr>stdbound(1) & serr<stdbound(2);
            %         MT_std = 1;
            MT_mean = mean(temp.MT)>meanbound(1) & mean(temp.MT)<meanbound(2);
        end
        A = addstruct(A ,temp);
    end
    h1 = figure('color' , 'white');
    subplot(211);[~,p,~] = lineplot(A.Horizon , A.MT , 'plotfcn' , 'nanmean' , 'errorfcn' , 'nanstd');
    title('subset')
    subplot(212);[~ , p0 , ~] = lineplot(ANA0.Horizon , ANA0.MT, 'plotfcn' , 'nanmean', 'errorfcn' , 'nanstd');
    hold on
    [~,p,~] = lineplot(A.Horizon , A.MT , 'plotfcn' , 'nanmean', 'errorfcn' , 'nanstd');
    %
    
    ANA = A;
    model = @(param,T,runNum , i) slm_optimSimTrial(param , T , runNum , i , parName , 'optim' , noise); % Model Function
    
    SeqLength = unique(ANA.seqlength);
    T.TN = ANA.TN;
    T.Horizon =repmat(ANA.Horizon , 1, SeqLength) .*(ones(length(ANA.TN),SeqLength));
    for tn = 1:length(ANA.TN)
        T.Horizon(tn , 1:ANA.Horizon(tn)) = NaN;
    end
    T.numPress = ANA.seqlength;
    T.stimTime = zeros(length(ANA.TN) , SeqLength);
    T.stimulus = ANA.AllPress;
    T.forcedPressTime = nan(length(ANA.TN) , SeqLength);
    
    % Set up the cost function
    x_desired = [ANA.AllPressTimes(:,1) - 1500 ANA.IPI(: , 1:3) mean(ANA.IPI(: ,4:10) , 2) ANA.IPI(: , 11:13)];
    if ~noise
        T = getrow(T , 1); % when the noise is off all the trials will turn out identical
        x_desired = nanmedian(x_desired);
    end
    % curate the initial parameters if told to do so
    if customizeInitParam
        org = nanmedian(x_desired); 
        org = (org - min(org))/max(org);
        rawinitParam = mapminmax(org , .3,.6);
        initParam = rawinitParam;
    end
    
    OLS = @(param) nansum(nansum((model(param,T,runNum , i) - x_desired).^2));
    
    opts = optimset('MaxIter', ItrNum ,'TolFun',1000,'Display','iter');
    [Param Fval] = fminsearchbnd(OLS,initParam,loBound,hiBound, opts);
    
end

