function R = slm_optimSimulate(Dall , par  , varargin)
samNum = 20;
tol = [0.5 0.03];
c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'parName'}
            % names of paramters
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
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end

Dall = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day , Day) & ismember(Dall.Horizon , Horizon));
if ~isempty(poolHorizons)
    Dall.Horizon(ismember(Dall.Horizon , poolHorizons)) = poolHorizons(1);
end
Horizon = unique(Dall.Horizon);
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
ANA = A;
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




R = slm_optimSimTrial(par , T , [] ,[] , parName ,'sim');
%% Horizon
clear A
A.singleH  = nanmean(T.Horizon, 2);
A.MT = R(:,end);
A.RT = R(:,1);
figure('color' , 'white')
subplot(211)
lineplot(A.singleH , A.MT , 'style_thickline')
title('MT')
grid on
set(gca , 'FontSize' , 16)
xlabel('Horizon')
subplot(212)
lineplot(A.singleH , A.RT,'style_thickline')
title('RT')
xlabel('Horizon')
grid on
set(gca , 'FontSize' , 16)


%% IPI
A.IPI = R(:,2:end-1);
A.singleH  = repmat(nanmean(T.Horizon, 2) , 1 , size(T.stimulus,2)-1);
A.ipiNum = repmat(1:size(T.stimulus,2)-1 , length(A.singleH) , 1);
figure('color' , 'white')
index = fliplr([reshape(A.ipiNum , numel(A.ipiNum) , 1) reshape(A.singleH, numel(A.singleH) , 1)]);
data  = reshape(A.IPI , numel(A.IPI) , 1);
lineplot(index , data , 'style_shade');
title('IPIs')
xlabel('IPIs number')
grid on
set(gca , 'FontSize' , 16)