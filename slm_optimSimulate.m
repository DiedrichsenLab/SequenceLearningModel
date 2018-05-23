function [R , IPI , MTRT] = slm_optimSimulate(Dall , par  , varargin)
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
        case {'noise'}
            % 1 or 0 -
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
if ~noise
    T = getrow(T , 1:10); % when the noise is off all the trials will turn out identical
end



R = slm_optimSimTrial(par , T , [] ,[] , parName ,'sim' , noise);
%% Horizon
All = [];
Act = [];
Fit = [];

Act.singleH = ANA.Horizon;
Fit.singleH  = nanmean(T.Horizon, 2);
Act.MT = ANA.MT;
Fit.MT = R(:,end);
Act.RT = ANA.AllPressTimes(:,1) - 1500;
Fit.RT = R(:,1);

Act.fitoract = ones(size(Act.MT));
Fit.fitoract = zeros(size(Fit.MT));

All  = addstruct(Fit , Act);

figure('color' , 'white')
subplot(211)
colorz = {[0.840000000000000,0.360000000000000,0.501176470588235],[0.360000000000000,0.456470588235294,0.760000000000000]};
lineplot(All.singleH , All.MT , 'plotfcn' , 'nanmean',...
    'split', All.fitoract  , 'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'});

title('MT')
grid on
set(gca , 'FontSize' , 16)
xlabel('Horizon')
subplot(212)
lineplot(All.singleH , All.RT , 'plotfcn' , 'nanmean',...
    'split', All.fitoract  , 'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'});
title('RT')
xlabel('Horizon')
grid on
set(gca , 'FontSize' , 16)

MTRT = All;
%% IPI
All = [];
Act = [];
Fit = [];
clear Fit Act
Fit.IPI = R(:,2:end-1);
Fit.IPI = reshape(Fit.IPI , numel(Fit.IPI) , 1);
Act.IPI = ANA.IPI;
Act.IPI = reshape(Act.IPI , numel(Act.IPI) , 1);

Fit.singleH  = repmat(nanmean(T.Horizon, 2) , 1 , size(T.stimulus,2)-1);
Fit.singleH  = reshape(Fit.singleH , numel(Fit.IPI) , 1);
Act.singleH  = repmat(nanmean(ANA.Horizon, 2) , 1 , size(ANA.AllPress,2)-1);
Act.singleH  = reshape(Act.singleH , numel(Act.IPI) , 1);


Fit.ipiNum = repmat(1:size(T.stimulus,2)-1 , size(T.stimulus,1) , 1);
Fit.ipiNum = reshape(Fit.ipiNum , numel(Fit.ipiNum) , 1);
Act.ipiNum = repmat(1:size(T.stimulus,2)-1 , length(ANA.AllPress) , 1);
Act.ipiNum = reshape(Act.ipiNum , numel(Act.ipiNum) , 1);

Act.fitoract = ones(size(Act.ipiNum));
Fit.fitoract = zeros(size(Fit.ipiNum));


All  = addstruct(Fit , Act);

figure('color' , 'white')

lineplot(All.ipiNum , All.IPI , 'plotfcn' , 'nanmean',...
    'split', All.fitoract  , 'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'});

title('IPIs')
xlabel('IPIs number')
grid on
set(gca , 'FontSize' , 16)
IPI = All;