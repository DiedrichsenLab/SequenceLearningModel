load('se2_SimData.mat')
%% prep
poolHorizons = [5:13];
Dall.RT = Dall.AllPressTimes(:,1)-1500;
A = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0]) & ~Dall.isError & ismember(Dall.Day ,[1 4 5]));
if ~isempty(poolHorizons)
    A.Horizon(ismember(A.Horizon , poolHorizons)) = poolHorizons(1);
end
A.Day(ismember(A.Day , [4 5])) = 2;
F = getrow(AllR , ~AllR.isError);
Fit.IPI = F.IPI;
Fit.IPI = reshape(Fit.IPI , numel(Fit.IPI) , 1);
Act.IPI = A.IPI;
Act.IPI = reshape(Act.IPI , numel(Act.IPI) , 1);

Fit.singleH  = repmat(nanmean(F.Horizon, 2), 1 , size(A.AllPress,2)-1);
Fit.singleH  = reshape(Fit.singleH , numel(Fit.IPI) , 1);
Act.singleH  = repmat(nanmean(A.Horizon, 2) , 1 , size(A.AllPress,2)-1);
Act.singleH  = reshape(Act.singleH , numel(Act.IPI) , 1);

Fit.Day  = repmat(nanmean(F.Day, 2), 1 , size(A.AllPress,2)-1);
Fit.Day  = reshape(Fit.Day , numel(Fit.IPI) , 1);
Act.Day  = repmat(nanmean(A.Day, 2) , 1 , size(A.AllPress,2)-1);
Act.Day  = reshape(Act.Day , numel(Act.IPI) , 1);


Fit.ipiNum = repmat(1:size(F.stimulus,2)-1 , size(F.stimulus,1) , 1);
Fit.ipiNum = reshape(Fit.ipiNum , numel(Fit.ipiNum) , 1);
Act.ipiNum = repmat(1:size(F.stimulus,2)-1 , length(A.AllPress) , 1);
Act.ipiNum = reshape(Act.ipiNum , numel(Act.ipiNum) , 1);

Act.fitoract = ones(size(Act.ipiNum));
Fit.fitoract = zeros(size(Fit.ipiNum));

All  = addstruct(Fit , Act);

%%



c1 = [255, 153, 179]/255; % Random Red Tones
ce = [153, 0, 51]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 6);
end
for i = 1:length(tempcol)
    colz{i,1} = tempcol(: , i)';
end

clear tempcol
c1 = [153, 194, 255]/255; % structured blue tones
ce = [0, 0, 102]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 6);
end
for i = 1:length(tempcol)
    colz{i,2} = tempcol(: , i)';
    avgCol{i} = mean([colz{i,2} ; colz{i,1}],1);
end
plot(R_seq.planFunc(1,:), '-o' ,'LineWidth' , 2 , 'MarkerSize' , 5)
%% MTs
figure('color' , 'white')
subplot(211)
colorz = colz([2,end],1); % red for fitted
dayz = [F.Day;A.Day];
F.fitted = ones(size(F.MT));
A.fitted = zeros(size(A.MT));
A = normData(A , {'MT'});
lineplot([[F.fitted;A.fitted] [F.singleH;A.Horizon]]  , [F.MT;A.normMT] ,  'plotfcn','nanmean' ,'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 1 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 5, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'} , 'split' , dayz)
title('MT')
set(gca,'FontSize' , 7,'GridAlpha' , .2 , 'Box' , 'off','YLim' , [3000 7500] , 'YTick' , [3000:1000: 7000],...
    'YTickLabel' , [3:7] ,'XTickLabel' ,repmat({'1' , '2' , '3' , '4' , '5 - 13'} ,1,2));
ylabel('Execution time [s]')
xlabel('Viewing window size (W)')


subplot(212)
colorz = colz([2,end],2); % red for fitted

lineplot([[F.fitted;A.fitted] [F.singleH;A.Horizon]]  , [F.MT;A.normMT] ,  'plotfcn','nanmean' ,'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 1 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 5, 'markercolor' , colorz , 'leg' , {'Fitted' , 'Actual'} , 'split' , dayz)

title('MT - different color')
set(gca,'FontSize' , 7,'GridAlpha' , .2 , 'Box' , 'off','YLim' , [3000 7500] , 'YTick' , [3000:1000: 7000],...
    'YTickLabel' , [3:7] ,'XTickLabel' ,repmat({'1' , '2' , '3' , '4' , '5 - 13'} ,1,2));
ylabel('Execution time [s]')
xlabel('Viewing window size (W)')

%% planning functions
T.planFunc = F.planFunc	;
T.planFunc = reshape(T.planFunc , numel(T.planFunc) , 1);

T.N = repmat([1:14] , size(F.stimulus,1) , 1);
T.N = reshape(T.N , numel(T.N) , 1);

T.Day  = repmat(nanmean(F.Day, 2), 1 , size(A.AllPress,2));
T.Day  = reshape(T.Day , numel(Fit.planFunc) , 1);


figure('color' , 'white')
colorz = colz(3,1:2); % red for fitted
lineplot(T.N  , T.planFunc ,  'plotfcn','nanmean' ,'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,'none',...
    'linewidth' , 1 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 5, 'markercolor' , colorz , 'leg' , {'Early learning stage' , 'Late learning stage'} , 'split' , T.Day)
title('Planning function')
set(gca,'FontSize' , 7,'GridAlpha' , .2 , 'Box' , 'off','XLim' , [0.5 14] , 'XTick' , [1:14],...
    'XTickLabel' , [0:13],'YLim' , [0 1] , 'YTick' , [0 1]);
ylabel('Planning function value [au]')
xlabel('Distance from the upcoming cue [digits]')
%% IPIs
figure('color' , 'white')
subplot(211)
colorz = colz(:,2);
lineplot([ All.Day All.ipiNum] , All.IPI , 'plotfcn' , 'nanmean',...
    'split', All.singleH  , 'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 1.5 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 5, 'markercolor' , colorz , 'leg' , {'W = 1' , 'W = 2' , 'W = 3' , 'W = 4' , 'W = 5-13'} , ...
    'subset' , All.fitoract ==1);
ylabel('Inter-press interval [s]')
xlabel('IPIs number')
set(gca , 'FontSize' , 7 , 'Box' , 'off' , 'YLim' , [150 650],'YTick' , [200:100: 600],...
    'YTickLabel' , [0.2:.1:0.6] )

subplot(212)
colorz = colz(:,1); % red for fitted
lineplot([ All.Day All.ipiNum] , All.IPI , 'plotfcn' , 'nanmean',...
    'split', All.singleH  , 'linecolor' , colorz,...
    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
    'linewidth' , 1.5 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
    'markersize' , 5, 'markercolor' , colorz , 'leg' , {'W = 1' , 'W = 2' , 'W = 3' , 'W = 4' , 'W = 5-13'} , ...
    'subset' , All.fitoract ==0);
ylabel('Inter-press interval [s]')
xlabel('IPIs number')
set(gca , 'FontSize' , 7 , 'Box' , 'off' , 'YLim' , [150 650],'YTick' , [200:100: 600],...
    'YTickLabel' , [0.2:.1:0.6] )
