function [varargout] = slm_findOptimalParams(varargin)
% function [varargout] = slm_findOptimalParams(varargin)
%
% init_params = [0.998, 0.0084, 0.0125]; % [Aintegrate, ThetaStim, SigmaNoise]
% init_params = [0.005]; % [Aintegrate, ThetaStim, SigmaNoise]
% example call: [optimal_params] = slm_findOptimalParams([0.005]);

% paths
pathToAnalyze = '/Users/gariani/Documents/data/SequenceRepetition/sr2/analyze';

% plot defaults
fs = 20; %default font size for all figures
lw = 4; %3; %default line width for all figures
ms = 12; %10; %default marker size for all figures

% colors
cbs_red = [213 94 0]/255;
cbs_blue = [0 114 178]/255;
cbs_yellow = [240 228 66]/255;
cbs_pink = [204 121 167]/255;
cbs_green = [0 158 115]/255;
blue = [49,130,189]/255;
lightblue = [158,202,225]/255;
red = [222,45,38]/255;
lightred = [252,146,114]/255;
green = [49,163,84]/255;
lightgreen = [161,217,155]/255;
orange = [253,141,60]/255;
yellow = [254,196,79]/255;
lightyellow = [255,237,160]/255;
purple = [117,107,177]/255;
lightpurple = [188,189,220]/255;
darkgray = [50,50,50]/255;
gray = [150,150,150]/255;
lightgray = [200,200,200]/255;
silver = [240,240,240]/255;
black = [0,0,0]/255;

% styles
style.reset;
style.custom({blue,lightblue,red,lightred,orange,yellow,lightyellow,purple,lightpurple,darkgray,gray,lightgray,green,lightgreen,black,silver,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink});
styred = style.custom({red}, 'markersize',ms, 'linewidth',lw);
styblue = style.custom({blue}, 'markersize',ms, 'linewidth',lw);


%% load actual data
D = load( fullfile(pathToAnalyze, 'sr2_training_all_data.mat') );

% create summary table for ACC
T = tapply(D, {'SN', 'prepTime'},...
    {(1-D.pressError) * 100, 'nanmean', 'name', 'ACC'}, ...
    'subset', D.dummy==0 & D.rtt==1 & D.timingError==0);

% normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
T = normData(T, {'ACC'}, 'sub');

figure; subplot(1,2,1);
[~, y] = plt.line(T.prepTime, T.normACC, 'errorbars','shade', 'style',styred);
title('Data'); xticklabels({'200','','300','','400','','500','','600',''}); xlim([180 670]); ylim([11 100]);
xlabel('Preparation time (ms)'); ylabel('Finger selection accuracy (%)');
axis square; set(gca,'fontsize', fs);


%% run simulation with initial parameters and fit simulated data to actual data (find best fit, search min loss, iteratively)
init_params = varargin{:};
[estim_params] = fitSim(y, init_params);


%% plot comparison between actual data and simulated data with optimal params
[R, ~, ~, ~, ~, ~, M] = slm_testModel('singleResp', ...
    'subj',1:20, 'block',[1,2], 'trial',1:55, 'plotSim',0, ...
    'Aintegrate',1, 'theta_stim',0.003548764584234, 'SigEps',estim_params(1));
R1 = tapply(R, {'SN', 'prepTime'}, {R.isError, 'nanmean', 'name', 'ER', 'subset',R.timingError==0});
subplot(1,2,1); hold on;
plt.line(R1.prepTime, (1-R1.ER)*100, 'errorbars','shade', 'style',styblue);
title( sprintf('Aint: %1.4f, Theta: %1.4f, Noise: %1.4f, Initial bound: %1.4f', M.Aintegrate, M.theta_stim, M.SigEps, M.Bound));
xlabel('Preparation time (ms)'); ylabel('Finger selection accuracy (%)');
axis square; set(gca,'fontsize', fs); hold off;


%% return the optimal parameters
varargout = {estim_params};

end

function [estim_params] = fitSim(y, init_params)
fcn          = @(params)bernoulli_loss(y, runSim(params));
estim_params = fminsearch(fcn, init_params);
end

function [y_hat] = runSim(params)
[R, ~, ~, ~, ~, ~, M] = slm_testModel('singleResp', ...
    'subj',1:20, 'block',[1,2], 'trial',1:55, 'plotSim',0, ...
    'Aintegrate',1, 'theta_stim',0.003548764584234, 'SigEps',params(1));
if all(R.isError==0) || all(R.isError==1)
    [~, y_hat] = plt.line(R.prepTime,(1-R.isError)*100,'errorbars','shade');
    title( sprintf('Aint: %1.4f, Theta: %1.4f, Noise: %1.4f, Initial bound: %1.4f', M.Aintegrate, M.theta_stim, M.SigEps, M.Bound));
    xlabel('Preparation time (ms)'); ylabel('Finger selection accuracy (%)'); xticklabels({'200','','300','','400','','500','','600',''}); xlim([180 670]); ylim([11 100]);
    axis square; set(gca,'fontsize', 20);
else
    R1 = tapply(R, {'SN', 'prepTime'}, {R.isError, 'nanmean', 'name', 'ER', 'subset',R.timingError==0});
    subplot(1,2,2); hold on;
    [~, y_hat] = plt.line(R1.prepTime, (1-R1.ER)*100, 'errorbars','shade');
    title( sprintf('Aint: %1.4f, Theta: %1.4f, Noise: %1.4f, Initial bound: %1.4f', M.Aintegrate, M.theta_stim, M.SigEps, M.Bound));
    xlabel('Preparation time (ms)'); ylabel('Finger selection accuracy (%)'); xticklabels({'200','','300','','400','','500','','600',''}); xlim([180 670]); ylim([11 100]);
    axis square; set(gca,'fontsize', 20);
end
end

function [loss] = bernoulli_loss(y, y_hat)
n = numel(y);
loss = zeros(n, 1);
for i = 1:n; loss(i, 1) = (y(i) * log(y_hat(i)) + (1-y(i)) * log(1-y_hat(i))); end
loss = -(nansum(loss));
end
