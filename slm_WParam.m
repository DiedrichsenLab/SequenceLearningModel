%% N  =13
% M.Bound      = 0.45;
M.numOptions = 5;
M.dT_visual  = 90;
M.Ainhibit   = 0;
M.Capacity   =5;% 7;
M.DecayParam   = 7;
M.dT_motor   = 150;
M.dtGrowth = 1;
M.TSDecayParam  = 7.75;%5;
M.TSmin = 0.03775;
M.Aintegrate  = 0.94173;
if~noise
    M.SigEps      = 0;%0.01;
else
    M.SigEps      = 0;
end
bAll = 0.6;
M.Bound = bAll.*ones(M.Capacity,size(T.stimulus , 2));
M.Bound(:,1) = [.60042 ;.601295 ;.606129; .6037 ;.602804];
M.B0 = 4;%3.85478167568902;
M.B_coef = 0.351509146252228;
origCap = M.Capacity; % to preserve the original M.Capacity in designs that the capacity/horizon keeps changing
% for pn = 1:length(parName)
%     eval(['M.' , parName{pn} , ' = par(pn);'] )
% end
%% modulating theta_stim with capacity

% ========= creating an exponential decay
% M.theta_stim  = M.TSmin*exp(-([M.Capacity:-1:1]-1)./M.TSDecayParam);%linspace(0.01,0.04,M.Capacity);

% ========= modulating with actual MTs
load([mainDir ,'SequenceLearningModel/MTRTday5.mat'])
% M.theta_stim = M.TSmin + (K.MT_norm).* M.TSmin;
M.theta_stim = [0.034989 0.04042 0.047 0.054 0.0622];


%% N  =15;


M.numOptions = 5;
M.dT_visual  = 90;
M.Ainhibit   = 0;
M.Capacity   =5;% 7;
M.DecayParam   = 7;
M.dT_motor   = 150;
M.dtGrowth = 1;
M.TSDecayParam  = 7.75;
M.TSmin = 0.03775;
M.Aintegrate  = 0.941727795;
if~noise
    M.SigEps      = 0;
else
    M.SigEps      = 0.02;
end
bAll = 0.6;
M.Bound = bAll.*ones(M.Capacity,size(T.stimulus , 2));
M.Bound(:,1) = [.60041;.60143  ;.606171; .605914 ;.603727];
M.B0 = 4;
M.B_coef = 0.351509146252228;
origCap = M.Capacity; % to preserve the original M.Capacity in designs that the capacity/horizon keeps changing

%% modulating theta_stim with capacity

% ========= creating an exponential decay
% M.theta_stim  = M.TSmin*exp(-([M.Capacity:-1:1]-1)./M.TSDecayParam);%linspace(0.01,0.04,M.Capacity);

% ========= modulating with actual MTs
% load([mainDir ,'SequenceLearningModel/MTRTday5.mat'])
% M.theta_stim = M.TSmin + (K.MT_norm).* M.TSmin;
% M.theta_stim = [0.034989 0.04042 0.047 0.054 0.0622]; N = 13
% M.theta_stim = [0.0349899858731604,0.0404294904253830,0.0474978731759359,0.0544961706444883,0.0622975992041242]; % N = 15
M.theta_stim =   [0.0349877979,   0.0404294904253830,   0.047005,          0.0542 ,0.0622975992041242]; % N = 15

