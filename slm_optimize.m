function [Param] = slm_optimize(param_init , varargin)
% initizlize and optimize
% Param(1) = M.capacity;
% Param(2) = M.theta_stim;
% Param(3) = M.Aintegrate;
% Param(4) = M.dtGrowth;
% Param(5) = M.SigEps;
% Param(6) = M.Ainhibit;

ANA = A;
model = @(param,T) slm_OptimSimTrial(param , T); % Model Function

SeqLength = unique(ANA.seqlength);
T.TN = ANA.TN;
T.Horizon =ANA.Horizon .*(ones(length(ANA.TN),SeqLength));
for tn = 1:length(ANA.TN)
    T.Horizon(tn , 1:ANA.Horizon(tn)) = NaN;
end
T.numPress = ANA.seqlength;
T.stimTime = zeros(length(ANA.TN) , SeqLength);
T.stimulus = ANA.AllPress;
T.forcedPressTime = nan(length(ANA.TN) , SeqLength);

% T = getrow(T  , [1:10]);
% ANA = getrow(ANA  , [1:10]);
% Set up the cost function
x_desired = [ANA.AllPressTimes(:,1) - 1500 , ANA.IPI];
OLS = @(param) nansum(nansum((model(param,T) - x_desired).^2));


param_init  = [4,0.0084,0.9785,0.6,0.01,3];

opts = optimset('MaxIter', 50,'TolFun',1e-5,'Display','iter');
[Param Fval] = fminsearch(OLS, param_init, opts);
% [Param Fval] = fminsearchbnd(OLS,param_init,[1 0.006 0.8 0.1 0.007 2],[5 0.02 0.98 10 0.01 5], opts);


[IPIs, Exam] = slm_diagModel( 'numSimulations' , 10,...
     'SigEps' , Param(5) ,'DecayParam' , ceil(Param(6)) ,'Aintegrate' , Param(3) , 'theta_stim' , Param(2) , 'Capacity' , ceil(Param(1)) ,...
     'SeqLength' , 14,'Horizons' , [1:2:14],'Ainhibit' , 0 , 'dtGrowth' , Param(4));