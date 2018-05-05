function R = slm_optimSimulate(par , T)
% theta_stim = Param(1);
% Aintegrate = Param(2);
% Ainhibit   = Param(3);
% dtGrowth   = Param(4);
% SigEps     = Param(5);
% dt_motor   = Param(6);
% [IPIs, Exam] = slm_diagModel( 'numSimulations' , 10,...
%     'SigEps' , SigEps ,'DecayParam' , 4 ,'Aintegrate' , Aintegrate , 'theta_stim' , theta_stim , 'Capacity' , 3 ,...
%     'SeqLength' , 14,'Horizons' , [1:2:14],'Ainhibit' , Ainhibit , 'dtGrowth' , dtGrowth , 'dt_motor' , dt_motor);
R = slm_optimSimTrial(par , T , [] ,[]);
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