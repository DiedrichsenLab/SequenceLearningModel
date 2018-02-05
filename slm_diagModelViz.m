function slm_diagModelViz(IPIs , what , varargin)

c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        
        case {'DecayParam_select'}
            % defines the parameters for the decay function
            % for 'exp' this would be the time constant (defaul = 1)
            % for 'linear' this would be a negative slope (default = -1/seqlength)
            % for 'boxcar' this would be the number of 1s in a row (default = 5)
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Aintegrate_select'}
            % Diagonal of A: determines the efect of the previous value of each horserace on it's next value
            % (the one-back autoregressive coefficinnt) Default :  0.98
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Ainhibit_select'}
            % off-Diagonal of A: determines the effect of the previous
            % values of other horseraces on eachother's next values aka Lateral inhibition
            % Lateral inhibition Default 0
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'theta_select'}
            % constant rate of information integration Default = 0.01
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

if ~exist('Aintegrate_select')
    Aintegrate_select = 0.98;
end
if ~exist('Ainhibit_select')
    Ainhibit_select = 0;
end
if ~exist('theta_select')
    theta_select = 0.01;
end
colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3],[.4 .4 .4] [.2 .5 .7]};
%% plot MTs and RTs for different decay constants and different Theta stims IPIs.Aintegrate == 0.98 for full horizon
theta_stim = unique(IPIs.theta_stim);
Aintegrate = unique(IPIs.Aintegrate);
DecayParam = unique(IPIs.DecayParam);
Horizon    = unique(IPIs.Horizon);

switch what
    case 'MT_RT_vs_theta'
        
        h = figure;
        for ts = 1:length(theta_stim)
            A = getrow(IPIs , IPIs.theta_stim == theta_stim(ts) & IPIs.Aintegrate == Aintegrate_select);
            temp = repmat(A.DecayParam , 1,length(A.MT{1}))';
            DecayLabel = reshape(temp , numel(temp) , 1);
            [xm{ts} , pm{ts} , em{ts}] = lineplot(DecayLabel , cell2mat(A.MT) , 'plotfcn','nanmean' ,'style_thickline');
            hold on
            [x{ts} , p{ts} , e{ts}] = lineplot(DecayLabel , cell2mat(A.RT) , 'plotfcn','nanmean' ,'style_thickline');
        end
        
        % *************************************************** MTs
        close(h)
        figure('color' , 'white')
        subplot(221)
        legenslabel = {};
        for ts = 1:length(theta_stim)
            legenslabel = [legenslabel , ['theta stim = ' , num2str(theta_stim(ts))]];
            errorbar(xm{ts} , pm{ts} , em{ts} , 'LineWidth' , 3);
            hold on
        end
        legend(legenslabel)
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title(['Movement Time for Different Exponential Decay Constants - Aintegrate = ' , num2str(Aintegrate_select)])
        % for better visibility of the other theta stims
        subplot(223)
        for ts = 2:length(theta_stim)
            errorbar(xm{ts} , pm{ts} , em{ts} , 'LineWidth' , 3, 'color' , colorz{dp});
            hold on
        end
        legend(legenslabel(2:end))
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title(['Movement Time for Different Exponential Decay Constants - Aintegrate = ' ,num2str(Aintegrate_select)])
        % histplot([IPIs.start{1} , IPIs.start{2}])
        
        
        % *************************************************** RTs
        
        subplot(222)
        legenslabel = {};
        for ts = 1:length(theta_stim)
            legenslabel = [legenslabel , ['theta stim = ' , num2str(theta_stim(ts))]];
            errorbar(x{ts} , p{ts} , e{ts} , 'LineWidth' , 3, 'color' , colorz{ts});
            hold on
        end
        legend(legenslabel)
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title(['Initial Reaction Time for Different Exponential Decay Constants - Aintegrate = ', num2str(Aintegrate_select)])
        % for better visibility of the other theta stims
        subplot(224)
        for ts = 2:length(theta_stim)
            errorbar(x{ts} , p{ts} , e{ts} , 'LineWidth' , 3, 'color' , colorz{ts});
            hold on
        end
        legend(legenslabel(2:end))
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title(['Initial Reaction Time for Different Exponential Decay Constants - Aintegrate = ', num2str(Aintegrate_select)])
    %%'MT_RT_vs_Aintegrate'
    case 'MT_RT_vs_Aintegrate'
        
        h = figure;
        hold on
        for ai = 1:length(Aintegrate)
            A = getrow(IPIs , IPIs.Aintegrate == Aintegrate(ai) & IPIs.theta_stim == theta_select);
            temp = repmat(A.DecayParam , 1,length(A.MT{1}))';
            DecayLabel = reshape(temp , numel(temp) , 1);
            [x{ai} , p{ai} , e{ai}] = lineplot(DecayLabel , cell2mat(A.MT) , 'plotfcn','nanmean' ,'style_thickline');
            [xr{ai} , pr{ai} , er{ai}] = lineplot(DecayLabel , cell2mat(A.RT) , 'plotfcn','nanmean' ,'style_thickline');
        end
        close(h)
        legenslabel = {};
        figure('color' , 'white')
        subplot(211)
        for ai = 1:length(Aintegrate)
            legenslabel = [legenslabel , ['A integrate = ' , num2str(Aintegrate(ai))]];
            errorbar(x{ai} , p{ai} , e{ai} , 'LineWidth' , 3, 'color' , colorz{ai});
            hold on
        end
        legend(legenslabel)
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title('Movement Time for Different Exponential Decay Constants - theta_stim = 0.01')
        
        subplot(212)
        for ai = 1:length(Aintegrate)
            errorbar(xr{ai} , pr{ai} , er{ai} , 'LineWidth' , 3);
            hold on
        end
        legend(legenslabel)
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title(['Initial Reaction Time for Different Exponential Decay Constants - theta_stim = ' , num2str(theta_select)])
        
    case 'MT_RT_vs_Horizon_constantAiTs'
        h1 = figure;
        for h = 1:length(Horizon)
            A = getrow(IPIs , IPIs.Horizon == Horizon(h) & IPIs.Aintegrate == Aintegrate_select  & IPIs.theta_stim == theta_select);
            temp = repmat(A.DecayParam , 1,length(A.MT{1}))';
            DecayLabel = reshape(temp , numel(temp) , 1);
            [xm{h} , pm{h} , em{h}] = lineplot(DecayLabel , cell2mat(A.MT) , 'plotfcn','nanmean' ,'style_thickline');
            hold on
            [x{h} , p{h} , e{h}] = lineplot(DecayLabel , cell2mat(A.RT) , 'plotfcn','nanmean' ,'style_thickline');
        end
        
        % *************************************************** MTs
        close(h1)
        figure('color' , 'white')
        subplot(211)
        legenslabel = {};
        for h = 1:length(Horizon)
            legenslabel = [legenslabel , ['Horizon = ' , num2str(Horizon(h))]];
            errorbar(xm{h} , pm{h} , em{h} , 'LineWidth' , 3, 'color' , colorz{h});
            hold on
        end
        legend(legenslabel)
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title(['Movement Time for Different Exponential Decay Constants - Aintegrate = ' num2str(Aintegrate_select),...
            'theta_stim = ' , num2str(theta_select)])
        
        
        % *************************************************** RTs
        
        subplot(212)
        legenslabel = {};
        for h = 1:length(Horizon)
            legenslabel = [legenslabel , ['Horizon = ' , num2str(Horizon(h))]];
            errorbar(x{h} , p{h} , e{h} , 'LineWidth' , 3, 'color' , colorz{h});
            hold on
        end
        legend(legenslabel)
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title(['Initial Reaction Time for Different Exponential Decay Constants - Aintegrate = ' num2str(Aintegrate_select),...
            'theta_stim = ' , num2str(theta_select)])
    case 'MT_RT_vs_Horizon'
        h1 = figure;
        for dp = 1:length(DecayParam)
            A = getrow(IPIs , IPIs.DecayParam == DecayParam(dp));
            temp = repmat(A.Horizon , 1,length(A.MT{1}))';
            HLabel = reshape(temp , numel(temp) , 1);
            [xm{dp} , pm{dp} , em{dp}] = lineplot(HLabel , cell2mat(A.MT) , 'plotfcn','nanmean' ,'style_thickline');
            hold on
            [x{dp} , p{dp} , e{dp}] = lineplot(HLabel , cell2mat(A.RT) , 'plotfcn','nanmean' ,'style_thickline');
        end
        
        % *************************************************** MTs
        close(h1)
        figure('color' , 'white')
        subplot(211)
        legenslabel = {};
        for dp = 1:length(DecayParam)
            legenslabel = [legenslabel , ['Decay = ' , num2str(DecayParam(dp))]];
            errorbar(xm{dp} , pm{dp} , em{dp} , 'LineWidth' , 3, 'color' , colorz{dp});
            hold on
        end
        legend(legenslabel)
        grid on
        ylabel('ms')
        xlabel('Horizon')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title(['Movement Time for Different Exponential Decay Constants'])
        
        
        % *************************************************** RTs
        
        subplot(212)
        legenslabel = {};
        for dp = 1:length(Horizon)
            legenslabel = [legenslabel , ['Decay = ' , num2str(DecayParam(dp))]];
            errorbar(x{dp} , p{dp} , e{dp} , 'LineWidth' , 3, 'color' , colorz{dp});
            hold on
        end
        legend(legenslabel)
        grid on
        ylabel('ms')
        xlabel('Horizon')
        set(gca, 'FontSize' , 20 ,'Box' , 'off')% , 'GridAlpha' , 1)
        title(['Initial Reaction Time for Different Exponential Decay Constants'])
        
        
        %% plot MTs and RTs for different decay constants and different Aintegrates IPIs.theta_stim = theta_select
        
    case 'IPIs_vs_Aintegrate'
        %% plot satrt and end IPIs with different decay constants and different Aintegrates IPIs.theta_stim = theta_select
        % ********************** varying decay
        colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3],[.4 .4 .4] [.2 .5 .7]};
        legenslabel = {};
        xlab = {'IPI 1,2' , 'IPI 3:9' };
        for ai = 1:length(Aintegrate)
            legenslabel = [legenslabel , ['A integrate = ' , num2str(Aintegrate(ai))]];
            A = getrow(IPIs , IPIs.Aintegrate == Aintegrate(ai) & IPIs.theta_stim == theta_select);
            ipiLabels = [0 1];
            x = 0;
            for dp = 1:length(DecayParam)
                ipis_e{ai,dp} = [A.E_start(dp) A.E_end(dp)];
                ipis_p{ai,dp} = [A.P_start(dp) A.P_end(dp)];
                X{ai,dp} = [x x+1];
                hold on
                x = x+2;
            end
        end
        
        figure('color' , 'white')
        for dp = 1:length(DecayParam)
            for ai = 1:length(Aintegrate)
                errorbar( X{ai,dp} , ipis_p{ai,dp} , ipis_e{ai,dp} , 'LineWidth' , 3 , 'color' , colorz{ai});
                hold on
            end
            legend(legenslabel)
        end
        xlab = repmat(xlab , 1 , dp);
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off' , 'XtickLabel' , xlab , 'XTick' , [0: 2*dp-1],...
            'XLim' ,  [-1 2*dp],'XTickLabelRotation' , 45)% , 'GridAlpha' , 1)
        title(['Movement Time for Different Exponential Decay Constants - thetastim = ' , num2str(theta_select)])
        
    case 'IPIs_vs_theta'
        %% plot satrt and end IPIs with different decay constants and different Theta stims IPIs.Aintegrate == 0.98
        % ********************** varying decay
        
        legenslabel = {};
        
        for ts = 1:length(Aintegrate)
            legenslabel = [legenslabel , ['theta stim = ' , num2str(theta_stim(ts))]];
            A = getrow(IPIs , IPIs.theta_stim == theta_stim(ts) & IPIs.Aintegrate == Aintegrate_select);
            ipiLabels = [0 1];
            x = 0;
            for dp = 1:length(DecayParam)
                ipis_e{ts,dp} = [A.E_start(dp) A.E_end(dp)];
                ipis_p{ts,dp} = [A.P_start(dp) A.P_end(dp)];
                X{ts,dp} = [x x+1];
                hold on
                x = x+2;
            end
        end
        xlab = {'IPI 1,2' , 'IPI 3:9' };
        figure('color' , 'white')
        subplot(211)
        for dp = 1:length(DecayParam)
            for ts = 1:length(Aintegrate)
                errorbar( X{ts,dp} , ipis_p{ts,dp} , ipis_e{ts,dp} , 'LineWidth' , 3 , 'color' , colorz{ts});
                hold on
            end
            legend(legenslabel)
        end
        xlab = repmat(xlab , 1 , dp);
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off' , 'XtickLabel' , xlab , 'XTick' , [0: 2*dp-1],...
            'XLim' ,  [-1 2*dp],'XTickLabelRotation' , 45)% , 'GridAlpha' , 1)
        title(['Movement Time for Different Exponential Decay Constants - Aintegrate = ', num2str(Aintegrate_select)])
        
        % for better visibility of the other theta stims
        subplot(212)
        for dp = 1:length(DecayParam)
            for ts = 2:length(Aintegrate)
                errorbar( X{ts,dp} , ipis_p{ts,dp} , ipis_e{ts,dp} , 'LineWidth' , 3 , 'color' , colorz{ts});
                hold on
            end
            legend(legenslabel(2:end))
        end
        xlab = repmat(xlab , 1 , dp);
        grid on
        ylabel('ms')
        xlabel('Exponential Decay Constant')
        set(gca, 'FontSize' , 20 ,'Box' , 'off' , 'XtickLabel' , xlab , 'XTick' , [0: 2*dp-1],...
            'XLim' ,  [-1 2*dp],'XTickLabelRotation' , 45)% , 'GridAlpha' , 1)
        title(['Movement Time for Different Exponential Decay Constants - Aintegrate = ', num2str(Aintegrate_select)])
        
end
