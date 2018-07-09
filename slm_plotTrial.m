function slm_plotTrial(what,SIM,T,R,M,varargin)
%% function slm_plotTrial(SIM,T,R,M,varargin);
% Plots the horse race on the current sequence trial

vararginoptions(varargin,{'style'});

switch what
    case 'TrialHorseRace'
        figure('color' , 'white')
        [~,~,numPresses] = size(SIM.X);
        
        numDec = unique(T.decisionTime);
        for i = 1:length(numDec)
            T.decNum(T.decisionTime == numDec(i)) = i;
        end
        if size(SIM.B , 2)<length(T.stimTime)
            SIM.B = repmat(SIM.B , 1,length(T.stimTime));
        end
        for i=1:numPresses
            subplot(numPresses,1,i);
            plot(SIM.t,SIM.X(:,:,i));
            hold on;
            plot(SIM.t,SIM.B(:,i),'k', 'LineWidth' , 1.5);
            ylim = get(gca , 'YLim');
            h1 = line([T.stimTime(i) T.stimTime(i)] , ylim ,'color','r','linestyle',':' , 'LineWidth' , 2);
            h2 = line([T.decisionTime(i) T.decisionTime(i)] , ylim,'color','r', 'LineWidth' , 2);
            h3 = line([T.pressTime(i) T.pressTime(i)],ylim,'color','k', 'LineWidth' , 2);
            
            % just legend the first one, the rest are the same
            if ~isnan(T.forcedPressTime(1,1))
                h4 = line([0 0],ylim,'color','b', 'LineWidth' , 1);
                line([800 800],ylim,'color','b', 'LineWidth' , 1);
                line([1600 1600],ylim,'color','b', 'LineWidth' , 1);
                line([2400 2400],ylim,'color','b', 'LineWidth' , 1);
                h5 = line([2300 2300],ylim,'color','g', 'LineWidth' , 1);
                line([2500 2500],ylim,'color','g', 'LineWidth' , 1);
                legend([h1 h2 h3 h4 h5],{'Stimulus came on' , 'Decision boundry reached' , 'Press executed',...
                    'Forced-RT tones','Response window'})
            else
                if i==1
                    legend([h1 h2 h3],{'Stimulus came on' , 'Decision boundry reached' , 'Press executed'})
                end
            end
 
%             title(['Decision No. ' ,num2str(T.decNum(i)), ', press No.' , num2str(i)])
            set(gca , 'Box' , 'off' , 'FontSize' , 10)
        end;
        
    case 'TrialHorseRace_pres'
%         clear tempcol
%         c1 = [153, 0, 0]/255; % structured blue tones
%         ce = [102, 204, 255]/255;
%         for rgb = 1:3
%             colz(:, rgb) = linspace(c1(rgb),ce(rgb)' , 5);
%         end
        colz = [83, 255, 26;26, 255, 255;255, 26, 255;255, 26, 26;255, 153, 0]/255;
%         colz = repmat(linspace(0, 200 , 5)' , 1 , 3)/255;
        figure('color' , 'white')
        [~,~,numPresses] = size(SIM.X);
        
        numDec = unique(T.decisionTime);
        for i = 1:length(numDec)
            T.decNum(T.decisionTime == numDec(i)) = i;
        end
        if size(SIM.B , 2)<length(T.stimTime)
            SIM.B = repmat(SIM.B , 1,length(T.stimTime));
        end
        % for better visualization of the noise free state so that the traces are not super-imposed
        addToFing = [0 -.05 -.1 -.15 -.2];
        addToFing = [0 0 0 0 0];
        B(1,:,:) = SIM.B; 
        SIM.B = repmat(B , 5,1,1);
        for opt = 1:5
            SIM.X(opt,:,:) =SIM.X(opt,:,:)+addToFing(opt);
            SIM.B(opt,:,:) = SIM.B(opt,:,:)+addToFing(opt);
        end
        for i=1:3
            subplot(3,1,i);
            for opts = 1:5
                plot(SIM.t,SIM.X(opts,:,i), 'LineWidth' , 1,'color' , colz(opts , :));
                hold on
            end
            plot(SIM.t,SIM.B(T.response(i),:,i),'color', [184 184 184]/255 ,'LineWidth' , 1);
            ylim = get(gca , 'YLim');
            h1 = line([T.stimTime(i) T.stimTime(i)] , ylim ,'color','r','linestyle',':' , 'LineWidth' , 1);
            h2 = line([T.decisionTime(i) T.decisionTime(i)] , ylim,'color','r', 'LineWidth' , 1);
            h3 = line([T.pressTime(i) T.pressTime(i)],ylim,'color','k', 'LineWidth' , 1);
            
            % just legend the first one, the rest are the same
            if ~isnan(T.forcedPressTime(1,1))
                h4 = line([0 0],ylim,'color','b', 'LineWidth' , 1);
                line([800 800],ylim,'color','b', 'LineWidth' , 1);
                line([1600 1600],ylim,'color','b', 'LineWidth' , 1);
                line([2400 2400],ylim,'color','b', 'LineWidth' , 1);
                h5 = line([2300 2300],ylim,'color','g', 'LineWidth' , 1);
                line([2500 2500],ylim,'color','g', 'LineWidth' , 1);
                legend([h1 h2 h3 h4 h5],{'Stimulus came on' , 'Decision boundry reached' , 'Press executed',...
                    'Forced-RT tones','Response window'})
            else
                if i==1
                    legend([h1 h2 h3],{'Stimulus came on' , 'Decision boundry reached' , 'Press executed'})
                end
            end
 
%             title(['Decision No. ' ,num2str(T.decNum(i)), ', press No.' , num2str(i)])
            set(gca , 'Box' , 'off' , 'FontSize' , 7 , 'Ylim' , [-.4 1.5],'ColorOrder',colz , 'Xlim' , [0 2000])
        end;
    case 'BlockMT'
        %colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3]};
        for i = 1:length(T.TN)
            T.H(i,1) = unique(unique(T.Horizon(i , ~isnan(T.Horizon(i,:)))));
        end
        S = tapply(T , {'H' , 'bufferSize'} , {'MT' , 'nanmean'});
        
        bs = unique(S.bufferSize);
        H = unique(S.H);
        figure('color' , 'white')
        subplot(211)
        if length(bs)>1 && length(H)>1
            [x, ~ , ~] = lineplot([T.H , T.bufferSize]  ,  T.MT , 'style_shade');
            count = 1;
            for i = 1 :length(x)/length(H) : length(x)
                text(x(i)+1 , max(T.MT) , ['Horizon = ' , num2str(H(count))])
                count = count +1;
            end
            xlabel('Buffer size')
        elseif length(bs)==1 && length(H)>1
            [~, ~ , ~] = lineplot([T.H]  ,  T.MT , 'style_shade');
            xlabel('Horizon size')
        elseif length(bs)>1 && length(H)==1
            [~, ~ , ~] = lineplot([T.bufferSize]  ,  T.MT , 'style_shade');
            xlabel('Buffer size')
        elseif length(bs)==1 && length(H)==1
            [~, ~ , ~] = lineplot([T.bufferSize]  ,  T.MT , 'style_shade');
            xlabel('Buffer size')
        end
        %P = reshape(p  , length(bs),length(H));
        %E = reshape(e , length(bs),length(H));
        
        grid on
        
        title('Movement Time in Correct Trials')
        set(gca , 'Box' , 'off' , 'FontSize', 20)
        subplot(212)
        leg = cell(1,length(bs));
        for i = 1:length(bs)
            %hi = errorbar([1:length(H)] , P(i,:) , E(i,:) , 'LineWidth' , 3 , 'color' , colorz{i})
            leg{i} = {'BufferSize = ' , num2str(bs(i))};
            hold on
        end
        grid on
        xlabel('Horizon size')
        ylabel('msec')
        legend(leg)
        set(gca , 'Box' , 'off' , 'FontSize', 20)
    case 'IPIHist'
        figure('color' , 'white')
        histogram(T.pressTime(~T.isError , :));
        hold on
        histogram(T.pressTime(logical(T.isError) , :));
        legend({'Correct Trials', 'Error Trials'})
        title('Distribution of Press Times')
    case 'IPIHorizon'
        colorz = {[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0],[.3 .3 .3],[.4 0 .4]};
        for i = 1:length(T.TN)
            T.H(i,1) = unique(unique(T.Horizon(i , ~isnan(T.Horizon(i,:)))));
        end
        T.IPI = diff(T.pressTime,2);
        S = tapply(T , {'H' , 'bufferSize'} , {'MT' , 'nanmean'} , {'IPI' , 'nanmedian'});
        
        bs = unique(S.bufferSize);
        H = unique(S.H);
        figure('color' , 'white')
        if length(bs)>1 && length(H)>1
            buff = input('which buffer size?');
            S = getrow(S , S.bufferSize == buff);
            leg=cell(1,length(H));
            for h = 1:length(H)
                plot(S.IPI(S.H==H(h), :) , 'Linewidth' , 3 , 'color' , colorz{H(h)});
                hold on
                leg{h} = {'Horizon = ' , num2str(H(h))};
            end
            xlabel('IPI number')
            title(['Inter-press intervals  - buffer size = ' , num2str(buff)])
            
        elseif length(bs)==1 && length(H)>1
            leg=cell(1,length(H));
            for h = 1:length(H)
                plot(S.IPI(S.H==H(h) , :) , 'Linewidth' , 3 , 'color' , colorz{H(h)});
                hold on
                leg{h} = {'Horizon = ' , num2str(H(h))};
            end
            xlabel('IPI number')
            title('Inter-press intervals')
        elseif length(bs)>=1 && length(H)==1
            leg=cell(1,length(bs));
            for buff = 1:length(bs)
                plot(S.IPI(S.bufferSize==bs(buff) , :) , 'Linewidth' , 3 , 'color' , colorz{bs(buff)});
                hold on
                leg{buff} = {'buffer size = ' , num2str(bs(buff))};
            end
            xlabel('IPI number')
            title('Inter-press intervals')
        end
        set(gca , 'Box' , 'off' , 'FontSize', 20)
        grid on
        legend(leg)
    case 'plotPrep'
        figure;
        if all(R.isError==0) || all(R.isError==1)
            plt.line(R.prepTime,(1-R.isError)*100,'errorbars','shade');
            xlabel('prep time'); ylabel('accuracy %');
        else
            R1=tapply(R,{'SN','prepTime'},{R.isError,'nanmean','name','ER','subset',R.timingError==0});
            plt.line(R1.prepTime,(1-R1.ER)*100,'errorbars','shade');
            xlabel('prep time'); ylabel('accuracy %'); axis square; title(sprintf('Aint: %1.3f, Theta: %1.4f, Noise: %1.3f, Initial bound: %1.2f',M.Aintegrate,M.theta_stim,M.SigEps,M.Bound/2));
            figure;
            plt.hist(R.pressTime(:,1)); 
            xlabel('RT (ms)'); ylabel('Distribution (n trials)'); axis square; title(sprintf('Aint: %1.4f, Theta: %1.4f, Noise: %1.4f, Initial bound: %1.4f',M.Aintegrate,M.theta_stim,M.SigEps,M.Bound/2));
        end
end