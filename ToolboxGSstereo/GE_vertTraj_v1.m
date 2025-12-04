function [IC, FC] = GE_vertTraj_v1(Data, headerTraj, time, FootMrk, side, fs, fig)

    % marker trajectories
    H = Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3);
    T = Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3);
    %
    minDist = 0.7;
    nMinDist = fs*minDist;
    [~,FC_locs] =...
        findpeaks((-T*1000),'MinPeakDistance',nMinDist,'MinPeakProminence',0.2,'Threshold',0.001);
    [~,IC_locs] = findpeaks((-H*1000),'MinPeakDistance',nMinDist,'MinPeakProminence',0.18);
    
    %% Remove extra events
    % HS before each peak of the heel mrk traj
    IC = [];
    [~,H_peaks] = findpeaks(H,'MinPeakHeight',.3*range(H),'MinPeakDistance',nMinDist);
    if length(H_peaks)<length(IC_locs) % extra events
        for i = 1:length(H_peaks)
            sol_now = find(IC_locs<H_peaks(i));
            if ~isempty(sol_now)
                IC = [IC;IC_locs(sol_now(end))];
            end
        end
    else
        IC = IC_locs;
    end
    
    % TO after each peak of the toe mrk traj
    FC = [];
    [~,T_peaks] = findpeaks(T,'MinPeakHeight',.2*range(T),'MinPeakDistance',nMinDist);
    if length(T_peaks)<length(FC_locs) % extra events
        for i = 1:length(T_peaks)
            if i <= length(T_peaks)
                sol_now = find(FC_locs<T_peaks(i));
                if ~isempty(sol_now)
                    FC = [FC;FC_locs(sol_now(end))];
                end
            end
        end
    else
        FC = FC_locs;
    end
    
    %% For each swing, check the following IC
    deltaSwingfToIC = nan(size(H_peaks));
    for i = 1:size(deltaSwingfToIC,1)
        posNow = find(IC(:,1)<H_peaks(i),1);
        if ~isempty(posNow)
            deltaSwingfToIC(i) = abs(IC(posNow(end),1)-H_peaks(i));
        end
    end
    Mean_deltaSwingfToIC = nanmean(deltaSwingfToIC);
    STD_deltaSwingfToIC = nanstd(deltaSwingfToIC);
    posToCheck = find(deltaSwingfToIC>Mean_deltaSwingfToIC+STD_deltaSwingfToIC | deltaSwingfToIC<Mean_deltaSwingfToIC-STD_deltaSwingfToIC);
    if ~isempty(posToCheck)
        for k = 1:size(posToCheck,1)
            if posToCheck(k) == 1
                startAOI = 1;
            else
                startAOI = H_peaks(posToCheck(k)-1);
            end
            sigNow = H(startAOI:H_peaks(posToCheck(k)));
            pos = find(islocalmin(sigNow)); % min(Hnow(posM:end));
            if length(pos) ==1
                pos = startAOI + pos;
                ICnowPos = find(IC(:,1)>startAOI,1);
                if pos ~= IC(ICnowPos)
                    IC(ICnowPos,1)=pos;                
                end
            else
                % check
                hAll = zeros(size(pos));
                for p = 1:length(pos)
                    hAll(p)= H(startAOI+pos(p));
                end                
                check = true;
                [~,Greatmin] = sort(hAll);
                m = 1;
                while (check)                    
                    if pos(Greatmin(m))< H_peaks(posToCheck(k))-(startAOI+pos(Greatmin(m)))
                        % OK!
                        posNow = startAOI + pos(Greatmin(m));
                        ICnowPos = find(IC(:,1)>startAOI,1);
                        if posNow ~= IC(ICnowPos)
                            IC(ICnowPos,1)=posNow;
                        end
                        check = false;
                    elseif m<length(Greatmin)
                        m = m+1;
                    else
                        check = false;
                    end
                end
            end            
        end
    end
    
     %% FIGURE
    if fig == 1
        figure('Name',strcat(side,'MrkTraj'),'NumberTitle','off')
        plot(time(1:end), T)
        hold on
        plot(time(1:end), H, 'r')

        scatter((FC_locs/fs),T(FC_locs),'v')
        text((FC_locs(1)/fs),T(FC_locs(1)),'TO')
        scatter((IC_locs/fs),H(IC_locs),'^')
        text((IC_locs(1)/fs),H(IC_locs(1)),'HS')
        ylabel(strcat(side, ' Vertical Mark Traj [m]'))
        xlabel('Time [s]')
        sz = 40;
        scatter((FC/fs),T(FC), sz,'MarkerEdgeColor',[0 .5 .5],...
            'MarkerFaceColor',[0 .7 .7],...
            'LineWidth',1.5)
        scatter((IC/fs),H(IC), sz,'MarkerEdgeColor',[0 .5 .5],...
            'MarkerFaceColor',[0 .7 .7],...
            'LineWidth',1.5)
        
        legend(FootMrk{1}, FootMrk{2})
        ylabel('Marker vetical displacement [m]')
        xlabel('Time [s]')
    end
    
    %% ADD A FLAG
    IC(:,2) = zeros(size(IC));
    FC(:,2) = zeros(size(FC));
 end
