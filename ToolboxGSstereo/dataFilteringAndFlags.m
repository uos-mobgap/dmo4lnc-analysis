function [Data, MrkGaps, TurnsNotDec, timeFlag, t, FlagCheck] = dataFilteringAndFlags(Data, headerTraj, markers, fs, time, N)
    %% TimeFLAG - sections where the DMOs can be properly identified
    timeFlag = ones(1,size(time,2));

    %% Filter - O'Connor et al. 2007 - GEs identification
    fc = 7;                             % Cutoff frequency
    Wn = fc/(fs/2);                     % Normalized cutoff frequency
    [b,a] = butter(4, Wn, 'low');       % Butterworth filter

    MrkGaps = [];
    for i = 1:length(markers)
        checkNan = find(isnan(Data.(headerTraj).(markers{i})(:,1)));
        if isempty(checkNan)
            Data.(headerTraj).(markers{i}) = Data.(headerTraj).(markers{i});
            Data.(headerTraj).(markers{i}) = filtfilt(b,a, Data.(headerTraj).(markers{i})/1000); % in m and filter
        else
            MrkGaps = [MrkGaps; markers(i)];
            Gaps.((markers{i})) = checkNan;
            startStopGaps = [];
            s = find(diff(checkNan)>1);
            if ~isempty(s)
                for j = 1:size(s,1)
                    if j == 1
                        startStopGaps = [startStopGaps; checkNan(1),checkNan(s(j))];
                    else
                        startStopGaps = [startStopGaps; checkNan(s(j-1)+1),checkNan(s(j))];
                    end
                end
                startStopGaps = [startStopGaps; checkNan(s(j)+1),checkNan(end)];
            else
                startStopGaps = [startStopGaps; checkNan(1),checkNan(end)];
            end
            Gaps.StartStop.((markers{i})) = startStopGaps;

            Data.(headerTraj).(markers{i}) = Data.(headerTraj).(markers{i})/1000;  % in m
            % filter those part of the signal without gaps
            for k = 1:size(startStopGaps,1)
                if k ==1
                    if size(1:startStopGaps(1,1)-1,2)>50
                        Data.(headerTraj).(markers{i})(1:startStopGaps(1,1)-1,:) = filtfilt(b,a,Data.(headerTraj).(markers{i})(1:startStopGaps(1,1)-1,:));
                    end
                else
                    if size(startStopGaps(k-1,2)+1:startStopGaps(k,1)-1,2)>50
                        Data.(headerTraj).(markers{i})(startStopGaps(k-1,2)+1:startStopGaps(k,1)-1,:) = ...
                            filtfilt(b,a,Data.(headerTraj).(markers{i})(startStopGaps(k-1,2)+1:startStopGaps(k,1)-1,:));
                    end
                end
                if k == size(startStopGaps,1)
                    if startStopGaps(k,2)~=N && size(startStopGaps(k,2)+1:N,2)>50
                        Data.(headerTraj).(markers{i})(startStopGaps(k,2)+1:end,:) = ...
                            filtfilt(b,a,Data.(headerTraj).(markers{i})(startStopGaps(k,2)+1:end,:));
                    end
                end
            end
        end
    end

    %% Check for gaps on mrk traj on the pelvic segment
    mrkNow = [];
    TurnsNotDec = [];
    GapsSelMrk = [];
    if ~isempty(MrkGaps)

        %% Flags for GAPS on feets and pelvis
        selMrk = markers(contains(markers,'BACK')|contains(markers,'HEEL')| contains(markers,'TOE'),1);
        for i = 1:length(selMrk)
            GapsSelMrk.(selMrk{i}) = find(isnan(Data.(headerTraj).(selMrk{i})(:,1)));
        end

        %% sections of signal where turns and GEs with Zeni's algo cannot be assessed
        mrkDyn = MrkGaps(contains(MrkGaps,'BACK'));
        if ~isempty(mrkDyn)
            mrkDyn = mrkDyn(~contains(mrkDyn,'REF'));
        end
        x1 = [];        
        if length(mrkDyn) > 1          
            for i = 1:length(mrkDyn)
                mrkDNow = mrkDyn;
                mrkDNow(i)=[];
                x1 = [x1; GapsSelMrk.(mrkDyn{i})];
                for j = 1:length(mrkDNow)
                    x1 = [x1;intersect(GapsSelMrk.(mrkDyn{i}),GapsSelMrk.(mrkDNow{j}))];
                end
            end
            
            x1 = unique(sort(x1));
            fprintf('Gaps on %.2f%% of pelvic marker trajectories.\n',(length(x1)/N)*100);
            
            if ~isempty(x1)
                s = find(diff(x1)>1);
                if ~isempty(s)
                    for i = 1:size(s,1)
                        if i == 1
                            TurnsNotDec = [TurnsNotDec; x1(1),x1(s(i))];
                        else
                            TurnsNotDec = [TurnsNotDec; x1(s(i-1)+1),x1(s(i))];
                        end
                    end
                    TurnsNotDec = [TurnsNotDec; x1(s(i)+1),x1(end)];
                else
                    TurnsNotDec = [TurnsNotDec; x1(1),x1(end)];
                end
            end
            if ~isempty(TurnsNotDec)
                SizeGaps = (TurnsNotDec(:,2)-TurnsNotDec(:,1))/fs;
                FlagTrial = find(SizeGaps>0.5);
                if ~isempty(FlagTrial) && length(selMrk) > 1
                    fprintf('This trial should be discarded. Gaps larger than 0.5 found!\n');
                end
            end
        elseif isempty(mrkDyn)
            fprintf('Gaps on %.2f%% of pelvic marker trajectories.\n',(length(x1)/N)*100);
        else
            fprintf('Gaps only on ONE mkr among the pelvic marker trajectories!\n');
            fprintf('Gaps on %.2f%% of pelvic marker trajectories.\n',(length(x1)/N)*100);
        end
    end

    % Add this info into the TimeFLAG vector
    if ~isempty(TurnsNotDec)
        for i = 1:size(TurnsNotDec,1)
            timeFlag(TurnsNotDec(i,1):TurnsNotDec(i,2)) = zeros(1,TurnsNotDec(i,2)-TurnsNotDec(i,1)+1);
        end
    end

    %% sections where GEs cannot be found - gaps
    MrkFoot = markers(contains(markers,'HEEL')| contains(markers,'TOE'),1);
    x2 = [];
    if ~isempty(GapsSelMrk)
        for i = 1:length(MrkFoot)
            MrkFootNow = MrkFoot;
            MrkFootNow(i)=[];
            x2 = [x2; GapsSelMrk.(MrkFoot{i})];
            for j = 1:length(MrkFootNow)
                x2 = [x2;intersect(GapsSelMrk.(MrkFoot{i}),GapsSelMrk.(MrkFootNow{j}))];
            end
        end
    end
    x2 = unique(sort(x2));
    GEsNotDec = [];
    if ~isempty(x2)
        s = find(diff(x2)>1);
        if ~isempty(s)
            for i = 1:size(s,1)
                if i == 1
                    GEsNotDec = [GEsNotDec; x2(1),x2(s(i))];
                else
                    GEsNotDec = [GEsNotDec; x2(s(i-1)+1),x2(s(i))];
                end
            end
            GEsNotDec = [GEsNotDec; x2(s(i)+1),x2(end)];
        else
            GEsNotDec = [GEsNotDec; x2(1),x2(end)];
        end
    end

    % Add this info into the TimeFLAG vector
    if ~isempty(GEsNotDec)
        for i = 1:size(GEsNotDec,1)
            timeFlag(GEsNotDec(i,1):GEsNotDec(i,2)) = zeros(1,GEsNotDec(i,2)-GEsNotDec(i,1)+1);
        end
    end
    
    [t] = find(timeFlag);
    fprintf('Gaps on %.2f%% of relevant marker trajectories for DMO estimation.\n',((N- length(t))/N)*100);
    
    %% Data quality
    OverallQuality = checkDataQuality(Data.(headerTraj), fs);  
    
    if OverallQuality.Quality_Pelvis.Mean > 50 && (OverallQuality.Quality_Lfoot.Mean > 50 || OverallQuality.Quality_Rfoot.Mean > 50) %(isempty(find(OverallQuality.LGaps>2, 1))||isempty(find(OverallQuality.RGaps>2, 1)))
        % ok data processing
        FlagCheck = true;
        
    elseif (OverallQuality.Quality_Lfoot.Mean > 25 && OverallQuality.Quality_Rfoot.Mean > 25)
        % try data processing
        FlagCheck = true;
        
    else
        FlagCheck = false;
    end
    
    
%     if ((N- length(t))/N)*100> 85
%         FlagCheck = false;
%     else
%         FlagCheck = true;
%     end
end