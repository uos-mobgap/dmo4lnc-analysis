function [Euler_Angles, TurnM, TurnDur, act_dopc]= identifyTurnsFromFeet_eScience(Data, headerTraj, fs, turnThres, side, N)
        
    %% sections of signal where turns from feet cannot be assessed
    markers = fieldnames(Data.(headerTraj));
    k = 1;
    for i = 1:length(markers)
        if contains(markers{i},side)
            if strcmp(side, 'l')
                if strcmp(markers{i},strcat(side,'_toe'))|| strcmp(markers{i},strcat(side,'_heel'))||...
                        strcmp(markers{i},strcat(side,'_ank'))
                    mrkFoot{k} = markers{i};
                    k = k+1;
                end
            else
                if  strcmp(markers{i},strcat(side,'_toe'))|| strcmp(markers{i},strcat(side,'_heel'))||...
                        strcmp(markers{i},strcat(side,'_ank'))
                    mrkFoot{k} = markers{i};
                    k = k+1;
                end
            end
        end
    end
    for i = 1:length(mrkFoot)
        GapsSelMrk.(mrkFoot{i}) = find(isnan(Data.(headerTraj).(mrkFoot{i})(:,1)));
    end
    x1 = [];
    for i = 1:length(mrkFoot)
        mrkDNow = mrkFoot;
        mrkDNow(i)=[];
        x1 = [x1; GapsSelMrk.(mrkFoot{i})];
        for j = 1:length(mrkFoot)
            x1 = [x1;intersect(GapsSelMrk.(mrkFoot{i}),GapsSelMrk.(mrkFoot{j}))];
        end
    end
    x1 = unique(sort(x1));
    fprintf('Gaps on %.2f%% of %s foot marker trajectories.\n',(length(x1)/N), side);
    TurnsNotDec = [];
    if ~isempty(x1)
        s = find(diff(x1)>1);
        if isempty(s)
            TurnsNotDec = [TurnsNotDec; x1(1),x1(end)];
        else
            for i = 1:size(s,1)
                if i == 1
                    TurnsNotDec = [TurnsNotDec; x1(1),x1(s(i))];
                else
                    TurnsNotDec = [TurnsNotDec; x1(s(i-1)+1),x1(s(i))];
                end
            end
        end
    end

    %% Mrks on the foot segment are filled but this information will not be used in the sections identified above
    for k = 1:length(mrkFoot)
        Data.(headerTraj).(mrkFoot{k}) = fillgaps(Data.(headerTraj).(mrkFoot{k}));
    end

    %% Define a LRF using the skin-markers placed on each foot
    for i = 1:N
        xAP_l(i,:) = (Data.(headerTraj).(strcat(side,'_toe'))(i,:) - Data.(headerTraj).(strcat(side,'_heel'))(i,:))/...
            norm(Data.(headerTraj).(strcat(side,'_toe'))(i,:) - Data.(headerTraj).(strcat(side,'_heel'))(i,:));               % Anterior-Posterior axis

        v(i,:) = (Data.(headerTraj).(strcat(side,'_ank'))(i,:) - Data.(headerTraj).(strcat(side,'_heel'))(i,:))/...
            norm(Data.(headerTraj).(strcat(side,'_ank'))(i,:) - Data.(headerTraj).(strcat(side,'_heel'))(i,:));

        zML_l(i,:) = cross(xAP_l(i,:),v(i,:))/norm(cross(xAP_l(i,:),v(i,:)));               % Medio-lateral axis

        yV_l(i,:) = cross(zML_l(i,:), xAP_l(i,:))/norm(cross(zML_l(i,:), xAP_l(i,:)));      % Vertical axis

        R(:,:,i) = [xAP_l(i,:)' yV_l(i,:)' zML_l(i,:)'];                                    % LRF Rotation Matix

        CosTheta(i) = dot(yV_l(i,:),[0 0 1])/(norm(yV_l(i,:))*norm([0 0 1]));               % Define inclination of the LRF with respect to the vertical (GRF)
        ThetaInDegrees(i) = acosd(CosTheta(i));
        Rm = rotz(ThetaInDegrees(i));
        R_align(:,:,i) = R(:,:,i)*Rm;                                                        % Alignment of the LRF with respect to the vertical direction
    end

    %% Turns
    % Tetha1: Flexion-Extension (about Z proximal SCS axis)
    Euler_Angles(:,1) = atan2(R(2,1,:),R(1,1,:))*180/pi;
    Euler_Angles(:,1) = Euler_Angles(:,1)- mean(Euler_Angles(1:5,1));     % remove offset

    % Remove singularities
    for i = 1:length(Euler_Angles)-1
        if abs(Euler_Angles(i+1)-Euler_Angles(i))>179
            if (Euler_Angles(i+1)-Euler_Angles(i))<0
                Euler_Angles(i+1:end) = Euler_Angles(i+1:end)+ abs(Euler_Angles(i+1)-Euler_Angles(i));
            else
                Euler_Angles(i+1:end) = Euler_Angles(i+1:end)- abs(Euler_Angles(i+1)-Euler_Angles(i));
            end
        end
    end

    % Mellone thesis && El-Gohary et al., 2014 - adapted for feet
    fc = 1;                             % Cutoff frequency
    Wn = fc/(fs/2);                     % Normalized cutoff frequency
    [b,a] = butter(6, Wn, 'low');       % Butterworth filter
    MaxVal = 15;
    CrossVal = 5;
    Euler_Angles = filtfilt(b,a,Euler_Angles);

    AngVel = diff(Euler_Angles)/(1/fs);
    max_dopc = 180;
    [~,posMaxDOPC] = findpeaks(abs(AngVel),'MinPeakHeight',max_dopc);

    %% remove max change of direction found in the region where the signal has been filled
    act_dopc = posMaxDOPC;
    if ~isempty(posMaxDOPC)
        for i = 1:size(TurnsNotDec,1)
            [~,d] = intersect(act_dopc,TurnsNotDec(i,1):1:TurnsNotDec(i,2));
            act_dopc(d)=[];
        end
    end

    [~,locMax] = findpeaks(abs(AngVel),'MinPeakHeight',MaxVal);
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    zx = zci(AngVel-CrossVal);
    Length0 = length(zx);
    while length(zx) < length(locMax)+2
        CrossVal = CrossVal + 1;
        zx_now = zci(abs(AngVel) - CrossVal);
        LengthNow = length(zx_now);
        if  LengthNow < Length0
            break;
        else
            zx = zx_now;
        end
    end

    x1 = [];
    x2 = [];
    Turns_m = [];
    Turns_d = [];
    Extremes = [];
    for i = 1:length(locMax)
        bef = zx<locMax(i);
        if isempty(find(bef,1))
            x1 = [x1;1];
        else
            bef = zx(bef);
            x1 = [x1;bef(end)];
        end
        aft = zx>locMax(i);
        if isempty(find(aft,1))
            x2 = [x2;length(AngVel)];
        else
            aft = zx(aft);
            x2 = [x2;aft(1)];
        end

        turnNow = Euler_Angles(x1(end))-Euler_Angles(x2(end));
        Turns_m = [Turns_m; turnNow];
        turn_d = (x2(end)- x1(end))/fs;
        Turns_d = [Turns_d;turn_d];
        Extremes = [Extremes; x1(end),x2(end)];
    end

    [Turns_m, i] = unique(Turns_m, 'stable');
    Turns_d = Turns_d(i);
    Extremes = Extremes(i,:);

    %%  Turn definition as described in the article
    allChanges = unique(sort([x1;x2]));

    TurnM = [];
    TurnDur = [];
    for i = 2:length(allChanges)
        turnNow = Euler_Angles(allChanges(i))-Euler_Angles(allChanges(i-1));
        turn_d = (allChanges(i)- allChanges(i-1))/fs;
        if (abs(turnNow) > turnThres) && (turn_d > 0.5)
            TurnM = [TurnM; turnNow];                               % turn magnitude - in degrees
            TurnDur = [TurnDur; allChanges(i-1), allChanges(i)];    % beginning and end of each turn (frames)                       
        elseif ((turn_d < 0.05) && (abs(turnNow) < turnThres)) && i<length(allChanges)
            j = i+1;
            turnNow = Euler_Angles(allChanges(j))-Euler_Angles(allChanges(i-1));
            turn_d = (allChanges(j)- allChanges(i-1))/fs;
            while (turn_d < 0.05) && (abs(turnNow) < turnThres)
                turnNow = Euler_Angles(allChanges(j))-Euler_Angles(allChanges(i-1));
                turn_d = (allChanges(j)- allChanges(i-1))/fs;
                j = j+1;
            end
            if (abs(turnNow) > turnThres) && (turn_d > 0.5)
                TurnM = [TurnM; turnNow];                             % turn magnitude - in degrees
                TurnDur = [TurnDur; allChanges(i-1), allChanges(i)];  % beginning and end of each turn (frames)
            end
        end
    end
end