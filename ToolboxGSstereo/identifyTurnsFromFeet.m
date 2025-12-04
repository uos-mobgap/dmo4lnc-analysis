function [Euler_Angles, TurnM, TurnDur, Euler_Angles_2]= identifyTurnsFromFeet(Data, headerTraj, fs, turnThres, side)
        
    %% sections of signal where turns from feet cannot be assessed
    markers = fieldnames(Data.(headerTraj));
    k = 1;
    for i = 1:length(markers)
        if contains(markers{i},side)
            if strcmp(side, 'L')
                if ~contains(markers{i},'WRI')&& ~strcmp(markers{i},'RHEEL')
                    mrkFoot{k} = markers{i};
                    k = k+1;
                end
            else
                if ~contains(markers{i},'WRI')&& ~contains(markers{i},'BACK')&& ~strcmp(markers{i},'LREF')
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
    X = sprintf('Gaps on %.2f%% of %s foot marker trajectories.',(length(x1)/length(Data.(headerTraj).LHEEL))*100, side);
    disp(X)
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
    
    %% POSSIBLE mrk FILL
    check = zeros(length(Data.(headerTraj).(strcat(side,'TOE'))),5);
    for i = 1:length(Data.(headerTraj).(strcat(side,'TOE')))
        for k = 1:length(mrkFoot)
            if ~isnan(Data.(headerTraj).(mrkFoot{k})(i,1))
                check(i,k) = 1;
            end
        end
        check(i,5) = sum(check(i,1:4));
    end
    
    % at least other 3 mrk trajectories have to be available
    x3 = find(check(:,5)==3);
    StartStopGaps3 = [];
    if ~isempty(x3)
        s = find(diff(x3)>1);
        if isempty(s)
            StartStopGaps3 = [StartStopGaps3; x3(1),x3(end)];
        else
            for i = 1:size(s,1)
                if i == 1
                    StartStopGaps3 = [StartStopGaps3; x3(1),x3(s(i))];
                else
                    StartStopGaps3 = [StartStopGaps3; x3(s(i-1)+1),x3(s(i))];
                end
            end
        end
    end
    for g = 1:size(StartStopGaps3,1)
        % Rigid body recontruction
        if g ==1
            fprintf('Only gaps on %d marker on the %s foot segment. RB reconstruction is possible!\n',1,side);
        end
        % find mrk with GAPs
        mrkNowPos = find(check(StartStopGaps3(g,1),:)==0);
        mrksPos = find(check(StartStopGaps3(g,1),:)==1);
        gap = find(isnan(Data.(headerTraj).(mrkFoot{mrkNowPos})(:,1)));
        newTraj = techFrame(Data.(headerTraj),{mrkFoot{mrksPos}}, mrkFoot{mrkNowPos}, gap);
        Data.(headerTraj).(mrkFoot{mrkNowPos})(StartStopGaps3(1):StartStopGaps3(2),:) = newTraj(StartStopGaps3(1):StartStopGaps3(2),:);
    end

    %% Mrks on the foot segment are filled but this information will not be used in the sections identified above
    for k = 1:length(mrkFoot)
        Data.(headerTraj).(mrkFoot{k}) = fillgaps(Data.(headerTraj).(mrkFoot{k}));
    end

    %% Define a LRF using the skin-markers placed on each foot
    for i = 1:length(Data.(headerTraj).(strcat(side,'TOE')))
        xAP_l(i,:) = (Data.(headerTraj).(strcat(side,'TOE'))(i,:) - Data.(headerTraj).(strcat(side,'HEEL'))(i,:))/...
            norm(Data.(headerTraj).(strcat(side,'TOE'))(i,:) - Data.(headerTraj).(strcat(side,'HEEL'))(i,:));               % Anterior-Posterior axis

        v(i,:) = (Data.(headerTraj).(strcat(side,'REF'))(i,:) - Data.(headerTraj).(strcat(side,'HEEL'))(i,:))/...
            norm(Data.(headerTraj).(strcat(side,'REF'))(i,:) - Data.(headerTraj).(strcat(side,'HEEL'))(i,:)); % changed to ref from INDIP

        zML_l(i,:) = cross(xAP_l(i,:),v(i,:))/norm(cross(xAP_l(i,:),v(i,:)));               % Medio-lateral axis

        yV_l(i,:) = cross(zML_l(i,:), xAP_l(i,:))/norm(cross(zML_l(i,:), xAP_l(i,:)));      % Vertical axis

        R(:,:,i) = [xAP_l(i,:)' yV_l(i,:)' zML_l(i,:)'];                                    % LRF Rotation Matix

        CosTheta(i) = dot(yV_l(i,:),[0 0 1])/(norm(yV_l(i,:))*norm([0 0 1]));               % Define inclination of the LRF with respect to the vertical (GRF)
        ThetaInDegrees(i) = acosd(CosTheta(i));
        
        if ~isnan(ThetaInDegrees(i))
            Rm = rotz(ThetaInDegrees(i));
        else
            %% if there are gaps, previous information is used!!!
            if ~isempty(find(~isnan(ThetaInDegrees), 1))
                posNow = find(~isnan(ThetaInDegrees)); posNow = posNow(end);
                Rm = rotz(ThetaInDegrees(posNow));    
            else
                Rm = rotz(0); 
            end
            %%
        end
        R_align(:,:,i) = R(:,:,i)*Rm;
    end

    %% Turns
    Euler_Angles(:,1) = atan2(R(2,1,:),R(1,1,:))*180/pi;
    Euler_Angles_2(:,1) = asin(R_align(3,1,:))*180/pi;
    Euler_Angles_2(:,2) = atan2(R_align(3,2,:),R_align(3,3,:)); 

    Euler_Angles(:,1) = Euler_Angles(:,1)- mean(Euler_Angles(1:5,1));     % remove offset

    % Remove singularities
    for i = 1:length(Euler_Angles)-1
        if abs(Euler_Angles(i+1)-Euler_Angles(i))>179
            if (Euler_Angles(i+1)-Euler_Angles(i))<0
                Euler_Angles(i+1:end) = Euler_Angles(i+1:end)+ abs(Euler_Angles(i+1)-Euler_Angles(i));
            else
                Euler_Angles(i+1:end) = Euler_Angles(i+1:end)- abs(Euler_Angles(i+1)-Euler_Angles(i));
            end
            if (Euler_Angles_2(i+1,1)-Euler_Angles_2(i,1))<0
                Euler_Angles_2(i+1:end,1) = Euler_Angles_2(i+1:end,1)+ abs(Euler_Angles_2(i+1,1)-Euler_Angles_2(i,1));
            else
                Euler_Angles_2(i+1:end,1) = Euler_Angles_2(i+1:end,1)- abs(Euler_Angles_2(i+1,1)-Euler_Angles_2(i,1));
            end
        end
    end
    
    % Mellone thesis && El-Gohary et al., 2014 - adapted for feet
    fc = 1;                             % Cutoff frequency 
    Wn = fc/(fs/2);                     % Normalized cutoff frequency        
    [b,a] = butter(4, Wn, 'low');       % Butterworth filter
    Euler_Angles = filtfilt(b,a,Euler_Angles); 
    
    AngVel = diff(Euler_Angles)/(1/fs);
    time = 1/fs:1/fs:length(AngVel)/fs;
    
    [Turns] = gyro_sgm_v1(Euler_Angles, AngVel, abs(AngVel), time, fs, turnThres);
    TurnM = Turns.Magnitude;
    TurnDur = Turns.TurnDur;  
end