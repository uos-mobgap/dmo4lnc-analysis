function [IC, FC, R_align, MrkL, foreFoot, checkOutputs, deltaEventHS, deltaEventTO] = GE_Zeni_v4_new(Data, headerTraj,time, FootMrk, side, fs, mrkDyn, fig, TurnsNotDec, MrkGaps, InclinationFootStatic, InclinationFootGait, GE_INDIP, twindow)

%% GE detection based on Zeni et al. (2008)
% Zeni Jr, J. A., Richards, J. G., & Higginson, J. S. (2008).
% Two simple methods for determining gait events during treadmill and
% overground walking using kinematic data.
% Gait & posture, 27(4), 710-714.
% https://doi.org/10.1016/j.gaitpost.2007.07.007
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% New version - 18.12.2020
% check if the foot-to-ground contact has been performed with the forefoot!
%
%%
%% New version - 21.12.2020
% Add info about foot inclination from gait performed during the
% personalisaztion trial
%
%%

%% Mrks on the pelvic segment are filled
if isempty(TurnsNotDec) && length(MrkGaps) == 1 && ~contains(MrkGaps,'REF') && ~contains(MrkGaps,'WRIST') && contains(MrkGaps,'DYN')
    check = false;
else
    check = true;
end

for k = 1:length(mrkDyn)
    Data.(headerTraj).(mrkDyn{k}) = fillgaps(Data.(headerTraj).(mrkDyn{k}));
end

%% Define a LRF using the skin-markers placed on the MIMU located on the Lower Back
if check
    for i = 1:length(Data.(headerTraj).DYNA0)
        zML(i,:) = (Data.(headerTraj).DYNAX(i,:) - Data.(headerTraj).DYNA0(i,:))/...
            norm(Data.(headerTraj).DYNAX(i,:) - Data.(headerTraj).DYNA0(i,:));          % Medio-lateral axis
        v(i,:) = (Data.(headerTraj).DYNAY(i,:) - Data.(headerTraj).DYNA0(i,:))/...
            norm(Data.(headerTraj).DYNAY(i,:) - Data.(headerTraj).DYNA0(i,:));
        xAP(i,:) = cross(v(i,:),zML(i,:))/norm(cross(v(i,:),zML(i,:)));                 % Anterior-Posterior axis
        
        yV(i,:) = cross(zML(i,:), xAP(i,:))/norm(cross(zML(i,:), xAP(i,:)));            % Vertical axis
        
        R(:,:,i) = [xAP(i,:)' yV(i,:)' zML(i,:)'];                                      % LRF Rotation Matix
        
        CosTheta(i) = dot(yV(i,:),[0 0 1])/(norm(yV(i,:))*norm([0 0 1]));               % Define inclination of the LRF with respect to the vertical (GRF)
        ThetaInDegrees(i) = acosd(CosTheta(i));
        
        if ~isnan(ThetaInDegrees(i))
            Rm = rotz(ThetaInDegrees(i));
        else
            %% if there are gaps, previous information is used!!!
            posNow = find(~isnan(ThetaInDegrees));
            if ~isempty(posNow)
                posNow = posNow(end);
                Rm = rotz(ThetaInDegrees(posNow));
            else
                Rm = rotz(90);
            end
            
            %%
        end
        R_align(:,:,i) = R(:,:,i)*Rm;                                                   % Alignment of the LRF with respect to the vertical direction
    end
    t = Data.(headerTraj).DYNA0;                                                        % Origin of the LRF
    
else
    if strcmp(MrkGaps, 'DYNAY')
        for i = 1:length(Data.(headerTraj).DYNA0)
            zML(i,:) = (Data.(headerTraj).DYNAX(i,:) - Data.(headerTraj).DYNA0(i,:))/...
                norm(Data.(headerTraj).DYNAX(i,:) - Data.(headerTraj).DYNA0(i,:));          % Medio-lateral axis
            v(i,:) = (Data.(headerTraj).DYNAREF(i,:) - Data.(headerTraj).DYNA0(i,:))/...
                norm(Data.(headerTraj).DYNAREF(i,:) - Data.(headerTraj).DYNA0(i,:));
            xAP(i,:) = cross(v(i,:),zML(i,:))/norm(cross(v(i,:),zML(i,:)));                 % Anterior-Posterior axis
            
            yV(i,:) = cross(zML(i,:), xAP(i,:))/norm(cross(zML(i,:), xAP(i,:)));            % Vertical axis
            
            R(:,:,i) = [xAP(i,:)' yV(i,:)' zML(i,:)'];                                      % LRF Rotation Matix
            
            CosTheta(i) = dot(yV(i,:),[0 0 1])/(norm(yV(i,:))*norm([0 0 1]));               % Define inclination of the LRF with respect to the vertical (GRF)
            ThetaInDegrees(i) = acosd(CosTheta(i));
            Rm = rotz(ThetaInDegrees(i));
            R_align(:,:,i) = R(:,:,i)*Rm;                                                   % Alignment of the LRF with respect to the vertical direction
        end
        t = Data.(headerTraj).DYNA0;                                                        % Origin of the LRF
        
    elseif strcmp(MrkGaps, 'DYNA0')
        warning('This condition has not been properly tested yet')
        %             for i = 1:length(Data.(headerTraj).DYNAREF)
        %                 zML(i,:) = (Data.(headerTraj).DYNAX(i,:) - Data.(headerTraj).DYNAREF(i,:))/...
        %                     norm(Data.(headerTraj).DYNAX(i,:) - Data.(headerTraj).DYNAREF(i,:));          % Medio-lateral axis
        %                 v(i,:) = (Data.(headerTraj).DYNAY(i,:) - Data.(headerTraj).DYNA0(i,:))/...
        %                     norm(Data.(headerTraj).DYNAY(i,:) - Data.(headerTraj).DYNA0(i,:));
        %                 xAP(i,:) = cross(v(i,:),zML(i,:))/norm(cross(v(i,:),zML(i,:)));                 % Anterior-Posterior axis
        %
        %                 yV(i,:) = cross(zML(i,:), xAP(i,:))/norm(cross(zML(i,:), xAP(i,:)));            % Vertical axis
        %
        %                 R(:,:,i) = [xAP(i,:)' yV(i,:)' zML(i,:)'];                                      % LRF Rotation Matix
        %
        %                 CosTheta(i) = dot(yV(i,:),[0 0 1])/(norm(yV(i,:))*norm([0 0 1]));               % Define inclination of the LRF with respect to the vertical (GRF)
        %                 ThetaInDegrees(i) = acosd(CosTheta(i));
        %                 Rm = rotz(ThetaInDegrees(i));
        %                 R_align(:,:,i) = R(:,:,i)*Rm;                                                   % Alignment of the LRF with respect to the vertical direction
        %             end
        %             t = Data.(headerTraj).DYNA0;                                                        % Origin of the LRF
    else % MrkGaps == DYNAX
        warning('This condition has not been properly tested yet')
    end
end

%% TOE and HEEL markers represented in the LRF aligned with the vertical
for i = 1:length(zML)
    HEEL_l(i,:) = (R_align(:,:,i)'*(Data.(headerTraj).(strcat(side,FootMrk{1}))(i,:) - t(i,:))')';   % from GRF to LRF
    TOE_l(i,:) = (R_align(:,:,i)'*(Data.(headerTraj).(strcat(side,FootMrk{2}))(i,:) - t(i,:))')';    % from GRF to LRF
end
MrkL.HEEL = HEEL_l;
MrkL.TOE = TOE_l;
MrkL.COM = t;
MrkL.R = R_align;

TrajH = HEEL_l(:,1);
TrajT = TOE_l(:,1);

%% check Local mrk traj velocity
vTrajH = diff(TrajH)*fs;
vTrajT = diff(TrajT)*fs;

IC = [];
FC = [];

%% IC identification
[~, HS_n] = findpeaks(TrajH,'MinPeakProminence',0.1);
v_Heel = diff(TrajH)/(1/fs);

if ~isempty(HS_n)
    for i = 1:length(HS_n)
        if (v_Heel(HS_n(i))<0) && (v_Heel(HS_n(i)-1)>0)
            IC = [IC; HS_n(i)];
        end
    end
else
    fprintf('IC events not found. Please use another method!\n');
end

%% TO identification
ws1 = 0.08*fs;
[~,FC_n] = findpeaks(-TrajT,'MinPeakDistance',ws1, 'MinPeakProminence',0.2);
v_Toe = diff(TrajT)/(1/fs);

if ~isempty(FC_n)
    for i = 1:length(FC_n)
        if (v_Toe(FC_n(i))>0) && (v_Toe(FC_n(i)-1)<0)
            FC = [FC; FC_n(i)];
        end
    end
else
    fprintf('FC events not found with Zeni algo. Please use another method!\n');
end

%% CHECK - plot
if fig == 1
    figure('Name',strcat(side,'_Zeni'),'NumberTitle','off')
    subplot(211)
    plot(time,TrajH, 'b')
    hold on
    plot(time, TrajT, 'r')
    if ~isempty(IC)
        scatter(IC/fs, TrajH(IC,1),'filled','o','MarkerFaceColor','cyan','MarkerEdgeColor','blue')
        text((IC(1)/fs),TrajH(IC(1,1)),'HS')
    end
    if  ~isempty(FC)
        scatter(FC/fs, TrajT(FC,1), 'filled','o','MarkerFaceColor','green')
        text((FC(1)/fs),TrajT(FC(1,1)),'TO')
    end
    
    minV = min([min(TrajH),min(TrajT)]);
    maxV = max([max(TrajH),max(TrajT)]);
    for i = 1:size(TurnsNotDec,1)
        X2=([TurnsNotDec(i,1)/fs, TurnsNotDec(i,2)/fs]);
        Y2=([minV,minV]);
        Y3=([maxV,maxV]);
        I = patch([X2 fliplr(X2)],[Y2 fliplr(Y3)], 'b', 'EdgeColor','none');
        alpha(0.1);
    end
    
    subplot(212)
    plot(time, Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3),'b')
    hold on
    plot(time, Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3),'r')
    
    if ~isempty(FC)
        scatter((FC/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC,3),'filled','o','MarkerFaceColor','green')
        text((FC(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC(1),3),'TO')
    end
    
    if ~isempty(IC)
        scatter((IC/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC,3),'filled','o',...
            'MarkerFaceColor','cyan','MarkerEdgeColor','blue')
        text((IC(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(1),3),'HS')
    end
    
    minV = min([min(Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3)),...
        min(Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3))]);
    maxV = max([max(Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3)),...
        max(Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3))]);
    
    for i = 1:size(TurnsNotDec,1)
        X2=([TurnsNotDec(i,1)/fs, TurnsNotDec(i,2)/fs]);
        Y2=([minV,minV]);
        Y3=([maxV,maxV]);
        I = patch([X2 fliplr(X2)],[Y2 fliplr(Y3)], 'b', 'EdgeColor','none');
        alpha(0.1);
    end
end

if fig == 1
    figure('Name',strcat(side,'_check mrk v'),'NumberTitle','off')
    plot(vTrajH,'b','LineWidth',2)
    hold on
    plot(vTrajT,'r','LineWidth',2)
end

if ~isempty(IC) && ~isempty(FC)
    checkOutputs = true;
    
    IC(:,2) = ones(size(IC));
    FC(:,2) = ones(size(FC));
    HeelG = Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3);
    
    %% check portion of the signal with spikes in foot velocity --> possible mrk swap!
    if range(vTrajH)>10 % there are spikes due to mrk swap!
        plot(5*((vTrajH<-2)| (vTrajT<-2)),':m','LineWidth',2)
        vIC = nan(size(IC));
        for k = 1:size(IC,1)
            vIC(k) = vTrajH(IC(k));
            scatter(IC(k),vTrajH(IC(k)),'filled','o','MarkerFaceColor','blue','MarkerEdgeColor','k')
        end
        vFC = nan(size(FC));
        for k = 1:size(FC,1)
            vFC(k) = vTrajT(FC(k));
            scatter(FC(k),vTrajT(FC(k)),'filled','o','MarkerFaceColor','red','MarkerEdgeColor','k')
        end
        
        x1 = find(vTrajH<-2);
        u =  find(vTrajT<-2);
        x1 = [x1;u]; x1 = unique(x1);
        PortionToCheck = [];
        if ~isempty(x1)
            s = find(diff(x1)>1);
            if isempty(s)
                PortionToCheck = [PortionToCheck; x1(1),x1(end)];
            else
                for i = 1:size(s,1)
                    if i == 1
                        PortionToCheck = [PortionToCheck; x1(1),x1(s(i))];
                    else
                        PortionToCheck = [PortionToCheck; x1(s(i-1)+1),x1(s(i))];
                    end
                end
            end
        end
        
        % check if this portions are in the parts of the signal where
        % there are gaps --> these would be removed!!!
        checkArea = nan(size(PortionToCheck,1),size(TurnsNotDec,1));
        for k = 1:size(PortionToCheck,1)
            for turns = 1:size(TurnsNotDec,1)
                checkArea(k,turns) = isempty(intersect(PortionToCheck(k,1):PortionToCheck(k,2),TurnsNotDec(turns,1):TurnsNotDec(turns,2)));
            end
        end
        if size(checkArea,2)>1
            checkArea(:,end+1) = zeros(size(checkArea,1),1);
            for k = 1:size(size(checkArea,1))
                if ~isempty(find(checkArea(k,:), 1))
                    checkArea(k,end)=1;
                else
                    checkArea(k,end)=0;
                end
            end
            checkArea = checkArea(:,end);
        end
        
        % area(s) of interest are those found outside the segments with
        % gaps
        pos = find(checkArea);
        for p = 1:length(pos)
            %check GEs in this area(s)
            posNowIC = find((IC(:,1)>(PortionToCheck(p,1)-0.2*fs)) & (IC(:,1)<(PortionToCheck(p,2)+0.2*fs)));
            if ~isempty(posNowIC)
                if size(posNowIC,1)==2
                    if max(abs(vTrajH(IC(posNowIC(1),1):IC(posNowIC(2),1))))<3
                        IC(posNowIC(1),1) = round(mean(IC(posNowIC,1)));
                        IC(posNowIC(2),2) = 0;
                    else
                        IC(posNowIC,2) = 0;
                    end
                else
                    IC(posNowIC,2) = 0;
                end
            else
                %check mrk traj
                posNowICdown = find(IC(:,1)<(PortionToCheck(p,2)+0.2*fs));
                if ~isempty(posNowICdown)
                    NowICdown = IC(posNowICdown(end),1);
                    if posNowICdown(end)>1
                        if posNowICdown(end)+1 <= size(IC,1)
                            sigNow = HeelG(IC(posNowICdown(end)-1,1):IC(posNowICdown(end)+1,1));
                        else
                            sigNow = HeelG(IC(posNowICdown(end)-1,1):end);
                        end
                        %find peaks --> swing
                        [val,locs] = findpeaks(sigNow);
                        [~,pks] = sort(val);
                        if length(pks)>1
                            if locs(pks(end-1)) < locs(pks(end))
                                startAOI = locs(pks(end-1));
                                endAOI = locs(pks(end));
                            else
                                startAOI = locs(pks(end));
                                endAOI = locs(pks(end-1));
                            end
                            [~,locMin] = min(sigNow(startAOI:endAOI));
                            hNow = sigNow(startAOI+locMin);
                            hZeni = HeelG(IC(posNowICdown(end),1));
                            deltaT = (startAOI+locMin+IC(posNowICdown(end)-1,1))-IC(posNowICdown(end),1);
                            if (hZeni-hNow)>0.03 || abs(deltaT)>10% > 3cm
                                IC(posNowICdown(end),2) = 0;
                            end
                        end
                    else
                        sigNow = HeelG(1:IC(posNowICdown(1)+1,1));
                        %find peaks --> swing
                        [val,locs] = findpeaks(sigNow);
                        [~,pks] = sort(val);
                        if locs(pks(end-1)) < locs(pks(end))
                            startAOI = locs(pks(end-1));
                            endAOI = locs(pks(end));
                        else
                            startAOI = locs(pks(end));
                            endAOI = locs(pks(end-1));
                        end
                        [~,locMin] = min(sigNow(startAOI:endAOI));
                        hNow = sigNow(startAOI+locMin);
                        hZeni = HeelG(IC(posNowICdown(1),1));
                        deltaT = startAOI+locMin-IC(posNowICdown(1),1);
                        if (hZeni-hNow)>0.03 || abs(deltaT)>10 % > 3cm
                            IC(posNowICdown(1),2) = 0;
                        end
                    end
                end
                posNowICup = find(IC(:,1)>(PortionToCheck(p,1)-0.2*fs));
                if ~isempty(posNowICup)
                    NowICup = IC(posNowICup(1),1);
                    if posNowICup(1)<length(IC) && posNowICup(1)>1
                        sigNow = HeelG(IC(posNowICup(1)-1,1):IC(posNowICup(1)+1,1));
                        %find peaks --> swing
                        [val,locs] = findpeaks(sigNow);
                        [~,pks] = sort(val);
                        if length(pks) ==2
                            if locs(pks(end-1)) < locs(pks(end))
                                startAOI = locs(pks(end-1));
                                endAOI = locs(pks(end));
                            else
                                startAOI = locs(pks(end));
                                endAOI = locs(pks(end-1));
                            end
                            [~,locMin] = min(sigNow(startAOI:endAOI));
                            hNow = sigNow(startAOI+locMin);
                            hZeni = HeelG(IC(posNowICup(1),1));
                            deltaT =(startAOI+locMin+IC(posNowICup(1)-1,1))-IC(posNowICup(1),1);
                            if (hZeni-hNow)>0.03 || abs(deltaT)>10 % > 3cm
                                IC(posNowICup(1),2) = 0;
                            end
                        end
                    elseif posNowICup(1)==1
                        if size(IC,1)> 1
                            % more than one IC has been found
                            sigNow = HeelG(1:IC(posNowICup(1)+1,1));
                            %find peaks --> swing
                            [val,locs] = findpeaks(sigNow);
                            [~,pks] = sort(val);
                            if length(pks)>1
                                if locs(pks(end-1)) < locs(pks(end))
                                    startAOI = locs(pks(end-1));
                                    endAOI = locs(pks(end));
                                else
                                    startAOI = locs(pks(end));
                                    endAOI = locs(pks(end-1));
                                end
                                [~,locMin] = min(sigNow(startAOI:endAOI));
                                hNow = sigNow(startAOI+locMin);
                                hZeni = HeelG(IC(posNowICup(1),1));
                                deltaT = startAOI+locMin-IC(posNowICup(1),1);
                                if (hZeni-hNow)>0.03 || abs(deltaT)>10 % > 3cm
                                    IC(posNowICup(1),2) = 0;
                                end
                            end
                        end
                    elseif posNowICup(1)==length(IC)
                        sigNow = HeelG(IC(posNowICup(1)-1,1):end);
                        %find peaks --> swing
                        [val,locs] = findpeaks(sigNow);
                        [~,pks] = sort(val);
                        if locs(pks(end-1)) < locs(pks(end))
                            startAOI = locs(pks(end-1));
                            endAOI = locs(pks(end));
                        else
                            startAOI = locs(pks(end));
                            endAOI = locs(pks(end-1));
                        end
                        [~,locMin] = min(sigNow(startAOI:endAOI));
                        hNow = sigNow(startAOI+locMin);
                        hZeni = HeelG(IC(posNowICup(1),1));
                        deltaT =(startAOI+locMin+IC(posNowICup(1)-1,1))-IC(posNowICup(1),1);
                        if (hZeni-hNow)>0.03 || abs(deltaT)>10 % > 3cm
                            IC(posNowICup(1),2) = 0;
                        end
                    else
                        warning('CHECK THIS CONDITION!');
                    end
                end
            end
            posNowFC = find((FC(:,1)>(PortionToCheck(p,1)-0.2*fs)) & (FC(:,1)<(PortionToCheck(p,2)+0.2*fs)));
            if ~isempty(posNowFC)
                FC(posNowFC,2) = 0;
            end
        end
    end % check portion of the signal with spikes in foot velocity --> possible mrk swap!
    %%
    
    ylim([-10 20])
    IC = IC(IC(:,2)==1,1);
    FC = FC(FC(:,2)==1,1);
    
    %% CHECK HS/TO
    DistInTO = diff(FC);
    PmHS = find(DistInTO> mean(DistInTO));   % check for possible missing HS
    xTO = FC(PmHS);                          % If there is a TO, there should be an HS afterwards
    HSd = [];
    for j = 1:length(xTO)
        checkHS = find(IC>xTO(j));
        if ~isempty(checkHS)
            if (IC(checkHS(1))-xTO(j))>150
                [~, addHS] = findpeaks(TrajH(xTO(j):xTO(j)+100,1),'MinPeakProminence',0.05);
                if isempty(addHS)
                else
                    HSd = [HSd; addHS+xTO(j)];
                    scatter((HSd(end)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(HSd(end),3),'b*')
                end
            end
        end
    end
    IC = sort([IC;HSd]);
    
    
    %% check spikes in the signal
    strideD = diff(IC);
    possibleFalseICs = find(strideD<25);
    IC(:,2)=ones(size(IC));
    if ~isempty(possibleFalseICs)
        for fIC = 1:size(possibleFalseICs,1)
            sigNow = HeelG(IC(possibleFalseICs(fIC)):IC(possibleFalseICs(fIC)+1));
            [~,pos]= min(sigNow);
            IC(possibleFalseICs(fIC),1) = pos+IC(possibleFalseICs(fIC));
            IC(possibleFalseICs(fIC)+1,2)=0;
        end
    end
    
    %% check first IC
    check = true;
    max_break = 3; %3 s
    while check
        if IC(1,1)-max_break*fs <= 0
            SigNow = HeelG(1:IC(1,1));
        else
            SigNow = HeelG(IC(1,1)-max_break*fs:IC(1,1));
        end
        if max(SigNow)< 0.45*max(HeelG(IC(1,1)+1:end)) % && (max(SigNow)-H(1))< 0.1*(max(H(IC(1,1)+1:end))-H(1))
            if (max(SigNow)-HeelG(1))<0.05 && IC(1,2)==0 % (max(SigNow)-mean(SigNow))<0.05
                IC(1,:) = [];
            else
                check = false;
            end
        else
            check = false;
        end
    end
    IC = IC(IC(:,2)==1,1);
    
    %% Refine GE positions & Check if an IC has been perfomed forefoot
    ToeG = Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3);
    ToeH_IC = zeros(size(IC));
    HeelH_IC = zeros(size(IC));
    for k = 1:size(IC,1)
        ToeH_IC(k) = ToeG(IC(k));
        HeelH_IC(k) = HeelG(IC(k));
    end
    minH_G = min(HeelG);
    medianStrideD = median(diff(IC));
    foreFoot = false(size(IC));
    
    h = 15;
    ICh = zeros(size(IC));
    hHeel_now = zeros(size(IC));
    PosHeel_3d = zeros(size(IC,1),3);
    PosToe_3d = zeros(size(IC,1),3);
    AngleNow = zeros(size(IC));
    
    % this would be empty for the personalisation trial
    % if ~isempty(InclinationFootGait)
    %     meanAngleSide = mean(InclinationFootGait.(side));
    %     minAngleSide = min(InclinationFootGait.(side));
    % end
    
    % identify the orientation of the foot with respect to the vertical
    for i = 1:size(IC,1)
        hHeel_now(i) = HeelG(IC(i))-minH_G;
        PosHeel_3d(i,:) = Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(i),:);
        PosToe_3d(i,:) = Data.(headerTraj).(strcat(side,FootMrk{2}))(IC(i),:);
        VectorFootNow = PosToe_3d(i,:)-PosHeel_3d(i,:);
        AngleNow(i) = dot(VectorFootNow,[0 0 1])/(norm(VectorFootNow));
        AngleNow(i) = 90 - rad2deg(acos(AngleNow(i)));
        % correct for the original inclination of the markers on the foot
        AngleNow(i) = AngleNow(i) - InclinationFootStatic.(side);
        if AngleNow(i) < 5
            if AngleNow(i) < 0
                %% POSSIBLE forefoot IC
                % identify the portion of signal to check of the Heel mrk
                if i ==1
                    startAOI = 1;
                else
                    startAOI = IC(i-1);
                end
                if IC(i)+h < size(HeelG,1)
                    endAOI = IC(i)+h;
                else
                    endAOI = size(HeelG,1);
                end
                signNowH = HeelG(startAOI:endAOI);
                
                % ~~
                % identify the mid-swing
                [PKS,LOCS]=findpeaks(signNowH);
                [~,MaxPeak] = max(PKS);
                
                % refined signal after the mid-swing - Heel vert Traj
                endAOI = LOCS(MaxPeak)+startAOI+round(medianStrideD);
                if endAOI < size(HeelG,1)
                    refinedSigH = HeelG(LOCS(MaxPeak)+startAOI:endAOI);
                else
                    refinedSigH = HeelG(LOCS(MaxPeak)+startAOI:end);
                end
                TF_H = ischange(refinedSigH,'linear','MaxNumChanges',1);
                TF_H = find(TF_H); TF_H = TF_H+LOCS(MaxPeak)+startAOI;
                
                % refined signal Toe vert Traj
                if (IC(i)-round(medianStrideD/2))< 1
                    startAOI = 1;
                else
                    startAOI = IC(i)-round((medianStrideD/2));
                end
                if (IC(i)+round(medianStrideD/2)) < size(HeelG,1)
                    endAOI = IC(i)+round((medianStrideD/2));
                else
                    endAOI = size(HeelG,1);
                end
                signNowT = ToeG(startAOI:endAOI);
                
                % identify the max
                [PKS,LOCS]=findpeaks(signNowT);
                [~,MaxPeak] = max(PKS);
                
                % refined signal
                refinedSigT = ToeG(LOCS(MaxPeak)+startAOI:endAOI);
                TF_T = ischange(refinedSigT,'linear','MaxNumChanges',1);
                TF_T = find(TF_T); TF_T = TF_T+LOCS(MaxPeak)+startAOI;
                
                if TF_T<TF_H
                    if range(refinedSigT) > (0.25*range(Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3)))
                        IC(i) = TF_T;
                    end
                    foreFoot(i) = true;
                else
                    %             warning('Check this!');
                    %             IC(i) = TF_T;
                    foreFoot(i) = true;
                end
                
            else % almost flat foot
                foreFoot(i) = true;
            end
            %%
            
        else
            if IC(i)+h > length(HeelG) || IC(i)-h <0
                if IC(i)+h > length(HeelG)
                    sigNow = HeelG(IC(i)-h:end);
                else
                    sigNow = HeelG(1:IC(i)+h);
                end
            else
                sigNow = HeelG(IC(i)-h:IC(i)+h);
            end
            pos = find(islocalmin(sigNow));
            if length(pos) ==1
                IC(i) = IC(i)-h + pos;
            elseif length(pos)>1
                hNow = nan(size(pos));
                for k = 1:length(pos)
                    hNow(k) = sigNow(pos(k));
                end
                [~,Greatmin] = sort(hNow);
                IC(i) = IC(i)-h + pos(Greatmin(1));
            end
            ICh(i) = HeelG(IC(i))-minH_G;
        end
    end
    
    if fig == 1
        h =  findobj('type','figure');
        n = length(h);
        figure(n-1);
        hold on
        subplot(212)
        for k = 1:size(IC,1)
            scatter(IC(k)/fs,HeelG(IC(k)),'filled','o','MarkerFaceColor','blue')
        end
    end
    
    %% ADD A FLAG - IC identified with the Zeni's ALGO
    if ~isempty(IC); IC(:,2) = ones(size(IC)); end
    if ~isempty(FC);FC(:,2) = ones(size(FC)); end
    
    if strcmp(side,'R')
        s = 'right';
    else
        s = 'left';
    end
    
    %% Match events and remove those outside the ones identified by the INDIP system
    [IC, FC, deltaEventHS, deltaEventTO] = matchEventsINDIPandDeltaT(IC, FC, GE_INDIP, fs, s, twindow);
    
%     deltaEventHS = nan(size(IC,1),3);
%     IC_r = round(GE_INDIP.HS.(s)*fs);
%     for e = 1:size(IC,1)
%         deltaEv_now = IC(e,1)-IC_r;
%         [~, posFCConsider] = min(abs(deltaEv_now));
%         if abs(deltaEv_now(posFCConsider))<twindow
%             deltaEventHS(e,1) = deltaEv_now(posFCConsider);
%             deltaEventHS(e,2) = IC(e,1);
%             deltaEventHS(e,3) = IC_r(posFCConsider);
%         end
%     end
%     if length(unique(deltaEventHS(~isnan(deltaEventHS(:,3)),3))) ~= length(deltaEventHS(~isnan(deltaEventHS(:,3)),3))
%         % an event has been associated twice!
%         for t = 1:size(deltaEventHS,1)-1
%             if (deltaEventHS(t+1,3)-deltaEventHS(t,3))== 0
%                 if abs(deltaEventHS(t+1,1))< abs(deltaEventHS(t,1))
%                     deltaEventHS(t,1) = nan;
%                     deltaEventHS(t,3) = nan;
%                 else
%                     deltaEventHS(t+1,1) = nan;
%                     deltaEventHS(t+1,3) = nan;
%                 end
% 
%             else
%             end
%         end
%     end
%     deltaEventHS = deltaEventHS(:,1);
% 
%     deltaEventTO = nan(size(FC,1),3);
%     FC_r = round(GE_INDIP.TO.(s)*fs);
%     for e = 1:size(FC,1)
%         deltaEv_now = FC(e,1)-FC_r;
%         [~, posFCConsider] = min(abs(deltaEv_now));
%         if abs(deltaEv_now(posFCConsider))<twindow
%             deltaEventTO(e,1) = deltaEv_now(posFCConsider);
%             deltaEventTO(e,2) = FC(e,1);
%             deltaEventTO(e,3) = FC_r(posFCConsider);
%         end
%     end
%     if length(unique(deltaEventTO(~isnan(deltaEventTO(:,3)),3))) ~= length(deltaEventTO(~isnan(deltaEventTO(:,3)),3))
%         % an event has been associated twice!
%         for t = 1:size(deltaEventTO,1)-1
%             if (deltaEventTO(t+1,3)-deltaEventTO(t,3))== 0
%                 if abs(deltaEventTO(t+1,1))< abs(deltaEventTO(t,1))
%                     deltaEventTO(t,1) = nan;
%                     deltaEventTO(t,3) = nan;
%                 else
%                     deltaEventTO(t+1,1) = nan;
%                     deltaEventTO(t+1,3) = nan;
%                 end
% 
%             else
%             end
%         end
%     end
%     deltaEventTO = deltaEventTO(:,1);
    
else
    checkOutputs = false;
    foreFoot = [];
    deltaEventHS = [];
    deltaEventTO = [];
end

end