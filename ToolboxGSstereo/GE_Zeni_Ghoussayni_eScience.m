function [IC, FC, R_align, MrkL, foreFoot, checkOutputs, deltaEventHS, deltaEventTO] = GE_Zeni_Ghoussayni_eScience(Data, headerTraj,time, FootMrk, side, fs, mrkDyn, fig, TurnsNotDec, MrkGaps, twindow, RWS)

%% GE detection based on Zeni et al. (2008)
% Zeni Jr, J. A., Richards, J. G., & Higginson, J. S. (2008).
% Two simple methods for determining gait events during treadmill and
% overground walking using kinematic data.
% Gait & posture, 27(4), 710-714.
% https://doi.org/10.1016/j.gaitpost.2007.07.007
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

% [GE_gh.HS, GE_gh.TO] = GE_Ghoussayni_500_new(Data, headerTraj, time, FootMrk, side, fs, GE_INDIP, twindow, fig);
[GE_gh.HS, GE_gh.TO] = GE_Ghoussayni_3D_variablethr_RWS_eScience(Data, headerTraj, time, FootMrk, side, fs, RWS, twindow, fig);

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

%% IC identification
v_Heel = diff(TrajH)/(1/fs);

IC = [];
[~, HS_n] = findpeaks(TrajH,'MinPeakProminence',0.1);
if ~isempty(HS_n)
    for i = 1:length(HS_n)
        if (v_Heel(HS_n(i))<0) && (v_Heel(HS_n(i)-1)>0)
            IC = [IC; HS_n(i)];
        end
    end
else
    fprintf('IC events not found. Please use another method!\n');
end
IC_Zeni = IC;

%% Adjust detected ICs
N = size(Data.(headerTraj).(strcat(side,FootMrk{1})),1);
HEEL = Data.(headerTraj).(strcat(side,FootMrk{1}));
TOE = Data.(headerTraj).(strcat(side,FootMrk{2}));
for t = 1:N-1
    % Heel markers velocity
    Hvelocity(t,1) = sqrt((HEEL(t+1,1)- HEEL(t,1))^2+ ...
        (HEEL(t+1,2)- HEEL(t,2))^2 +...
        (HEEL(t+1,3)- HEEL(t,3))^2)/ ...
        (1/fs);
    
    % Toe markers velocity
    Tvelocity(t,1) = sqrt((TOE(t+1,1)- TOE(t,1)).^2+ ...
        (TOE(t+1,2)- TOE(t,2)).^2 +...
        (TOE(t+1,3)- TOE(t,3)).^2)/ ...
        (1/fs);
end
IC_G = GE_gh.HS;
% adjust values IC values based on the Ghoussayni_3D method
if ~isempty(IC)
    for i = 1:length(IC)
        % match the relevant FC
        [val, posToMatch] = min(abs(IC(i)-IC_G(:,1)));
        if IC_G(posToMatch,1)>IC(i) && val<2*twindow
            IC(i) = IC_G(posToMatch,1);
        else
            
        end
    end
else
    
end

%% TO identification
FC = [];
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

%% Refine with the heel mrk velocity
FC_G = GE_gh.TO;
FC_G(:,end+1)=zeros(size(FC_G));
if strcmp(FootMrk{1},'HEEL')
    if ~isempty(FC_G)
        for i = 1:size(FC_G,1)
            if FC_G(i,1)-(twindow/5) > 0 && FC_G(i,1)+(twindow/5) <length(Hvelocity)
                sigNow = Hvelocity(FC_G(i,1)-(twindow/5):FC_G(i,1)+(twindow/5));
                [~,posNow] = findpeaks(sigNow);
                if ~isempty(posNow)
                    if length(posNow)==1
                        FC_G(i,1) = FC_G(i,1)+posNow-(twindow/5)+1;
                        FC_G(i,2)=1;
                    end
                else
                    
                end
            end
        end
    else
        
    end
end

% adjust values
if ~isempty(FC)
    FC(:,2) = nan(size(FC,1),1);
    for i = 1:size(FC,1)
        % match the relevant FC
        [~, posToMatch] = min(abs(FC(i,1)-FC_G(:,1)));
        if abs(FC(i,1)-FC_G(posToMatch,1))<twindow
            if FC_G(posToMatch,2)==1 && FC_G(posToMatch,1)>FC(i,1)
                %                 FC(i,1) = FC_G(posToMatch,1);
                FC(i,2) = FC_G(posToMatch,1);
            else
                
            end
        end
    end
else
end
if ~isempty(FC)
    if size(unique(FC(~isnan(FC(:,2)),2)),1) == size(FC(~isnan(FC(:,2)),2),1)
        % OK
        for k = 1:size(FC,1)
            if ~isnan(FC(k,2))
                FC(k,1) = FC(k,2);
            end
        end
    else
        % Check if there are repetitions
        mm = 1;
    end
    FC = FC(:,1);
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

if isempty(IC)
    if ~isempty(IC_G)
        IC = IC_G;
    end
end
if isempty(FC)
    if ~isempty(FC_G)
        FC = FC_G(:,1);
    end
end

%% Check if an IC has been found which the heel velicity is decreasing
if strcmp(FootMrk{1},'HEEL')
    IC(:,2) = ones(size(IC,1),1);
    for k = 1:size(IC,1)
        if IC(k,1)-2>0 && IC(k,1)+2<length(Hvelocity)
            if Hvelocity(IC(k,1)-2)<Hvelocity(IC(k,1)+2)
                IC(k,2) = 0;
            end
        end
    end
    IC = IC(IC(:,2)==1,1);
end

%% Check if an FC has been found which the TOE velicity is decreasing
if ~isempty(FC)
    FC(:,2) = ones(size(FC,1),1);
    for k = 1:size(FC,1)
        if Tvelocity(FC(k,1)-1)> Tvelocity(FC(k,1)+1)
            FC(k,2) = 0;
        end
    end
    FC = FC(FC(:,2)==1,1);
end

if ~isempty(IC) || ~isempty(FC)
    checkOutputs = true;
    %% check portion of the signal with spikes in foot velocity --> possible mrk swap!
    if range(vTrajH)>10 % there are spikes due to mrk swap!
        figure('Name',strcat(side,'_check mrk v'),'NumberTitle','off')
        plot(vTrajH,'b','LineWidth',2)
        hold on
        plot(vTrajT,'r','LineWidth',2)
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
        if ~isempty(checkArea)
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
        end
        
        % area(s) of interest are those found outside the segments with gaps
        pos_N = find(checkArea);
        h = 20;
        for p = 1:length(pos_N)
            %check GEs in this area(s)
            IC(:,2) = ones(size(IC,1),1);
            for k = 1:size(IC,1)
                if ~isempty(intersect(IC(k,1),PortionToCheck(p,1)-h:PortionToCheck(p,2)+h))
                    GE_new = GE_gh.HS;
                    [~,pos] = min(abs(IC(k,1)-GE_new));
                    if abs(IC(k,1)-GE_new(pos))<twindow+h
                        IC(k,1) = GE_new(pos);
                    else
                        IC(k,2) = 0;
                    end
                end
            end
            IC = IC(IC(:,2)==1);
            
            FC(:,2) = ones(size(FC,1),1);
            for k = 1:size(FC,1)
                if ~isempty(intersect(FC(k),PortionToCheck(p,1)-h:PortionToCheck(p,2)+h))
                    GE_new = GE_gh.TO;
                    [~,pos] = min(abs(FC(k)-GE_new));
                    if abs(FC(k,1)-GE_new(pos))<twindow+h
                        FC(k,1) = GE_new(pos);
                    else
                        FC(k,2) = 0;
                    end
                end
            end
            FC = FC(FC(:,2)==1);
        end
    end % check portion of the signal with spikes in foot velocity --> possible mrk swap!
    %%
    
    %% Remove events found where turns cannot be assessed
    h = 50;
    if ~isempty(TurnsNotDec)
        for tr = 1:size(TurnsNotDec,1)
            IC(:,2) = ones(size(IC,1),1);
            for k = 1:size(IC,1)
                if ~isempty(intersect(IC(k,1),TurnsNotDec(tr,1)-h:TurnsNotDec(tr,2)+h))
                    GE_new = GE_gh.HS;
                    [~,pos] = min(abs(IC(k,1)-GE_new));
                    if abs(IC(k,1)-GE_new(pos))<twindow+h
                        IC(k,1) = GE_new(pos);
                    else
                        IC(k,2) = 0;
                    end
                end
            end
            IC = IC(IC(:,2)==1);
            IC = unique(IC);
            
            FC(:,2) = ones(size(FC,1),1);
            for k = 1:size(FC,1)
                if ~isempty(intersect(FC(k),TurnsNotDec(tr,1)-h:TurnsNotDec(tr,2)+h))
                    GE_new = GE_gh.TO;
                    [~,pos] = min(abs(FC(k)-GE_new));
                    if abs(FC(k,1)-GE_new(pos))<twindow+h
                        FC(k,1) = GE_new(pos);
                    else
                        FC(k,2) = 0;
                    end
                end
            end
            FC = FC(FC(:,2)==1);
            FC = unique(FC);
        end
    end
    
    %% Merge other ICs
    IC_temp = IC;
    IC_temp(:,2) = ones(size(IC));
    for k = 1:size(IC,1)
        posNow = find(GE_gh.HS<IC(k));
        if ~isempty(posNow)
            candidateICs = GE_gh.HS(posNow);
            if k>1
                refPos = find((IC(k)-candidateICs)>twindow & abs(IC(k-1)-candidateICs)>twindow);
                if ~isempty(refPos)
                    IC_temp = [IC_temp; GE_gh.HS(posNow(refPos)) zeros(size( GE_gh.HS(posNow(refPos))))];
                    [~,p] = unique(IC_temp(:,1));
                    IC_temp = IC_temp(p,:);
                end
            else
                refPos = find((IC(k)-candidateICs)>twindow);
                if ~isempty(refPos)
                    IC_temp = [IC_temp; GE_gh.HS(posNow(refPos)) zeros(size( GE_gh.HS(posNow(refPos))))];
                    [~,p] = unique(IC_temp(:,1));
                    IC_temp = IC_temp(p,:);
                end
            end
        end
        if k == size(IC,1)
            posNow = find(GE_gh.HS>IC(k));
            if ~isempty(posNow)
                candidateICs = GE_gh.HS(posNow);
                refPos = find(abs(IC(k)-candidateICs)>twindow);
                if ~isempty(refPos)
                    IC_temp = [IC_temp; GE_gh.HS(posNow(refPos)) zeros(size( GE_gh.HS(posNow(refPos))))];
                    [~,p] = unique(IC_temp(:,1));
                    IC_temp = IC_temp(p,:);
                end
            end
        end
    end
    % check introduced values
    IC_temp(:,3) = ones(size(IC_temp,1),1);
    for k = 1:size(IC_temp,1)-1
        if IC_temp(k,2)==0
            if IC_temp(k+1,1)-IC_temp(k,1)<twindow
                if IC_temp(k+1,2)~=0
                    IC_temp(k,3)=0;
                end
            end
        end
    end
    IC = IC_temp(IC_temp(:,3)==1);
    clear ICtemp;
    
    %% CHECK HS/TO
    % for each FC there should be an IC before
    if ~isempty(IC)
        for k = 1:size(FC,1)
            if k ==1
                IC_now = intersect(IC(:,1),1:FC(k,1));
            else
                IC_now = intersect(IC(:,1),FC(k-1,1):FC(k,1));
            end
            if length(IC_now)>1 % too many ICs found
                for m = 1:length(IC_now)
                    if isempty(find(abs(IC_now(m) - IC)<3*twindow, 1))
                        IC = [IC; IC_now];
                    end
                end
            elseif isempty(IC_now) % FC NOT found
                GE_new = GE_gh.HS;
                if ~isempty(GE_new)
                    if k ==1
                        IC_now = intersect(GE_new(:,1),1:FC(k,1));
                    else
                        IC_now = intersect(GE_new(:,1),FC(k-1,1):FC(k,1));
                    end
                    if length(IC_now)==1
                        if isempty(find(abs(IC_now - IC)<3*twindow, 1))
                            IC = [IC; IC_now];
                        end
                    end
                end
            end
        end
    end
    IC = sort(IC);
    
    deltaICs = diff(IC);
    PmHS = find(deltaICs> mean(deltaICs)); % possible missing ICs
    for k = 1:size(PmHS)
        GE_new = GE_gh.HS;
        IC_now = intersect(GE_new(:,1),IC(PmHS(k),1)+twindow:IC(PmHS(k)+1,1)-twindow);
        if length(IC_now)==1
            if isempty(find(abs(IC_now - IC)<3*twindow, 1))
                IC = [IC; IC_now];
            end
        elseif length(IC_now)>1
            %% too many event have been found!
            for m = 1:length(IC_now)
                if isempty(find(abs(IC_now(m) - IC)<3*twindow, 1))
                    IC = [IC; IC_now];
                end
            end
        end
    end
    IC = sort(IC);
    
    if strcmp(FootMrk{1},'HEEL')
        HeelG = Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3);
    else
        HeelG = Data.(headerTraj).(strcat(side,'HEEL'))(:,3);
    end
    
    ToeG = Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3);
    AUC_vToe = zeros(size(IC,1)-1,1);
    strideDur = diff(IC(:,1));
    for k = 1:size(IC,1)-1
        vToeNow = v_Toe(IC(k,1):IC(k+1,1));
        x = 1/fs:1/fs:size(vToeNow,1)/fs;
        myInt = cumtrapz(x,abs(vToeNow));
        myIntv = @(a,b) max(myInt(x<=b)) - min(myInt(x>=a));
        AUC_vToe(k) = myIntv(x(1), x(end));
    end
    % for each IC there should be a FC afterwards
    IC(:,2)=ones(size(IC,1),1);
    for k = 1:size(IC,1)-1
        FC_now = intersect(FC(:,1),IC(k,1):IC(k+1,1));
        if length(FC_now)>1 % too many FCs found
            for m = 1:length(FC_now)
                if isempty(find(abs(FC_now(m) - FC)<3*twindow, 1))
                    FC = [FC; FC_now];
                end
            end
        elseif isempty(FC_now) % FC NOT found
            GE_new = GE_gh.TO;
            FC_now = intersect(GE_new(:,1),IC(k,1):IC(k+1,1));
            if length(FC_now)==1
                if isempty(find(abs(FC_now - FC)<3*twindow, 1))
                    FC = [FC; FC_now];
                end
            end
            if isempty(FC_now)
                % check that there are no gaps or a few
                if sum(isnan(ToeG(IC(k,1):IC(k+1,1))))<0.5*size(ToeG(IC(k,1):IC(k+1,1)),1)
                    % check heel signal
                    deltaH_heel = max(HeelG(IC(k,1):IC(k+1,1))) - min(HeelG(IC(k,1):IC(k+1,1)));
                    h1 = HeelG(IC(k,1)); h2 = HeelG(IC(k+1,1));
                    if deltaH_heel<0.01
                        %                     if h1>0.7*deltaH_heel && h2>0.7*deltaH_heel % && deltaH_heel<0.02
                        vToeNow = v_Toe(IC(k,1):IC(k+1,1));
                        x = 1/fs:1/fs:size(vToeNow,1)/fs;
                        myInt = cumtrapz(x,abs(vToeNow));
                        myIntv = @(a,b) max(myInt(x<=b)) - min(myInt(x>=a));
                        AUC = myIntv(x(1), x(end));
                        % check if both GEs have been found with the Zeni
                        % algorithm
                        if isempty(find(IC_Zeni == IC(k,1), 1)) || isempty(find(IC_Zeni == IC(k+1,1), 1))
                            f1 = figure('Name',strcat(side,'Check DATA'),'NumberTitle','off');
                            subplot(211)
                            plot(time(1:end-1), Hvelocity,'r')
                            hold on
                            plot(time(1:end-1), Tvelocity,'g')
                            
                            if ~isempty(FC)
                                scatter((FC/fs),Tvelocity(FC),'filled','o','MarkerFaceColor','green')
                                text((FC(1)/fs),Tvelocity(FC(1),1),'FC')
                            end
                            if ~isempty(IC)
                                scatter((IC(:,1)/fs),Hvelocity(IC(:,1),1),'filled','o','MarkerFaceColor','blue')
                                text((IC(:,1)/fs),Hvelocity(IC(:,1),1),'IC')
                                
                                scatter((IC(k,1)/fs),Hvelocity(IC(k,1),1),'filled','o','MarkerFaceColor','red')
                                scatter((IC(k+1,1)/fs),Hvelocity(IC(k+1,1),1),'filled','o','MarkerFaceColor','red')
                            end
                            
                            legend(FootMrk{1}, FootMrk{2})
                            ylabel('Marker velocity [m/s]')
                            xlabel('Time [s]')
                            if IC(k,1)-500 >0 && IC(k+1,1)+500<size(Hvelocity,1)
                                xlim([(IC(k,1)-500)/fs, (IC(k+1,1)+500)/fs])
                            elseif IC(k,1)-500 >0
                                xlim([(IC(k,1)-500)/fs, (size(Hvelocity,1))/fs])
                            else
                                xlim([1/fs, (IC(k+1,1)+500)/fs])
                            end
                            ylim([0 4])
                            
                            subplot(212)
                            hold on
                            if strcmp(FootMrk{1},'HEEL')
                                plot(time, Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3),'r')
                            else
                                plot(time, Data.(headerTraj).(strcat(side,'HEEL'))(:,3),'r')
                            end
                            plot(time, Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3),'g')
                            if ~isempty(FC)
                                scatter((FC/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC,3),'filled','o','MarkerFaceColor','green')
                                text((FC(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC(1),3),'FC')
                            end
                            if ~isempty(IC)
                                if strcmp(FootMrk{1},'HEEL')
                                    scatter((IC(:,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(:,1),3),'filled','o','MarkerFaceColor','blue')
                                    
                                    scatter((IC(k,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(k,1),3),'filled','o','MarkerFaceColor','red')
                                    scatter((IC(k+1,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(k+1,1),3),'filled','o','MarkerFaceColor','red')
                                else
                                    scatter((IC(:,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(:,1),3),'filled','o','MarkerFaceColor','blue')
                                    
                                    scatter((IC(k,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(k,1),3),'filled','o','MarkerFaceColor','red')
                                    scatter((IC(k+1,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(k+1,1),3),'filled','o','MarkerFaceColor','red')
                                end
                            end
                            if IC(k,1)-500 >0 && IC(k+1,1)+500<size(Hvelocity,1)
                                xlim([(IC(k,1)-500)/fs, (IC(k+1,1)+500)/fs])
                            elseif IC(k,1)-500 >0
                                xlim([(IC(k,1)-500)/fs, (size(Hvelocity,1))/fs])
                            else
                                xlim([1/fs, (IC(k+1,1)+500)/fs])
                            end
                            ylim([0 0.35])
                            hold off
                            
                            fprintf('AUC toe v %4.2f m/s, h1: %4.1f m and h2: %4.1f m  \n', AUC, h1, h2)
                            prompt = 'Which ICs have to be REMOVED? F(first)/S(second)/B(both)/N(none): ';
                            str = input(prompt,'s');
                            if strcmpi(str,'F')
                                IC(k,2)=0;
                            elseif strcmpi(str,'S')
                                IC(k+1,2)=0;
                            elseif strcmpi(str,'B')
                                IC(k,2)=0;
                                IC(k+1,2)=0;
                            end
                            
                            close(f1)
                            
                        else % Both found with the Zeni Algo
                            IC(k,2)=0;
                            IC(k+1,2)=0;
                        end
                    else
                        if isempty(find(IC_Zeni == IC(k,1), 1)) || isempty(find(IC_Zeni == IC(k+1,1), 1))
                            % both found with the G. meth
                            vToeNow = v_Toe(IC(k,1):IC(k+1,1));
                            x = 1/fs:1/fs:size(vToeNow,1)/fs;
                            myInt = cumtrapz(x,abs(vToeNow));
                            myIntv = @(a,b) max(myInt(x<=b)) - min(myInt(x>=a));
                            AUC = myIntv(x(1), x(end));
                            
                            f1 = figure('Name',strcat(side,'Check DATA'),'NumberTitle','off');
                            subplot(211)
                            plot(time(1:end-1), Hvelocity,'r')
                            hold on
                            plot(time(1:end-1), Tvelocity,'g')
                            
                            if ~isempty(FC)
                                scatter((FC/fs),Tvelocity(FC),'filled','o','MarkerFaceColor','green')
                                text((FC(1)/fs),Tvelocity(FC(1),1),'FC')
                            end
                            if ~isempty(IC)
                                scatter((IC(:,1)/fs),Hvelocity(IC(:,1),1),'filled','o','MarkerFaceColor','blue')
                                text((IC(:,1)/fs),Hvelocity(IC(:,1),1),'IC')
                                
                                scatter((IC(k,1)/fs),Hvelocity(IC(k,1),1),'filled','o','MarkerFaceColor','red')
                                scatter((IC(k+1,1)/fs),Hvelocity(IC(k+1,1),1),'filled','o','MarkerFaceColor','red')
                            end
                            
                            legend(FootMrk{1}, FootMrk{2})
                            ylabel('Marker velocity [m/s]')
                            xlabel('Time [s]')
                            if IC(k,1)-500 >0 && IC(k+1,1)+500<size(Hvelocity,1)
                                xlim([(IC(k,1)-500)/fs, (IC(k+1,1)+500)/fs])
                            elseif IC(k,1)-500 >0
                                xlim([(IC(k,1)-500)/fs, (size(Hvelocity,1))/fs])
                            else
                                xlim([1/fs, (IC(k+1,1)+500)/fs])
                            end
                            ylim([0 4])
                            
                            subplot(212)
                            hold on
                            if strcmp(FootMrk{1},'HEEL')
                                plot(time, Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3),'r')
                            else
                                plot(time, Data.(headerTraj).(strcat(side,'HEEL'))(:,3),'r')
                            end
                            plot(time, Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3),'g')
                            if ~isempty(FC)
                                scatter((FC/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC,3),'filled','o','MarkerFaceColor','green')
                                text((FC(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC(1),3),'FC')
                            end
                            if ~isempty(IC)
                                if strcmp(FootMrk{1},'HEEL')
                                    scatter((IC(:,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(:,1),3),'filled','o','MarkerFaceColor','blue')
                                    
                                    scatter((IC(k,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(k,1),3),'filled','o','MarkerFaceColor','red')
                                    scatter((IC(k+1,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(k+1,1),3),'filled','o','MarkerFaceColor','red')
                                else
                                    scatter((IC(:,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(:,1),3),'filled','o','MarkerFaceColor','blue')
                                    
                                    scatter((IC(k,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(k,1),3),'filled','o','MarkerFaceColor','red')
                                    scatter((IC(k+1,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(k+1,1),3),'filled','o','MarkerFaceColor','red')
                                end
                            end
                            if IC(k,1)-500 >0 && IC(k+1,1)+500<size(Hvelocity,1)
                                xlim([(IC(k,1)-500)/fs, (IC(k+1,1)+500)/fs])
                            elseif IC(k,1)-500 >0
                                xlim([(IC(k,1)-500)/fs, (size(Hvelocity,1))/fs])
                            else
                                xlim([1/fs, (IC(k+1,1)+500)/fs])
                            end
                            ylim([0 0.35])
                            hold off
                            
                            fprintf('AUC toe v %4.2f m/s, h1: %4.1f m and h2: %4.1f m  \n', AUC, h1, h2)
                            prompt = 'Which ICs have to be REMOVED? F(first)/S(second)/B(both)/N(none): ';
                            str = input(prompt,'s');
                            if strcmpi(str,'F')
                                IC(k,2)=0;
                            elseif strcmpi(str,'S')
                                IC(k+1,2)=0;
                            elseif strcmpi(str,'B')
                                IC(k,2)=0;
                                IC(k+1,2)=0;
                            end
                            
                            close(f1)
                        
                        else % Both found with the Zeni Algo and not refined
                            IC(k,2)=0;
                            IC(k+1,2)=0;
                        end
                        
                    end
                end % check that there are no GAPs
            end % Check for those strides where an FC has not been detected
        end
        
        if strideDur(k)<2*twindow && ~isempty(FC_now)
            % check for spurious ICs
            h1 = HeelG(IC(k,1)); h2 = HeelG(IC(k+1,1));
            vToeNow = v_Toe(IC(k,1):IC(k+1,1));
            x = 1/fs:1/fs:size(vToeNow,1)/fs;
            myInt = cumtrapz(x,abs(vToeNow));
            myIntv = @(a,b) max(myInt(x<=b)) - min(myInt(x>=a));
            AUC = myIntv(x(1), x(end));
            
            % check if both GEs have been found with the Zeni
            % algorithm
            if isempty(find(IC_Zeni == IC(k,1), 1)) || isempty(find(IC_Zeni == IC(k+1,1), 1))
                f1 = figure('Name',strcat(side,'Check DATA'),'NumberTitle','off');
                subplot(211)
                plot(time(1:end-1), Hvelocity,'r')
                hold on
                plot(time(1:end-1), Tvelocity,'g')
                
                if ~isempty(FC)
                    scatter((FC/fs),Tvelocity(FC),'filled','o','MarkerFaceColor','green')
                    text((FC(1)/fs),Tvelocity(FC(1),1),'FC')
                end
                if ~isempty(IC)
                    scatter((IC(:,1)/fs),Hvelocity(IC(:,1),1),'filled','o','MarkerFaceColor','blue')
                    text((IC(:,1)/fs),Hvelocity(IC(:,1),1),'IC')
                    
                    scatter((IC(k,1)/fs),Hvelocity(IC(k,1),1),'filled','o','MarkerFaceColor','red')
                    scatter((IC(k+1,1)/fs),Hvelocity(IC(k+1,1),1),'filled','o','MarkerFaceColor','red')
                end
                if IC(k,1)-500 >0 && IC(k+1,1)+500<size(Hvelocity,1)
                    xlim([(IC(k,1)-500)/fs, (IC(k+1,1)+500)/fs])
                elseif IC(k,1)-500 >0
                    xlim([(IC(k,1)-500)/fs, (size(Hvelocity,1))/fs])
                else
                    xlim([1/fs, (IC(k+1,1)+500)/fs])
                end
                ylim([0 4])
                
                legend(FootMrk{1}, FootMrk{2})
                ylabel('Marker velocity [m/s]')
                xlabel('Time [s]')
                
                subplot(212)
                hold on
                if strcmp(FootMrk{1},'HEEL')
                    plot(time, Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3),'r')
                else
                    plot(time, Data.(headerTraj).(strcat(side,'HEEL'))(:,3),'r')
                end
                plot(time, Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3),'g')
                if ~isempty(FC)
                    scatter((FC/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC,3),'filled','o','MarkerFaceColor','green')
                    text((FC(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC(1),3),'FC')
                end
                if ~isempty(IC)
                    if strcmp(FootMrk{1},'HEEL')
                        scatter((IC(:,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(:,1),3),'filled','o','MarkerFaceColor','blue')
                        
                        scatter((IC(k,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(k,1),3),'filled','o','MarkerFaceColor','red')
                        scatter((IC(k+1,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(k+1,1),3),'filled','o','MarkerFaceColor','red')
                    else
                        scatter((IC(:,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(:,1),3),'filled','o','MarkerFaceColor','blue')
                        
                        scatter((IC(k,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(k,1),3),'filled','o','MarkerFaceColor','red')
                        scatter((IC(k+1,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(k+1,1),3),'filled','o','MarkerFaceColor','red')
                    end
                end
                if IC(k,1)-500 >0 && IC(k+1,1)+500<size(Hvelocity,1)
                    xlim([(IC(k,1)-500)/fs, (IC(k+1,1)+500)/fs])
                elseif IC(k,1)-500 >0
                    xlim([(IC(k,1)-500)/fs, (size(Hvelocity,1))/fs])
                else
                    xlim([1/fs, (IC(k+1,1)+500)/fs])
                end
                ylim([0 0.35])
                hold off
                
                fprintf('AUC toe v %4.2f m/s, h1: %4.1f m and h2: %4.1f m  \n', AUC, h1, h2)
                prompt = 'Which ICs have to be REMOVED? F(first)/S(second)/B(both)/N(none): ';
                str = input(prompt,'s');
                if strcmpi(str,'F')
                    IC(k,2)=0;
                elseif strcmpi(str,'S')
                    IC(k+1,2)=0;
                elseif strcmpi(str,'B')
                    IC(k,2)=0;
                    IC(k+1,2)=0;
                end
                
                close(f1)
            else % Both found with the Zeni Algo
                IC(k,2)=0;
                IC(k+1,2)=0;
            end
            
        end % check for spurious strides
    end
    IC = IC(IC(:,2)==1);
    FC = sort(FC);
    
    %% check spikes in the signal
    IC(:,2)=ones(size(IC,1),1);
    
    %% check first IC
    check = true;
    max_break = 3; %3 s
    if ~isempty(IC)
        while check
            if IC(1,1)-max_break*fs <= 0
                if strcmp(FootMrk{1},'HEEL')
                    SigNow = HeelG(1:IC(1,1));
                else
                    SigNow = Data.(headerTraj).(strcat(side,FootMrk{1}))(1:IC(1,1),3);
                end
            else
                if strcmp(FootMrk{1},'HEEL')
                    SigNow = HeelG(IC(1,1)-max_break*fs:IC(1,1));
                else
                    SigNow = Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(1,1)-max_break*fs:IC(1,1),3);
                end
            end
            if strcmp(FootMrk{1},'HEEL')
                if max(SigNow)< 0.45*max(HeelG(IC(1,1)+1:end)) % && (max(SigNow)-H(1))< 0.1*(max(H(IC(1,1)+1:end))-H(1))
                    if (max(SigNow)-HeelG(IC(1,1)))<0.01 % 0.02 && IC(1,2)==0 % (max(SigNow)-mean(SigNow))<0.05
                        if sum(isnan(SigNow))<0.6*length(SigNow)
                            IC(1,:) = [];
                        elseif sum(~isnan(SigNow))<5
                            IC(1,:) = [];
                        else
                            check = false;
                        end
                    else
                        check = false;
                    end
                elseif sum(isnan(SigNow))>0.95*length(SigNow)
                    IC(1,:) = [];
                elseif length(SigNow)< 10
                    IC(1,:) = [];
                else
                    check = false;
                end                
            else
                if range(SigNow)< 0.45*range(Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(1,1)+1:end,3)) 
                    if (max(abs(SigNow))-Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(1,1),3))<0.01 % 0.02 && IC(1,2)==0 % (max(SigNow)-mean(SigNow))<0.05
                        if sum(isnan(SigNow))<0.6*length(SigNow)
                            IC(1,:) = [];
                        elseif sum(~isnan(SigNow))<5
                            IC(1,:) = [];
                        else
                            check = false;
                        end
                    else
                        check = false;
                    end
                elseif sum(isnan(SigNow))>0.95*length(SigNow)
                    IC(1,:) = [];
                elseif length(SigNow)< 10
                    IC(1,:) = [];
                else
                    check = false;
                end
            end
        end
    end
    IC = IC(IC(:,2)==1,1);
    
    %% check vertical mrk traj
    IC(:,2) = ones(size(IC,1),1);
    h_all = nan(size(IC,1),1);
    for k = 1:size(IC,1)
        h_all(k) = HeelG(IC(k,1));
    end
    h_mean = nanmean(h_all);
    h_std = nanstd(h_all);
    thr = h_mean + 2*h_std;
    thrFixed = min(HeelG)+0.03; % 5 cm
    if ~isempty(find(h_all>thr,1)) || ~isempty(find(h_all>thrFixed,1))% && size(IC,1)>20
        posToCheck = [find(h_all>thr);find(h_all>thrFixed)]; posToCheck = unique(posToCheck);
        for k = 1:size(posToCheck,1)
            if posToCheck(k) ==1
                if IC(posToCheck(k)+1,1)-IC(posToCheck(k),1)>2*twindow
                    % duration is OK
                else
                    IC(posToCheck(k),2)=0;
                end
            elseif posToCheck(k) ==size(IC,1)
                if IC(posToCheck(k),1)-IC(posToCheck(k)-1,1)>2*twindow
                    % duration is OK
                else
                    IC(posToCheck(k),2)=0;
                end
            else
                if IC(posToCheck(k),1)-IC(posToCheck(k)-1,1)>2*twindow &&  IC(posToCheck(k)+1,1)-IC(posToCheck(k),1)>2*twindow
                    % duration is OK
                else
                    IC(posToCheck(k),2)=0;
                end
            end
        end
        
        f1 = figure('Name',strcat(side,'Check DATA - vertical components'),'NumberTitle','off');
        subplot(211)
        plot(time(1:end-1), Hvelocity,'r')
        hold on
        plot(time(1:end-1), Tvelocity,'g')
        
        if ~isempty(FC)
            scatter((FC/fs),Tvelocity(FC),'filled','o','MarkerFaceColor','green')
            text((FC(1)/fs),Tvelocity(FC(1),1),'FC')
        end
        if ~isempty(IC)
            scatter((IC(:,1)/fs),Hvelocity(IC(:,1),1),'filled','o','MarkerFaceColor','blue')
            text((IC(:,1)/fs),Hvelocity(IC(:,1),1),'IC')
            
            scatter((IC(posToCheck,1)/fs),Hvelocity(IC(posToCheck,1),1),'filled','o','MarkerFaceColor','red')
        end
        
        legend(FootMrk{1}, FootMrk{2})
        ylabel('Marker velocity [m/s]')
        xlabel('Time [s]')
        
        subplot(212)
        hold on
        if strcmp(FootMrk{1},'HEEL')
            plot(time, Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3),'r')
        else
            plot(time, Data.(headerTraj).(strcat(side,'HEEL'))(:,3),'r')
            plot(time, Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3),'k')
        end
        plot(time, Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3),'g')
        if ~isempty(FC)
            scatter((FC/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC,3),'filled','o','MarkerFaceColor','green')
            text((FC(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC(1),3),'FC')
        end
        if ~isempty(IC)
            if strcmp(FootMrk{1},'HEEL')
                scatter((IC(:,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(:,1),3),'filled','o','MarkerFaceColor','blue')
                
                for k = 1:size(posToCheck,1)
                    scatter((IC(posToCheck(k),1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(posToCheck(k),1),3),'filled','o','MarkerFaceColor','red')
                    text((IC(posToCheck(k),1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(posToCheck(k),1),3),num2str(k))
                end
            else
                scatter((IC(:,1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(:,1),3),'filled','o','MarkerFaceColor','blue')
                for k = 1:size(posToCheck,1)
                    scatter((IC(posToCheck(k),1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(posToCheck(k),1),3),'filled','o','MarkerFaceColor','red')
%                     text((IC(posToCheck(k),1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(posToCheck(k),1),3),num2str(k))
                end
                
                scatter((IC(:,1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(:,1),3),'filled','o','MarkerFaceColor','blue')
                
                for k = 1:size(posToCheck,1)
                    scatter((IC(posToCheck(k),1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(posToCheck(k),1),3),'filled','o','MarkerFaceColor','red')
                    text((IC(posToCheck(k),1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(posToCheck(k),1),3),num2str(k))
                end
            end
        end
        hold off
        
        fprintf('Positions to check %4.0f \n', size(posToCheck,1))
        prompt = 'Which ICs have to be REMOVED? Insert number - for more events use squared brackets (e.g., [1 2]): ';
        str_n = input(prompt);
        IC(posToCheck(str_n),2)=0;
        
        close(f1)
    end
    
    hMax_all = max(HeelG);
    for k = 1:size(IC,1)
        if IC(k,1)-(2*twindow)>0 && IC(k,1)+(2*twindow)<length(HeelG)
            SigNow = HeelG(IC(k,1)-(2*twindow):IC(k,1)+(2*twindow));
        elseif IC(k,1)-(2*twindow)>0
            SigNow = HeelG(IC(k,1)-(2*twindow):end);
        else
            SigNow = HeelG(1:IC(k,1)+(2*twindow));
        end
        hIC = HeelG(IC(k,1));
        hMaxNow = max(SigNow);
        if abs(hMaxNow-hIC)<0.03*hMaxNow
            if k>1 && k<size(IC,1)
                if IC(k,1)-IC(k-1,1)>2*twindow &&  IC(k+1,1)-IC(k,1)>2*twindow
                    % duration is OK
                else
                    IC(k,2)=0;
                end
            else
                IC(k,2)=0;
            end
        end
    end
    if ~isempty(find(IC(:,2)==0,1))
        IC_OK = IC(IC(:,2)==1,1);
        IC_v = IC(IC(:,2)==0,1);
        figure('Name',strcat(side,'Zeni_Ghoussayni'),'NumberTitle','off')
        subplot(211)
        plot(time(1:end-1), Hvelocity,'r')
        hold on
        plot(time(1:end-1), Tvelocity,'g')
        
        if ~isempty(FC)
            scatter((FC/fs),Tvelocity(FC),'filled','o','MarkerFaceColor','green')
            text((FC(1)/fs),Tvelocity(FC(1),1),'FC')
        end
        if ~isempty(IC_OK)
            scatter((IC_OK/fs),Hvelocity(IC_OK,1),'filled','o','MarkerFaceColor','blue')
            text((IC_OK(1)/fs),Hvelocity(IC_OK(1),1),'IC')
        end
        if ~isempty(IC_v)
            scatter((IC_v/fs),Hvelocity(IC_v,1),'filled','o','MarkerFaceColor','red','MarkerEdgeColor','black')
            text((IC_v(1)/fs),Hvelocity(IC_v(1),1),'IC')
        end
        
        legend(FootMrk{1}, FootMrk{2})
        ylabel('Marker velocity [m/s]')
        xlabel('Time [s]')
        
        subplot(212)
        hold on
        if strcmp(FootMrk{1},'HEEL')
            plot(time, Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3),'r')
        else
            plot(time, Data.(headerTraj).(strcat(side,'HEEL'))(:,3),'r')
        end
        plot(time, Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3),'g')
        if ~isempty(FC)
            scatter((FC/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC,3),'filled','o','MarkerFaceColor','green')
            text((FC(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC(1),3),'FC')
        end
        if ~isempty(IC_OK)
            if strcmp(FootMrk{1},'HEEL')
                scatter((IC_OK/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC_OK,3),'filled','o','MarkerFaceColor','blue')
                text((IC_OK(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC_OK(1),3),'IC')
                
                scatter((IC_v/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC_v,3),'filled','o','MarkerFaceColor','red','MarkerEdgeColor','black')
            else
                scatter((IC_OK/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC_OK,3),'filled','o','MarkerFaceColor','blue')
                text((IC_OK(1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC_OK(1),3),'IC')
                
                scatter((IC_v/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC_v,3),'filled','o','MarkerFaceColor','red','MarkerEdgeColor','black')
            end
        end
        hold off
    end
    IC = IC(IC(:,2)==1,1);
    
    figure('Name',strcat(side,'Zeni_Ghoussayni'),'NumberTitle','off')
    subplot(211)
    plot(time(1:end-1), Hvelocity,'r')
    hold on
    plot(time(1:end-1), Tvelocity,'g')
    
    if ~isempty(FC)
        scatter((FC/fs),Tvelocity(FC),'filled','o','MarkerFaceColor','green')
        text((FC(1)/fs),Tvelocity(FC(1),1),'FC')
    end
    if ~isempty(IC)
        scatter((IC/fs),Hvelocity(IC,1),'filled','o','MarkerFaceColor','blue')
        text((IC(1)/fs),Hvelocity(IC(1),1),'IC')
    end
    
    legend(FootMrk{1}, FootMrk{2})
    ylabel('Marker velocity [m/s]')
    xlabel('Time [s]')
    
    subplot(212)
    hold on
    if strcmp(FootMrk{1},'HEEL')
        plot(time, Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3),'r')
    else
        plot(time, Data.(headerTraj).(strcat(side,'HEEL'))(:,3),'r')
    end
    plot(time, Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3),'g')
    if ~isempty(FC)
        scatter((FC/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC,3),'filled','o','MarkerFaceColor','green')
        text((FC(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{2}))(FC(1),3),'FC')
    end
    if ~isempty(IC)
        if strcmp(FootMrk{1},'HEEL')
            scatter((IC/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC,3),'filled','o','MarkerFaceColor','blue')
            text((IC(1)/fs),Data.(headerTraj).(strcat(side,FootMrk{1}))(IC(1),3),'IC')
        else
            scatter((IC/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC,3),'filled','o','MarkerFaceColor','blue')
            text((IC(1)/fs),Data.(headerTraj).(strcat(side,'HEEL'))(IC(1),3),'IC')
        end
    end
    hold off
    
    %% ADD A FLAG - IC identified with the Zeni's ALGO
    [~,pos] = sort(IC(:,1)); IC = IC(pos,:);
    [~,pos] = sort(FC(:,1)); FC = FC(pos,:);
    if ~isempty(IC); IC(:,2) = ones(size(IC)); end
    if ~isempty(FC);FC(:,2) = ones(size(FC)); end
    
    foreFoot = [];
else
    checkOutputs = false;
    foreFoot = [];
    deltaEventHS = [];
    deltaEventTO = [];
end

end