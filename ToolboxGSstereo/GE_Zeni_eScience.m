function [IC, FC, R_align, MrkL] = GE_Zeni_eScience(Data, headerTraj,time, FootMrk, side, fs, PelvicMrk, fig,  TurnsNotDec)

    %% Mrks on the pelvic segment are filled
    for k = 1:length(PelvicMrk)
        Data.(headerTraj).(PelvicMrk{k}) = fillgaps(Data.(headerTraj).(PelvicMrk{k}));
    end

    %% Define a LRF using the skin-markers placed on the pelvis
    % Pelvic Anatomical Reference Frame (ARF), as defined in Cappozzo et al., 1995
    midASIS = (Data.(headerTraj).r_asis+Data.(headerTraj).l_asis)/2; % midpoint between the anterior superior iliac spines
    midPSIS = (Data.(headerTraj).r_psis+Data.(headerTraj).l_psis)/2; % midpoint between the posterior superior iliac spines

    for i = 1:length(Data.(headerTraj).l_asis)
        zML(i,:) = (Data.(headerTraj).r_asis(i,:) - Data.(headerTraj).l_asis(i,:))/...
            norm(Data.(headerTraj).r_asis(i,:) - Data.(headerTraj).l_asis(i,:));          % Medio-lateral axis
        v(i,:) = (midASIS(i,:) - midPSIS(i,:))/...
            norm(midASIS(i,:) - midPSIS(i,:));
        yV(i,:) = cross(zML(i,:), v(i,:))/norm(cross(zML(i,:), v(i,:)));                % Vertical axis
        xAP(i,:) = cross(yV(i,:),zML(i,:))/norm(cross(yV(i,:),zML(i,:)));               % Anterior-Posterior axis

        R(:,:,i) = [xAP(i,:)' yV(i,:)' zML(i,:)'];                                      % LRF Rotation Matix

        CosTheta(i) = dot(yV(i,:),[0 0 1])/(norm(yV(i,:))*norm([0 0 1]));               % Define inclination of the LRF with respect to the vertical (GRF)
        ThetaInDegrees(i) = acosd(CosTheta(i));
        Rm = rotz(ThetaInDegrees(i));
        R_align(:,:,i) = R(:,:,i)*Rm;                                                   % Alignment of the LRF with respect to the vertical direction
    end
    t = midPSIS;                    % Origin of the LRF

    %% TOE and HEEL markers represented in the LRF aligned with the vertical
    for i = 1:length(zML)
        HEEL_l(i,:) = (R_align(:,:,i)'*(Data.(headerTraj).(strcat(side,'_',FootMrk{1}))(i,:) - t(i,:))')';   % from GRF to LRF
        TOE_l(i,:) = (R_align(:,:,i)'*(Data.(headerTraj).(strcat(side,'_',FootMrk{2}))(i,:) - t(i,:))')';    % from GRF to LRF
    end
    MrkL.HEEL = HEEL_l;
    MrkL.TOE = TOE_l;
    MrkL.COM = t;
    MrkL.R = R_align;
    
    sideW = 0;
    
%     if range(HEEL_l(:,1)) > range(HEEL_l(:,3))
        TrajH = HEEL_l(:,1);
        TrajT = TOE_l(:,1);
%     else
% %         TrajH = HEEL_l(:,3);
% %         TrajT = TOE_l(:,3);
%         sideW = 1;
%         fprintf('Sidewalk identified - Zenis algo cannot be used!\n');
%     end
    
    IC = []; 
    FC = [];
    if sideW == 0
        %% IC identification
        ws1 = 0.8*fs;
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

        IC_s = [];
%         if sideW == 1
%             for k = 1:size(IC,1)-1
%                 k1 = find(FC>IC(k),1);
%                 if ~isempty(k1)
%                     SigN = TrajH(IC(k):FC(k1));
%                     [~, ICn] = findpeaks(-SigN);
%                     IC_s = [IC_s; ICn+IC(k)];
%                 end
%             end
%         end

        %% CHECK - plot
        if fig == 1
            figure('Name',strcat(side,'_Zeni'),'NumberTitle','off')
            subplot(211)
            plot(time,TrajH, 'm')
            hold on
            plot(time, TrajT, 'b')
            if ~isempty(IC)
                scatter(IC/fs, TrajH(IC,1),'c^')
%                 scatter(IC_s/fs, TrajH(IC_s,1),'c*')
                text((IC(1)/fs),TrajH(IC(1,1)),'HS')
            end
            if  ~isempty(FC)
                scatter(FC/fs, TrajT(FC,1), 'gv')
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
            plot(time, Data.(headerTraj).(strcat(side,'_',FootMrk{1}))(:,3),'r')
            hold on
            plot(time, Data.(headerTraj).(strcat(side,'_',FootMrk{2}))(:,3),'g')

            if ~isempty(FC)
                scatter((FC/fs),Data.(headerTraj).(strcat(side,'_',FootMrk{2}))(FC,3),'rv')
                text((FC(1)/fs),Data.(headerTraj).(strcat(side,'_',FootMrk{2}))(FC(1),3),'TO')
            end

            if ~isempty(IC)
                scatter((IC/fs),Data.(headerTraj).(strcat(side,'_',FootMrk{1}))(IC,3),'b^')
%                 scatter((IC_s/fs),Data.(headerTraj).(strcat(side,'_',FootMrk{1}))(IC_s,3),'b*')
                text((IC(1)/fs),Data.(headerTraj).(strcat(side,'_',FootMrk{1}))(IC(1),3),'HS')
            end
            
            minV = min([min(Data.(headerTraj).(strcat(side,'_',FootMrk{1}))(:,3)),...
                min(Data.(headerTraj).(strcat(side,'_',FootMrk{2}))(:,3))]);
            maxV = max([max(Data.(headerTraj).(strcat(side,'_',FootMrk{1}))(:,3)),...
                max(Data.(headerTraj).(strcat(side,'_',FootMrk{2}))(:,3))]);
            
            for i = 1:size(TurnsNotDec,1)
                X2=([TurnsNotDec(i,1)/fs, TurnsNotDec(i,2)/fs]);
                Y2=([minV,minV]);
                Y3=([maxV,maxV]);
                I = patch([X2 fliplr(X2)],[Y2 fliplr(Y3)], 'b', 'EdgeColor','none');
                alpha(0.1);
            end
        end

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
                        scatter((HSd(end)/fs),Data.(headerTraj).(strcat(side,'_',FootMrk{1}))(HSd(end),3),'b*')
                    end
                end
            end
        end

        IC = sort([IC;HSd]);
    end
    
    %% ADD A FLAG - IC identified with the Zeni's ALGO
    if ~isempty(IC); IC(:,2) = ones(size(IC)); end
    if ~isempty(FC);FC(:,2) = ones(size(FC)); end
end