function [data, Stereo, InclinationFootGait, InclinationFootStatic] = extractDMOsFromData_diffMrkForGaps_v3(data,TM, tm, Test,testNow, strideParam, WB_Param,...
    TurnParam, CWP_Param, Fig_Settings, thrs_ALL, InclinationFootStatic, InclinationFootGait, Pers, versResult, folder, sbjNow, originalZeni, saveFileDiffFolder, TVS)

%% Extract parameters from the relevant structures
max_st = strideParam.max_st;
min_st = strideParam.min_st;
min_sl = strideParam.min_sl;
max_h = strideParam.max_h;

n_min = WB_Param.n_min;
max_break = WB_Param.max_break;
n_min_Strides = WB_Param.n_strides;

turnThres = TurnParam.turnThres;
maxTurn = TurnParam.maxTurn;
averVel = TurnParam.averVel;

max_h_CWP = CWP_Param.max_h_CWP;
maxTurn_CWP = CWP_Param.maxTurn_CWP;
averVel_CWP = CWP_Param.averVel_CWP;

fig = Fig_Settings.fig;
figV = Fig_Settings.figV;
figAngle = Fig_Settings.figAngle;
saveFigure = Fig_Settings.saveFigure;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trials = fieldnames(data.(TM{tm}).(Test{testNow}));
for tl = 1:size(Trials,1)
    standardsAll = fieldnames(data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards);
    pos = find(strcmp(standardsAll,'Stereophoto_raw'));
    if ~isempty(pos)
        close all;
        %% Skin-markers data Headers
        Data = data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards.(standardsAll{pos});
        headerTraj = 'Mrks'; %'mrks'
        markers = fieldnames(Data.(headerTraj));
        fs = Data.Fs; %Data.(headerTraj).fs
        N = length(Data.(headerTraj).(markers{1}));
        time = 1/fs:1/fs:N/fs;
        
        %% Filter data and TimeFLAG
        [Data, MrkGaps, TurnsNotDec, timeFlag, t, FlagCheck] = dataFilteringAndFlags(Data, headerTraj, markers, fs, time, N);
        
        if FlagCheck %|| ~FlagCheck
            
            %% Identify when the subject is moving
            mrkDyn = markers(contains(markers, 'DYN'));
            [Act, PercAct] = actRec_v3(Data, headerTraj, fs, mrkDyn, fig);
            fprintf('Detected activity Percentage:  %.2f%% \n',PercAct);
            
            %% Identify how a static portion
            if Pers == 1 && isempty(InclinationFootStatic)
                posNowStatic = find(Act,1)-20;
                % select a static portion of the signal
                if ~isnan(Data.(headerTraj).LHEEL(posNowStatic,1))&&~isnan(Data.(headerTraj).LTOE(posNowStatic,1))&&...
                        ~isnan(Data.(headerTraj).RHEEL(posNowStatic,1))&&~isnan(Data.(headerTraj).RTOE(posNowStatic,1))
                    sides = {'R','L'};
                    for s = 1:length(sides)
                        PosHeel_3d = eval(['Data.(headerTraj).' sides{s} 'HEEL(posNowStatic,:)']);
                        PosToe_3d = eval(['Data.(headerTraj).' sides{s} 'TOE(posNowStatic,:)']);
                        VectorFootNow = PosToe_3d-PosHeel_3d;
                        AngleNow = dot(VectorFootNow,[0 0 1])/(norm(VectorFootNow));
                        InclinationFootStatic.(sides{s}) = 90 - rad2deg(acos(AngleNow));
                    end
                else
                    warning('check for gaps in the static!')
                end
            end
            
            %% GE detection
            FootMrk = {'HEEL','TOE'};
            % Ghoussayni 3D - GEs
            twindow = 25;
            figN = 0;
            [GEs4.HS.right, GEs4.TO.right] = GE_Ghoussayni_3D_500(Data, headerTraj, time, FootMrk, 'R', fs, twindow, figN);
            [GEs4.HS.left, GEs4.TO.left] = GE_Ghoussayni_3D_500(Data, headerTraj, time, FootMrk, 'L', fs, twindow, figN);
            
            mrk_to_use = 'HEEL';
            RWS = strideDet_eScience_mrk_GEeffect_v1(GEs4,max_st, min_st, min_sl, max_h_CWP, Data, headerTraj, fs,mrk_to_use, []);
            
            % Zeni- Ghoussayni - GEs identification
            if ~originalZeni
                [GEs.HS.right, GEs.TO.right, Rl, rMrkL, foreFootICs.right, checkEvents.right] = GE_Zeni_Ghoussayni_eScience(Data, headerTraj,time, FootMrk, 'R', fs, mrkDyn, fig, TurnsNotDec, MrkGaps, twindow, RWS);
                [GEs.HS.left, GEs.TO.left, ~, lMrkL, foreFootICs.left, checkEvents.left] = GE_Zeni_Ghoussayni_eScience(Data, headerTraj,time, FootMrk, 'L', fs, mrkDyn, fig, TurnsNotDec, MrkGaps, twindow, RWS);
            else
                [GEs.HS.right, GEs.TO.right, Rl, rMrkL, foreFootICs.right, checkEvents.right] = GE_Zeni_original(Data, headerTraj, time, FootMrk, 'R', fs, mrkDyn, fig,TurnsNotDec, MrkGaps, InclinationFootStatic, InclinationFootGait);
                [GEs.HS.left, GEs.TO.left, ~, lMrkL, foreFootICs.left, checkEvents.left] = GE_Zeni_original(Data, headerTraj, time, FootMrk, 'L', fs, mrkDyn, fig, TurnsNotDec, MrkGaps, InclinationFootStatic, InclinationFootGait);
            end
            
            if checkEvents.right && checkEvents.left
                
                %% Check if there are still additional events or missed events
                %                 if ~originalZeni
                %                     [GEs.HS.right, GEs.TO.right, misingIC.right] = GE_checkEvents_v4(Data, headerTraj, time, FootMrk, 'R', fs, GEs, GEs_v, max_break, foreFootICs, fig);
                %                     [GEs.HS.left, GEs.TO.left, misingIC.left] = GE_checkEvents_v4(Data, headerTraj, time, FootMrk, 'L', fs, GEs, GEs_v, max_break, foreFootICs, fig);
                %                 else
                %                     [GEs.HS.right, GEs.TO.right, misingIC.right] = GE_checkEvents_original(Data, headerTraj, time, FootMrk, 'R', fs, GEs, GEs_v, max_break, foreFootICs, fig);
                %                     [GEs.HS.left, GEs.TO.left, misingIC.left] = GE_checkEvents_original(Data, headerTraj, time, FootMrk, 'L', fs, GEs, GEs_v, max_break, foreFootICs, fig);
                %                 end
                
                %% Identify turns where possible
                [Euler_Angles, ~, ~] = identifyTurns_v1(Rl, fs, turnThres);
                [Euler_Angles_LF, ~, ~, otherAngleFoot.L] = identifyTurnsFromFeet(Data, headerTraj, fs, turnThres, 'L');
                [Euler_Angles_RF, ~, ~, otherAngleFoot.R] = identifyTurnsFromFeet(Data, headerTraj, fs, turnThres, 'R');
                [Euler_Angles_Corr, TurnM, TurnDur] = mergePelvisAndFootAngles(Euler_Angles, Euler_Angles_LF, Euler_Angles_RF, fs, turnThres);
                
                if figAngle == 1
                    f = figure('Name','Angles and turns','NumberTitle','off');
                    subplot(121)
                    plot(Euler_Angles_Corr,'m','LineWidth',2)
                    hold on
                    plot(Euler_Angles_LF,'--g')
                    plot(Euler_Angles_RF,'--r')
                    plot(Euler_Angles,'--k')
                    plotEvents(TurnM, TurnDur,max(max([Euler_Angles_Corr,Euler_Angles_LF,Euler_Angles_RF,Euler_Angles])),...
                        min(min([Euler_Angles_Corr,Euler_Angles_LF,Euler_Angles_RF,Euler_Angles])))
                    legend('PelvisCorr','L_foot','R_foot','Pelvis')
                    subplot(122)
                    plot(Data.(headerTraj).DYNA0(:,1),Data.(headerTraj).DYNA0(:,2),'c','LineWidth',2)
                    f.GraphicsSmoothing = 'off';
                end
                
                %% Check for correct stride detection
                [stride_list] = strideDet_eScience_mrk(GEs, max_st, min_st, min_sl, max_h, Data, headerTraj, fs, 'HEEL');
                
                %% Select only those strides identified when activity is detected
                [stride_listA] = strideDetAct_v1(stride_list, Act);
                if ~isempty(stride_listA)
                    if size(fieldnames(stride_listA),1) == 1
                        stride_listA = stride_list;
                    end
                end
                
                %% "CWP" identification
                %% Check for correct stride detection - strides with a change of elevation have to be included
                [stride_list_CWP] = strideDet_eScience_mrk(GEs, max_st, min_st, min_sl, max_h_CWP, Data, headerTraj, fs, 'HEEL');
                
                %% Select only those strides identified when activity is detected
                [stride_listA_CWP] = strideDetAct_v1(stride_list_CWP, Act);
                
                if ~isempty(stride_listA_CWP)
                    if size(fieldnames(stride_listA_CWP),1) == 1
                        stride_listA_CWP = stride_list;
                    end
                end
                
                if ~isempty(stride_list_CWP.R) || ~isempty(stride_list_CWP.L)
                    th = 6; thrs = thrs_ALL(th);
                    if Pers == 1
                        % only CWP is used to identify mean and  STD
                        % inclination of the feet
                        %% CWP identification
                        if~isempty(stride_listA_CWP)
                            [Standards_CWP, ~] = CWPdet_eScience(Data, stride_listA_CWP, n_min, n_min_Strides, max_break, rMrkL, lMrkL, headerTraj, fs, ...
                                TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn_CWP, averVel_CWP, maxTurn, averVel, mrkDyn, otherAngleFoot, thrs);
                            %                             (Data, stride_listA_CWP, n_min, max_break, rMrkL, lMrkL, headerTraj, fs, ...
                            %                                 TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn_CWP, averVel_CWP, maxTurn, averVel, mrkDyn, otherAngleFoot, thrs);
                        else
                            X = sprintf('CWPs not found!');
                            disp(X)
                            Standards_CWP = [];
                        end
                    else
                        %% CWP identification
                        if~isempty(stride_listA_CWP)
                            [Standards_CWP, ~] = CWPdet_eScience(Data, stride_listA_CWP, n_min, n_min_Strides, max_break, rMrkL, lMrkL, headerTraj, fs, ...
                                TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn_CWP, averVel_CWP, maxTurn, averVel, mrkDyn, otherAngleFoot, thrs);
                            %                             (Data, stride_listA_CWP, n_min, max_break, rMrkL, lMrkL, headerTraj, fs, ...
                            %                                 TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn_CWP, averVel_CWP, maxTurn, averVel, mrkDyn, otherAngleFoot, thrs);
                        else
                            X = sprintf('CWPs not found!');
                            disp(X)
                            Standards_CWP = [];
                        end
                        
                        %% WB identification
                        if~isempty(stride_listA)
                            [Standards_MicroWB, ~] = WBdet_eScience(Data, stride_listA, n_min, n_min_Strides, max_break, rMrkL, lMrkL, headerTraj, fs,...
                                TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn, averVel, mrkDyn, otherAngleFoot, Standards_CWP, thrs);
                            %                             (Data, stride_listA, n_min, max_break, rMrkL, lMrkL, headerTraj, fs,...
                            %                                 TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn, averVel, mrkDyn, otherAngleFoot, Standards_CWP, thrs);
                        else
                            X = sprintf('WBs not found!');
                            disp(X)
                            Standards_MicroWB = [];
                        end
                    end
                else
                    Standards_CWP = [];
                    Standards_MicroWB = [];
                end
                
            else
                % check event flag --> Bad data quality
                Standards_MicroWB = [];
                Standards_CWP = [];
            end
            
        else
            % FlagCheck == False --> Bad data quality
            Standards_MicroWB = [];
            Standards_CWP = [];
        end % if flagCheck
        
        %% Identify how a person walks
        if Pers == 1
            if ~isempty(Standards_CWP)
                ICs = Standards_CWP.InitialContact_Event;
                ICsLR = Standards_CWP.InitialContact_LeftRight;
                ICsNow.R = ICs(strcmp(ICsLR,'Right'));
                ICsNow.L = ICs(strcmp(ICsLR,'Left'));
                % Evaluate how the foot is usally positioned during gait
                side = {'R','L'};
                for s = 1:2
                    for k = 1:length(ICsNow.(side{s}))
                        if ~isnan(ICsNow.(side{s})(k))
                            PosHeel_3d = eval(['Data.(headerTraj).' side{s} 'HEEL(' num2str(ICsNow.(side{s})(k)*fs) ',:)']);
                            PosToe_3d = eval(['Data.(headerTraj).' side{s} 'TOE(' num2str(ICsNow.(side{s})(k)*fs) ',:)']);
                            VectorFootNow = PosToe_3d-PosHeel_3d;
                            AngleNow = dot(VectorFootNow,[0 0 1])/(norm(VectorFootNow));
                            AngleNow = 90 - rad2deg(acos(AngleNow));
                            InclinationFootGait.(side{s})(k) = AngleNow - InclinationFootStatic.(side{s});
                        end
                    end
                end
            end
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        minV = -.03;
        maxV = 1.1*max([max(Data.(headerTraj).LHEEL(:,3)),max(Data.(headerTraj).RHEEL(:,3))]);
        figure('Name','WBs, CWP and Flag','NumberTitle','off');
        %                         f.GraphicsSmoothing = 'off';
        subplot(3,1,[1,2])
        plot(time, Data.(headerTraj).LHEEL(:,3),'g','LineWidth',2)
        hold on
        plot(time, Data.(headerTraj).RHEEL(:,3),'r','LineWidth',2)
        if Pers ~= 1
            WB_flag = zeros(size(timeFlag));
            n_WB = size(Standards_MicroWB,2);
            if ~isempty(n_WB)
                for i = 1:n_WB
                    n_now = size(round(Standards_MicroWB(i).Start*fs):round(Standards_MicroWB(i).End*fs),2);
                    WB_flag(round(Standards_MicroWB(i).Start*fs):round(Standards_MicroWB(i).End*fs))=ones(1,n_now);
                end
            end
            plot(time, WB_flag*0.25,'m','LineWidth',2)
        end
        
        CWP_flag = zeros(size(timeFlag));
        n_CWP = size(Standards_CWP,2);
        if ~isempty(n_CWP)
            for i = 1:n_CWP
                n_now = size(round(Standards_CWP(i).Start*fs):round(Standards_CWP(i).End*fs),2);
                CWP_flag(round(Standards_CWP(i).Start*fs):round(Standards_CWP(i).End*fs))=ones(1,n_now);
            end
        end
        plot(time, CWP_flag*0.22,'c','LineWidth',1.5)
        
        ylim([minV maxV])
        ylabel('Vertical marker trajectories [m]');xlabel('Time [s]')
        
        %% Annotations Field and Standards
        Standards.Stereophoto.ContinuousWalkingPeriod = Standards_CWP;
        if Pers ~= 1
            Standards.Stereophoto.MicroWB = Standards_MicroWB;
            Stereo.MicroWB = Standards_MicroWB;
        end
        Standards.Stereophoto.Fs = fs;
        
        Stereo.ContinuousWalkingPeriod = Standards_CWP;
        Stereo.Fs = fs;
        
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        dev = ischange(timeFlag,'Threshold',0.1);
        changeSig = find(dev);
        if timeFlag(1) == 1
            oddN = 1:2:length(changeSig);
            evenN = 2:2:length(changeSig);
        else
            evenN = 1:2:length(changeSig);
            oddN = 2:2:length(changeSig);
        end
        GS_Info.MissingInfoStart = changeSig(oddN)/fs;
        GS_Info.MissingInfoStop = changeSig(evenN)/fs;
        if timeFlag(end) == 0
            GS_Info.MissingInfoStop(1,end+1) = length(timeFlag)/fs;
        end
        if isempty(GS_Info.MissingInfoStop)
            GS_Info.MissingInfoStop = size(length(timeFlag));
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Annotations.Stereophoto.Flag = timeFlag;
        Annotations.Stereophoto.GapsPercentage = ((N- length(t))/N)*100;
        % Annotations.Stereophoto.GS_Info = GS_Info;
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for i = 1:size(GS_Info.MissingInfoStart,2)
            X2=([GS_Info.MissingInfoStart(i), GS_Info.MissingInfoStop(i)]);
            Y2=([minV,minV]);
            Y3=([maxV,maxV]);
            I = patch([X2 fliplr(X2)],[Y2 fliplr(Y3)], 'b', 'EdgeColor','none');
            alpha(0.1);
        end
        
        if ~isempty(n_CWP)
            for i = 1:n_CWP
                if ~isempty(Standards_CWP(i).Turn_Start)
                    for k = 1:length(Standards_CWP(i).Turn_Start)
                        line([Standards_CWP(i).Turn_Start(k), Standards_CWP(i).Turn_Start(k)],...
                            [minV,maxV],'Color','black','LineStyle','--')
                        line([Standards_CWP(i).Turn_End(k), Standards_CWP(i).Turn_End(k)],...
                            [minV,maxV],'Color','black','LineStyle','--')
                        text(Standards_CWP(i).Turn_Start(k)+0.1,...
                            maxV,num2str(Standards_CWP(i).Turn_Angle(k)))
                    end
                end
                sz = 30;
                for k = 1:size(Standards_CWP(i).InitialContact_Event,2)
                    if ~isnan(Standards_CWP(i).InitialContact_Event(k))
                        if strcmp(Standards_CWP(i).InitialContact_LeftRight{k},'Left')
                            scatter(Standards_CWP(i).InitialContact_Event(k),Data.(headerTraj).LHEEL(round(Standards_CWP(i).InitialContact_Event(k)*fs),3), sz,...
                                'filled','o', 'MarkerFaceColor','blue')
                        else
                            scatter(Standards_CWP(i).InitialContact_Event(k),Data.(headerTraj).RHEEL(round(Standards_CWP(i).InitialContact_Event(k)*fs),3), sz,...
                                'filled','o','MarkerFaceColor','blue')
                        end
                    end
                end
            end
        end
        
        %% Check GEs
        for k = 1:size(Standards_CWP,2)
            if size(Standards_CWP(i).InitialContact_Event,2)~=size(Standards_CWP(i).FinalContact_Event,2)+2
                warning(strcat('Check ICs/FCs for WB ',k))
            end
            if size(Standards_CWP(i).InitialContact_Event,2)~= Standards_CWP(i).NumberStrides+2
                warning(strcat('Check ICs/nStrides',k))
            end
        end
        
        %% add info in the ORIGINAL data structure
        data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards.Stereophoto = Stereo;
        data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards.Annotations = Annotations;
        
        fprintf('Trial:  %s \n',strcat(Test{testNow},'_',Trials{tl}));
        fprintf('n CWP:  %d \n',size(Stereo.ContinuousWalkingPeriod,2));
        if ~isempty(Stereo.ContinuousWalkingPeriod)
            fprintf('CWP n strides:  %d \n',Stereo.ContinuousWalkingPeriod.NumberStrides);
        end
        if Pers ~= 1
            fprintf('\nn WB:  %d \n',size(Stereo.MicroWB,2));
            if ~isempty(Stereo.MicroWB)
                fprintf('WB n strides:  %d \n',Stereo.MicroWB.NumberStrides);
            end
        end
        
        pos2 = find(strcmp(standardsAll,'INDIP'));
        if ~isempty(pos2)
            INDIP = data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards.(standardsAll{pos2});
            fprintf('INDIP n CWP:  %d \n',size(INDIP.ContinuousWalkingPeriod,2));
            CWP_flag = zeros(size(timeFlag));
            if ~isempty(INDIP.ContinuousWalkingPeriod)
                fprintf('INDIP CWP n strides:  %d \n',INDIP.ContinuousWalkingPeriod.NumberStrides);
                n_CWP = size(INDIP.ContinuousWalkingPeriod,2);
                if ~isempty(n_CWP)
                    for i = 1:n_CWP
                        n_now = size(round(INDIP.ContinuousWalkingPeriod(i).Start*fs):round(INDIP.ContinuousWalkingPeriod(i).End*fs),2);
                        CWP_flag(round(INDIP.ContinuousWalkingPeriod(i).Start*fs):round(INDIP.ContinuousWalkingPeriod(i).End*fs))=ones(1,n_now);
                        sz = 30;
                        for k = 1:size(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event,2)
                            if ~isnan(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k))
                                if strcmp(INDIP.ContinuousWalkingPeriod(i).InitialContact_LeftRight{k},'Left')
                                    scatter(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k),Data.(headerTraj).LHEEL(round(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k)*fs),3), sz,...
                                        'o', 'MarkerEdgeColor','magenta','LineWidth',2)
                                else
                                    scatter(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k),Data.(headerTraj).RHEEL(round(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k)*fs),3), sz,...
                                        'o','MarkerEdgeColor','magenta','LineWidth',2)
                                end
                            end
                        end
                    end
                end
            end
            subplot(313)
            hold on
            plot(time, data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.RightFoot.Gyr(:,2),'r','LineWidth',2)
            plot(time, data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.LeftFoot.Gyr(:,2),'g','LineWidth',2)
            plot(time, WB_flag*250,'m','LineWidth',2)
            plot(time, CWP_flag*400,'b','LineWidth',1.5)
            fprintf('\nn INDIP WB:  %d \n',size(INDIP.MicroWB,2));
            if ~isempty(INDIP.MicroWB)
                fprintf('WB INDIP n strides:  %d \n',INDIP.MicroWB.NumberStrides);
            end
            sz = 30;
            n_CWP = size(Stereo.ContinuousWalkingPeriod,2);
            for i = 1:n_CWP
                for k = 1:size(Standards_CWP(i).InitialContact_Event,2)
                    if ~isnan(Standards_CWP(i).InitialContact_Event(k))
                        if strcmp(Standards_CWP(i).InitialContact_LeftRight{k},'Left')
                            scatter(Standards_CWP(i).InitialContact_Event(k),...
                                data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.LeftFoot.Gyr(round(Standards_CWP(i).InitialContact_Event(k)*fs),2), sz,...
                                'filled','o', 'MarkerFaceColor','blue')
                        else
                            scatter(Standards_CWP(i).InitialContact_Event(k),...
                                data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.RightFoot.Gyr(round(Standards_CWP(i).InitialContact_Event(k)*fs),2), sz,...
                                'filled','o','MarkerFaceColor','blue')
                        end
                    end
                end
            end
            if ~isempty(INDIP.ContinuousWalkingPeriod)
                n_CWP = size(INDIP.ContinuousWalkingPeriod,2);
                if ~isempty(n_CWP)
                    for i = 1:n_CWP
                        sz = 30;
                        for k = 1:size(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event,2)
                            if ~isnan(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k))
                                if strcmp(INDIP.ContinuousWalkingPeriod(i).InitialContact_LeftRight{k},'Left')
                                    scatter(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k),...
                                        data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.LeftFoot.Gyr(round(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k)*fs),2), sz,...
                                        'o', 'MarkerEdgeColor','magenta','LineWidth',2)
                                else
                                    scatter(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k),...
                                        data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.RightFoot.Gyr(round(INDIP.ContinuousWalkingPeriod(i).InitialContact_Event(k)*fs),2), sz,...
                                        'o','MarkerEdgeColor','magenta','LineWidth',2)
                                end
                            end
                        end
                    end
                end
            end
        else
            fprintf('INDIP n CWP:  %d \n',0);
            subplot(313)
            hold on
            IMU = fieldnames(data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP);
            if ~isempty(find(contains(IMU, 'RightFoot'), 1))
                plot(time, data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.RightFoot.Gyr(:,2),'r','LineWidth',2)
            end
            if ~isempty(find(contains(IMU, 'LeftFoot'), 1))
                plot(time, data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.LeftFoot.Gyr(:,2),'g','LineWidth',2)
            end
            %             plot(time, WB_flag*0.25,'m','LineWidth',2)
            plot(time, CWP_flag*400,'b','LineWidth',1.5)
            sz = 30;
            if ~isempty(Standards_CWP)
                n_CWP = size(Standards_CWP,2);
                for k = 1:size(Standards_CWP(i).InitialContact_Event,2)
                    if ~isnan(Standards_CWP(i).InitialContact_Event(k))
                        if strcmp(Standards_CWP(i).InitialContact_LeftRight{k},'Left')
                            if ~isempty(find(contains(IMU, 'LeftFoot'), 1))
                                scatter(Standards_CWP(i).InitialContact_Event(k),...
                                    data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.LeftFoot.Gyr(round(Standards_CWP(i).InitialContact_Event(k)*fs),2), sz,...
                                    'filled','o', 'MarkerFaceColor','blue')
                            end
                        else
                            if ~isempty(find(contains(IMU, 'RightFoot'), 1))
                                scatter(Standards_CWP(i).InitialContact_Event(k),...
                                    data.(TM{tm}).(Test{testNow}).(Trials{tl}).SU_INDIP.RightFoot.Gyr(round(Standards_CWP(i).InitialContact_Event(k)*fs),2), sz,...
                                    'filled','o','MarkerFaceColor','blue')
                            end
                        end
                    end
                end
            end
        end
        
        if saveFigure == 1 && Pers ~=1
            if TVS == 0
                h =  findobj('type','figure');
                n = length(h);
                if contains(folder,'Lab')
                    FolderSP = strcat('G:\My Drive\MOBILISE-D (UNISS-USFD)\_Validation_\',sbjNow,'\In Lab\Results');
                else
                    FolderSP = strcat('G:\My Drive\MOBILISE-D (UNISS-USFD)\_Validation_\',sbjNow,'\Results');
                end
                if ~saveFileDiffFolder
                    newSubFolder = sprintf('%s/AllResults_%s/FigureResultsSPvsINDIP', FolderSP, num2str(versResult));
                    currentFolderFig = strcat(FolderSP,'\AllResults_',num2str(versResult),'\FigureResultsSPvsINDIP');
                else
                    newSubFolder = sprintf('%s/AllResults_%s/FigureResultsSPvsINDIP', FolderSP, num2str(versResult)+1);
                    currentFolderFig = strcat(FolderSP,'\AllResults_',num2str(versResult)+1,'\FigureResultsSPvsINDIP');
                end
                if ~exist(newSubFolder, 'dir')
                    mkdir(newSubFolder);
                end
                saveas(figure(n),fullfile(currentFolderFig,['Subj', sbjNow, Test{testNow}, Trials{tl}, '.fig']))
            else
                h =  findobj('type','figure');
                n = length(h);
                if ~saveFileDiffFolder
                    newSubFolder = sprintf('%s/AllResults_%s/FigureResultsSPvsINDIP', folder, num2str(versResult));
                    currentFolderFig = strcat(folder,'\AllResults_',num2str(versResult),'\FigureResultsSPvsINDIP');
                else
                    newSubFolder = sprintf('%s/AllResults_%s/FigureResultsSPvsINDIP', folder, num2str(versResult)+1);
                    currentFolderFig = strcat(folder,'\AllResults_',num2str(versResult)+1,'\FigureResultsSPvsINDIP');
                end
                if ~exist(newSubFolder, 'dir')
                    mkdir(newSubFolder);
                end
                saveas(figure(n),fullfile(currentFolderFig,['Subj', sbjNow, Test{testNow}, Trials{tl}, '.fig']))
            end
        end
        %%
        
        if Pers ~=1
            %% evaluate errors
            thr = 5;
            pos2 = find(strcmp(standardsAll,'INDIP'));
            if ~isempty(pos2)
                INDIP = data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards.(standardsAll{pos2});
            else
                INDIP = '';
            end
            [checkErrors, flag] = evaluateSPvsINDIP(Stereo,INDIP,'ContinuousWalkingPeriod', fs, thr);
            fprintf('Average Stride Speed error: %.1f %% \n',checkErrors.AverageStrideSpeed)
            fprintf('Average Step Cadence error: %.1f %% \n',checkErrors.AverageStepCadence)
            
            %% check if there are gaps at the edges of each WBs
            Stereo_raw = data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards.(standardsAll{pos}).(headerTraj);
            CheckOutputs = WB_GAPS_edges_check(Stereo, Stereo_raw, INDIP, 'ContinuousWalkingPeriod',fs);
            
            if flag
                for k = 1:size(CheckOutputs.Quality_Pelvis,2)
                    fprintf('Overall Pelvis data quality: %.1f (%.1f) %% \n',CheckOutputs.Quality_Pelvis(k).Mean,CheckOutputs.Quality_Pelvis(k).STD)
                    fprintf('Overall R foot data quality: %.1f (%.1f) %% \n',CheckOutputs.Quality_Rfoot(k).Mean,CheckOutputs.Quality_Rfoot(k).STD)
                    fprintf('Overall L foot data quality: %.1f (%.1f) %% \n\n',CheckOutputs.Quality_Rfoot(k).Mean,CheckOutputs.Quality_Lfoot(k).STD)
                end
            end
            
            if isstruct(CheckOutputs.HeelMrkFlag) || isempty([checkErrors.Start])
                countWBremoved = 0;
                elementToRem = [];
                for n = 1:size(Stereo.ContinuousWalkingPeriod,2)
                    if isstruct(CheckOutputs.HeelMrkFlag)
                        if ~isempty(CheckOutputs.HeelMrkFlag(n).Start)&&~isempty(CheckOutputs.HeelMrkFlag(n).End)
                            if CheckOutputs.HeelMrkFlag(n).Start || CheckOutputs.HeelMrkFlag(n).End
                                elementToRem = [elementToRem; n];
                                countWBremoved = countWBremoved+1;
                            end
                        end
                    end
                end
                fprintf('Outputs should be removed for %d WB/CWP(s) for ID %s due to gaps on heel mkr traj at the edges of the WB/CWP \n GS outputs obtained from %s\n\n',countWBremoved, sbjNow, CheckOutputs.RefSystem)
                
                %% Calculate GEs and events using a different marker
                %% Zeni et al. 2008 - GEs identification
                FootMrk = {'INDIP','TOE'};
                if ~originalZeni
                    [GEs_i.HS.right, GEs_i.TO.right, ~, ~, ~, checkEvents_i.right] = GE_Zeni_Ghoussayni_eScience(Data, headerTraj,time, FootMrk, 'R', fs, mrkDyn, fig, TurnsNotDec, MrkGaps, twindow, RWS);
                    [GEs_i.HS.left, GEs_i.TO.left, ~, ~, ~, checkEvents_i.left] = GE_Zeni_Ghoussayni_eScience(Data, headerTraj,time, FootMrk, 'L', fs, mrkDyn, fig, TurnsNotDec, MrkGaps, twindow, RWS);
                    
                    %                     [GEs_i.HS.right, GEs_i.TO.right, ~, ~, foreFootICs_i.right, checkEvents_i.right] = GE_Zeni_v4_INDIP(Data, headerTraj, time, FootMrk, 'R', fs, mrkDyn, fig,TurnsNotDec, MrkGaps, InclinationFootStatic, InclinationFootGait);
                    %                     [GEs_i.HS.left, GEs_i.TO.left, ~, ~, foreFootICs_i.left, checkEvents_i.left] = GE_Zeni_v4_INDIP(Data, headerTraj, time, FootMrk, 'L', fs, mrkDyn, fig, TurnsNotDec, MrkGaps, InclinationFootStatic, InclinationFootGait);
                else
                    [GEs_i.HS.right, GEs_i.TO.right, ~, ~, foreFootICs_i.right, checkEvents_i.right] = GE_Zeni_original(Data, headerTraj, time, FootMrk, 'R', fs, mrkDyn, fig,TurnsNotDec, MrkGaps, InclinationFootStatic, InclinationFootGait);
                    [GEs_i.HS.left, GEs_i.TO.left, ~, ~, foreFootICs_i.left, checkEvents_i.left] = GE_Zeni_original(Data, headerTraj, time, FootMrk, 'L', fs, mrkDyn, fig, TurnsNotDec, MrkGaps, InclinationFootStatic, InclinationFootGait);
                end
                
                if checkEvents_i.right || checkEvents_i.left
                    %% Add missing GEs
                    if ~isempty(CheckOutputs.HeelMrkFlag)% there are gaps at the edges
                        if isstruct(CheckOutputs.HeelMrkFlag) 
                            if ~isempty(find([CheckOutputs.HeelMrkFlag.Start], 1))
                                side = {'right','left'};
                                GEs_names = {'HS','TO'};
                                for s = 1:length(side)
                                    for g = 1:size(GEs_names)
                                        if ~isempty(find(GEs_i.(GEs_names{g}).(side{s})(:,1)<GEs.(GEs_names{g}).(side{s})(1,1),1))
                                            GEs.(GEs_names{g}).(side{s}) = [GEs.(GEs_names{g}).(side{s});...
                                                GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)<GEs.(GEs_names{g}).(side{s})(1,1),1),...
                                                zeros(size(GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)<GEs.(GEs_names{g}).(side{s})(1,1),1)))];
                                            [~,pos_n] = sort(GEs.(GEs_names{g}).(side{s})(:,1));
                                            GEs.(GEs_names{g}).(side{s}) = GEs.(GEs_names{g}).(side{s})(pos_n,:);
                                        end
                                    end
                                end
                            end
                            
                            if ~isempty(find([CheckOutputs.HeelMrkFlag.End], 1))% CheckOutputs.HeelMrkFlag.End == 1
                                side = {'right','left'};
                                GEs_names = {'HS','TO'};
                                for s = 1:length(side)
                                    for g = 1:size(GEs_names)
                                        if ~isempty(find(GEs_i.(GEs_names{g}).(side{s})(:,1)>GEs.(GEs_names{g}).(side{s})(end,1),1))
                                            GEs.(GEs_names{g}).(side{s}) = [GEs.(GEs_names{g}).(side{s});...
                                                GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)>GEs.(GEs_names{g}).(side{s})(end,1),1),...
                                                zeros(size(GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)>GEs.(GEs_names{g}).(side{s})(end,1),1)))];
                                            [~,pos_n] = sort(GEs.(GEs_names{g}).(side{s})(:,1));
                                            GEs.(GEs_names{g}).(side{s}) = GEs.(GEs_names{g}).(side{s})(pos_n,:);
                                        end
                                    end
                                end
                            end
                        else
                            side = {'right','left'};
                            GEs_names = {'HS','TO'};
                            for s = 1:length(side)
                                for g = 1:size(GEs_names)
                                    if ~isempty(GEs_i.(GEs_names{g}).(side{s})) && ~isempty(GEs.(GEs_names{g}).(side{s}))
                                        if ~isempty(find(GEs_i.(GEs_names{g}).(side{s})(:,1)<GEs.(GEs_names{g}).(side{s})(1,1),1))
                                            GEs.(GEs_names{g}).(side{s}) = [GEs.(GEs_names{g}).(side{s});...
                                                GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)<GEs.(GEs_names{g}).(side{s})(1,1),1),...
                                                zeros(size(GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)<GEs.(GEs_names{g}).(side{s})(1,1),1)))];
                                            [~,pos_n] = sort(GEs.(GEs_names{g}).(side{s})(:,1));
                                            GEs.(GEs_names{g}).(side{s}) = GEs.(GEs_names{g}).(side{s})(pos_n,:);
                                        end
                                    elseif ~isempty(GEs_i.(GEs_names{g}).(side{s})) && isempty(GEs.(GEs_names{g}).(side{s}))
                                        GEs.(GEs_names{g}).(side{s}) = [GEs_i.(GEs_names{g}).(side{s})(:,1),...
                                            zeros(size(GEs_i.(GEs_names{g}).(side{s})(:,1)))];
                                    end
                                end                                
                            end
                            side = {'right','left'};
                            GEs_names = {'HS','TO'};
                            for s = 1:length(side)
                                for g = 1:size(GEs_names)
                                    if ~isempty(GEs_i.(GEs_names{g}).(side{s})) && ~isempty(GEs.(GEs_names{g}).(side{s}))
                                        if ~isempty(find(GEs_i.(GEs_names{g}).(side{s})(:,1)>GEs.(GEs_names{g}).(side{s})(end,1),1))
                                            GEs.(GEs_names{g}).(side{s}) = [GEs.(GEs_names{g}).(side{s});...
                                                GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)>GEs.(GEs_names{g}).(side{s})(end,1),1),...
                                                zeros(size(GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)>GEs.(GEs_names{g}).(side{s})(end,1),1)))];
                                            [~,pos_n] = sort(GEs.(GEs_names{g}).(side{s})(:,1));
                                            GEs.(GEs_names{g}).(side{s}) = GEs.(GEs_names{g}).(side{s})(pos_n,:);
                                        end
                                    elseif ~isempty(GEs_i.(GEs_names{g}).(side{s})) && isempty(GEs.(GEs_names{g}).(side{s}))
                                         GEs.(GEs_names{g}).(side{s}) = [GEs_i.(GEs_names{g}).(side{s})(:,1),...
                                             zeros(size(GEs_i.(GEs_names{g}).(side{s})(:,1)))];
                                    end
                                end
                            end
                        end
                    else
                        side = {'right','left'};
                        GEs_names = {'HS','TO'};
                        for s = 1:length(side)
                            for g = 1:size(GEs_names)
                                if ~isempty(find(GEs_i.(GEs_names{g}).(side{s})(:,1)<GEs.(GEs_names{g}).(side{s})(1,1),1))
                                    GEs.(GEs_names{g}).(side{s}) = [GEs.(GEs_names{g}).(side{s});...
                                        GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)<GEs.(GEs_names{g}).(side{s})(1,1),1),...
                                        zeros(size(GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)<GEs.(GEs_names{g}).(side{s})(1,1),1)))];
                                    [~,pos_n] = sort(GEs.(GEs_names{g}).(side{s})(:,1));
                                    GEs.(GEs_names{g}).(side{s}) = GEs.(GEs_names{g}).(side{s})(pos_n,:);
                                end
                            end
                            
                        end
                        side = {'right','left'};
                        GEs_names = {'HS','TO'};
                        for s = 1:length(side)
                            for g = 1:size(GEs_names)
                                if ~isempty(find(GEs_i.(GEs_names{g}).(side{s})(:,1)>GEs.(GEs_names{g}).(side{s})(end,1),1))
                                    GEs.(GEs_names{g}).(side{s}) = [GEs.(GEs_names{g}).(side{s});...
                                        GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)>GEs.(GEs_names{g}).(side{s})(end,1),1),...
                                        zeros(size(GEs_i.(GEs_names{g}).(side{s})(GEs_i.(GEs_names{g}).(side{s})(:,1)>GEs.(GEs_names{g}).(side{s})(end,1),1)))];
                                    [~,pos_n] = sort(GEs.(GEs_names{g}).(side{s})(:,1));
                                    GEs.(GEs_names{g}).(side{s}) = GEs.(GEs_names{g}).(side{s})(pos_n,:);
                                end
                            end
                        end
                    end
                    
                    %% Check added strides
                    for s = 1:length(side)
                        for g = 1:size(GEs_names)
                            eventsNow = GEs.(GEs_names{g}).(side{s});
                            deltaStrides = diff(eventsNow(:,1));
                            if ~isempty(find(deltaStrides<2*twindow,1))
                                % at least an extra event has been added
                                posToRem = find(deltaStrides<2*twindow);
                                for k = 1:length(posToRem)
                                    if GEs.(GEs_names{g}).(side{s})(posToRem(k),2)==1
                                       posToRem(k) = posToRem(k)+1;
                                    end
                                end
                                GEs.(GEs_names{g}).(side{s})(posToRem,:)=[];
                            end
                        end
                    end
                    
                    %% Check for correct stride detection
                    [stride_list] = strideDet_eScience_mrk(GEs, max_st, min_st, min_sl, max_h, Data, headerTraj, fs, 'INDIP');
                    
                    %% Select only those strides identified when activity is detected
                    [stride_listA] = strideDetAct_v1(stride_list, Act);
                    if ~isempty(stride_listA)
                        if size(fieldnames(stride_listA),1) == 1
                            stride_listA = stride_list;
                        end
                    end
                    
                    %% "CWP" identification
                    %% Check for correct stride detection - strides with a change of elevation have to be included
                    [stride_list_CWP] = strideDet_eScience_mrk(GEs, max_st, min_st, min_sl, max_h_CWP, Data, headerTraj, fs, 'INDIP');
                    
                    %% Select only those strides identified when activity is detected
                    [stride_listA_CWP] = strideDetAct_v1(stride_list_CWP, Act);
                    
                    if ~isempty(stride_listA_CWP)
                        if size(fieldnames(stride_listA_CWP),1) == 1
                            stride_listA_CWP = stride_list;
                        end
                    end
                    
                    if ~isempty(stride_list_CWP.R) || ~isempty(stride_list_CWP.L)
                        th = 6; thrs = thrs_ALL(th);
                        if Pers == 1
                            % only CWP is used to identify mean and  STD
                            % inclination of the feet
                            %% CWP identification
                            if~isempty(stride_listA_CWP)
                                [Standards_CWP_2, ~] = CWPdet_eScience(Data, stride_listA_CWP, n_min, n_min_Strides, max_break, rMrkL, lMrkL, headerTraj, fs, ...
                                    TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn_CWP, averVel_CWP, maxTurn, averVel, mrkDyn, otherAngleFoot, thrs);
                                %                             (Data, stride_listA_CWP, n_min, max_break, rMrkL, lMrkL, headerTraj, fs, ...
                                %                                     TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn_CWP, averVel_CWP, maxTurn, averVel, mrkDyn, otherAngleFoot, thrs);
                            else
                                X = sprintf('CWPs not found!');
                                disp(X)
                                Standards_CWP_2 = [];
                            end
                        else
                            %% CWP identification
                            if~isempty(stride_listA_CWP)
                                [Standards_CWP_2, ~] = CWPdet_eScience(Data, stride_listA_CWP, n_min, n_min_Strides, max_break, rMrkL, lMrkL, headerTraj, fs, ...
                                    TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn_CWP, averVel_CWP, maxTurn, averVel, mrkDyn, otherAngleFoot, thrs);
                                %                             (Data, stride_listA_CWP, n_min, max_break, rMrkL, lMrkL, headerTraj, fs, ...
                                %                                     TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn_CWP, averVel_CWP, maxTurn, averVel, mrkDyn, otherAngleFoot, thrs);
                            else
                                X = sprintf('CWPs not found!');
                                disp(X)
                                Standards_CWP_2 = [];
                            end
                            
                            %% WB identification
                            if~isempty(stride_listA)
                                [Standards_MicroWB_2, ~] = WBdet_eScience(Data, stride_listA, n_min, n_min_Strides, max_break, rMrkL, lMrkL, headerTraj, fs,...
                                    TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn, averVel, mrkDyn, otherAngleFoot, Standards_CWP_2, thrs);
                                %                             (Data, stride_listA, n_min, max_break, rMrkL, lMrkL, headerTraj, fs,...
                                %                                     TurnM, TurnDur, TurnsNotDec, Euler_Angles_Corr, maxTurn, averVel, mrkDyn, otherAngleFoot, Standards_CWP, thrs);
                            else
                                X = sprintf('WBs not found!');
                                disp(X)
                                Standards_MicroWB_2 = [];
                            end
                            
                            
                        end
                    else
                        Standards_CWP_2 = [];
                        Standards_MicroWB_2 = [];
                    end
                    
                    [Standards, Stereo, data] = plotAndStandard_WB_CWP(Data, headerTraj, data, TM, tm, Test, testNow, Trials, tl, standardsAll, time, fs, Standards_MicroWB_2, Standards_CWP_2, Pers, timeFlag, N, t);
                    
                    [checkErrors, ~] = evaluateSPvsINDIP(Stereo,INDIP,'ContinuousWalkingPeriod', fs, thr);
                    fprintf('Average Stride Speed error: %.1f %% \n',checkErrors.AverageStrideSpeed)
                    fprintf('Average Step Cadence error: %.1f %% \n',checkErrors.AverageStepCadence)
                    
                    %% check if the start/end are now matching
                    Stereo_raw = data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards.(standardsAll{pos}).(headerTraj);
                    CheckOutputs = WB_GAPS_edges_check(Stereo, Stereo_raw, INDIP, 'ContinuousWalkingPeriod',fs);
                    
                    if isstruct(CheckOutputs.HeelMrkFlag)
                        countWBremoved = 0;
                        elementToRem = [];
                        for n = 1:size(Stereo.ContinuousWalkingPeriod,2)
                            if ~isempty(CheckOutputs.HeelMrkFlag(n).Start)&&~isempty(CheckOutputs.HeelMrkFlag(n).End)
                                if CheckOutputs.HeelMrkFlag(n).Start || CheckOutputs.HeelMrkFlag(n).End
                                    elementToRem = [elementToRem; n];
                                    countWBremoved = countWBremoved+1;
                                end
                            end
                        end
                        fprintf('Outputs cannot be recovered %s - %s\n',(Test{testNow}),(Trials{tl}))
                        fprintf('Outputs should be removed for %d WB/CWP(s) for ID %s due to gaps on heel mkr traj at the edges of the WB/CWP \n GS outputs obtained from %s\n\n',countWBremoved, sbjNow, CheckOutputs.RefSystem)
                    else
                        fprintf('Outputs have been recovered %s - %s\n\n',(Test{testNow}),(Trials{tl}))
                    end
                end
                
            else
                %% OK! No gaps at the edges
            end
        end % if Stereo Raw data available
    end % Trials
end