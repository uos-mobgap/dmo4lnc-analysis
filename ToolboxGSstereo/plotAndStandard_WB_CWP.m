function [Standards, Stereo, data] = plotAndStandard_WB_CWP(Data, headerTraj, data, TM, tm, Test, testNow, Trials, tl, standardsAll, time, fs, Standards_MicroWB, Standards_CWP, Pers, timeFlag, N, t)

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

end