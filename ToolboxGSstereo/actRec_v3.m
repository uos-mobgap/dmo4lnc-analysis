function [Act, PercAct] = actRec_v3(Data, headerTraj, fs, PelvicMrk, fig)
    
%% Pre-Processing - ACTIVITY RECOGNITION
    checkNan = zeros(1,4);
    for i = 1:size(PelvicMrk,2)
        checkNan(1,i) = size(find(isnan(Data.(headerTraj).(PelvicMrk{i})(:,1))),1);
    end
    [~,pos] = min(checkNan);
    MrkToUse = PelvicMrk{pos}; 
    COM = Data.(headerTraj).(MrkToUse);
    
    N = length(COM);
    
    % velocity (m/s)
    vP = diff(COM)/(1/fs);
    vRF = diff(Data.(headerTraj).RHEEL)/(1/fs);
    vLF = diff(Data.(headerTraj).LHEEL)/(1/fs);
    
    % acceleration (g)
    aP = (diff(vP)/(1/fs))/9.81;
    aRF = (diff(vRF)/(1/fs))/9.81;
    aLF = (diff(vLF)/(1/fs))/9.81;
    
    [ActP,~, Vp] = getStationaryPeriodsFromMrks(aP, fs,'wlen',1,'actThresh', 0.0005);
    [ActRF,~, VrF] = getStationaryPeriodsFromMrks(aRF, fs,'wlen',1,'actThresh', 0.004);
    [ActLF,~, VlF] = getStationaryPeriodsFromMrks(aLF, fs,'wlen',1,'actThresh', 0.004);
    
    %% Overall Activity
    Act = zeros(N,1);
    for i = 1:length(ActP)
        if ActP(i)==1 && (ActRF(i)==1 || ActLF(i)==1)
            Act(i)=1;
        end
    end
    
    %% Check activity on-off
    onOffAct = [];
    if ~isempty(Act)
        s1 = find(diff(Act)==1);
        s2 = find(diff(Act)==-1);
        if isempty(s1)
            % OK
        else
            for i = 1:size(s1,1)                
                onOffAct = [onOffAct; s1(i),s2(i)];
            end
        end
    end
    onOffAct_temp = onOffAct;
    onOffAct_temp(:,end+1) = ones(size(onOffAct,1),1);
    for i = 1:size(onOffAct,1)
        if onOffAct(i,2)-onOffAct(i,1)<fs
            % activity less than a second
            if i<size(onOffAct,1)
                if onOffAct(i+1,1)-onOffAct(i,2)<fs
                    onOffAct_temp(i,3) = 2;
                else
                    onOffAct_temp(i,3) = 0;
                end
            else
                onOffAct_temp(i,3) = 0;
            end
        end
    end
    if ~isempty(onOffAct_temp)
        if ~isempty(find(onOffAct_temp(:,3)==2,1))
            % portions to merge
            pos = find(onOffAct_temp(:,3)==2);
            startStopNow = [onOffAct_temp(pos,2),onOffAct_temp(pos+1,1)];
            for i = 1:size(startStopNow,1)
                Act(startStopNow(i,1):startStopNow(i,2)) = ones(size(startStopNow(i,1):startStopNow(i,2),2),1);
            end
        end
        if ~isempty(find(onOffAct_temp(:,3)==0,1))
            % portion to remove
            startStopNow = onOffAct_temp(onOffAct_temp(:,3)==0,1:2);
            for i = 1:size(startStopNow,1)
                Act(startStopNow(i,1):startStopNow(i,2)) = zeros(size(startStopNow(i,1):startStopNow(i,2),2),1);
            end
        end
    end
    if fig == 1
        figure('Name','Detected Activity','NumberTitle','off')
        subplot(2,2,[1,3])
        plot(COM(:,3), 'k')
        hold on
        plot(Data.(headerTraj).RHEEL(:,3), 'r')
        plot(Data.(headerTraj).LHEEL(:,3), 'g')
        plot(Act*1.2*max(COM(:,3)), 'b','LineWidth',2) 
        title('Vertical Mrk Traj [m]')
        legend('COM','R foot','L foot','Activity')
        
        subplot(2,2,2)
        plot(aP, 'k')
        hold on
        plot(aRF, 'r')
        plot(aLF, 'g')
        plot(Act, 'b','LineWidth',2) 
        title('Acceleration values [g]')
        
        subplot(2,2,4)
        plot(Vp, 'k')
        hold on
        plot(VrF, 'r')
        plot(VlF, 'g')
        plot(Act, 'b','LineWidth',2)
        title('AMVD [g]')
    end   
    
    PercAct = (length(find(Act))/N)*100;   
end  