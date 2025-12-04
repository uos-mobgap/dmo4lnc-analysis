function [IC,FC] = GE_Ghoussayni_3D_variablethr_RWS_eScience(Data, headerTraj, time, FootMrk, side, fs, RWS, twindow, fig)

%% GE detection based on Ghoussayni et al. (2004)
% Ghoussayni, S., Stevens, C., Durham, S., & Ewins, D. (2004).
% Assessment and validation of a simple automated method for
% the detection of gait events and intervals.
% Gait & Posture, 20(3), 266-272.
% https://doi.org/10.1016/j.gaitpost.2003.10.001

%% Modified as suggested in
% Bruening, D. A., & Ridge, S. T. (2014).
% Automated event detection algorithms in pathological gait.
% Gait & posture, 39(1), 472-477.
% https://doi.org/10.1016/j.gaitpost.2013.08.023
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

StrideL = [RWS.L.Velocity]'; StrideR = [RWS.R.Velocity]';
v_now = nanmean([StrideL; StrideR]);
if strcmp(FootMrk{1},'HEEL')
    thr_IC = .5*v_now;
else
    thr_IC = .8*v_now;
end
thr_FC = .8*v_now;

detectChange = find(Hvelocity>thr_IC);
StartStop = [];
if ~isempty(detectChange)
    s = find(diff(detectChange)>1);
    if ~isempty(s)
        for i = 1:size(s,1)
            if i == 1
                StartStop = [StartStop; detectChange(1),detectChange(s(i))];
            else
                StartStop = [StartStop; detectChange(s(i-1)+1),detectChange(s(i))];
            end
        end
        StartStop = [StartStop; detectChange(s(i)+1),detectChange(end)];
    else
        StartStop = [StartStop; detectChange(1),detectChange(end)];
    end
end

if ~isempty(StartStop)
    duration = StartStop(:,2)-StartStop(:,1);
    StartStop(duration<10,:)=[];
    IC = StartStop(:,2)+1;
else
    IC = [];
end

if ~isempty(IC)
    if isnan(Hvelocity(1,1))
        newIC = find(Hvelocity<thr_IC,1);
        if IC(1,1)-newIC > 2*twindow
            IC = [IC; newIC]; IC = unique(IC);
        end
    end
end

if ~isempty(IC)
    IC(:,2) = ones(size(IC,1),1);
    for i = 1:size(IC,1)
        if IC(i,1)<size(Hvelocity,1)
            if isnan(Hvelocity(IC(i,1)))
                IC(i,2)=0;
            end
        else
            IC(i,2)=0;
        end
    end
    IC = IC(IC(:,2)==1,1);
end
if ~isempty(find((diff(IC))<twindow,1))
    posToCheck = find((diff(IC))<twindow);
    IC(:,2) = ones(size(IC,1),1);
    % check v between these events
    for k = 1:length(posToCheck)
        v1 = Hvelocity(IC(posToCheck(k),1));
        v2 = Hvelocity(IC(posToCheck(k)+1,1));
        minV = nanmin(Hvelocity(IC(posToCheck(k),1):IC(posToCheck(k)+1,1)));
        maxV = nanmax(Hvelocity(IC(posToCheck(k),1):IC(posToCheck(k)+1,1)));
        if maxV>v2 && maxV>v1
            % there is a drop in the velocity curve
            [val,npos] = findpeaks(Hvelocity(IC(posToCheck(k),1):IC(posToCheck(k)+1,1)),'MinPeakHeight',minV);
            [~,pos2]=max(val);
            if max(val) == maxV && npos(pos2)+IC(posToCheck(k),1)< IC(posToCheck(k)+1,1)
                IC(posToCheck(k),2)=0;
            end
        end        
    end
    IC = IC(IC(:,2)==1);
end
detectChange = find(Tvelocity>thr_FC);
StartStop = [];
if ~isempty(detectChange)
    s = find(diff(detectChange)>1);
    if ~isempty(s)
        for i = 1:size(s,1)
            if i == 1
                StartStop = [StartStop; detectChange(1),detectChange(s(i))];
            else
                StartStop = [StartStop; detectChange(s(i-1)+1),detectChange(s(i))];
            end
        end
        StartStop = [StartStop; detectChange(s(i)+1),detectChange(end)];
    else
        StartStop = [StartStop; detectChange(1),detectChange(end)];
    end
end

if ~isempty(StartStop)
    duration = StartStop(:,2)-StartStop(:,1);
    StartStop(duration<10,:)=[];
    FC = StartStop(:,1);
    IC_2 = StartStop(:,2)+1;
else
    FC =  [];
end
if ~isempty(FC)
    FC(:,2) = ones(size(FC,1),1);
    for i = 1:size(FC,1)
        if isnan(Tvelocity(FC(i,1)))
            FC(i,2)=0;
        end
    end
    FC = FC(FC(:,2)==1,1);
end

%% Refine with the heel mrk velocity
if strcmp(FootMrk{1},'HEEL')
    if ~isempty(FC)
        FC(:,2)=ones(size(FC,1),1);
        for i = 1:size(FC,1)
            if FC(i,1)-(twindow/5) > 0 && FC(i,1)+(twindow/5) <length(Hvelocity)
                sigNow = Hvelocity(FC(i,1)-(twindow/5):FC(i,1)+(twindow/5));
                [~,posNow] = findpeaks(sigNow);
                if ~isempty(posNow)
                    if length(posNow)==1
                        FC(i,1) = FC(i,1)+posNow-(twindow/5)+1;      
                    end
                end
            end
            if FC(i,1)-(twindow/5) > 0 && FC(i,1)+(twindow*4) <length(Hvelocity)
                sigNow = Tvelocity(FC(i,1)-(twindow/5):FC(i,1)+(twindow*2));
                mVnow = max(sigNow);
            else
                sigNow = Tvelocity(FC(i,1)-(twindow/5):end);
                mVnow = max(sigNow);
            end
            if Tvelocity(FC(i,1))>.7*mVnow % > 2 m/s
                FC(i,2)=0;
            end
        end
        FC = FC(FC(:,2)==1,1);
    else
        fprintf('FC events not found. Please use another method!\n');
    end    
end

%% refine ICs
if ~isempty(IC)    
    for i = 1:size(IC,1)
        % find matching ICs - TOE marker
        [~,pos] = min(abs(IC(i)-IC_2));
        if abs(IC(i,1)-IC_2(pos))<twindow
            % check if this is a possible IC in forefoot
            if IC_2(pos,1)<IC(i,1) 
                % check if this is a possible IC in forefoot using vertical
                % mrk traj
                if HEEL(IC_2(pos,1),3)>TOE(IC_2(pos,1),3)
                    IC(i,1)=IC_2(pos,1);
                end
            end
        else            
        end
    end
end

if fig == 1
    figure('Name',strcat(side,'_Ghoussayni'),'NumberTitle','off')
    plot(time(1:end-1), Hvelocity,'r')
    hold on
    plot(time(1:end-1), Tvelocity,'g')
    
    if ~isempty(FC)
        scatter((FC/fs),Tvelocity(FC),'v')
        text((FC(1)/fs),Tvelocity(FC(1),1),'FC')
    end
    if ~isempty(IC)
        scatter((IC/fs),Hvelocity(IC,1),'^')
        text((IC(1)/fs),Hvelocity(IC(1),1),'IC')
    end
    hold off
    
    legend(FootMrk{1}, FootMrk{2})
    ylabel('Marker vetical displacement [m]')
    xlabel('Time [s]')
end
end