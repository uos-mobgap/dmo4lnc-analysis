function [IC,FC] = GE_Ghoussayni_3D_500(Data, headerTraj, time, FootMrk, side, fs, twindow, fig)

%% GE detection based on Ghoussayni et al. (2004)
% Ghoussayni, S., Stevens, C., Durham, S., & Ewins, D. (2004).
% Assessment and validation of a simple automated method for
% the detection of gait events and intervals.
% Gait & Posture, 20(3), 266-272.
% https://doi.org/10.1016/j.gaitpost.2003.10.001
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

thr = 0.5; % 500 mm/s  --> 0.5 m/s

detectChange = find(Hvelocity>thr);
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
    IC = StartStop(:,2)+1;
else
    IC = [];
end

thr2 = 1.0; % 500 mm/s  --> 0.5 m/s
detectChange = find(Tvelocity>thr2);
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
    FC = StartStop(:,1);
else
    FC =  [];
end

%% Refine with the heel mrk velocity
if ~isempty(FC)
    for i = 1:length(FC)
        if FC(i)-(twindow/5) > 0 && FC(i)+(twindow/5) <length(Hvelocity)
            sigNow = Hvelocity(FC(i)-(twindow/5):FC(i)+(twindow/5));
            [~,posNow] = findpeaks(sigNow);
            if ~isempty(posNow)
                if length(posNow)==1
                    FC(i) = FC(i)+posNow-(twindow/5)+1;      
                end
            end
        end        
    end
else
    fprintf('FC events not found. Please use another method!\n');
end

if fig ==1
    figure('Name',strcat(side,'_Ghoussayni'),'NumberTitle','off')
    plot(time(1:end-1), Hvelocity,'r')
    hold on
    plot(time(1:end-1), Tvelocity,'g')
    
    if~isempty(FC)
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