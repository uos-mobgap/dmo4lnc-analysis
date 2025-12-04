function [Turns] = gyro_sgm(yaw_angle, wz_filt, wz_filtabs, time, fs, turnThres)

MaxVal = 15; % Find maxima of absolute value
[hpeaks,locpeak]=findpeaks(wz_filtabs,'MinPeakHeight',MaxVal);
n = size(wz_filtabs,1);

figure
subplot(2,1,1)
title('Filtered Angular Velocity (deg/s)');
plot(time, wz_filtabs,'Color',[0.7 0.7 0.7],'DisplayName','wz filtered [deg/s]')
hold on
plot(locpeak/fs,hpeaks,'*r');

CrossVal = 5 ;
zci = @(v) find(v(:).*circshift(v(:),[-1 0]) <= 0);
zx = zci(wz_filt - CrossVal);
Length0 = length(zx);
while length(zx) < length(locpeak)+2
    LengthNow = length(zx);
    CrossVal = CrossVal + 1;
    zx = zci(wz_filt - CrossVal);
    if  LengthNow < Length0
        break;
    end
end

% identify 5deg/s crossings preceding and following each maximum 
x1 = [];
x2 = [];
for i = 1:length(locpeak)
    bef = zx<locpeak(i);
    if isempty(find(bef,1))
        x1 = [x1;1];
    else
        bef = zx(bef);
        x1 = [x1;bef(end)];
    end
    aft = zx>locpeak(i);
    if isempty(find(aft,1))
        x2 = [x2;length(wz_filt)];
    else
        aft = zx(aft);
        x2 = [x2;aft(1)];
    end    
end

allChanges = unique(sort([x1;x2]));

subplot(2,1,2)
plot(time,yaw_angle(1:n))
hold on

% Turn definitions
TurnM = [];
TurnDur = [];
for i = 2:length(allChanges)
    turnNow = yaw_angle(allChanges(i))-yaw_angle(allChanges(i-1));
    turn_d = (allChanges(i)- allChanges(i-1))/fs;
    if (abs(turnNow) > turnThres) && (turn_d > 0.05)
        TurnM = [TurnM; turnNow]; % turn magnitude - in degrees
        TurnDur = [TurnDur; allChanges(i-1), allChanges(i)]; % beginning and end of each turn (frames)
        line([allChanges(i-1)/fs allChanges(i-1)/fs],get(gca,'ylim'),'Color',[0 0 1]);
        line([allChanges(i)/fs allChanges(i)/fs],get(gca,'ylim'),'Color',[0 0 1]);
    elseif ((turn_d < 0.05) && (abs(turnNow) < turnThres)) && i<length(allChanges)
        j = i+1;
        turnNow = yaw_angle(allChanges(j))-yaw_angle(allChanges(i-1));
        turn_d = (allChanges(j)- allChanges(i-1))/fs;
        while (turn_d < 0.05) && (abs(turnNow) < turnThres)
            turnNow = yaw_angle(allChanges(j))-yaw_angle(allChanges(i-1));
            turn_d = (allChanges(j)- allChanges(i-1))/128;
            j = j+1;
        end
        TurnM = [TurnM; turnNow];
        TurnDur = [TurnDur; allChanges(i-1), allChanges(i)];
        line([allChanges(i-1)/fs allChanges(i-1)/fs],get(gca,'ylim'),'Color',[0 0 1]);
        line([allChanges(i)/fs allChanges(i)/fs],get(gca,'ylim'),'Color',[0 0 1]);
    end
end

Turns.Magnitude = TurnM;
Turns.TurnDur = TurnDur;
end