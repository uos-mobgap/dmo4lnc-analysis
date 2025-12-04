function [Euler_Angles, TurnM, TurnDur]= identifyTurns_v1(R, fs, turnThres)
    
    %% Turns
    Euler_Angles(:,1) = atan2(R(2,1,:),R(1,1,:))*180/pi;
    Euler_Angles(:,1) = Euler_Angles(:,1)- mean(Euler_Angles(1:5,1));     % remove offset

    % Remove singularities
    for i = 1:length(Euler_Angles)-1
        if abs(Euler_Angles(i+1)-Euler_Angles(i))>179
            if (Euler_Angles(i+1)-Euler_Angles(i))<0
                Euler_Angles(i+1:end) = Euler_Angles(i+1:end)+ abs(Euler_Angles(i+1)-Euler_Angles(i));
            else
                Euler_Angles(i+1:end) = Euler_Angles(i+1:end)- abs(Euler_Angles(i+1)-Euler_Angles(i));
            end
        end
    end

    % Mellone thesis && El-Gohary et al., 2014
    fc = 1.5;                           % Cutoff frequency 
    Wn = fc/(fs/2);                     % Normalized cutoff frequency        
    [b,a] = butter(4, Wn, 'low');       % Butterworth filter   
    
    if sum(isnan(Euler_Angles)) == size(Euler_Angles,1)       
        TurnM = [];
        TurnDur = [];
    else
        Euler_Angles = filtfilt(b,a,Euler_Angles);
        AngVel = diff(Euler_Angles)/(1/fs);
        time = 1/fs:1/fs:length(AngVel)/fs;

        [Turns] = gyro_sgm_v1(Euler_Angles, AngVel, abs(AngVel), time, fs, turnThres);
        TurnM = Turns.Magnitude;
        TurnDur = Turns.TurnDur;
    end
end