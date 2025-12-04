function [AngleCorr, TurnM, TurnDur] = mergePelvisAndFootAngles(Euler_Angles, Euler_Angles_LF, Euler_Angles_RF, fs, turnThres)
    
    angleDiff = 20;
    angleDiffF = 40;
    n = size(Euler_Angles_LF,1);
    
    a = abs(Euler_Angles-Euler_Angles_RF)<angleDiff;
    b = abs(Euler_Angles-Euler_Angles_LF)<angleDiff;
    c = abs(Euler_Angles_RF-Euler_Angles_LF)<angleDiff;
    
    AngleCorr = nan(n,1);
    for i = 1:n
        if a(i)&& b(i)&& c(i)
            %% agreement among the different signals
            AngleCorr(i) = Euler_Angles(i);
            
        elseif (a(i)&& b(i)) ||(a(i)&& c(i))
            %% agreement between pelvis segment and a foot
            AngleCorr(i) = Euler_Angles(i);
            
        elseif (b(i)&& c(i)) 
            %% agreement between feet
            AngleCorr(i) = Euler_Angles_LF(i);
        end
    end
    
    %% Sections of the signal where there is no agreement
    checkNan = find(isnan(AngleCorr));
    AngleNotAg = [];
    if ~isempty(checkNan)
        s = find(diff(checkNan)>1);
        if ~isempty(s)
            for i = 1:size(s,1)
                if i == 1
                    AngleNotAg = [AngleNotAg; checkNan(1),checkNan(s(i))];
                else
                    AngleNotAg = [AngleNotAg; checkNan(s(i-1)+1),checkNan(s(i))];
                end
            end
            AngleNotAg = [AngleNotAg; checkNan(s(i)+1),checkNan(end)];
        else
            AngleNotAg = [AngleNotAg; checkNan(1),checkNan(end)];
        end
    end
    
    d = 30;
    t = 1;
    if ~isempty(AngleNotAg)
        for i = 1:size(AngleNotAg,1)
            if AngleNotAg(i,1)-d >0
            start = AngleNotAg(i,1)-d;
            else
                start = 1;
                t = 0;
            end
            l = size(start:AngleNotAg(i,2),2);
            l = (l/n)*100;
            if l<10
                AngleNowP = Euler_Angles(start:AngleNotAg(i,2));
                AngleNowRF = Euler_Angles_RF(start:AngleNotAg(i,2));
                AngleNowLF = Euler_Angles_LF(start:AngleNotAg(i,2));
                R1 = corrcoef(AngleNowP,AngleNowRF);
                R2 = corrcoef(AngleNowP,AngleNowLF);
                R3 = corrcoef(AngleNowRF,AngleNowLF);
                allR = [R1(1,2), R2(1,2), R3(1,2)];
                [~, pos] = max(allR);
                
                if pos == 1
                    shapeOfAngle = mean([AngleNowP, AngleNowRF],2);                  
                elseif pos == 2
                    shapeOfAngle = mean([AngleNowP, AngleNowLF],2);
                elseif pos == 3
                    shapeOfAngle = mean([AngleNowRF, AngleNowLF],2);
                end
                
                if t ==1
                    if start == 1
                        delta = AngleCorr(start)-shapeOfAngle(1);
                    else                        
                        delta = AngleCorr(start-1)-shapeOfAngle(1);
                    end
                else
                    delta = 0;
                end
                AngleCorr(start:AngleNotAg(i,2)) = shapeOfAngle+delta;
            end
        end
    end
    
    %% Sections of the signal where there is STILL no agreement
    checkNan = find(isnan(AngleCorr));
    AngleNotAg = [];
    if ~isempty(checkNan)
        s = find(diff(checkNan)>1);
        for i = 1:size(s,1)
            if i == 1
                AngleNotAg = [AngleNotAg; checkNan(1),checkNan(s(i))];
            else
                AngleNotAg = [AngleNotAg; checkNan(s(i-1)+1),checkNan(s(i))];
            end
        end
        if ~isempty(s)
            AngleNotAg = [AngleNotAg; checkNan(s(i)+1),checkNan(end)];
        else
            AngleNotAg = [AngleNotAg; checkNan(1),checkNan(end)];
        end
    end
       
    AngleFoot = nan(n,1);
    for i = 1:n
        if abs(Euler_Angles_LF(i)-Euler_Angles_RF(i))<angleDiffF
            AngleFoot(i,1) = mean([Euler_Angles_LF(i),Euler_Angles_RF(i)]);
        end
    end
        
    %% Sections of the signal where there is no agreement
    checkNan = find(isnan(AngleFoot));
    AngleNotAgFoot = [];
    if ~isempty(checkNan)
        s = find(diff(checkNan)>1);
        if ~isempty(s)
            for i = 1:size(s,1)
                if i == 1
                    AngleNotAgFoot = [AngleNotAgFoot; checkNan(1),checkNan(s(i))];
                else
                    AngleNotAgFoot = [AngleNotAgFoot; checkNan(s(i-1)+1),checkNan(s(i))];
                end
            end
            AngleNotAgFoot = [AngleNotAgFoot; checkNan(s(i)+1),checkNan(end)];
        else
            AngleNotAgFoot = [AngleNotAgFoot; checkNan(1),checkNan(end)];
        end
    end
    
    %% Remove parts that have been already filled
    for i = 1:size(AngleNotAg,1)
        for j = 1:size(AngleNotAgFoot,1)
            if (AngleNotAgFoot(j,1)<=AngleNotAg(i,1)) && (AngleNotAg(i,2)>=AngleNotAgFoot(j,2))
                AngleNotAgFoot(j,:)= [0, 0];
            end
        end
    end
    if ~isempty(AngleNotAgFoot)
        AngleNotAgFoot = [AngleNotAgFoot(AngleNotAgFoot(:,1)>0,1),AngleNotAgFoot(AngleNotAgFoot(:,2)>0,2)];
    end
    for j = 1:size(AngleNotAgFoot,1)
        if AngleNotAgFoot(j,1)-d >0
            start = AngleNotAgFoot(j,1)-d;
        else
            start = 1;
        end
        stop = AngleNotAgFoot(j,2);        
        shapeOfAngle = AngleFoot(start:stop);
        delta = AngleCorr(start)-shapeOfAngle(1);
        AngleCorr(start:stop) = shapeOfAngle+delta;        
    end
    
    checkNan = find(isnan(AngleCorr));
    AngleNotAg = [];
    if ~isempty(checkNan)
        s = find(diff(checkNan)>1);
        if ~isempty(s)
            for i = 1:size(s,1)
                if i == 1
                    AngleNotAg = [AngleNotAg; checkNan(1),checkNan(s(i))];
                else
                    AngleNotAg = [AngleNotAg; checkNan(s(i-1)+1),checkNan(s(i))];
                end
            end
            AngleNotAg = [AngleNotAg; checkNan(s(i)+1),checkNan(end)];
        else
            AngleNotAg = [AngleNotAg; checkNan(1),checkNan(end)];
        end
    end
    
    t = 1;
    if ~isempty(AngleNotAg)
        for i = 1:size(AngleNotAg,1) 
            shapeOfAngle = [];
            if AngleNotAg(i,1)-d >0
                ind1 = AngleNotAg(i,1)-d;
            else
                ind1 = 1;
                t = 0;
            end
            AngleNowP = Euler_Angles(ind1:AngleNotAg(i,2));
            AngleNowRF = Euler_Angles_RF(ind1:AngleNotAg(i,2));
            AngleNowLF = Euler_Angles_LF(ind1:AngleNotAg(i,2));
            R1 = corrcoef(AngleNowP,AngleNowRF);
            R2 = corrcoef(AngleNowP,AngleNowLF);
            R3 = corrcoef(AngleNowRF,AngleNowLF);
            allR = [R1(1,2), R2(1,2), R3(1,2)];
            [~, pos] = max(allR);
            
            shapeOfAngle(:,1) = mean([AngleNowP, AngleNowRF],2);
            shapeOfAngle(:,2) = mean([AngleNowP, AngleNowLF],2);
            shapeOfAngle(:,3) = mean([AngleNowLF, AngleNowRF],2);
            
            if range(diff(shapeOfAngle(:,pos))/fs)<0.02
                check = [];
                m = 50;
                for j = 1:size(shapeOfAngle(:,pos),1)-m
                    sigNow = shapeOfAngle(i:i+m);
                    check = [check; range(sigNow)<150];
                end
                if size(find(check),1)==size(check,1)
                    if t ==1 && ind1 >1
                        delta = AngleCorr(ind1-1)-shapeOfAngle(1,pos);
                    else
                        delta = 0;
                    end
                    AngleCorr(ind1:AngleNotAg(i,2)) = shapeOfAngle(:,pos)+delta;
                else
                    allR(pos) = 0;
                    [~, pos] = max(allR);
                    if t ==1
                        delta = AngleCorr(ind1-1)-shapeOfAngle(1,pos);
                    else
                        delta = 0;
                    end
                    AngleCorr(ind1:AngleNotAg(i,2)) = shapeOfAngle(:,pos)+delta;
                end
            else
                allR(pos) = 0;
                [~, pos] = max(allR);
                if t ==1 && ind1 >2
                    delta = AngleCorr(ind1-1)-shapeOfAngle(1,pos);
                else
                    delta = 0;
                end
                AngleCorr(ind1:AngleNotAg(i,2)) = shapeOfAngle(:,pos)+delta;
            end  
        end
    end
    
    %%   
    fc = 1.5;                           % Cutoff frequency 
    Wn = fc/(fs/2);                     % Normalized cutoff frequency        
    [b,a] = butter(4, Wn, 'low');       % Butterworth filter
    AngleCorr = filtfilt(b,a,AngleCorr); 
    
    AngVel = diff(AngleCorr)/(1/fs);
    time = 1/fs:1/fs:length(AngVel)/fs;
    
    [Turns] = gyro_sgm_v1(AngleCorr, AngVel, abs(AngVel), time, fs, turnThres);
    TurnM = Turns.Magnitude;
    TurnDur = Turns.TurnDur;
end