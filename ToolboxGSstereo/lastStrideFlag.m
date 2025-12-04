% the last stride has to be removed - final IC is not a
% real IC for another stride
function stride_list = lastStrideFlag(stride_list, otherAngleFoot, fs, max_break, thrs, Data, headerTraj, mrk)
fc = 20;                            % Cutoff frequency
Wn = fc/(fs/2);                     % Normalized cutoff frequency
[b,a] = butter(6, Wn, 'low');       % Butterworth filter

TrMrkV = 0.3;
deltaN = 0.05; % 5 cm

if ~isempty(find(isnan(otherAngleFoot.R(:,1)), 1))
    % LARGE GAPS in the SIGNALs
    posToFill = find(isnan(otherAngleFoot.R(:,1)));
    otherAngleFoot.R(:,1) = sigToFill(otherAngleFoot.R(:,1),posToFill);
end
if ~isempty(find(isnan(otherAngleFoot.L(:,1)), 1))
    % LARGE GAPS in the SIGNALs
    posToFill = find(isnan(otherAngleFoot.L(:,1)));
    otherAngleFoot.L(:,1) = sigToFill(otherAngleFoot.L(:,1),posToFill);
end

if isempty(find(isnan(otherAngleFoot.R(:,1)), 1))
    AngVel_RF = diff(filtfilt(b,a,otherAngleFoot.R(:,1)))/(1/fs);
else
    % large gaps could NOT be filled! 
    AngVel_RF = zeros(size(otherAngleFoot.R,1)-1,1);
end
if isempty(find(isnan(otherAngleFoot.L(:,1)), 1))
    AngVel_LF = diff(filtfilt(b,a,otherAngleFoot.L(:,1)))/(1/fs);
else
    % large gaps could NOT be filled!
    AngVel_LF = zeros(size(otherAngleFoot.L,1)-1,1);
end

Max_AngVel_RF = max(AngVel_RF);
Sw_ph_find_R =(AngVel_RF>(Max_AngVel_RF*thrs))*Max_AngVel_RF;

Max_AngVel_LF = max(AngVel_LF);
Sw_ph_find_L =(AngVel_LF>(Max_AngVel_LF*thrs))*Max_AngVel_LF;

Rmrk = strcat('R',mrk); Lmrk = strcat('L',mrk);

figure
subplot(211)
plot(AngVel_RF,'g')
hold on
plot(Sw_ph_find_R,'r');
if ~isempty(stride_list.R)
    IC_R = [stride_list.R.ICEvents]; IC_R = unique(IC_R);plot(IC_R(~isnan(IC_R)), AngVel_RF(IC_R(~isnan(IC_R))),'*k');
end

subplot(212)
hold on
plot(AngVel_LF,'b')
plot(Sw_ph_find_L,'m');
if ~isempty(stride_list.L)
    IC_L = [stride_list.L.ICEvents]; IC_L = unique(IC_L);    plot(IC_L(~isnan(IC_L)), AngVel_LF(IC_L(~isnan(IC_L))),'*k');
end

%% ADD last stride flag
if ~isempty(stride_list.R)
    rStridesD = [stride_list.R.Duration]; posToKeep = ~isnan([stride_list.R.Duration]); rStridesD = rStridesD(posToKeep);
    rTraj = Data.(headerTraj).(Rmrk)(:,3)-min(Data.(headerTraj).(Rmrk)(:,3)); maxRtraj = max(rTraj);
    sw_dR = [];
    for i = 1:size(stride_list.R,2)
        IC_now = stride_list.R(i).ICEvents(1,2)+1;
        sw_dR = [sw_dR; size(find(Sw_ph_find_R(stride_list.R(i).ICEvents(1,1):stride_list.R(i).ICEvents(1,2))),1)/fs];
        if IC_now + max_break*fs < length(Sw_ph_find_R)
            check = find(Sw_ph_find_R(IC_now:IC_now + max_break*fs));
            trajNow = Data.(headerTraj).(Rmrk)(IC_now:IC_now + max_break*fs,3)-min(Data.(headerTraj).(Rmrk)(:,3));
            trajNowOriginal = Data.(headerTraj).(Rmrk)(IC_now:IC_now + max_break*fs,3);
        else
            check = find(Sw_ph_find_R(IC_now:end));
            trajNow = Data.(headerTraj).(Rmrk)(IC_now:end,3)-min(Data.(headerTraj).(Rmrk)(:,3));
            trajNowOriginal = Data.(headerTraj).(Rmrk)(IC_now:end,3);
        end
        if isempty(check)
            if abs(trajNowOriginal(1)-max(trajNowOriginal))<0.03
                stride_list.R(i).LastStrideFlag = true(1);
            else
                if  i == size(stride_list.R,2) && abs(trajNowOriginal(1)-max(trajNowOriginal))<0.04 % this is the last stride that has been identified
                    stride_list.R(i).LastStrideFlag = true(1);
                else
                    stride_list.R(i).LastStrideFlag = false(1);
                end
            end
        elseif length(check)<10
            check2 = find(check>10, 1);
            if isempty(check2)
                if abs(trajNowOriginal(1)-max(trajNowOriginal))<0.03
                    stride_list.R(i).LastStrideFlag = true(1);
                else
                    stride_list.R(i).LastStrideFlag = false(1);
                end
            else
                if abs(trajNowOriginal(1)-max(trajNowOriginal))<0.03
                    stride_list.R(i).LastStrideFlag = true(1);
                else
                    stride_list.R(i).LastStrideFlag = false(1);
                end
            end
        else
            if  i == size(stride_list.R,2) % this is the last stride that has been identified
                sw_Now = length(check)/fs;
                if max(trajNow)>TrMrkV*maxRtraj
                    stride_list.R(i).LastStrideFlag = false(1);
                else
                    stride_list.R(i).LastStrideFlag = true(1);
                end
            else
                stride_list.R(i).LastStrideFlag = false(1);
            end
        end
    end
end

if ~isempty(stride_list.L)
    lStridesD = [stride_list.L.Duration]; posToKeep = ~isnan([stride_list.L.Duration]); lStridesD = lStridesD(posToKeep);
    lTraj = Data.(headerTraj).(Lmrk)(:,3)-min(Data.(headerTraj).(Lmrk)(:,3)); maxLtraj = max(lTraj);
    sw_dL = [];
    for i = 1:size(stride_list.L,2)
        IC_now = stride_list.L(i).ICEvents(1,2)+1;
        sw_dL = [sw_dL; size(find(Sw_ph_find_L(stride_list.L(i).ICEvents(1,1):stride_list.L(i).ICEvents(1,2))),1)/fs];
        if IC_now + max_break*fs < length(Sw_ph_find_L)
            check = find(Sw_ph_find_L(IC_now:IC_now + max_break*fs));
            trajNow = Data.(headerTraj).(Lmrk)(IC_now:IC_now + max_break*fs,3)-min(Data.(headerTraj).(Lmrk)(:,3));
            trajNowOriginal = Data.(headerTraj).(Lmrk)(IC_now:IC_now + max_break*fs,3);
        else
            check = find(Sw_ph_find_L(IC_now:end));
            trajNow = Data.(headerTraj).(Lmrk)(IC_now:end,3)-min(Data.(headerTraj).(Lmrk)(:,3));
            trajNowOriginal = Data.(headerTraj).(Lmrk)(IC_now:end,3);
        end
        if isempty(check)
            if abs(trajNowOriginal(1)-max(trajNowOriginal))<0.03
                stride_list.L(i).LastStrideFlag = true(1);
            else
                if  i == size(stride_list.L,2) && abs(trajNowOriginal(1)-max(trajNowOriginal))<0.04 % this is the last stride that has been identified
                    stride_list.L(i).LastStrideFlag = true(1);
                else
                    stride_list.L(i).LastStrideFlag = false(1);
                end
            end
        elseif length(check)<10
            check2 = find(check>10, 1);
            if isempty(check2)
                if abs(trajNowOriginal(1)-max(trajNowOriginal))<0.03
                    stride_list.L(i).LastStrideFlag = true(1);
                else
                    stride_list.L(i).LastStrideFlag = false(1);
                end
            else
                if abs(trajNowOriginal(1)-max(trajNowOriginal))<0.03
                    stride_list.L(i).LastStrideFlag = true(1);
                else
                    stride_list.L(i).LastStrideFlag = false(1);
                end
            end
        else
            if  i == size(stride_list.L,2) % this is the last stride that has been identified
                sw_Now = length(check)/fs;
                if max(trajNow)>TrMrkV*maxLtraj
                    stride_list.L(i).LastStrideFlag = false(1);
                else
                    stride_list.L(i).LastStrideFlag = true(1);
                end
            else
                stride_list.L(i).LastStrideFlag = false(1);
            end
        end
    end
end
end