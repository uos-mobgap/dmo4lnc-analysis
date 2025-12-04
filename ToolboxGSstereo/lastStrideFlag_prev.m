% the last stride has to be removed - final IC is not a
% real IC for another stride
function stride_list = lastStrideFlag_prev(stride_list, otherAngleFoot, fs, max_break, thrs, Data, headerTraj, mrk)    
    fc = 20;                            % Cutoff frequency
    Wn = fc/(fs/2);                     % Normalized cutoff frequency
    [b,a] = butter(6, Wn, 'low');       % Butterworth filter    
    
    AngVel_RF = diff(filtfilt(b,a,otherAngleFoot.R(:,1)))/(1/fs);
    AngVel_LF = diff(filtfilt(b,a,otherAngleFoot.L(:,1)))/(1/fs);    
    
    Max_AngVel_RF = max(AngVel_RF);
    Sw_ph_find_R =(AngVel_RF>(Max_AngVel_RF*thrs))*Max_AngVel_RF; 
    
    Max_AngVel_LF = max(AngVel_LF);
    Sw_ph_find_L =(AngVel_LF>(Max_AngVel_LF*thrs))*Max_AngVel_LF; 
    
    Rmrk = strcat('R',mrk); Lmrk = strcat('L',mrk);
    
%     N = 3;
%     ZVU_R = abs(AngVel_RF);
%     for i=1:size(AngVel_RF,1)-N
%         ZVU_R(i) = sum(ZVU_R(i:i+N))/N;
%     end
%     
%     thrs1 = 0.05;
%     Max_ZVU_R = max(ZVU_R);
%     Sw_ph_find_R2 =(ZVU_R>(Max_ZVU_R*thrs1))*Max_AngVel_RF;   
    
%     ZVU_L = abs(AngVel_LF);
%     for i=1:size(AngVel_LF,1)-N
%         ZVU_L(i)=sum(ZVU_L(i:i+N))/N;
%     end
%     
%     Max_ZVU_L = max(ZVU_L);
%     Sw_ph_find_L2 =(ZVU_L>(Max_ZVU_L*thrs1))*Max_AngVel_LF;
    
    figure
    subplot(211)
    plot(AngVel_RF,'g')
    hold on    
    plot(Sw_ph_find_R,'r');
%     plot(Sw_ph_find_R2,'--r')
    IC_R = [stride_list.R.ICEvents]; IC_R = unique(IC_R);
    plot(IC_R, AngVel_RF(IC_R),'*k');
    
    subplot(212)
    hold on
    plot(AngVel_LF,'b')
    plot(Sw_ph_find_L,'m');
%     plot(Sw_ph_find_L2,'--m')
    IC_L = [stride_list.L.ICEvents]; IC_L = unique(IC_L);
    plot(IC_L, AngVel_LF(IC_L),'*k');
    
    %% ADD last stride flag
    rTraj = Data.(headerTraj).(Rmrk)(:,3); maxRtraj = max(rTraj);
    sw_dR = [];
    for i = 1:size(stride_list.R,2)
        IC_now = stride_list.R(i).ICEvents(1,2)+1;
        sw_dR = [sw_dR; size(find(Sw_ph_find_R(stride_list.R(i).ICEvents(1,1):stride_list.R(i).ICEvents(1,2))),1)/fs];
        if IC_now + max_break*fs + fs < length(Sw_ph_find_R)
            check = find(Sw_ph_find_R(IC_now:IC_now + max_break*fs + fs));
            trajNow = Data.(headerTraj).(Rmrk)(IC_now:IC_now + max_break*fs + fs,3);
%             check2 = find(Sw_ph_find_R2(IC_now:IC_now + max_break*fs + fs),1);
        else
            check = find(Sw_ph_find_R(IC_now:end));
            trajNow = Data.(headerTraj).(Rmrk)(IC_now:end,3);
%             check2 = find(Sw_ph_find_R2(IC_now:end),1);
        end        
        if isempty(check) 
            if max(trajNow)>0.5*maxRtraj                
                stride_list.R(i).LastStrideFlag = false(1);
            else
                stride_list.R(i).LastStrideFlag = true(1);
            end
        else
%             if isempty(check2)
%                 stride_list.R(i).LastStrideFlag = true(1);
%             else
                stride_list.R(i).LastStrideFlag = false(1);
%             end
        end
    end
    
    lTraj = Data.(headerTraj).(Lmrk)(:,3); maxLtraj = max(lTraj);
    sw_dL = [];
    for i = 1:size(stride_list.L,2)
        IC_now = stride_list.L(i).ICEvents(1,2)+1;
        sw_dL = [sw_dL; size(find(Sw_ph_find_L(stride_list.L(i).ICEvents(1,1):stride_list.L(i).ICEvents(1,2))),1)/fs];
        if IC_now + max_break*fs + fs < length(Sw_ph_find_L)
            check = find(Sw_ph_find_L(IC_now:IC_now + max_break*fs + fs));
            trajNow = Data.(headerTraj).(Lmrk)(IC_now:IC_now + max_break*fs + fs,3);
%             check2 = find(Sw_ph_find_L2(IC_now:IC_now + max_break*fs + fs),1);
        else
            check = find(Sw_ph_find_L(IC_now:end));
            trajNow = Data.(headerTraj).(Lmrk)(IC_now:end,3);
%             check2 = find(Sw_ph_find_L2(IC_now:end));
        end 
        if isempty(check) 
            if max(trajNow)>0.5*maxLtraj
                stride_list.L(i).LastStrideFlag = false(1);
            else
                stride_list.L(i).LastStrideFlag = true(1);
            end
        else
%             if isempty(check2)
%                 stride_list.L(i).LastStrideFlag = true(1);
%             else
                if  i == size(stride_list.L,2) % this is the last stride that has been identified
                    sw_Now = length(check)/fs;
                    if sw_Now < mean(sw_dL)-2*std(sw_dL)
                        stride_list.L(i).LastStrideFlag = true(1);
                    else
                        stride_list.L(i).LastStrideFlag = false(1);
                    end
                    
                else
                    stride_list.L(i).LastStrideFlag = false(1);
                end
                
%             end
        end
    end    
end