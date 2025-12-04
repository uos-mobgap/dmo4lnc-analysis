function [stride_list] = strideDet_eScience_mrk_GEeffect_v1 (GEs, max_st, min_st, min_sl,  max_h, Data, headerTraj, fs, mrkToUse, Turn)

%% INPUTs
% GEs (HS and TO for L and R - i.e., rHS, lHS, rTO, lTO) [n x 1]
% max_st = 3;         % Maximum stride duration accepted for each stride to be considered a correct stride [s]
% min_st = 0.2;       % Minimum stride duration accepted for each stride to be considered a correct stride [s]
% min_sl = 0.15;      % Minimum stride length accepted for each stride to be considered a correct stride [m]
% max_h = 0.21;       % Max height difference between initial and final stride instants [m]

%% OUTPUTs
% rStrides             start, stop, and length and H for the R limb; TO
% lStrides             start, stop, and length and H for the L limb; TO

%%
rHS = round(GEs.HS.right);
lHS = round(GEs.HS.left);
rTO = round(GEs.TO.right);
lTO = round(GEs.TO.left);

if (~isempty(rHS)||size(rHS,1)>1) && (~isempty(lHS)||size(lHS,1)>1)
    %% Strides that last longer than max_st and less than min_st - 0
    RstrideD = diff(rHS(:,1))/fs;
    RstrideD_ToRem = find(RstrideD > max_st | RstrideD < min_st);
    if size(RstrideD_ToRem,2)==0
        RstrideD_ToRem = double.empty(0,1);
    end
    
    LstrideD = diff(lHS(:,1))/fs ;
    LstrideD_ToRem = find(LstrideD > max_st | LstrideD < min_st);
    if size(LstrideD_ToRem,2)==0
        LstrideD_ToRem = double.empty(0,1);
    end

    %% Stride Length - 1
    rL = [];
    rCL = []; % Foot Vertical Clearance
    for  i = 1:size(rHS,1)-1
        rL = [rL; sqrt((Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i+1,1),1)-Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i,1),1))^2+...
                        (Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i+1,1),2)-Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i,1),2))^2+...
                        (Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i+1,1),3)-Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i,1),3))^2)];
        rCL = [rCL; min(Data.(headerTraj).RTOE(rHS(i:i+1,1),3))];
    end
    rL_ToRem = find(rL < min_sl);
    if size(rL_ToRem,2)==0
        rL_ToRem = double.empty(0,1);
    end

    lL = [];
    lCL = []; % Foot Vertical Clearance
    for  i = 1:size(lHS,1)-1
        lL = [lL; sqrt((Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i+1,1),1)-Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i,1),1))^2+...
                            (Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i+1,1),2)-Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i,1),2))^2+...
                            (Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i+1,1),3)-Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i,1),3))^2)];
        lCL = [lCL; min(Data.(headerTraj).LTOE(lHS(i:i+1,1),3))];
    end
    lL_ToRem = find(lL < min_sl);
    if size(lL_ToRem,2)==0
        lL_ToRem = double.empty(0,1);
    end

    %% Stride Height - 2
    rH = [];
    for  i = 1:size(rHS,1)-1
        rH = [rH; (Data.(headerTraj).RHEEL(rHS(i+1,1),3)-Data.(headerTraj).RHEEL(rHS(i,1),3))];
    end
    rH_ToRem = find(rH > max_h);
    if size(rH_ToRem,2)==0
        rH_ToRem = double.empty(0,1);
    end

    lH = [];
    for  i = 1:size(lHS,1)-1
        lH = [lH; (Data.(headerTraj).LHEEL(lHS(i+1,1),3)-Data.(headerTraj).LHEEL(lHS(i,1),3))];
    end
    lH_ToRem = find(lH > max_h);
    if size(lH_ToRem,2)==0
        lH_ToRem = double.empty(0,1);
    end

    %% Select events to be removed
    rSide = unique([RstrideD_ToRem(:,1); rL_ToRem(:,1); rH_ToRem(:,1)]);
    rSide = rSide(rSide>0);
    lSide = unique([LstrideD_ToRem(:,1); lL_ToRem(:,1); lH_ToRem(:,1)]);
    lSide = lSide(lSide>0);    

    %% For each stride, start, stop, and length
    if ~isempty(rHS) rHS = rHS(:,1); end
    if~isempty(rTO) rTO = rTO(:,1); end
    [~,selStrideR] = strideSelection_GEmeth(rHS, rSide, rL, rH, rTO, rCL, RstrideD, max_h, min_st, max_st, min_sl,'R');
    % add stance/swing duration
    for m = 1:size(selStrideR,2)
        if size(selStrideR(m).FCEvents,1)==1
            if selStrideR(m).FCEvents ~= 0
                selStrideR(m).StanceDuration = selStrideR(m).FCEvents - selStrideR(m).ICEvents(1,1);
                selStrideR(m).SwingDuration = selStrideR(m).ICEvents(1,2)- selStrideR(m).FCEvents;
            else
                selStrideR(m).StanceDuration = [];
                selStrideR(m).SwingDuration = [];
            end
        else
            selStrideR(m).StanceDuration = [];
            selStrideR(m).SwingDuration = [];
        end
    end
    
    if ~isempty(lHS) lHS = lHS(:,1); end
    if ~isempty(lTO) lTO = lTO(:,1); end
    [~,selStrideL] = strideSelection_GEmeth(lHS, lSide, lL, lH, lTO, lCL, LstrideD, max_h, min_st, max_st, min_sl,'L');
    % add stance/swing duration
    for m = 1:size(selStrideL,2)
        if size(selStrideL(m).FCEvents,1)==1
            if selStrideL(m).FCEvents ~= 0
                selStrideL(m).StanceDuration = selStrideL(m).FCEvents - selStrideL(m).ICEvents(1,1);
                selStrideL(m).SwingDuration = selStrideL(m).ICEvents(1,2)- selStrideL(m).FCEvents;
            else
                selStrideL(m).StanceDuration = [];
                selStrideL(m).SwingDuration = [];
            end
        else
            selStrideL(m).StanceDuration = [];
            selStrideL(m).SwingDuration = [];
        end
    end
    
    stride_list.R =  selStrideR;
    stride_list.L =  selStrideL;
    
    %% Stride on a Turn/Step/Other
    side = fieldnames(stride_list);
    for s = 1:length(side)
        nStrides = size(stride_list.(side{s}),2);
        for k = 1:nStrides
            
            % check stride on a step 
            if abs(stride_list.(side{s})(k).VerticalDispl)>0.10
                stride_list.(side{s})(k).StrideStep = true;
            else
                stride_list.(side{s})(k).StrideStep = false;
            end
            
            % check stride in a turn
            eventNow(:,1) = stride_list.(side{s})(k).ICEvents';
            eventNow(:,2) = zeros(size(eventNow,1),1);
            for turn = 1:size(Turn,1)
                for mm = 1:size(eventNow,1)
                    if ~isempty(intersect(eventNow(mm,1),Turn(turn,1):Turn(turn,2)))
                        eventNow(mm,2)=1;
                    end
                end
            end
            if ~isempty(find(eventNow(:,2), 1))
                stride_list.(side{s})(k).StrideTurn = true;
            else
                stride_list.(side{s})(k).StrideTurn = false;
            end
            
            % other strides
            if ~stride_list.(side{s})(k).StrideTurn && ~stride_list.(side{s})(k).StrideStep
                stride_list.(side{s})(k).StrideOther = true;
            else
                stride_list.(side{s})(k).StrideOther = false;
            end
        end
    end
    
else
    stride_list.R =  [];
    stride_list.L =  [];
    fprintf('Not enough GEs for defining strides!\n');
end
end