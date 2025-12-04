function [stride_list] = stride_INDIP_GEs_INDIPmrk (GEsINDIP, GEsIN_LR, max_st, min_st, min_sl,  max_h, Data, headerTraj, fs, mrkToUse)

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
rHS = round(GEsINDIP.IC(strcmp(GEsIN_LR.IC,'Right'))*fs);
rHS = rHS'; rHS(:,2) = ones(length(rHS),1);
lHS = round(GEsINDIP.IC(strcmp(GEsIN_LR.IC,'Left'))*fs);
lHS = lHS'; lHS(:,2) = ones(length(lHS),1);
rTO = round(GEsINDIP.FC(strcmp(GEsIN_LR.FC,'Right'))*fs);
rTO = rTO'; rTO(:,2) = ones(length(rTO),1);
lTO = round(GEsINDIP.FC(strcmp(GEsIN_LR.FC,'Left'))*fs);
lTO = lTO'; lTO(:,2) = ones(length(lTO),1);

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
        if ~isnan(rHS(i,1)) && ~isnan(rHS(i+1,1))
        rL = [rL; sqrt((Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i+1,1),1)-Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i,1),1))^2+...
                        (Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i+1,1),2)-Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i,1),2))^2+...
                        (Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i+1,1),3)-Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i,1),3))^2)];
        rCL = [rCL; min(Data.(headerTraj).RTOE(rHS(i:i+1,1),3))];
        else
            rL = [rL; 0];
            rCL = [rCL; 0];
        end
    end
    rL_ToRem = find(rL < min_sl);
    if size(rL_ToRem,2)==0
        rL_ToRem = double.empty(0,1);
    end

    lL = [];
    lCL = []; % Foot Vertical Clearance
    for  i = 1:size(lHS,1)-1
        if ~isnan(lHS(i,1)) && ~isnan(lHS(i+1,1))
            lL = [lL; sqrt((Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i+1,1),1)-Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i,1),1))^2+...
                            (Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i+1,1),2)-Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i,1),2))^2+...
                            (Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i+1,1),3)-Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i,1),3))^2)];
            lCL = [lCL; min(Data.(headerTraj).LTOE(lHS(i:i+1,1),3))];
        else
            lL = [lL; 0];
            lCL = [lCL; 0];
        end
    end
    lL_ToRem = find(lL < min_sl);
    if size(lL_ToRem,2)==0
        lL_ToRem = double.empty(0,1);
    end

    %% Stride Height - 2
    rH = [];
    for  i = 1:size(rHS,1)-1
        if ~isnan(rHS(i,1)) && ~isnan(rHS(i+1,1))
            rH = [rH; (Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i+1,1),3)-Data.(headerTraj).(strcat('R', mrkToUse))(rHS(i,1),3))];
        else
            rH = [rH; 0];
        end
    end
    rH_ToRem = find(rH > max_h);
    if size(rH_ToRem,2)==0
        rH_ToRem = double.empty(0,1);
    end

    lH = [];
    for  i = 1:size(lHS,1)-1
        if ~isnan(lHS(i,1)) && ~isnan(lHS(i+1,1))
            lH = [lH; (Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i+1,1),3)-Data.(headerTraj).(strcat('L', mrkToUse))(lHS(i,1),3))];
        else
            lH = [lH; 0];
        end
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

    if ~isempty(rSide)
        for i = 1:length(rSide)
            d = find(RstrideD_ToRem(:,1) == rSide(i));
            if ~isempty(d) a = d; else a = nan; end
            l = find(rL_ToRem(:,1) == rSide(i));
            if ~isempty(l) b = l; else b = nan; end
            h = find(rH_ToRem(:,1) == rSide(i));
            if ~isempty(h) c = h; else c = nan; end
            reasonR(i,:) = [a, b, c];
        end
    else
        reasonR(1,:) = [nan, nan, nan];
    end
    if ~isempty(lSide)
        for i = 1:length(lSide)
            d = find(LstrideD_ToRem(:,1) == lSide(i));
            if ~isempty(d) a = d; else a = nan; end
            l = find(lL_ToRem(:,1) == lSide(i));
            if ~isempty(l) b = l; else b = nan; end
            h = find(lH_ToRem(:,1) == lSide(i));
            if ~isempty(h) c = h; else c = nan; end
            reasonL(i,:) = [a, b, c];
        end
    else
        reasonL(1,:) = [nan, nan, nan];
    end

    %% For each stride, start, stop, and length
    [~,selStrideR] = strideSelection(rHS, rSide, rL, rH, rTO, reasonR, rCL, RstrideD, max_h, min_st, max_st, min_sl,'R');
    [~,selStrideL] = strideSelection(lHS, lSide, lL, lH, lTO, reasonL, lCL, LstrideD, max_h, min_st, max_st, min_sl,'L');
    ICs = [selStrideR.ICEvents]; ICs = ICs(1:2:end);
    if ~isempty(find(isnan(ICs), 1))
        stride_list.R =  selStrideR;
        pos = find(isnan(ICs));
        if pos>1 && length(pos)==1
            selStrideR(pos-1).ICEvents = [selStrideR(pos-1).ICEvents(1,1), selStrideR(pos).ICEvents(1,2)];
            selStrideR(pos) = [];
        else
            warning('this case has not been addressed yet')
        end
    end
    stride_list.R =  selStrideR;
    
    ICs = [selStrideL.ICEvents]; ICs = ICs(1:2:end);
    if ~isempty(find(isnan(ICs), 1))
        pos = find(isnan(ICs));
        if pos>1 && length(pos)==1
            selStrideL(pos-1).ICEvents = [selStrideL(pos-1).ICEvents(1,1), selStrideL(pos).ICEvents(1,2)];
            selStrideL(pos) = [];
        else
            warning('this case has not been addressed yet')
        end
    end
    stride_list.L =  selStrideL;
    
    
    %% CHECK if there are NaNs in the GEs
    if ~isreal(lL) || ~isreal(rL)
        checkthis = 1;
    end
else
    stride_list.R =  [];
    stride_list.L =  [];
    fprintf('Not enough GEs for defining strides!\n');
end
end