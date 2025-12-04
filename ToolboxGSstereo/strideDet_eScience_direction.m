function [stride_list] = strideDet_eScience_direction(GEs, max_st, min_st, min_sl,  max_h, Data, headerTraj, fs, rMrkL, lMrkL)

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
rHS = GEs.HS.right;
lHS = GEs.HS.left;
rTO = GEs.TO.right;
lTO = GEs.TO.left;

%% Strides that last longer than max_st and less than min_st - 0
RstrideD = diff(rHS(:,1))/fs;
RstrideD_ToRem = find(RstrideD > max_st | RstrideD < min_st);
LstrideD = diff(lHS(:,1))/fs ;
LstrideD_ToRem = find(LstrideD > max_st | LstrideD < min_st);

%% Stride Length - 1
rL = [];
rCL = []; % Foot Vertical Clearance
for  i = 1:size(rHS,1)-1
    rL = [rL; sqrt((Data.(headerTraj).RHEEL(rHS(i+1,1),1)-Data.(headerTraj).RHEEL(rHS(i,1),1))^2+...
                    (Data.(headerTraj).RHEEL(rHS(i+1,1),2)-Data.(headerTraj).RHEEL(rHS(i,1),2))^2)];
    rCL = [rCL; min(Data.(headerTraj).RTOE(rHS(i:i+1,1),3))];
end
rL_ToRem = find(rL < min_sl);

lL = [];
lCL = []; % Foot Vertical Clearance
for  i = 1:size(lHS,1)-1
    lL = [lL; sqrt((Data.(headerTraj).LHEEL(lHS(i+1,1),1)-Data.(headerTraj).LHEEL(lHS(i,1),1))^2+...
                    (Data.(headerTraj).LHEEL(lHS(i+1,1),2)-Data.(headerTraj).LHEEL(lHS(i,1),2))^2)];
    lCL = [lCL; min(Data.(headerTraj).LTOE(lHS(i:i+1,1),3))];
end
lL_ToRem = find(lL < min_sl);

%% Stride Height - 2
rH = [];
for  i = 1:size(rHS,1)-1
    rH = [rH; (Data.(headerTraj).RHEEL(rHS(i+1,1),3)-Data.(headerTraj).RHEEL(rHS(i,1),3))];
end
rH_ToRem = find(rH > max_h);

lH = [];
for  i = 1:size(lHS,1)-1
    lH = [lH; (Data.(headerTraj).LHEEL(lHS(i+1,1),3)-Data.(headerTraj).LHEEL(lHS(i,1),3))];
end
lH_ToRem = find(lH > max_h);

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

%% For each stride, identify the direction: forwards, lateral, backwards
selStrideR = stride_direction(selStrideR, rMrkL, 'R', Data);
selStrideL = stride_direction(selStrideL, lMrkL, 'L', Data);

stride_list.R =  selStrideR;
stride_list.L =  selStrideL;
end