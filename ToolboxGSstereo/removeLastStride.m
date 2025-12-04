% the last stride has to be removed - final IC is not a
% real IC for another stride
function [allStrides, WB_e, StrNowR, StrNowL] = removeLastStride(allStrides, StrNowR, StrNowL, stride_list)

%% 1) find the relevant strides to remove
Ls = allStrides(allStrides(:,7)==0,:);
allIC_l = [stride_list.L.ICEvents]; allIC_l = allIC_l(1:2:end);
if ~isempty(Ls)&& ((size(allIC_l,2)==1)&&isnan(allIC_l(1,1)))
    CheckL = false(1);
    LastLs = [];
    
elseif ~isempty(Ls)
    LastLs = Ls(end,:);
    if ~isnan(LastLs(1,1))
        pos_L = find(LastLs(1,1)==allIC_l);
        CheckL = stride_list.L(pos_L).LastStrideFlag;
    else
        CheckL = true(1);
    end
else
    LastLs = [];
end

Rs = allStrides(allStrides(:,7)==1,:);
allIC_r = [stride_list.R.ICEvents]; allIC_r = allIC_r(1:2:end);
if ~isempty(Rs) && ((size(allIC_r,2)==1)&&isnan(allIC_r(1,1)))
    CheckR = false(1);
    LastRs = [];
    
elseif ~isempty(Rs)
    LastRs = Rs(end,:);
    if ~isnan(LastRs(1,1))
        pos_R = find(LastRs(1,1)==allIC_r);
        CheckR = stride_list.R(pos_R).LastStrideFlag;
    else
        CheckR = true(1);
    end
else
    LastRs = [];
end

% 2) Remove them from all the strides
if ~isempty(LastRs) && CheckR
    pos = find(allStrides(:,1)== LastRs(1,1));
    if ~isempty(StrNowR)
        ALL_IC1r = [StrNowR.ICEvents]; ALL_IC1r = ALL_IC1r(1:2:end);
        pos2 = find(ALL_IC1r == LastRs(1,1));
        if ~isempty(pos)
            allStrides(pos,:)=[];
            StrNowR(pos2) = [];
            stride_list.R(pos2) = [];
        else % there is a row of NaN
            %             pos = find(allStrides(:,6)== 1, 6);
            %             allStrides(pos(end),:)=[];
            n = allStrides(:,6);
            if allStrides(n(end),7) == 1
                posToRem = n(end);
            else
                posToRem = n(end-1);
            end
            allStrides(posToRem,:)=[];
        end
    end
end

if ~isempty(LastLs) && CheckL
    pos = find(allStrides(:,1)== LastLs(1,1));
    if ~isempty(StrNowL)
        ALL_IC1l = [StrNowL.ICEvents]; ALL_IC1l = ALL_IC1l(1:2:end);
        pos2 = find(ALL_IC1l == LastLs(1,1));
        if ~isempty(pos)
            allStrides(pos,:)=[];
            StrNowL(pos2) = [];
            stride_list.L(pos2) = [];
        else % there is a row of NaN
            n = allStrides(:,6);
            if allStrides(n(end),7) == 0
                posToRem = n(end);
            else
                posToRem = n(end-1);
            end
            allStrides(posToRem,:)=[];
        end
    end
end

%% check if there are other "last strides" to remove
posToRemR = [];
if ((size(allIC_r,2)>1) && ~isnan(allIC_r(1,1)))
    lastStrideFlagR = find([StrNowR.LastStrideFlag]);
    if ~isempty(lastStrideFlagR)
        for k = 1:size(lastStrideFlagR,2)
            posToRemR = [posToRemR; find(allStrides(:,1)== StrNowR(lastStrideFlagR(k)).ICEvents(1,1))];
        end
    end
end

posToRemL = [];
if ((size(allIC_l,2)>1) && ~isnan(allIC_l(1,1)))
    lastStrideFlagL = find([StrNowL.LastStrideFlag]);
    if ~isempty(lastStrideFlagL)
        for k = 1:size(lastStrideFlagL,2)
            posToRemL = [posToRemL; find(allStrides(:,1) == StrNowL(lastStrideFlagL(k)).ICEvents(1,1))];
        end
    end
end

if ~isempty(posToRemR) && ~isempty(posToRemL)
    allpos = sort([posToRemR;posToRemL]);
    StrNowR(lastStrideFlagR) = [];
    StrNowL(lastStrideFlagL) = [];
    allStrides(allpos,:)=[];
    
elseif ~isempty(posToRemR)
    if min(posToRemR) == size(allStrides,1)
        allStrides(posToRemR,:)=[];
        StrNowR(lastStrideFlagR) = [];
    else
        if ~isempty(StrNowL)
            if StrNowL(end).LastStrideFlag

            elseif StrNowL(end).Length < 0.2* mean([StrNowL.Length])
                % also this stride should be removed
                posToRemL = [posToRemL; find(allStrides(:,1) == StrNowL(end).ICEvents(1,1))];
                allpos = sort([posToRemR;posToRemL]);

                StrNowR(lastStrideFlagR) = [];
                StrNowL(end) = [];
                allStrides(allpos,:)=[];
            end
        end
    end
    
elseif ~isempty(posToRemL)
    if min(posToRemL) == size(allStrides,1)
        allStrides(posToRemL,:)=[];
        StrNowL(lastStrideFlagL) = [];
    else
        if ~isempty(StrNowR)
            if StrNowR(end).LastStrideFlag

            elseif StrNowR(end).Length < 0.2* mean([StrNowR.Length])
                % also this stride should be removed
                posToRemR = [posToRemR; find(allStrides(:,1) == StrNowR(end).ICEvents(1,1))];
                allpos = sort([posToRemR;posToRemL]);

                StrNowR(end) = [];
                StrNowL(lastStrideFlagL) = [];
                allStrides(allpos,:)=[];
            end
        end
    end
end

%% Check if the last stride is a correct one
if ~isempty(allStrides)
    while allStrides(end, end) == 0
        ALL_IC1r = [StrNowR.ICEvents]; ALL_IC1r = ALL_IC1r(1:2:end);
        ALL_IC1l = [StrNowL.ICEvents]; ALL_IC1l = ALL_IC1l(1:2:end);
        pos_toRem = find(ALL_IC1r == allStrides(end, 1),1);
        if isempty(pos_toRem)
            pos_toRem = find(ALL_IC1l == allStrides(end, 1),1);
            StrNowL(pos_toRem) = [];
        else
            StrNowR(pos_toRem) = [];
        end
        allStrides(end, :) = [];
    end
end

%% Check if there are at the end two consecutive strides from the same side
if ~isempty(allStrides) && size(allStrides,1)>1
    if allStrides(end, end-1) == allStrides(end-1, end-1)
        % the last stride has to be removed!
        if allStrides(end, end-1) == 1
            StrNowR(end) = [];
        else
            StrNowL(end) = [];
        end
        allStrides(end, :) = [];
    end
end

if ~isempty(allStrides)
    while (isnan(allStrides(end,1)))&& size(allStrides,1)>1
        allStrides(end,:) = [];
    end
    if (isnan(allStrides(end,1)))&& size(allStrides,1)==1
        allStrides(end,:) = [];
    end
end

%% 3) update WB termination
if size(allStrides,1)>2
    WB_e = max(allStrides(end-1:end,2));
    while isnan(WB_e)
        allStrides(end,:) = [];
        if ~isempty(allStrides)
            WB_e = max(allStrides(end-1:end,2));
        else
            WB_e = [];
        end
    end
else
    if ~isempty(allStrides)
        WB_e = allStrides(end,2);
    else
        WB_e = [];
    end
end

allStrides(:,6) = 1:size(allStrides,1);
end