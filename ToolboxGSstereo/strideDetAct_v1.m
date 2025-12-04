function [stride_listA] = strideDetAct_v1 (stride_list,Act)

if ~isempty(stride_list.R) && ~isempty(stride_list.L)
    stride_listA = stride_list;    
% Activity recognition - Beginning and end
    beg = [];
    fin = [];
    if (Act(1)==1)
        beg = [beg; 1];
    end
    for i = 2:length(Act)
        if Act(i)>Act(i-1)
            beg = [beg; i];
        elseif Act(i)<Act(i-1)
            fin = [fin; i];
        end
    end
    if length(fin)<length(beg)
        fin = [fin; length(Act)];
    end
    
    % Left Limb
    for i = 1:size(stride_list.L,2)
%         if size(beg,1)<=size(stride_list.L,2)
            for j = 1:length(beg)                
                Flagstride_listA.L(i,j) = (stride_list.L(i).ICEvents(1,1)>beg(j)||stride_list.L(i).ICEvents(1,2)>beg(j)) && stride_list.L(i).ICEvents(1,2)<fin(j);
            end
%         end
    end
    
    % Identify strides that belong to actvity regions
    for k = 1:size(Flagstride_listA.L,1)
        stride_listA.L(k).Activity = ~isempty(find(Flagstride_listA.L(k,:),1));
        if ~stride_listA.L(k).Activity && stride_listA.L(k).TrueStrideFlag
            stride_listA.L(k).TrueStrideFlag = false;
        end
    end
    
    % Right Limb
    for i = 1:size(stride_list.R,2)
%         if size(beg,1)<=size(stride_list.R,2)
            for j = 1:length(beg)
                Flagstride_listA.R(i,j) = (stride_list.R(i).ICEvents(1,1)>beg(j)||stride_list.R(i).ICEvents(1,2)>beg(j)) && stride_list.R(i).ICEvents(1,2)<fin(j);
            end
%         end
    end  
    
    % Identify strides that belong to actvity regions
    for k = 1:size(Flagstride_listA.R,1)
        stride_listA.R(k).Activity = ~isempty(find(Flagstride_listA.R(k,:),1));
        if ~stride_listA.R(k).Activity && stride_listA.R(k).TrueStrideFlag
            stride_listA.R(k).TrueStrideFlag = false;
        end
    end    
else
    stride_listA.R = stride_list.R;
    stride_listA.L = stride_list.L;
end
end
