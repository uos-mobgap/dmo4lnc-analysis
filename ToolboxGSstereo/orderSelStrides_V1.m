function [allStrides, TurnStrides, allStridesStruct] = orderSelStrides_V1(SelStrR,SelStrL)
% function that puts in order the different strides
% 1 is assigned to R strides
% 0 is assigned to L strides

%% 
% OUTPUT: for each row a stride that has been identified - the following
% info are reported
% (1) IC1, (2) IC2, (3) Length, (4) VerticalDispl, (5) FC, (6)
% TrueStrideFlag, (7) n, (8) L(==0) or R(==1)
%%

    i = 1;
    n = 1;
    allStrides = [];
    TurnStrides = [];

    if ~isempty(SelStrR) && ~isempty(SelStrL)
        while i <= size(SelStrR,2) && i <= size(SelStrL,2)
            if SelStrR(i).ICEvents(1,1)<SelStrL(i).ICEvents(1,1)
                if isempty(SelStrR(i).FCEvents); SelStrR(i).FCEvents = 0; end
                allStrides = [allStrides; SelStrR(i).ICEvents, SelStrR(i).Length, ...
                    SelStrR(i).VerticalDispl, SelStrR(i).FCEvents, SelStrR(i).TrueStrideFlag,n,1];
                allStridesStruct(n) = SelStrR(i); n = n+1;
                TurnStrides = [TurnStrides; SelStrR(i).TurnFlag];
                
                if isempty(SelStrL(i).FCEvents); SelStrL(i).FCEvents = 0; end
                allStrides = [allStrides; SelStrL(i).ICEvents, SelStrL(i).Length, ...
                    SelStrL(i).VerticalDispl, SelStrL(i).FCEvents, SelStrL(i).TrueStrideFlag,n,0];
                allStridesStruct(n) = SelStrL(i); n = n+1;
                TurnStrides = [TurnStrides; SelStrL(i).TurnFlag];
            else
                if isempty(SelStrL(i).FCEvents); SelStrL(i).FCEvents = 0; end
                allStrides = [allStrides; SelStrL(i).ICEvents, SelStrL(i).Length, ...
                    SelStrL(i).VerticalDispl, SelStrL(i).FCEvents, SelStrL(i).TrueStrideFlag,n,0]; 
                allStridesStruct(n) = SelStrL(i); n = n+1;
                TurnStrides = [TurnStrides; SelStrL(i).TurnFlag];
                
                if isempty(SelStrR(i).FCEvents); SelStrR(i).FCEvents = 0; end
                allStrides = [allStrides; SelStrR(i).ICEvents, SelStrR(i).Length, ...
                    SelStrR(i).VerticalDispl, SelStrR(i).FCEvents, SelStrR(i).TrueStrideFlag,n,1];
                allStridesStruct(n) = SelStrR(i); n = n+1;
                TurnStrides = [TurnStrides; SelStrR(i).TurnFlag];
            end
            i = i+1;
        end

        if size(allStrides,1) ~= size(SelStrR,2)+size(SelStrL,2)
            i = size(allStrides(allStrides(:,end)==1),1)+1;
            if size(SelStrR,2)>size(SelStrL,2)                
                while i<=size(SelStrR,2)
                    if isempty(SelStrR(i).FCEvents); SelStrR(i).FCEvents = 0; end
                    allStrides = [allStrides; SelStrR(i).ICEvents, SelStrR(i).Length, ...
                        SelStrR(i).VerticalDispl, SelStrR(i).FCEvents, SelStrR(i).TrueStrideFlag,n,1];
                    TurnStrides = [TurnStrides; SelStrR(i).TurnFlag];
                    k_now = size(allStridesStruct,2);
                    allStridesStruct(k_now+1) = SelStrR(i); 
                    i = i+1;                
                end
            else
                
                while i<=size(SelStrL,2)
                    if isempty(SelStrL(i).FCEvents); SelStrL(i).FCEvents = 0; end
                    allStrides = [allStrides; SelStrL(i).ICEvents, SelStrL(i).Length, ...
                        SelStrL(i).VerticalDispl, SelStrL(i).FCEvents, SelStrL(i).TrueStrideFlag,n,0];
                    TurnStrides = [TurnStrides; SelStrL(i).TurnFlag];
                    k_now = size(allStridesStruct,2);
                    allStridesStruct(k_now+1) = SelStrL(i); 
                    i = i+1;
                end
            end
        end

    elseif ~isempty(SelStrR)
        while n<=size(SelStrR,2)
            if isempty(SelStrR(i).FCEvents); SelStrR(i).FCEvents = 0; end
            allStrides = [allStrides; SelStrR(i).ICEvents, SelStrR(i).Length, ...
                        SelStrR(i).VerticalDispl, SelStrR(i).FCEvents, SelStrR(i).TrueStrideFlag,n,1]; n = n+1;
            TurnStrides = [TurnStrides; SelStrR(i).TurnFlag];
            k_now = size(allStrides,1);
            allStridesStruct(k_now) = SelStrR(i);
            i = i+1;
        end
    else
        while i<=size(SelStrL,2)
           if isempty(SelStrL(i).FCEvents); SelStrL(i).FCEvents = 0; end
           allStrides = [allStrides; SelStrL(i).ICEvents, SelStrL(i).Length, ...
                        SelStrL(i).VerticalDispl, SelStrL(i).FCEvents, SelStrL(i).TrueStrideFlag,n,0]; n = n+1;
           TurnStrides = [TurnStrides; SelStrL(i).TurnFlag];
           k_now = size(allStrides,1);
           allStridesStruct(k_now) = SelStrL(i);
           i = i+1;
        end    
    end
    [~,pos] = sort(allStrides(:,1));
    allStrides = allStrides(pos,:);
    TurnStrides = TurnStrides(pos);
    allStridesStruct = allStridesStruct(pos);
    allStrides(:,7) = 1:size(allStrides,1)';
    
    % check L/R alternance
    durationOverlap = zeros(size(allStrides,1)-1,1);
    durationOverlap2 = zeros(size(allStrides,1)-1,1);
    for i = 1:size(allStrides,1)-1
        durationOverlap(i)= ((length(intersect(allStrides(i,1):allStrides(i,2),allStrides(i+1,1):allStrides(i+1,2))))/...
            length(allStrides(i,1):allStrides(i,2)))*100;
        durationOverlap2(i)= ((length(intersect(allStrides(i,1):allStrides(i,2),allStrides(i+1,1):allStrides(i+1,2))))/...
            length(allStrides(i+1,1):allStrides(i+1,2)))*100;
    end
end