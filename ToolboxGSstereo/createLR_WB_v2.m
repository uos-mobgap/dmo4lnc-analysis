function [WB, stride_list, TerminationReason] = createLR_WB_v2(stride_list, TurnDur, posTurn, quick_Turn_pos, n_min)
    WB = [];
    TerminationReason = [];
    
    %% ADD a SHARP turn FLAG
    QuickTurns = TurnDur(posTurn(quick_Turn_pos),:); 
    QuickTurns_temp = QuickTurns;
    k = 1;
    while size(QuickTurns_temp,1)>=1
        for i = 1:size(stride_list,2)
            d = intersect(stride_list(i).ICEvents(1,1):1:stride_list(i).ICEvents(1,2),QuickTurns_temp(1,1):1:QuickTurns_temp(1,2));
            stride_list(i).SharpTurnFlag(1,k) = ~isempty(d);
        end
        QuickTurns_temp(1,:) = [];
        k = k+1;
    end
    if isempty(QuickTurns)
        for i = 1:size(stride_list,2)
            stride_list(i).SharpTurnFlag(1,k) = k<0;
        end
    end
    if size(stride_list(i).SharpTurnFlag,2)>1
        for i = 1:size(stride_list,2)
            check = find(stride_list(i).SharpTurnFlag==1,1);
            stride_list(i).SharpTurnFlag = ~isempty(check);
        end 
    end
    
    %% ADD a turn FLAG
    Turns_temp = TurnDur;
    k = 1;
    while size(Turns_temp,1)>=1
        for i = 1:size(stride_list,2)
            d = intersect(stride_list(i).ICEvents(1,1):1:stride_list(i).ICEvents(1,2),Turns_temp(1,1):1:Turns_temp(1,2));
            stride_list(i).TurnFlag(1,k) = ~isempty(d);
        end
        Turns_temp(1,:) = [];
        k = k+1;
    end
    if isempty(Turns_temp)
        for i = 1:size(stride_list,2)
            stride_list(i).TurnFlag(1,k) = k<0;
        end
    end
    if size(stride_list(i).TurnFlag,2)>1
        for i = 1:size(stride_list,2)
            check = find(stride_list(i).TurnFlag==1,1);
            stride_list(i).TurnFlag = ~isempty(check);
        end 
    end
    
    %% WB ELIGIBLE STRIDES
    for i = 1:size(stride_list,2)
        if stride_list(i).TrueStrideFlag && ~stride_list(i).SharpTurnFlag && ~stride_list(i).MaxHeightFlag
            stride_list(i).WB_EligibleStrides = true(1);
        else
            stride_list(i).WB_EligibleStrides = false(1);
        end
    end
    
    %% WB ELIGIBLE STRIDES - also with shorter strides
    for i = 1:size(stride_list,2)
        if ~stride_list(i).SharpTurnFlag && ~stride_list(i).MaxHeightFlag
            stride_list(i).WB_allEligibleStrides = true(1);
        else
            stride_list(i).WB_allEligibleStrides = false(1);
        end
    end
    
    %% WB assembly    
    all_Prev = [stride_list.WB_EligibleStrides];
    for i = 1:size(stride_list,2)
        if ~stride_list(i).TrueStrideFlag || stride_list(i).SharpTurnFlag || stride_list(i).MaxHeightFlag
            if stride_list(i).MaxTimeFlag    
                now = i-1;                
                if i == 1
                    prev = [];                    
                else                    
                    prev_sel = find(all_Prev==0);
                     if isempty(prev_sel)                         
                         prev = 1;   % all correct strides
                     else
                         sel = find(prev_sel<i);
                         if isempty(sel)
                             prev = [];
                         else
                             prev = prev_sel(sel(end))+1;
                         end
                     end
                end
                if (i-prev) >= n_min
                    WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(now).ICEvents(1,2), size(prev:now,2)];
                    TerminationReason(size(WB,1)).WB = 'Pause';                    
                end
            end
            
            if stride_list(i).SharpTurnFlag                
                now = i-1;                
                if i == 1
                    prev = [];                    
                else                    
                    prev_sel = find(all_Prev==0);
                     if isempty(prev_sel)
                         prev = 1; % all correct strides  
                     else
                         sel = find(prev_sel<i);
                         if isempty(sel)
                             prev = [];
                         else
                             prev = prev_sel(sel(end))+1;
                         end
                     end
                end
                if (size(prev:now,2)) >= n_min
                    WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(now).ICEvents(1,2), size(prev:now,2)];
                    TerminationReason(size(WB,1)).WB = 'SharpTurn';
                end
               
            end
            if stride_list(i).MaxHeightFlag  
                if i< size(stride_list,2)
                    if stride_list(i+1).MaxHeightFlag %% two continuous strides with a change of H can split WB
                        now = i-1;
                        if i == 1
                            prev = [];
                        else
                            prev_sel = find(all_Prev==0);
                            if isempty(prev_sel)
                                prev = 1; % all correct strides  
                            else
                                sel = find(prev_sel<i);
                                if isempty(sel)
                                    prev = [];
                                else
                                    prev = prev_sel(sel(end))+1;
                                end
                            end
                        end
                        if (i-prev) >= n_min
                            WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(now).ICEvents(1,2), size(prev:now,2)];
                            TerminationReason(size(WB,1)).WB = 'InclineWalking';
                        end
                    end                   
                end
            end
        end
        if i == size(stride_list,2)    
            if stride_list(i).WB_EligibleStrides
                now = i;
            else
                now = i-1;
            end
            prev_sel = find(all_Prev==0);
            if isempty(prev_sel)
                prev = 1;                
            else
                sel = find(prev_sel<i);
                if isempty(sel)
                    prev = 1;
                else
                    prev = prev_sel(sel(end))+1;
                    if prev == now && (stride_list(prev_sel(sel(end))).MinLengthFlag || stride_list(prev_sel(sel(end))).MaxHeightFlag ...
                            || stride_list(prev_sel(sel(end))).MinTimeFlag)
                        sel = find(prev_sel<prev_sel(sel(end)));
                        if isempty(sel)
                            prev = 1;
                        else
                            prev = prev_sel(sel(end))+1;
                        end
                    end
                end
            end
            
            if (size(prev:now,2)) >= n_min
                WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(now).ICEvents(1,2), size(prev:now,2)];
                TerminationReason(size(WB,1)).WB = 'Pause';
            end
        end
    end
    
    WB_temp = WB;
    WB_temp(:,end+1) = ones(size(WB,1),1);
    for i = 1:size(WB,1)
        if WB(i,3) < n_min
            WB_temp(i,end) = 0;
        end
    end
    WB = WB(WB_temp(:,end)==1,:);
end