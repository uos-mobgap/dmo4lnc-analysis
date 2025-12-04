function [WB, stride_list, TerminationReason] = createLR_WB_v1(stride_list, TurnDur, posTurn, quick_Turn_pos, n_min)
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
    
    %% WB assembly
    k = 1;
    for i = 1:size(stride_list,2)
        if ~stride_list(i).TrueStrideFlag || stride_list(i).SharpTurnFlag || stride_list(i).MaxHeightFlag
            if stride_list(i).MaxTimeFlag                
                k = [k; i];
                if (k(end)-k(end-1)-1) >= n_min
                    if (k(end-1)==1)
                        WB = [WB; stride_list(k(end-1)).ICEvents(1,1), stride_list(k(end)-1).ICEvents(1,2), k(end)-k(end-1)-1];
                        TerminationReason(size(WB,1)).WB = 'Pause';
                    else
                        WB = [WB; stride_list(k(end-1)+1).ICEvents(1,1), stride_list(k(end)-1).ICEvents(1,2), k(end)-k(end-1)-1];
                        TerminationReason(size(WB,1)).WB = 'Pause';
                    end
                    
                end
            end
            if stride_list(i).SharpTurnFlag                
                k = [k; i];
                if (k(end)-k(end-1)-1) >= n_min
                    if (k(end-1)==1)
                        WB = [WB; stride_list(k(end-1)).ICEvents(1,1), stride_list(k(end)-1).ICEvents(1,2), k(end)-k(end-1)-1];
                        TerminationReason(size(WB,1)).WB = 'SharpTurn';
                    else
                        WB = [WB; stride_list(k(end-1)+1).ICEvents(1,1), stride_list(k(end)-1).ICEvents(1,2), k(end)-k(end-1)-1];
                        TerminationReason(size(WB,1)).WB = 'SharpTurn';
                    end
                end
            end
            if stride_list(i).MaxHeightFlag  
                if stride_list(i+1).MaxHeightFlag
                    k = [k; i];
                    if (k(end)-k(end-1)-1) >= n_min
                        if (k(end-1)==1)
                            WB = [WB; stride_list(k(end-1)).ICEvents(1,1), stride_list(k(end)-1).ICEvents(1,2), k(end)-k(end-1)-1];
                            TerminationReason(size(WB,1)).WB = 'InclineWalking';
                        else
                            WB = [WB; stride_list(k(end-1)+1).ICEvents(1,1), stride_list(k(end)-1).ICEvents(1,2), k(end)-k(end-1)-1];
                            TerminationReason(size(WB,1)).WB = 'InclineWalking';
                        end
                    end
                end
            end
        end
        if i == size(stride_list,2)             
            if stride_list(i).TrueStrideFlag
                k = [k; i];
            else
                k = [k; i-1];
            end
            if (k(end)-k(end-1)) >= n_min
                if (k(end-1)==1)
                    WB = [WB; stride_list(k(end-1)).ICEvents(1,1), stride_list(k(end)).ICEvents(1,2), k(end)-k(end-1)];
                    TerminationReason(size(WB,1)).WB = 'Pause';
                else
                    WB = [WB; stride_list(k(end-1)+1).ICEvents(1,1), stride_list(k(end)).ICEvents(1,2), k(end)-k(end-1)];
                    TerminationReason(size(WB,1)).WB = 'Pause';
                end
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