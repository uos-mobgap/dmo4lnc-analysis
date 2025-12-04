function [WB, stride_list, TerminationReason] = createLR_WB_v3(stride_list, TurnDur, quick_Turn_pos, n_min, max_break, fs)
WB = [];
TerminationReason = [];

if ~isempty(stride_list)
    if size(fieldnames(stride_list),1)>1
        %% ADD a SHARP turn FLAG
        QuickTurns = TurnDur(quick_Turn_pos,:);
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
        
        %% ADD a delta FLAG
        stride_list(1).deltaTime = true(1);
        for i = 2:size(stride_list,2)
            if (stride_list(i).ICEvents(1,1) - stride_list(i-1).ICEvents(1,2))/fs < max_break
                stride_list(i).deltaTime = true(1);
            else
                stride_list(i).deltaTime = false(1);
            end
        end
        
        %% WB assembly
        i = 1;
        while i <= size(stride_list,2)
            all_Prev = [stride_list.WB_EligibleStrides];
            if ~stride_list(i).WB_EligibleStrides
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
                                prev = 1;
                            else
                                prev = prev_sel(sel(end))+1;
                            end
                        end
                    end
                    if (i-prev) >= n_min
                        allEvents = [stride_list.ICEvents];
                        ICnow = allEvents(1:2:end);
                        IC2now = allEvents(2:2:end);
                        diffNow = [IC2now(prev+1:now)-IC2now(prev:now-1)]/fs;
                        if isempty(find(diffNow>max_break,1))
                            WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(now).ICEvents(1,2), size(prev:now,2)];
                            TerminationReason(size(WB,1)).WB = 'Pause';
                        else
                            checkAll = [stride_list.deltaTime];
                            checkNow = checkAll(prev:now);
                            stridesNow = prev:now;
                            stridesToKeep = stridesNow(checkNow);
                            if size(stridesToKeep,2)>= 2 && isempty(find(diff(stridesToKeep)>1, 1))
                                WB = [WB; stride_list(min(stridesToKeep)).ICEvents(1,1), stride_list(max(stridesToKeep)).ICEvents(1,2), size(stridesToKeep,2)];
                                TerminationReason(size(WB,1)).WB = 'Pause';
                            elseif size(stridesToKeep,2)>= 2
                                posToDivide = find(diff(stridesToKeep)>1);
                                for k = 1:length(posToDivide)
                                    if k == 1
                                        stridesToKeepNow = stridesNow(1:posToDivide(k)+1);
                                    else
                                        stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                    end
                                    % first portion
                                    if size(stridesToKeepNow,2)>= 2
                                        WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                        TerminationReason(size(WB,1)).WB = 'Pause';
                                    end
                                    % second portion
                                    if length(posToDivide) ==1
                                        stridesToKeepNow = stridesNow(posToDivide+2:end);
                                    else
                                        stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                    end
                                    if size(stridesToKeepNow,2)>= 2
                                        WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                        TerminationReason(size(WB,1)).WB = 'Pause';
                                    end
                                end
                            end
                        end
                    end
                    
                elseif stride_list(i).SharpTurnFlag
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
                                prev = 1;
                            else
                                prev = prev_sel(sel(end))+1;
                            end
                        end
                    end
                    if (size(prev:now,2)) >= n_min
                        allEvents = [stride_list.ICEvents];
                        ICnow = allEvents(1:2:end);
                        IC2now = allEvents(2:2:end);
                        diffNow = [IC2now(prev+1:now)-IC2now(prev:now-1)]/fs;
                        if isempty(find(diffNow>max_break,1))
                            WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(now).ICEvents(1,2), size(prev:now,2)];
                            TerminationReason(size(WB,1)).WB = 'SharpTurn';
                        else
                            checkAll = [stride_list.deltaTime];
                            checkNow = checkAll(prev:now);
                            stridesNow = prev:now;
                            stridesToKeep = stridesNow(checkNow);
                            if size(stridesToKeep,2)>= 2 && isempty(find(diff(stridesToKeep)>1, 1))
                                WB = [WB; stride_list(min(stridesToKeep)).ICEvents(1,1), stride_list(max(stridesToKeep)).ICEvents(1,2), size(stridesToKeep,2)];
                                TerminationReason(size(WB,1)).WB = 'SharpTurn';
                                
                            elseif size(stridesToKeep,2)>= 2
                                posToDivide = find(diff(stridesToKeep)>1);
                                for k = 1:length(posToDivide)
                                    if k == 1
                                        stridesToKeepNow = stridesNow(1:posToDivide(k)+1);
                                    else
                                        stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                    end
                                    % first portion
                                    if size(stridesToKeepNow,2)>= 2
                                        WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                        TerminationReason(size(WB,1)).WB = 'SharpTurn';
                                    end
                                    % second portion
                                    if length(posToDivide) ==1
                                        stridesToKeepNow = stridesNow(posToDivide+2:end);
                                    else
                                        stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                    end
                                    if size(stridesToKeepNow,2)>= 2
                                        WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                        TerminationReason(size(WB,1)).WB = 'SharpTurn';
                                    end
                                end
                            end
                        end
                    end
                    
                elseif stride_list(i).MaxHeightFlag
                    if i< size(stride_list,2)
                        if stride_list(i+1).MaxHeightFlag %% two continuous strides with a change of H can split WB
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
                                        prev = 1;
                                    else
                                        prev = prev_sel(sel(end))+1;
                                    end
                                end
                            end
                            if (i-prev) >= n_min
                                allEvents = [stride_list.ICEvents];
                                ICnow = allEvents(1:2:end);
                                IC2now = allEvents(2:2:end);
                                diffNow = [IC2now(prev+1:now)-IC2now(prev:now-1)]/fs;
                                if isempty(find(diffNow>max_break,1))
                                    WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(now).ICEvents(1,2), size(prev:now,2)];
                                    TerminationReason(size(WB,1)).WB = 'InclineWalking';
                                else
                                    checkAll = [stride_list.deltaTime];
                                    checkNow = checkAll(prev:now);
                                    stridesNow = prev:now;
                                    stridesToKeep = stridesNow(checkNow);
                                    if size(stridesToKeep,2)>= 2 && isempty(find(diff(stridesToKeep)>1, 1))
                                        WB = [WB; stride_list(min(stridesToKeep)).ICEvents(1,1), stride_list(max(stridesToKeep)).ICEvents(1,2), size(stridesToKeep,2)];
                                        TerminationReason(size(WB,1)).WB = 'InclineWalking';
                                        
                                    elseif size(stridesToKeep,2)>= 2
                                        posToDivide = find(diff(stridesToKeep)>1);
                                        for k = 1:length(posToDivide)
                                            if k == 1
                                                stridesToKeepNow = stridesNow(1:posToDivide(k)+1);
                                            else
                                                stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                            end
                                            % first portion
                                            if size(stridesToKeepNow,2)>= 2
                                                WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                                TerminationReason(size(WB,1)).WB = 'InclineWalking';
                                            end
                                            % second portion
                                            if length(posToDivide) ==1
                                                stridesToKeepNow = stridesNow(posToDivide+2:end);
                                            else
                                                stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                            end
                                            if size(stridesToKeepNow,2)>= 2
                                                WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                                TerminationReason(size(WB,1)).WB = 'InclineWalking';
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            stride_list(i).WB_EligibleStrides = true(1);
                        end
                    end
                else
                    % A stride that do not satisfy the conditions - too short (either
                    % in duration or in distance)
                    if i<size(stride_list,2) && i>1
                        % find previous correct stride
                        prev_s = find(all_Prev(1:i)==1);
                        if isempty(prev_s)
                            prev_s = [];
                        else
                            prev_s = prev_s(end);
                        end
                        % find next correct stride
                        next_s = find(all_Prev(i:end)==1);
                        if isempty(next_s)
                            next_s = [];
                        else
                            next_s = next_s(1)+i-1;
                        end
                        if ~isempty(next_s) && ~isempty(prev_s)
                            if (stride_list(next_s).ICEvents(1,1)-stride_list(prev_s).ICEvents(1,2))< fs*max_break
                                % OK!
                                stride_list(i).WB_EligibleStrides = true(1);
                            else
                                %% WB has to be terminated
                                % find previous correct stride
                                prev_sel = find(all_Prev==0);
                                if isempty(prev_sel)
                                    prev = 1; % all correct strides
                                else
                                    sel = find(prev_sel<i);
                                    if isempty(sel)
                                        prev = 1;
                                    else
                                        prev = prev_sel(sel(end))+1;
                                    end
                                end
                                if (size(prev:prev_s,2)) >= n_min
                                    allEvents = [stride_list.ICEvents];
                                    ICnow = allEvents(1:2:end);
                                    IC2now = allEvents(2:2:end);
                                    diffNow = [IC2now(prev+1:prev_s)-IC2now(prev:prev_s-1)]/fs;
                                    if isempty(find(diffNow>max_break,1))
                                        WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(prev_s).ICEvents(1,2), size(prev:prev_s,2)];
                                        TerminationReason(size(WB,1)).WB = 'Pause';
                                    else
                                        checkAll = [stride_list.deltaTime];
                                        checkNow = checkAll(prev:prev_s);
                                        stridesNow = prev:prev_s;
                                        stridesToKeep = stridesNow(checkNow);
                                        if size(stridesToKeep,2)>= 2 && isempty(find(diff(stridesToKeep)>1, 1))
                                            WB = [WB; stride_list(min(stridesToKeep)).ICEvents(1,1), stride_list(max(stridesToKeep)).ICEvents(1,2), size(stridesToKeep,2)];
                                            TerminationReason(size(WB,1)).WB = 'Pause';
                                            
                                        elseif size(stridesToKeep,2)>= 2
                                            posToDivide = find(diff(stridesToKeep)>1);
                                            for k = 1:length(posToDivide)
                                                if k == 1
                                                    stridesToKeepNow = stridesNow(1:posToDivide(k)+1);
                                                else
                                                    stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                                end
                                                % first portion
                                                if size(stridesToKeepNow,2)>= 2
                                                    WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                                    TerminationReason(size(WB,1)).WB = 'Pause';
                                                end
                                                % second portion
                                                if length(posToDivide) ==1
                                                    stridesToKeepNow = stridesNow(posToDivide+2:end);
                                                else
                                                    stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                                end
                                                if size(stridesToKeepNow,2)>= 2
                                                    WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                                    TerminationReason(size(WB,1)).WB = 'Pause';
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        elseif ~isempty(prev_s)
                            prev_sel = find(all_Prev(1:prev_s)==0);
                            if isempty(prev_sel)
                                prev = 1; % all correct strides
                            else
                                sel = find(prev_sel<i);
                                if isempty(sel)
                                    prev = 1;
                                else
                                    prev = prev_sel(sel(end))+1;
                                end
                            end
                            allEvents = [stride_list.ICEvents];
                            ICnow = allEvents(1:2:end);
                            IC2now = allEvents(2:2:end);
                            diffNow = [IC2now(prev+1:prev_s)-IC2now(prev:prev_s-1)]/fs;
                            if isempty(find(diffNow>max_break,1))
                                WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(prev_s).ICEvents(1,2), size(prev:prev_s,2)];
                                TerminationReason(size(WB,1)).WB = 'Pause';
                            else
                                checkAll = [stride_list.deltaTime];
                                checkNow = checkAll(prev:prev_s);
                                stridesNow = prev:prev_s;
                                stridesToKeep = stridesNow(checkNow);
                                if size(stridesToKeep,2)>= 2 && isempty(find(diff(stridesToKeep)>1, 1))
                                    WB = [WB; stride_list(min(stridesToKeep)).ICEvents(1,1), stride_list(max(stridesToKeep)).ICEvents(1,2), size(stridesToKeep,2)];
                                    TerminationReason(size(WB,1)).WB = 'Pause';
                                    
                                elseif size(stridesToKeep,2)>= 2
                                    posToDivide = find(diff(stridesToKeep)>1);
                                    for k = 1:length(posToDivide)
                                        if k == 1
                                            stridesToKeepNow = stridesNow(1:posToDivide(k)+1);
                                        else
                                            stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                        end
                                        % first portion
                                        if size(stridesToKeepNow,2)>= 2
                                            WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                            TerminationReason(size(WB,1)).WB = 'Pause';
                                        end
                                        % second portion
                                        if length(posToDivide) ==1
                                            stridesToKeepNow = stridesNow(posToDivide+2:end);
                                        else
                                            stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                        end
                                        if size(stridesToKeepNow,2)>= 2
                                            WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                            TerminationReason(size(WB,1)).WB = 'Pause';
                                        end
                                    end
                                end
                            end
                        end
                    else
                        
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
                    allEvents = [stride_list.ICEvents];
                    ICnow = allEvents(1:2:end);
                    IC2now = allEvents(2:2:end);
                    diffNow = [IC2now(prev+1:now)-IC2now(prev:now-1)]/fs;
                    if isempty(find(diffNow>max_break,1))
                        WB = [WB; stride_list(prev).ICEvents(1,1), stride_list(now).ICEvents(1,2), size(prev:now,2)];
                        if size(WB,1)>1
                            if WB(end,1)==WB(end-1,1)
                                WB(end,:)=[];
                            else
                                TerminationReason(size(WB,1)).WB = 'Pause';
                            end
                        else
                            TerminationReason(size(WB,1)).WB = 'Pause';
                        end
                    else
                        checkAll = [stride_list.deltaTime];
                        checkNow = checkAll(prev:now);
                        stridesNow = prev:now;
                        stridesToKeep = stridesNow(checkNow);
                        if size(stridesToKeep,2)>= 2 && isempty(find(diff(stridesToKeep)>1, 1))
                            WB = [WB; stride_list(min(stridesToKeep)).ICEvents(1,1), stride_list(max(stridesToKeep)).ICEvents(1,2), size(stridesToKeep,2)];
                            TerminationReason(size(WB,1)).WB = 'Pause';
                            
                        elseif size(stridesToKeep,2)>= 2
                            posToDivide = find(diff(stridesToKeep)>1);
                            for k = 1:length(posToDivide)
                                if k == 1
                                    stridesToKeepNow = stridesNow(1:posToDivide(k)+1);
                                else
                                    stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                end
                                % first portion
                                if size(stridesToKeepNow,2)>= 2
                                    WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                    TerminationReason(size(WB,1)).WB = 'Pause';
                                end
                                % second portion
                                if length(posToDivide) ==1
                                    stridesToKeepNow = stridesNow(posToDivide+2:end);
                                else
                                    stridesToKeepNow = stridesNow(posToDivide(k)+2:posToDivide(k)+3);
                                end
                                if size(stridesToKeepNow,2)>= 2
                                    WB = [WB; stride_list(min(stridesToKeepNow)).ICEvents(1,1), stride_list(max(stridesToKeepNow)).ICEvents(1,2), size(stridesToKeepNow,2)];
                                    TerminationReason(size(WB,1)).WB = 'Pause';
                                end
                            end
                        end
                    end
                end
            end
            i = i+1;
        end
        
        % check if there are duplicated rows
        WB_temp = WB;
        WB_temp(:,end+1) = ones(size(WB,1),1);
        for k = 1:size(WB,1)-1
            if isequal(WB(k,1:3),WB(k+1,1:3))
                WB_temp(k+1,end) = 0;
            end
        end
        WB = WB(WB_temp(:,end)==1,:);
        
        WB_temp = WB;
        WB_temp(:,end+1) = ones(size(WB,1),1);
        for i = 1:size(WB,1)
            if WB(i,3) < n_min
                WB_temp(i,end) = 0;
            end
        end
        WB = WB(WB_temp(:,end)==1,:);
    else
        
    end
else
end
end