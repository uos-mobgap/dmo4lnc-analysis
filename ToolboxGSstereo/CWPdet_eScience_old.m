function [Standards, WB_list] = CWPdet_eScience(Data, stride_list, n_min, max_break,  rMrkL, lMrkL, headerTraj, fs, ...
    TurnM, TurnDur, TurnsNotDec, Euler_Angles, maxTurn, minDurTurn, averVel)
    Standards = [];
    WB_list = [];
    AngVel = diff(Euler_Angles)/(1/fs);

    %% Idenitfy quick change of directions
    posTurn = find(abs(TurnM)> maxTurn);
    if ~isempty(posTurn)
        if ~isempty(minDurTurn)
            quick_Turn_pos  = find((TurnDur(posTurn,2) - TurnDur(posTurn,1))/fs <minDurTurn);
        else
            meanAngV = zeros(size(posTurn));
            for i = 1:size(posTurn,1)
                meanAngV(i) = mean(AngVel(TurnDur(posTurn(i),1):TurnDur(posTurn(i),2)));
            end
            quick_Turn_pos  = find(abs(meanAngV)> averVel);
        end
    else
        quick_Turn_pos = [];
    end

    %% WB based on correct strides
    [WB_R, stride_list.R, TerminationReason_R] = createLR_WB_v1(stride_list.R, TurnDur, posTurn, quick_Turn_pos, n_min);
    [WB_L, stride_list.L, TerminationReason_L] = createLR_WB_v1(stride_list.L, TurnDur, posTurn, quick_Turn_pos, n_min);

     %% Find intersections between WB_R and WB_L (and also WB_R with no intersection)
     FlagL = zeros(1,size(WB_L,1));
     if ~isempty(WB_R)
         for i = 1:size(WB_R,1)
             CandidateWBs(i).IndexR = i;
             CandidateWBs(i).IndexL = [];
         end
         for i = 1:size(WB_R,1)
             for j = 1:size(WB_L,1)
                 if (~isempty(intersect(WB_R(i,1):1:WB_R(i,2),WB_L(j,1):1:WB_L(j,2))))
                     CandidateWBs(i).IndexL = [CandidateWBs(i).IndexL j];
                     FlagL(j)=1;
                 else
                 end
             end
         end
     else
         CandidateWBs=[];
     end
 
     %% Save also WB_L with no intersections referring to FlagL
     if ~isempty(WB_L)
         k=length(CandidateWBs)+1;
         for i=1:length(FlagL)
             if FlagL(i)==0
                 %this is a WB_L with zero intersections
                 CandidateWBs(k).IndexR= [];
                 CandidateWBs(k).IndexL= i;
                 k = k+1;
             else
             end
         end
     end
 
     %% Check if there are overlapping WBs
     if ~isempty(CandidateWBs)
         CandidateWBs_temp = CandidateWBs;
         if ~isempty([CandidateWBs.IndexL])
             if size (unique(unique([CandidateWBs.IndexL],'rows')),2) ~= size([CandidateWBs.IndexL],2)  % check for the L Index   
                 WBs_int = [CandidateWBs_temp.IndexL];
                 pos_Nochanges = find(diff(WBs_int)==0);
                 rowTorem = min(pos_Nochanges):1:max(pos_Nochanges)+1;
                 CandidateWBs_temp(min(pos_Nochanges)).IndexR = [CandidateWBs(rowTorem).IndexR];
                 CandidateWBs_temp(min(pos_Nochanges)+1:max(pos_Nochanges)+1)=[];
             end
         end
        CandidateWBs = CandidateWBs_temp; clear ('CandidateWBs_temp');
     end

    %% WB (R-L-R-L-.... or L-R-L-R-....)
    for i = 1:size(CandidateWBs,2) % Number of possible WBs
        breaks = 0;
        %Select the relevant strides from the LIST
        if ~isempty(CandidateWBs(i).IndexR) && ~isempty(CandidateWBs(i).IndexL)
            caseWB = 1;
            
            % R side
            r = CandidateWBs(i).IndexR;
            begR = 1;
            while isempty(find(stride_list.R(1,begR).ICEvents(1,1) == WB_R(min(r),1), 1))
                begR = begR+1;
            end            
            finR = 1;
            while isempty(find(stride_list.R(1,finR).ICEvents(1,2) == WB_R(max(r),2), 1))
                finR = finR+1;
            end            
            StrNowR = stride_list.R(1,begR:finR);
            
            % L side
            l = CandidateWBs(i).IndexL;
            begL = 1;
            while isempty(find(stride_list.L(1,begL).ICEvents(1,1) == WB_L(min(l),1), 1))
                begL = begL+1;
            end            
            finL = 1;
            while isempty(find(stride_list.L(1,finL).ICEvents(1,2) == WB_L(max(l),2), 1))
                finL = finL+1;
            end            
            StrNowL = stride_list.L(1,begL:finL);

            [allStrides, ~, allStridesStruct] = orderSelStrides_V1(StrNowR,StrNowL);

        elseif ~isempty(CandidateWBs(i).IndexR) % only R WB detected
            caseWB = 2;
            r = CandidateWBs(i).IndexR;
            begR = find(stride_list.R(:,1) == WB_R(min(r),1));
            finR = find(stride_list.R(:,2) == WB_R(max(r),2));
            StrNowR = stride_list.R(begR:finR,:);

            SelStrR = StrNowR;
            SelStrL = [];
            [allStrides, ~, allStridesStruct] = orderSelStrides_V1(SelStrR,SelStrL);

        else % only L WB detected
            caseWB = 3;
            l = CandidateWBs(i).IndexL;
            begL = find(stride_list.L(:,1) == WB_L(min(l),1));
            finL = find(stride_list.L(:,2) == WB_L(max(l),2));
            StrNowL = stride_list.L(begL:finL,:);

            SelStrR = [];
            SelStrL = StrNowL;
            [allStrides, ~, allStridesStruct] = orderSelStrides_V1(SelStrR,SelStrL);
        end

        WB_b = allStrides(1,1);
        WB_e = max(allStrides(end-1:end,2));

        %% L and R strides have to be consecutive
        allStrides_temp = allStrides;
        allStridesStruct_temp = allStridesStruct;
        c = 0;
        for k = 2:size(allStrides)
            if allStrides(k-1,end) == allStrides(k,end) % a stride on the other side has been either missed or removed
                % check that both strides are correct strides
                if allStrides(k-1,6) == 1 && allStrides(k,6) == 1
                    if allStrides(k,end) == 1
                        side = 0; side_str = 'L';
                    else
                        side = 1; side_str = 'R';
                    end
                    allStrides_temp(k+c,:) = [nan(1,7),side];
                    fields = fieldnames(allStridesStruct_temp);
                    for fn = 1:size(fields,1)
                        if strcmp(fields{fn},'ICEvents')
                            allStridesStruct_temp(k+c).(fields{fn})= [nan, nan];
                        elseif strcmp(fields{fn},'side')
                            allStridesStruct_temp(k+c).(fields{fn})= side_str;
                        else
                            allStridesStruct_temp(k+c).(fields{fn})= nan;
                        end
                    end
                    allStrides_temp(k+1+c:k+c+size(allStrides(k:end,:),1),:) = allStrides(k:end,:);
                    for fn = k:size(allStridesStruct,2)
                        allStridesStruct_temp(1+c+fn) = allStridesStruct(fn);
                    end
                    c = c+1;
                end
            end
        end
        allStrides = allStrides_temp; clear 'allStrides_temp';
        allStridesStruct = allStridesStruct_temp; clear 'allStridesStruct_temp';
        temp = allStrides(:,6);
        allStrides(:,6) = [];
        allStrides(:,6) = 1:size(allStrides,1);
        allStrides(:,end+1) = temp; clear temp; % now last column identify if whether a stride was a correct one or not

        %% PAUSES-MISSING EVENTS
        r_strides = allStrides(allStrides(:,7)==1 & allStrides(:,end)==1,:);
        l_strides = allStrides(allStrides(:,7)==0 & allStrides(:,end)==1,:);

        if size(r_strides,1)>=2 && size(l_strides,1)>=2        
            %% identification of pauses
            r_pause = []; l_pause = [];
            for k = 2:size(r_strides,1)
                if r_strides(k-1,2)~=r_strides(k,1)
                    if  (r_strides(k,1)- r_strides(k-1,2))/fs<= max_break
                        r_pause = [r_pause; r_strides(k-1,2), r_strides(k,1)];
                    end
                end
            end

            for k = 2:size(l_strides,1)
                if l_strides(k-1,2)~=l_strides(k,1)
                    if  (l_strides(k,1)- l_strides(k-1,2))/fs <= max_break
                        l_pause = [l_pause; l_strides(k-1,2), l_strides(k,1)];
                    end
                end
            end

            % find overlaping pauses - overlapping pauses = 2, R pauses = 1, L pauses = 0
            rl_pause = []; %start, stop, side (0,1,2)
            for k = 1:size(r_pause,1)
                for j = 1:size(l_pause,1)
                    if ~isempty(intersect(r_pause(k,1):1:r_pause(k,2),l_pause(j,1):1:l_pause(j,2)))
                        rl_pause = [rl_pause; max([r_pause(k,1),l_pause(j,1)]), min([r_pause(k,2),l_pause(j,2)]), 2];
                    else
                        rl_pause = [rl_pause; r_pause(k,1) r_pause(k,2) 1];
                        rl_pause = [rl_pause; l_pause(j,1) l_pause(j,2) 0];
                    end
                end
                if isempty(l_pause)
                    rl_pause = [rl_pause; r_pause(k,1) r_pause(k,2) 1];
                end
            end
            if isempty(r_pause) && ~isempty(l_pause)
                rl_pause = [rl_pause; l_pause(:,1) l_pause(:,2) ones(size(l_pause,1),1)];
            end

            if ~isempty(rl_pause)
                pause_begins = rl_pause(:,1);
                [~,pos] = sort(pause_begins);
                rl_pause = rl_pause(pos,:);
            end

        %% Check reason why WB ends
        r = CandidateWBs(i).IndexR;
        l = CandidateWBs(i).IndexR;
        if strcmp(TerminationReason_L(l).WB,TerminationReason_R(r).WB)
            if strcmp(TerminationReason_L(l).WB,'SharpTurn') || strcmp(TerminationReason_L(r).WB,'InclineWalking')
                % the last ICs do NOT have to be removed
                if strcmp(TerminationReason_L(i).WB,'SharpTurn')
                    WB_list.(strcat('WB',num2str(i))).ReasonEnd = 'Turn';
                else
                    WB_list.(strcat('WB',num2str(i))).ReasonEnd = 'InclineWalking';
                end
            elseif strcmp(TerminationReason_L(i).WB,'Pause') 
                WB_list.(strcat('WB',num2str(i))).ReasonEnd = 'Pause';
                % the last stride has to be removed - final IC is not a
                % real IC for another stride   
                
                % THIS IS NOT REMOVED FOR CWP
%                 [allStrides, WB_e] = removeLastStride(allStrides);  
            end
        else
            t = 1
        end

        if size(allStrides,1)>=4          
            if caseWB == 1
                HS_r = [StrNowR.ICEvents]';
                HS_r = unique(HS_r); HS_r = HS_r(HS_r<=WB_e); HS_r = HS_r(HS_r>=WB_b);
                HS_l = [StrNowL.ICEvents]';
                HS_l = unique(HS_l); HS_l = HS_l(HS_l<=WB_e); HS_l = HS_l(HS_l>=WB_b);

            elseif caseWB == 2
                HS_r = [StrNowR.ICEvents]';
                HS_r = unique(HS_r); HS_r = HS_r(HS_r<=WB_e); HS_r = HS_r(HS_r>=WB_b);

            else
                HS_l = [StrNowL.ICEvents]';
                HS_l = unique(HS_l); HS_l = HS_l(HS_l<=WB_e); HS_l = HS_l(HS_l>=WB_b);
            end

            % FIGURE for the specific WB
            if (WB_e+100) < length(Data.(headerTraj).l_heel) && (WB_b -100)>0
                timeNow = 1/fs:1/fs:length(Data.(headerTraj).l_heel(WB_b -100:WB_e+100,3))/fs;
            elseif (WB_b -100)<0 && (WB_e+100) < length(Data.(headerTraj).l_heel)
                timeNow = 1/fs:1/fs:length(Data.(headerTraj).l_heel(1:WB_e+100,3))/fs;
            elseif (WB_b -100)<0
                timeNow = 1/fs:1/fs:length(Data.(headerTraj).l_heel(1:end,3))/fs;
            else
                timeNow = 1/fs:1/fs:length(Data.(headerTraj).l_heel(WB_b -100:end,3))/fs;
            end

            %% Centre of mass
            COM = (Data.(headerTraj).r_psis + Data.(headerTraj).l_psis)/2;

            figure('Name',strcat('CWP_', num2str(i)),'NumberTitle','off')
            subplot(211)
            plot(COM(:,1), COM(:,2))
            hold on
            plot(COM(WB_b:WB_e,1), COM(WB_b:WB_e,2),'Color','m','LineWidth',3)            

            %% SAVE WBs details
            WB_list.(strcat('WB',num2str(i))).Details = [WB_b, WB_e];
            WB_list.(strcat('WB',num2str(i))).Duration = (WB_e - WB_b)/fs;
            WB_list.(strcat('WB',num2str(i))).nPauses = size(rl_pause,1);

            if caseWB == 1
                if size(allStrides(allStrides(:,7)==1),1) == size(allStrides(allStrides(:,7)==0),1)
                    if find(isnan(allStrides(:,3)))
                        OverallDistance = [];
                    else
                        stridesLNow = sum(allStrides(allStrides(:,7)==0,3));
                        stridesRNow = sum(allStrides(allStrides(:,7)==1,3));
                        OverallDistance =  mean([stridesLNow,stridesRNow]);
                    end
                elseif size(allStrides(allStrides(:,7)==0),1)> size(allStrides(allStrides(:,7)==1),1)
                    % more L strides
                    stridesLNow = allStrides(allStrides(:,7)==0,3);
                    if find(isnan(stridesLNow))
                        OverallDistance = [];
                    else
                        OverallDistance =  sum(stridesLNow);
                    end
                else
                    % more R strides
                    stridesRNow = allStrides(allStrides(:,7)==1,3);
                    if find(isnan(stridesRNow))
                        OverallDistance = [];
                    else
                        OverallDistance =  sum(stridesRNow);
                    end
                end
                OverallDistance_COM = travDist(sort([HS_r; HS_l]), COM);

            elseif caseWB == 2
                OverallDistance =  sum(SelStrR(:,3));
                OverallDistance_COM = travDist(sort(HS_r), COM);
            else % caseWB == 3
                OverallDistance =  sum(SelStrL(:,3));
                OverallDistance_COM = travDist(sort(HS_l), COM);
            end
            v = OverallDistance/WB_list.(strcat('WB',num2str(i))).Duration;

            WB_list.(strcat('WB',num2str(i))).TravDist = OverallDistance;
            WB_list.(strcat('WB',num2str(i))).TravDist_COM = OverallDistance_COM;
            WB_list.(strcat('WB',num2str(i))).v = v;
            WB_list.(strcat('WB',num2str(i))).v_COM = OverallDistance_COM/WB_list.(strcat('WB',num2str(i))).Duration;

            % STRIDES dx and sx
            WB_list.(strcat('WB',num2str(i))).nStrides = size(allStrides,1);
            WB_list.(strcat('WB',num2str(i))).nRStrides = size(allStrides(allStrides(:,7)==1),1);
            WB_list.(strcat('WB',num2str(i))).nLStrides = size(allStrides(allStrides(:,7)==0),1);
            WB_list.(strcat('WB',num2str(i))).RStrides = allStrides(allStrides(:,7)==1,1:5); % [beginning, end, length, H, TO]
            WB_list.(strcat('WB',num2str(i))).LStrides = allStrides(allStrides(:,7)==0,1:5); % [beginning, end, length, H, TO]
            WB_list.(strcat('WB',num2str(i))).allStrides = allStrides;

            All_strides_dur = ((allStrides(:,2)-allStrides(:,1))/fs)';
            pos = find(allStrides(:,end)==0);
            if ~isempty(pos)
                All_strides_dur(pos) = nan;
            end
            All_strides_L = allStrides(:,3)';
            All_strides_H = allStrides(:,4)';
            All_strides_v = All_strides_L./All_strides_dur;

            [lv_foot, lvCOM, lSTRelevChange, orderL] = strideCounter(allStrides(allStrides(:,7)==0,:), COM, fs);
            [rv_foot, rvCOM, rSTRelevChange, orderR] = strideCounter(allStrides(allStrides(:,7)==1,:), COM, fs);

            % order ALL events
            all_v_foot = [lv_foot,orderL;rv_foot,orderR]; [~,pos] = sort(all_v_foot(:,2));
            all_v_foot = all_v_foot(pos);
            all_v_COM = [lvCOM;rvCOM]; all_v_COM = all_v_COM(pos);
            all_elevChange_COM = [lSTRelevChange;rSTRelevChange]; all_elevChange_COM = all_elevChange_COM(pos);

            % save info
            WB_list.(strcat('WB',num2str(i))).vRStrides = rv_foot;
            WB_list.(strcat('WB',num2str(i))).vLStrides = lv_foot;
            WB_list.(strcat('WB',num2str(i))).vAllStrides = all_v_foot;

            WB_list.(strcat('WB',num2str(i))).vR_COM_Strides = rvCOM;
            WB_list.(strcat('WB',num2str(i))).vL_COM_Strides = lvCOM;
            WB_list.(strcat('WB',num2str(i))).vAll_COM_Strides = all_v_COM;

            WB_list.(strcat('WB',num2str(i))).COM_rElevChange = rSTRelevChange;
            WB_list.(strcat('WB',num2str(i))).COM_lElevChange = lSTRelevChange;
            WB_list.(strcat('WB',num2str(i))).COM_AllElevChange = all_elevChange_COM;
            clear ('lv_foot', 'rv', 'lvCOM', 'rvCOM', 'orderL', 'orderR',...
                'pos','all_v_COM','all_v_COM','all_elevChange_COM');

            clear('StrNowR', 'StrNowL', 'SelStrR','SelStrL');

            %% Left/Right Steps
            % Left Step: From rHS to lHS
            % Right Step: From lHS to rHS
            [HS,IA] = unique([allStrides(:,1); allStrides(:,2)]);
            check = 1;
            if find(isnan(HS),1) % there are NaNs
                check = 0;
                HS_temp1 = [allStrides(:,1), allStrides(:,7); allStrides(end-1:end,2), allStrides(end-1:end,7)];     
                HS_temp1(:,3) = ones(size(HS_temp1,1),1);
                pos = find(isnan(HS_temp1(:,1)));
                pos_s = find(max(pos)>size(allStrides,1),1);
                if isempty(pos_s)
                    for s = 1:size(pos,1)
                        sideNow = allStrides(pos(s),7);
                        allStridesSideNow = allStrides((allStrides(1:pos(s),7)==sideNow),2:6);
                        if size(allStridesSideNow,1)>1
                            HS_temp1(pos(s),1) = allStridesSideNow(end-1,1);
                            HS_temp1(pos(s),3) = 0;
                        else
                            HS_temp1(pos(s),3) = 0;
                        end
                    end
                else
                    for s = 1:size(pos,1)
                        if pos(s)<= size(allStrides,1)
                            sideNow = allStrides(pos(s),7);
                            allStridesSideNow = allStrides((allStrides(1:pos(s),7)==sideNow),2:6);
                            if size(allStridesSideNow,1)>1
                                HS_temp1(pos(s),1) = allStridesSideNow(end-1,1);
                                HS_temp1(pos(s),3) = 0;
                            else
                                HS_temp1(pos(s),3) = 0;
                            end
                        else
                            pos_corr = size(allStrides,1)- pos(s)+size(allStrides,1);
                            sideNow = allStrides(pos_corr,7);
                            allStridesSideNow = allStrides((allStrides(1:pos_corr,7)==sideNow),2:6);
                            if size(allStridesSideNow,1)>1
                                HS_temp1(pos(s),1) = allStridesSideNow(end-1,1);
                                HS_temp1(pos(s),3) = 0;
                            else
                                HS_temp1(pos(s),3) = 0;
                            end
                        end
                    end
                end

                %% check
                pos_zeros = find(HS_temp1(:,3)==0);
                for z = 1:size(pos_zeros,1)
                    if HS_temp1(pos_zeros(z),1)>HS_temp1(pos_zeros(z)-1,1) && HS_temp1(pos_zeros(z),1)<HS_temp1(pos_zeros(z)+1,1)
                    %OK
                    else
                        HS_temp1(pos_zeros(z),1) = nan;
                    end
                end

                %% CHECK ORDER!
                for o = 2:size(HS_temp1,1)
                    if ~isnan(HS_temp1(o,1)) && ~isnan(HS_temp1(o-1,1))
                        if HS_temp1(o-1,1)> HS_temp1(o,1)
                            temp = HS_temp1(o-1,:);
                            HS_temp1(o-1,:) = HS_temp1(o,:);
                            HS_temp1(o,:) = temp;
                        end
                    end
                end
                HS = HS_temp1(:,1:2); clear HS_temp1;
            else
                temp = [allStrides(:,7); allStrides(:,7)];
                HS(:,2) = temp(IA);

                %% CHECK ORDER!
                for o = 2:size(HS,1)
                    if ~isnan(HS(o,1)) && ~isnan(HS(o-1,1))
                        if HS(o-1,1)> HS(o,1)
                            temp = HS(o-1,:);
                            HS(o-1,:) = HS(o,:);
                            HS(o,:) = temp;
                        end
                    end
                end
                clear temp;
            end

            %% correct step lenght can be assessed only where the markers on the pelvis are without gaps - add a FLAG
            HS(:,3) = ones(size(HS,1),1);
            for k = 1:size(TurnsNotDec,1)
                [~,d] = intersect(HS(:,1), TurnsNotDec(k,1):1:TurnsNotDec(k,2));
                HS(d,3) = zeros(1,size(d,1));
            end

            %% ADD a FLAG for those HS that belongs to real strides
            if check ==1
                temp = [allStrides(:,end); allStrides(:,end)];
                HS(:,4) = temp(IA); clear temp;
            else
                for k = 1:size(allStrides,1)
                    if ~isnan(HS(k,1))
                        pos = find(allStrides(:,1)== HS(k,1),1);
                        if ~isempty(pos)
                            HS(k,4) = allStrides(pos,end);
                        else
                            pos = find(allStrides(:,2)== HS(k,1),1);
                            HS(k,4) = allStrides(pos,end);
                        end
                    else
                        HS(k,4) = 0;
                    end
                end
                if ~isnan(allStrides(end-1:end,end))
                    HS(end-1:end,4) = allStrides(end-1:end,end);
                else
                    pos = find(isnan(allStrides(end-1:end,end)));
                    n_pos = size(allStrides,1);
                    allStrides(n_pos-2+pos,end)=zeros(size(pos));
                    HS(end-1:end,4) = allStrides(n_pos-2+pos,end);
                end
            end

            TO = [allStrides(:,5),allStrides(:,7)];
            TO(TO(:,1) == 0) = nan;

            if caseWB ==1 % both L and R GEs are needed

                RH = Data.(headerTraj).r_heel;
                LH = Data.(headerTraj).l_heel;
                [n_steps, dur_steps, Leng_steps, L_COM_steps, v_COM_steps, elevChangeCOM_steps] = stepIdentification(HS, COM, fs, lMrkL, rMrkL, LH, RH);

                WB_list.(strcat('WB',num2str(i))).Steps.n.Left = length(find(dur_steps(:,2)==0));
                WB_list.(strcat('WB',num2str(i))).Steps.n.Right = length(find(dur_steps(:,2)==1));
                WB_list.(strcat('WB',num2str(i))).Steps.n.Overall = n_steps;
                WB_list.(strcat('WB',num2str(i))).Steps.dur.Left = dur_steps(dur_steps(:,2)==0,1);
                WB_list.(strcat('WB',num2str(i))).Steps.dur.Right = dur_steps(dur_steps(:,2)==1,1);
                WB_list.(strcat('WB',num2str(i))).Steps.d.Left = Leng_steps(Leng_steps(:,2)==0,1);
                WB_list.(strcat('WB',num2str(i))).Steps.d.Right = Leng_steps(Leng_steps(:,2)==1,1);
                WB_list.(strcat('WB',num2str(i))).Steps.dCOM.Left = L_COM_steps(L_COM_steps(:,2)==0,1);
                WB_list.(strcat('WB',num2str(i))).Steps.dCOM.Right = L_COM_steps(L_COM_steps(:,2)==1,1);
                WB_list.(strcat('WB',num2str(i))).Steps.elevCOM.Left = elevChangeCOM_steps(elevChangeCOM_steps(:,2)==0,1);
                WB_list.(strcat('WB',num2str(i))).Steps.elevCOM.Right = elevChangeCOM_steps(elevChangeCOM_steps(:,2)==1,1);
                WB_list.(strcat('WB',num2str(i))).Steps.v_COM.Left = v_COM_steps(dur_steps(:,2)==0,1);
                WB_list.(strcat('WB',num2str(i))).Steps.v_COM.Right = v_COM_steps(dur_steps(:,2)==1,1);

                WB_list.(strcat('WB',num2str(i))).Cadence = (60*WB_list.(strcat('WB',num2str(i))).Steps.n.Overall)/WB_list.(strcat('WB',num2str(i))).Duration;
                clear ('ln','rn','ldur','rdur','lLeng','rLeng','lv','rv',...
                    'relevChange','lelevChange','relevChangeCOM','lelevChangeCOM');
            else

                WB_list.(strcat('WB',num2str(i))).Steps = [];
                WB_list.(strcat('WB',num2str(i))).Cadence = [];
            end

            %% Turns
            Turn = [];
            ManTurn = [];
            ManTurnWB = [];
            TurnDur(:,3) = ones(size(TurnDur,1),1);

            %% detect turns identified where the signal had gaps - to be removed - add FLAG 1 == turn correct, otherwise FLAG == 0
            for k = 1:size(TurnsNotDec,1)
                for l = 1:size(TurnDur)
                    [~,d] = intersect(TurnDur(l,1):1:TurnDur(l,2), TurnsNotDec(k,1):1:TurnsNotDec(k,2));
                    if ~isempty(d)
                        TurnDur(l,3)= 0;
                    else
                        TurnDur(l,3)= 1;
                    end
                end
            end

            if ~isempty(TurnDur)
                posB = find(TurnDur(:,1)<WB_e);
                posE = find(TurnDur(:,2)>WB_b);
                [~, kk, q] = unique([posB; posE],'first');
                indexToDupes = find(not(ismember(1:numel([posB; posE]),kk)));
                SelRow = q(indexToDupes);
                if ~isempty(SelRow)
                    SelTurn = TurnDur(SelRow(1):SelRow(end),:);
                    SelTurnM = TurnM(SelRow(1):SelRow(end));

                    for j = 1:size(SelTurn,1)
                        if SelTurn(j,1)>WB_b && SelTurn(j,2)<WB_e && SelTurn(j,3) == 1
                            %% Turn within the WB
                            Turn = [Turn; SelTurn(j,1), SelTurn(j,2), (SelTurn(j,2)-SelTurn(j,1))/fs, SelTurnM(j)];
                            plot(COM(SelTurn(j,1):SelTurn(j,2),1), COM(SelTurn(j,1):SelTurn(j,2),2),'b:','LineWidth',2)
                            %% Check for Maneuvering turns (turns not completed in a WB) at the beginning of the WB
                        elseif SelTurn(j,1) <WB_b  && SelTurn(j,3) == 1
                            ManTurn = [ManTurn; SelTurn(j,1), SelTurn(j,2), (SelTurn(j,2)-SelTurn(j,1))/fs, SelTurnM(j)];
                            MagNow = Euler_Angles(SelTurn(j,2))-Euler_Angles(WB_b);
                            ManTurnWB = [ManTurnWB; WB_b, SelTurn(j,2), (SelTurn(j,2)-WB_b)/fs, MagNow];
                            plot(COM(SelTurn(j,1):SelTurn(j,2),1), COM(SelTurn(j,1):SelTurn(j,2),2),'c:','LineWidth',2)
                            %% Check for Maneuvering turns (turns not completed in a WB) at the end of the WB
                        elseif SelTurn(j,3) == 1
                            ManTurn = [ManTurn; SelTurn(j,1), SelTurn(j,2), (SelTurn(j,2)-SelTurn(j,1))/fs, SelTurnM(j)];
                            MagNow = Euler_Angles(WB_b)-Euler_Angles(SelTurn(j,1));
                            ManTurnWB = [ManTurnWB; SelTurn(j,1), WB_b, (WB_b-SelTurn(j,1))/fs, MagNow];
                            plot(COM(SelTurn(j,1):SelTurn(j,2),1), COM(SelTurn(j,1):SelTurn(j,2),2),'c:','LineWidth',2)
                        end
                    end
                end
            end

            % save info
            WB_list.(strcat('WB',num2str(i))).Turn = Turn;
            WB_list.(strcat('WB',num2str(i))).MTurn.Turn = ManTurn;
            WB_list.(strcat('WB',num2str(i))).MTurn.TurnWB = ManTurnWB;
            
            legend({'CoM',strcat('CoM - CWP', num2str(i))},'Location','northwest')
            legend('boxoff'); xlabel('CoM x traj [m]'); ylabel('CoM y traj [m]')            
            clear Turn; clear ManTurn; clear ManTurnWB;

            %% PLOT
            subplot(212)

            if (WB_e+100) < length(Data.(headerTraj).l_heel) && (WB_b -100)>0
                plot(timeNow, Data.(headerTraj).l_heel(WB_b -100:WB_e+100,3), 'g')
            elseif (WB_b -100)<0 && (WB_e+100) < length(Data.(headerTraj).l_heel)
                plot(timeNow, Data.(headerTraj).l_heel(1:WB_e+100,3), 'g')
            elseif (WB_b -100)<0
                plot(timeNow, Data.(headerTraj).l_heel(1:end,3), 'g')
            else
                plot(timeNow, Data.(headerTraj).l_heel(WB_b -100:end,3), 'g')
            end

            hold on
            if (WB_e+100) < length(Data.(headerTraj).l_heel)&& (WB_b -100)>0
                plot(timeNow, Data.(headerTraj).r_heel(WB_b -100:WB_e+100,3), 'r')
            elseif (WB_b -100)<0 && (WB_e+100) < length(Data.(headerTraj).l_heel)
                plot(timeNow, Data.(headerTraj).r_heel(1:WB_e+100,3), 'r')
            elseif (WB_b -100)<0
                plot(timeNow, Data.(headerTraj).r_heel(1:end,3), 'r')
            else
                plot(timeNow, Data.(headerTraj).r_heel(WB_b -100:end,3), 'r')
            end
            xlabel('Time [s]'); ylabel('[m]')

            if caseWB ==1
                if (WB_e+100) < length(Data.(headerTraj).l_heel)&& (WB_b -100)>0
                    for j = 1:length(HS_r)
                        line([(HS_r(j)-WB_b+100)/fs, (HS_r(j)-WB_b+100)/fs],[0 0.3],'Color','r','LineWidth',1, 'LineStyle','--')
                    end
                    for j = 1:length(HS_l)
                        line([(HS_l(j)-WB_b+100)/fs, (HS_l(j)-WB_b+100)/fs],[0 0.3],'Color','g','LineWidth',1, 'LineStyle','--')
                    end
                elseif (WB_b -100)<0
                    for j = 1:length(HS_r)
                        line([(HS_r(j))/fs, (HS_r(j))/fs],[0 0.3],'Color','r','LineWidth',1, 'LineStyle','--')
                    end
                    for j = 1:length(HS_l)
                        line([(HS_l(j))/fs, (HS_l(j))/fs],[0 0.3],'Color','g','LineWidth',1, 'LineStyle','--')
                    end
                else
                    for j = 1:length(HS_r)
                        line([(HS_r(j)-WB_b+100)/fs, (HS_r(j)-WB_b+100)/fs],[0 0.3],'Color','r','LineWidth',1, 'LineStyle','--')
                    end
                    for j = 1:length(HS_l)
                        line([(HS_l(j)-WB_b+100)/fs, (HS_l(j)-WB_b+100)/fs],[0 0.3],'Color','g','LineWidth',1, 'LineStyle','--')
                    end
                end
            elseif caseWB == 2
                if (WB_e+100) < length(Data.(headerTraj).l_heel)&& (WB_b -100)>0
                    for j = 1:length(HS_r)
                        line([(HS_r(j)-WB_b+100)/fs, (HS_r(j)-WB_b+100)/fs],[0 0.3],'Color','r','LineWidth',1, 'LineStyle','--')
                    end
                elseif (WB_b -100)<0
                    for j = 1:length(HS_r)
                        line([(HS_r(j))/fs, (HS_r(j))/fs],[0 0.3],'Color','r','LineWidth',1, 'LineStyle','--')
                    end
                else
                    for j = 1:length(HS_r)
                        line([(HS_r(j)-WB_b+100)/fs, (HS_r(j)-WB_b+100)/fs],[0 0.3],'Color','r','LineWidth',1, 'LineStyle','--')
                    end
                end
            else
                if (WB_e+100) < length(Data.(headerTraj).l_heel)&& (WB_b -100)>0
                    for j = 1:length(HS_l)
                        line([(HS_l(j)-WB_b+100)/fs, (HS_l(j)-WB_b+100)/fs],[0 0.3],'Color','g','LineWidth',1, 'LineStyle','--')
                    end
                elseif (WB_b -100)<0

                    for j = 1:length(HS_l)
                        line([(HS_l(j))/fs, (HS_l(j))/fs],[0 0.3],'Color','g','LineWidth',1, 'LineStyle','--')
                    end
                else
                    for j = 1:length(HS_l)
                        line([(HS_l(j)-WB_b+100)/fs, (HS_l(j)-WB_b+100)/fs],[0 0.3],'Color','g','LineWidth',1, 'LineStyle','--')
                    end
                end
            end
            
            legend({'L heel','R heel'},'Location','northwest')
            legend('boxoff')

            %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %% save relevant outputs - STANDARD
            %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Standards(i).Start = WB_b/fs;
            Standards(i).End  = WB_e/fs;

            Standards(i).Cadence = WB_list.(strcat('WB',num2str(i))).Cadence;
            Standards(i).Duration = WB_list.(strcat('WB',num2str(i))).Duration;
            if isnan( WB_list.(strcat('WB',num2str(i))).TravDist_COM)
                Standards(i).Length = [];
            else
                Standards(i).Length = WB_list.(strcat('WB',num2str(i))).TravDist_COM;
            end
            Standards(i).Length_Feet = WB_list.(strcat('WB',num2str(i))).TravDist;
            if isnan(WB_list.(strcat('WB',num2str(i))).v_COM)
                Standards(i).Speed = [];
            else
                Standards(i).Speed = WB_list.(strcat('WB',num2str(i))).v_COM;
            end
            
            Standards(i).NumberStrides = WB_list.(strcat('WB',num2str(i))).nStrides;
            Standards(i).TerminationReason = WB_list.(strcat('WB',num2str(i))).ReasonEnd;
            
            %% only pauses recognised on both sides are considered
            RL_pause = rl_pause;
            if ~isempty(RL_pause)
                pos = find(RL_pause(:,3)~=2);
                RL_pause(pos,:) = [];
                Standards(i).StartBreaks = RL_pause(:,1)./fs;
                Standards(i).EndBreaks = RL_pause(:,2)./fs;
                Standards(i).NumberBreaks = size(RL_pause,1);
                Standards(i).DurationBreaks = ((RL_pause(:,2)-RL_pause(:,1))./fs)';
            else
                Standards(i).StartBreaks = [];
                Standards(i).EndBreaks = [];
                Standards(i).NumberBreaks = [];
                Standards(i).DurationBreaks = [];
            end

            if ~isempty(WB_list.(strcat('WB',num2str(i))).Turn)
                Standards(i).Turning_Flag = 1;
                Standards(i).Turning_Start = WB_list.(strcat('WB',num2str(i))).Turn(:,1)/fs';
                Standards(i).Turning_End = WB_list.(strcat('WB',num2str(i))).Turn(:,2)/fs';
                Standards(i).Turning_NumberTurns = size(WB_list.(strcat('WB',num2str(i))).Turn,1);
                Standards(i).Turning_AngleTurns = WB_list.(strcat('WB',num2str(i))).Turn(:,4)';
                Standards(i).Turning_DurationTurns = WB_list.(strcat('WB',num2str(i))).Turn(:,3)';
                AngVel_Peak_all = zeros(1,Standards(i).Turning_NumberTurns);
                AngVel_Mean_all = zeros(1,Standards(i).Turning_NumberTurns);
                AngVel_all = nan(Standards(i).Turning_NumberTurns,round(max(Standards(i).Turning_DurationTurns)*fs));
                n_strides_Turn = zeros(Standards(i).Turning_NumberTurns,1);
                for j = 1:Standards(i).Turning_NumberTurns
                    n_now = round((WB_list.(strcat('WB',num2str(i))).Turn(j,3)*fs)+1);
                    AngVel_all(j,1:n_now) = AngVel(WB_list.(strcat('WB',num2str(i))).Turn(j,1):...
                        WB_list.(strcat('WB',num2str(i))).Turn(j,2));
                    AngVel_Peak_all(1,j) = max(AngVel_all(j,1:n_now));
                    AngVel_Mean_all(1,j) = mean(AngVel_all(j,1:n_now));
                    nIC1 = find(allStrides(:,1)> WB_list.(strcat('WB',num2str(i))).Turn(j,1),1);
                    nIC2 = find(allStrides(:,2)> WB_list.(strcat('WB',num2str(i))).Turn(j,1),1);
                    if ~isempty(nIC1) && ~isempty(nIC2)
                        n_min = min(allStrides(nIC1,6), allStrides(nIC2,6));
                    elseif ~isempty(nIC1)
                        n_min = allStrides(nIC1,6);
                    else
                        n_min = allStrides(nIC2,6);
                    end
                    
                    nIC1 = find(allStrides(:,1)< WB_list.(strcat('WB',num2str(i))).Turn(j,2)); nIC1 = nIC1(end);
                    nIC2 = find(allStrides(:,2)< WB_list.(strcat('WB',num2str(i))).Turn(j,2)); nIC2 = nIC2(end);
                    if ~isempty(nIC1) && ~isempty(nIC2)
                        n_max = max(allStrides(nIC1,6), allStrides(nIC2,6));
                    elseif ~isempty(nIC1)
                        n_max = allStrides(nIC1,6);
                    else
                        n_max = allStrides(nIC2,6);
                    end
                    n_strides_Turn(j,1) = n_max - n_min+1;
                end
                Standards(i).Turning_NumberStrides = n_strides_Turn;
                Standards(i).Turning_AngularVelocity = AngVel_all;
                Standards(i).Turning_PeakAngularVelocity = AngVel_Peak_all;
                Standards(i).Turning_MeanAngularVelocity = AngVel_Mean_all;
                Standards(i).Turning_SharpTurn_Flag = zeros(size(Standards(i).Turning_Start));
            else
                Standards(i).Turning_Flag = 0;
                Standards(i).Turning_Start = [];
                Standards(i).Turning_End = [];
                Standards(i).Turning_DurationTurns = [];
                Standards(i).Turning_NumberTurns = [];
                Standards(i).Turning_AngleTurns = [];   
                Standards(i).Turning_NumberStrides = [];
                Standards(i).Turning_AngularVelocity = [];
                Standards(i).Turning_PeakAngularVelocity = [];
                Standards(i).Turning_MeanAngularVelocity = [];
                Standards(i).Turning_SharpTurn_Flag = [];
            end
            
            %% INCLINE WALKING
            InclineWalkingOutputs = InclineWalkingDet(allStridesStruct, fs);
            Standards(i).Incline_Start = InclineWalkingOutputs.Incline_Start;
            Standards(i).Incline_End = InclineWalkingOutputs.Incline_End;
            Standards(i).Incline_Duration = InclineWalkingOutputs.Incline_Duration;
            Standards(i).Incline_Number = InclineWalkingOutputs.Incline_Number;
            Standards(i).Incline_NumberStrides = InclineWalkingOutputs.Incline_NumberStrides;
            Standards(i).Incline_PositiveElevation = InclineWalkingOutputs.Incline_PositiveElevation;
            Standards(i).Incline_NegativeElevation = InclineWalkingOutputs.Incline_NegativeElevation;

            %% STRIDE INFO
            Standards(i).Stride_Duration = All_strides_dur;
            Standards(i).Stride_Length = All_strides_L;
            Standards(i).Stride_Height = All_strides_H;
            Standards(i).Stride_Velocity = All_strides_v;
            Standards(i).Stride_TrunkElevationChange = WB_list.(strcat('WB',num2str(i))).COM_AllElevChange';

            Standards(i).Stride_TrunkRotation = [];

            [ST_d, SW_d, SS_d, DS_d, ST_l, SW_l, ST_v, SW_v] = calculateStrideStanceSwing(allStrides,Data.(headerTraj),fs);
            Standards(i).Stance_Duration = ST_d;
            Standards(i).Swing_Duration = SW_d;
            Standards(i).SingleSupport_Duration = SS_d;
            Standards(i).DoubleSupport_Duration = DS_d;
            Standards(i).Stance_Length = ST_l;
            Standards(i).Swing_Length = SW_l;
            Standards(i).Stance_Velocity = ST_v;
            Standards(i).Swing_Velocity = SW_v;

            Standards(i).SingleSupport_Length = [];
            Standards(i).DoubleSupport_Length = [];

            Standards(i).SingleSupport_Velocity =[];
            Standards(i).DoubleSupport_Velocity = [];

            Standards(i).InitialContacts = allStrides(:,1:2)/fs;
            Standards(i).InitialContact_Event = HS(:,1)'/fs;

            for ev = 1:size(HS,1)
                if HS(ev,2) == 0
                    Standards(i).InitialContact_LeftRight(1,ev) ={'Left'};
                elseif HS(ev,2) == 1
                    Standards(i).InitialContact_LeftRight(1,ev) ={'Right'};
                else
                    Standards(i).InitialContact_LeftRight(1,ev) ={'NaN'};
                end
            end

            Standards(i).FinalContact_Event = TO(:,1)'/fs;
            for ev = 1:size(TO,1)
                if TO(ev,2) == 0
                    Standards(i).FinalContact_LeftRight(1,ev) ={'Left'};
                elseif TO(ev,2) == 1
                    Standards(i).FinalContact_LeftRight(1,ev) ={'Right'};
                else
                    Standards(i).FinalContact_LeftRight(1,ev) ={'NaN'};
                end
            end

            Standards(i).Step_TrunkRotation = [];
            if caseWB == 1
                Standards(i).Step_Duration = dur_steps(:,1)';
                Standards(i).Step_Length = L_COM_steps(:,1)';
                Standards(i).Step_Length_Feet = Leng_steps(:,1)';
                Standards(i).Step_Velocity = v_COM_steps';
                Standards(i).Step_TrunkElevationChange = elevChangeCOM_steps(:,1)';

            else % caseWB == 2 || caseWB == 3
                Standards(i).Step_Duration = [];
                Standards(i).Step_Length = [];
                Standards(i).Step_Velocity =  [];
                Standards(i).Step_TrunkElevationChange = [];
            end
        else
            fprintf('Strides have been removed during sharp turns/WB end conditions - not enough strides for satisfying the WB conditions!\n');
        end

        %% Remove WBs where there are less than 3 strides when only R or L events have been identifies
        if caseWB == 2 || caseWB == 3
            if size(Standards(i).Stride_Length,2)<3
                Standards(i) = [];
            end
        end
        else
            %         Standards(i) = [];
            %         WB_list = [];
            fprintf('Strides have been removed during sharp turns - not enough strides for satisfying the WB conditions!\n');
        end
    end
    if isempty(CandidateWBs)
        %% add empty fields to standards
        fieldsStandards = {'Start';'End';'Cadence';'Duration';...
            'Length';'Length_Feet';'Speed';'StartBreaks';'EndBreaks';...
            'NumberBreaks';'DurationBreaks';'NumberStrides';'Turning_Flag';'TerminationReason';...
            'Turning_Start';'Turning_End';'Turning_NumberTurns';'Turning_AngleTurns';...
            'Turning_DurationTurns';'Turning_NumberStrides';'Turning_AngularVelocity';'Turning_PeakAngularVelocity';...
            'Turning_MeanAngularVelocity';'Turning_SharpTurn_Flag';...
            'Incline_Start';'Incline_End';'Incline_Duration';'Incline_Number';...
            'Incline_NumberStrides'; 'Incline_PositiveElevation'; 'Incline_NegativeElevation';...
            'Stride_Duration';'Stride_Length';'Stride_Height';'Stride_Velocity';'Stride_TrunkElevationChange';...
            'Stride_TrunkRotation';'Stride_InitialContacts';'Stance_Duration';'Swing_Duration';...
            'SingleSupport_Duration';'DoubleSupport_Duration';'Stance_Length';'Swing_Length';'Stance_Velocity';'Swing_Velocity';...
            'SingleSupport_Length';'DoubleSupport_Length';'SingleSupport_Velocity';'DoubleSupport_Velocity';'InitialContact_Event';...
            'InitialContact_LeftRight';'FinalContact_Event';'FinalContact_LeftRight';'Step_TrunkRotation';'Step_Duration';'Step_Length';...
            'Step_Length_Feet';'Step_Velocity';'Step_TrunkElevationChange'};
        for k = 1:size(fieldsStandards,1)
            Standards(1).(fieldsStandards{k}) = [];
        end
        WB_list = [];
        fprintf('WBs not detected, WB conditions not satisfied!\n');
    end

    %% CHECK order of the WBs
    if ~isempty(Standards)
        %% remove empty rows
        check = 0; i = 1;
        while check == 0
            if isempty(Standards(i).Start)
                Standards(i) = [];
            else
                i = i+1;
            end
            if i > size(Standards,2)
                check = 1;
            end
        end
        WB_begins = [Standards.Start];
        [~,pos] = sort(WB_begins);
        Standards = Standards(pos);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%% STEPS
function [n, dur, Leng, L_COM, v, elevChangeCOM] = stepIdentification(HS, COMTraj, fs, lMrkL, rMrkL, LH, RH)

%% INPUTS:
% HS [nx4], n = identified HS, 1st column: L(0)/R(1); 2nd column: step L
% celculation possible(1)/notPossible(0); 3rd column: HS belongs to a
% "proper stride" true(1)/false(1)
%
% COMTraj [Nx3], 3D COM mrk traj in the GRF
%
% fs [Hz], sampling frequency 
%
% l/r MrlL, structure with mrks traj (Heel, Toe, COM and R) represented in the pelvic LRF
%
% L/RH, heel mrk traj in the GRF
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dur = [];
    Leng = [];
    L_COM = [];
    elevChangeCOM = [];

    for j = 2:size(HS,1)       
        if HS(j,2)<HS(j-1,2) && HS(j,4)==1 && HS(j-1,4)==1 % left step
            dur = [dur; (HS(j,1) - HS(j-1,1))/fs 0];
            if HS(j,3)== 1 && HS(j,3)== HS(j-1,3) && (~isnan(HS(j,1))&&~isnan(HS(j-1,1)))           
                % check if sideway walking occurs
                if range(lMrkL.HEEL(HS(j-1,1):HS(j,1),1)) > range(lMrkL.HEEL(HS(j-1,1):HS(j,1),3)) % normal walking                    
                    Rnow = rMrkL.R(:,1,HS(j-1,1):HS(j,1)); Rnow = squeeze(Rnow); 
                    dirwalk = mean(Rnow,2);
                    
                    stepG = [(LH(HS(j,1),1) - RH(HS(j-1,1),1)), (LH(HS(j,1),2) - RH(HS(j-1,1),2)), ...
                        (LH(HS(j,1),3) - RH(HS(j-1,1),3))];
                    
                    step_COM = [(COMTraj(HS(j,1),1) - COMTraj(HS(j-1,1),1)), (COMTraj(HS(j,1),2) - COMTraj(HS(j-1,1),2)),...
                        (COMTraj(HS(j,1),3) - COMTraj(HS(j-1,1),3))];
                    
                    L_COM = [L_COM; dot(step_COM,dirwalk) 0];
                    Leng = [Leng; dot(stepG,dirwalk) 0];
                    
                else % side walking
                    Rnow = rMrkL.R(:,3,HS(j-1,1):HS(j,1)); Rnow = squeeze(Rnow); 
                    dirwalk = mean(Rnow,2);
                    
                    stepG = [(LH(HS(j,1),1) - RH(HS(j-1,1),1)), (LH(HS(j,1),2) - RH(HS(j-1,1),2)), ...
                        (LH(HS(j,1),3) - RH(HS(j-1,1),3))];
                    
                    step_COM = [(COMTraj(HS(j,1),1) - COMTraj(HS(j-1,1),1)), (COMTraj(HS(j,1),2) - COMTraj(HS(j-1,1),2)),...
                        (COMTraj(HS(j,1),3) - COMTraj(HS(j-1,1),3))];
                    
                    L_COM = [L_COM; dot(step_COM,dirwalk) 0];
                    Leng = [Leng; dot(stepG,dirwalk) 0]; 
                end
                         
                elevChangeCOM = [elevChangeCOM; abs(COMTraj(HS(j,1),3) - COMTraj(HS(j-1,1),3)) 0];
            else
                Leng = [Leng; nan 0];
                L_COM = [L_COM; nan 0];
                elevChangeCOM = [elevChangeCOM; nan 0];
            end
            
        elseif HS(j,2)> HS(j-1,2) && HS(j,4)==1 && HS(j-1,4)==1 % right step
            dur = [dur; (HS(j,1) - HS(j-1,1))/fs 1];
            if HS(j,3)== 1 && HS(j,3)== HS(j-1,3) && (~isnan(HS(j,1))&&~isnan(HS(j-1,1)))                
                % check walking direction
                if range(rMrkL.HEEL(HS(j-1,1):HS(j,1),1)) > range(rMrkL.HEEL(HS(j-1,1):HS(j,1),3)) % normal walking 
                    Rnow = rMrkL.R(:,1,HS(j-1,1):HS(j,1)); Rnow = squeeze(Rnow); 
                    dirwalk = mean(Rnow,2);
                    
                    stepG = [(RH(HS(j,1),1) - LH(HS(j-1,1),1)), (RH(HS(j,1),2) - LH(HS(j-1,1),2)), ...
                        (RH(HS(j,1),3) - LH(HS(j-1,1),3))];
                    
                    step_COM = [(COMTraj(HS(j,1),1) - COMTraj(HS(j-1,1),1)), (COMTraj(HS(j,1),2) - COMTraj(HS(j-1,1),2)),...
                        (COMTraj(HS(j,1),3) - COMTraj(HS(j-1,1),3))];
                    
                    L_COM = [L_COM; dot(step_COM,dirwalk) 1];
                    Leng = [Leng; dot(stepG,dirwalk) 1];
                else % side walking
                    Rnow = rMrkL.R(:,3,HS(j-1,1):HS(j,1)); Rnow = squeeze(Rnow); 
                    dirwalk = mean(Rnow,2);
                    
                    stepG = [(RH(HS(j,1),1) - LH(HS(j-1,1),1)), (RH(HS(j,1),2) - LH(HS(j-1,1),2)), ...
                        (RH(HS(j,1),3) - LH(HS(j-1,1),3))];
                    
                    step_COM = [(COMTraj(HS(j,1),1) - COMTraj(HS(j-1,1),1)), (COMTraj(HS(j,1),2) - COMTraj(HS(j-1,1),2)),...
                        (COMTraj(HS(j,1),3) - COMTraj(HS(j-1,1),3))];
                    
                    L_COM = [L_COM; dot(step_COM,dirwalk) 1];
                    Leng = [Leng; dot(stepG,dirwalk) 1];  
                end
                elevChangeCOM = [elevChangeCOM; abs(COMTraj(HS(j,1),3) - COMTraj(HS(j-1,1),3)) 1];
                
            else
                Leng = [Leng; nan 1];
                L_COM = [L_COM; nan 1];
                elevChangeCOM = [elevChangeCOM; nan 1];
            end
            
        else % IC from the same foot or stride that does not fulfill the definition criteria
             dur = [dur; nan 2];
             Leng = [Leng; nan 2];
             L_COM = [L_COM; nan 2];
             elevChangeCOM = [elevChangeCOM; nan 2];
        end     
    end
    n = size(dur,1);
    v = L_COM(:,1)./dur(:,1);
end

%% Stride Counter
function [vFoot, vCOM, elevChange, order] = strideCounter(Stride, COMTraj, fs)
    duration = [];
    LengCOM = [];
    LFoot = [];
    elevChange = [];
    vFoot = [];
    order = [];

    for j = 1:size(Stride,1)
        order = [order; Stride(j,5)];
        if (Stride(j,7)==1)
            duration = [duration; (Stride(j,2)-Stride(j,1))/fs];
            LFoot = [LFoot; Stride(j,3)];
        else
            duration = [duration; nan];
            LFoot = [LFoot; nan];
        end
        
        if ~isnan(Stride(j,1)) && (Stride(j,7)==1)
            LengCOM = [LengCOM; sqrt((COMTraj(Stride(j,2),1) - COMTraj(Stride(j,1),1))^2+...
                (COMTraj(Stride(j,2),2) - COMTraj(Stride(j,1),2))^2)];
            elevChange = [elevChange; abs(COMTraj(Stride(j,2),3) - COMTraj(Stride(j,1),3))];
        else
            LengCOM = [LengCOM; nan];
            elevChange = [elevChange; nan];
        end
    end
    vCOM = LengCOM./duration;
    vFoot = LFoot./duration;
end

%% Traversed distance
function OverallDistance = travDist(HS, COM)
    OverallDistance = 0;
    for i = 2:size(HS)
        OverallDistance = OverallDistance + sqrt((COM(HS(i),1) - COM(HS(i-1),1))^2+...
            (COM(HS(i),2) - COM(HS(i-1),2))^2);
    end
end

%% Stance/Swing
function [ST_d, SW_d, SS_d, DS_d, ST_l, SW_l, ST_v, SW_v] = calculateStrideStanceSwing(Strides, headerTraj, fs)
    ST_d = []; % stance duration
    SW_d = []; % swing duration
    SS_d = []; % single support duration
    DS_d = []; % double support duration
    ST_l = []; % stance lenght
    SW_l = []; % swing lenght
    ST_v = []; % stance velocity
    SW_v = []; % swing velocity

    for k = 1:size(Strides,1)
        if Strides(k,5)>0  && Strides(k,end) == 1%% FC has been identified stance/swing phase can be identified
            % stance/swing duration
            ST_d = [ST_d, (Strides(k,5) - Strides(k,1))/fs];
            SW_d = [SW_d, (Strides(k,2) - Strides(k,5))/fs];

            % stance/swing length
            if Strides(k,7) == 1
                side = 'r';
            else
                side = 'l';
            end
            mrk_now = strcat(side,'_heel');
            ST_l = [ST_l, sqrt((headerTraj.(mrk_now)(Strides(k,5),1) - headerTraj.(mrk_now)(Strides(k,1),1))^2+...
                (headerTraj.(mrk_now)(Strides(k,5),2) - headerTraj.(mrk_now)(Strides(k,1),2))^2)];
            SW_l = [SW_l, sqrt((headerTraj.(mrk_now)(Strides(k,2),1) - headerTraj.(mrk_now)(Strides(k,5),1))^2+...
                (headerTraj.(mrk_now)(Strides(k,2),2) - headerTraj.(mrk_now)(Strides(k,5),2))^2)];

            % stance/swing length, velocity
            ST_v = [ST_v, ST_l(end)/ST_d(end)];
            SW_v = [SW_v, SW_l(end)/SW_d(end)];

            % single/double support - need info from the previous FC
            if k ==1
                SS_d = [SS_d, nan];
                DS_d = [DS_d, nan];
            else
                if Strides(k-1,5)> 0 && Strides(k-1,7)~=Strides(k,7)%previous FC is detected and L and R strides
                    IDS = (Strides(k-1,5) - Strides(k,1))/fs; %initial double support
                    TDS = (Strides(k,5) - Strides(k-1,2))/fs; %terminal double support
                    DS_d =  [DS_d, IDS+TDS];
                    GCD = (Strides(k,2) - Strides(k,1))/fs; % Gait cycle duration
                    SS_d = [SS_d, GCD-DS_d(end)];
                else
                    SS_d = [SS_d, nan];
                    DS_d = [DS_d, nan];
                end
            end

        else %% FC has NOT been identified
            ST_d = [ST_d, nan];
            SW_d = [SW_d, nan];

            SS_d = [SS_d, nan];
            DS_d = [DS_d, nan];

            ST_l = [ST_l, nan];
            SW_l = [SW_l, nan];

            ST_v = [ST_v, nan];
            SW_v = [SW_v, nan];
        end
    end
end