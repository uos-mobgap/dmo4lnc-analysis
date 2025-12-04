function [selStride, TableStride] = strideSelection(IC, strideToRem, L, H, FC, reason, CL, strideD_all, max_h, min_st, max_st, min_sl, side)
    selStride = [];
    TableStride = [];
    r = 1;
    
    for i = 1:size(IC,1)-1     
        
        %% save TableStride info
        TableStride(i).ICEvents = [IC(i,1), IC(i+1,1)];
        TableStride(i).side = side;
        if IC(i,2)==IC(i+1,2) && IC(i,2) == 1
            TableStride(i).ICEvents_Meth = 'Zeni';
        elseif IC(i,2)==IC(i+1,2) && IC(i,2) == 0
            TableStride(i).ICEvents_Meth = 'V_Tr';
        else
            TableStride(i).ICEvents_Meth = 'Both';
        end
                
        TableStride(i).Clearance = CL(i);                
        
        %% ~~~~~~
        if find(strideToRem == i)
            TableStride(i).Length = nan;
            TableStride(i).VerticalDispl = nan;
            TableStride(i).Duration = nan;
            DurNow = strideD_all(i);
            HNow = H(i);
            LNow = L(i);
            TableStride(i).Velocity = TableStride(i).Length/TableStride(i).Duration;
            TableStride(i).MaxHeightFlag = abs(HNow) > max_h;
            TableStride(i).MinLengthFlag = LNow < min_sl;
            TableStride(i).MaxTimeFlag = DurNow > max_st; 
            TableStride(i).MinTimeFlag = DurNow < min_st;
            TableStride(i).TrueStrideFlag = ~find(strideToRem == i);
            
%             if isnan(reason(r,1))
%                 r = r+1;
                % strides that should have been removed - indicated with a 0 in
                % the 6th column and NaNs in the lenght and H                                  
                possibleTO = find(FC(:,1)-IC(i,1)>0 & FC(:,1)-IC(i+1,1)<0); % Find relevant TO
                if isempty(possibleTO)
                    selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, 0,0, CL(i)]; % TO not identified                 
                    
                    TableStride(i).FCEvents = 0;
                    TableStride(i).FCEvents_Meth = '';                  
                    %% ~~
                else
                    if length(possibleTO)>1
                        pos = find(FC(possibleTO,2));
                        if size(pos,1) == 1
                            selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, FC(possibleTO(pos),1),0, CL(i)];
                            TableStride(i).FCEvents = FC(possibleTO(pos),1);
                            TableStride(i).FCEvents_Meth = 'Zeni';
                            
                        elseif size(pos,1) > 1 % more than one FC identified with the Zeni's algo
                            strideD = IC(i+1,1)-IC(i,1);
                            dur = []; stance = [];
                            for k = 1:size(pos,1)
                                dur = [dur; FC(possibleTO(pos(k)),1)-IC(i,1)];
                                stance = [stance; (dur(end)*100)/strideD];
                            end
                            posN = find(stance>50);
                            if length(posN)==1 
                                selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, FC(possibleTO(posN),1),0, CL(i)];
                                
                                TableStride(i).FCEvents = FC(possibleTO(posN),1);
                                TableStride(i).FCEvents_Meth = 'Zeni';
                                
                            elseif isempty(posN)
                                [~,posN] = max(stance);
                                selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, FC(possibleTO(posN),1),0, CL(i)];
                                TableStride(i).FCEvents = FC(possibleTO(posN),1);
                                TableStride(i).FCEvents_Meth = 'Zeni';
                            else
                                selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, 0, 0, CL(i)]; % not sure which is the correct FC!
                                
                                TableStride(i).FCEvents = 0;
                                TableStride(i).FCEvents_Meth = 'Zeni';
                            end

                        elseif isempty(pos) % no FC identified with the Zeni's algo
                            strideD = IC(i+1,1)-IC(i,1);
                            dur = []; stance = [];
                            for k = 1:size(possibleTO,1)
                                dur = [dur; FC(possibleTO(k),1)-IC(i,1)];
                                stance = [stance; (dur(end)*100)/strideD];
                            end
                            posN = find(stance>50);
                            if length(posN)==1 
                                selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, FC(possibleTO(posN),1),0, CL(i)];
                                
                                TableStride(i).FCEvents = FC(possibleTO(posN),1);
                                TableStride(i).FCEvents_Meth = 'V_Tr';
                            elseif isempty(posN)
                                [~,posN] = max(stance);
                                selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, FC(possibleTO(posN),1),0, CL(i)];
                                TableStride(i).FCEvents = FC(possibleTO(posN),1);
                                TableStride(i).FCEvents_Meth = 'V_Tr';
                            else
                                selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, 0,0, CL(i)]; % not sure which is the correct FC!
                                TableStride(i).FCEvents = 0;
                                TableStride(i).FCEvents_Meth = 'V_Tr';
                            end
                        end
                    else
                        selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, FC(possibleTO,1),0, CL(i)];
                        TableStride(i).FCEvents = FC(possibleTO,1);
                        TableStride(i).FCEvents_Meth = 'Zeni';
                    end
                end
%             else
%                 r = r+1;
%             end
        %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
        else    
            % Correct strides
            TableStride(i).Length = L(i);
            TableStride(i).VerticalDispl = H(i);
            TableStride(i).Duration = strideD_all(i);
            TableStride(i).Velocity = TableStride(i).Length/TableStride(i).Duration;
            
            TableStride(i).MinLengthFlag = TableStride(i).Length < min_sl;
            TableStride(i).MaxHeightFlag = abs(TableStride(i).VerticalDispl) > max_h;            
            TableStride(i).TrueStrideFlag = TableStride(i).Length > min_sl;
            TableStride(i).MaxTimeFlag = TableStride(i).Duration > max_st; 
            TableStride(i).MinTimeFlag = TableStride(i).Duration < min_st;
                       
            possibleTO = find(FC(:,1)-IC(i,1)>0 & FC(:,1)-IC(i+1,1)<0); % Find relevant TO
            if isempty(possibleTO)
                selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), 0, 1, CL(i)]; % TO not identified
                
                TableStride(i).FCEvents = 0;
                TableStride(i).FCEvents_Meth = '';
            else
                if length(possibleTO)>1
                    pos = find(FC(possibleTO,2));
                    if size(pos,1) == 1
                        selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), FC(possibleTO(pos),1), 1, CL(i)];
                        
                        TableStride(i).FCEvents = FC(possibleTO(pos),1);
                        TableStride(i).FCEvents_Meth = 'Zeni';
                        
                    elseif size(pos,1) > 1 % more than one FC identified with the Zeni's algo
                        strideD = IC(i+1,1)-IC(i,1);
                        dur = []; stance = [];
                        for k = 1:size(pos,1)
                            dur = [dur; FC(possibleTO(pos(k)),1)-IC(i,1)];
                            stance = [stance; (dur(end)*100)/strideD];
                        end
                        posN = find(stance>50);
                        if length(posN)==1 
                            selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), FC(possibleTO(posN),1), 1, CL(i)];
                            
                            TableStride(i).FCEvents = FC(possibleTO(posN),1);
                            TableStride(i).FCEvents_Meth = 'Zeni';
                        elseif isempty(posN)
                            [~,posN] = max(stance);
                            selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), FC(possibleTO(posN),1), 1, CL(i)];
                            
                            TableStride(i).FCEvents = FC(possibleTO(posN),1);
                            TableStride(i).FCEvents_Meth = 'Zeni';
                        else
                            selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), 0, 1, CL(i)]; % not sure which is the correct FC!
                            
                            TableStride(i).FCEvents = 0;
                            TableStride(i).FCEvents_Meth = 'Zeni';
                        end
                        
                    elseif isempty(pos) % no FC identified with the Zeni's algo
                        strideD = IC(i+1,1)-IC(i,1);
                        dur = []; stance = [];
                        for k = 1:size(possibleTO,1)
                            dur = [dur; FC(possibleTO(k),1)-IC(i,1)];
                            stance = [stance; (dur(end)*100)/strideD];
                        end
                        posN = find(stance>50);
                        if length(posN)==1 
                            selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), FC(possibleTO(posN),1), 1, CL(i)];
                            
                            TableStride(i).FCEvents = FC(possibleTO(posN),1);
                            TableStride(i).FCEvents_Meth = 'V_Tr';
                        elseif isempty(posN)
                            [~,posN] = max(stance);
                            selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), FC(possibleTO(posN),1), 1, CL(i)];
                            TableStride(i).FCEvents = FC(possibleTO(posN),1);
                            TableStride(i).FCEvents_Meth = 'V_Tr';
                        else
                            selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), 0, 1, CL(i)]; % not sure which is the correct FC!
                            TableStride(i).FCEvents = 0;
                            TableStride(i).FCEvents_Meth = 'V_Tr';
                        end
                    end
                else
                    selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), FC(possibleTO,1), 1, CL(i)];
                    TableStride(i).FCEvents = FC(possibleTO,1);
                    TableStride(i).FCEvents_Meth = 'Zeni';
                end
            end
        end
    end
end