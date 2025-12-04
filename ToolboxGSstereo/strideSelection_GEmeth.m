function [selStride, TableStride] = strideSelection_GEmeth(IC, strideToRem, L, H, FC, CL, strideD_all, max_h, min_st, max_st, min_sl, side)
selStride = [];
TableStride = [];

for i = 1:size(IC,1)-1
    
    %% save TableStride info
    TableStride(i).ICEvents = [IC(i,1), IC(i+1,1)];
    TableStride(i).side = side;
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
        TableStride(i).MaxHeightFlag = HNow > max_h;
        TableStride(i).MinLengthFlag = LNow < min_sl;
        TableStride(i).MaxTimeFlag = DurNow > max_st;
        TableStride(i).MinTimeFlag = DurNow < min_st;
        TableStride(i).TrueStrideFlag = ~find(strideToRem == i);
        
        if ~isempty(FC)
            possibleTO = find(FC(:,1)-IC(i,1)>0 & FC(:,1)-IC(i+1,1)<0); % Find relevant TO
            if isempty(possibleTO)
                selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, 0,0, CL(i)]; % TO not identified
                
                TableStride(i).FCEvents = 0;
                
            else
                if length(possibleTO)>1
                    TableStride(i).FCEvents = FC(possibleTO);
                else
                    selStride = [selStride; IC(i,1), IC(i+1,1), nan, nan, FC(possibleTO,1),0, CL(i)];
                    TableStride(i).FCEvents = FC(possibleTO,1);
                end
            end
            
        else
            TableStride(i).FCEvents = 0;
        end
    else
        % Correct strides
        TableStride(i).Length = L(i);
        TableStride(i).VerticalDispl = H(i);
        TableStride(i).Duration = strideD_all(i);
        TableStride(i).Velocity = TableStride(i).Length/TableStride(i).Duration;
        
        TableStride(i).MinLengthFlag = TableStride(i).Length < min_sl;
        TableStride(i).MaxHeightFlag = TableStride(i).VerticalDispl > max_h;
        TableStride(i).TrueStrideFlag = TableStride(i).Length > min_sl;
        TableStride(i).MaxTimeFlag = TableStride(i).Duration > max_st;
        TableStride(i).MinTimeFlag = TableStride(i).Duration < min_st;
        
        if ~isempty(FC)
            possibleTO = find(FC(:,1)-IC(i,1)>0 & FC(:,1)-IC(i+1,1)<0); % Find relevant TO
            if isempty(possibleTO)
                selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), 0, 1, CL(i)]; % TO not identified                
                TableStride(i).FCEvents = 0;
            else
                if length(possibleTO)>1
                    TableStride(i).FCEvents = FC(possibleTO);
                else
                    selStride = [selStride; IC(i,1), IC(i+1,1), L(i,1), H(i,1), FC(possibleTO,1), 1, CL(i)];
                    TableStride(i).FCEvents = FC(possibleTO,1);
                end
            end
        else
            TableStride(i).FCEvents = 0;
        end
    end
end
end