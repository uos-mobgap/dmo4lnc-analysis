function [allStrides, WB_b, WB_e, StrNowR, StrNowL] = CheckSpuriousStride(COM, fs, allStrides, StrNowR, StrNowL)

if ~isempty(allStrides)
    trh = 0.3;
    v_COM = diff(COM)/(1/fs);
    for i = 1:size(v_COM,1)
        vPelvis(i,1) = norm(COM(i,:));
    end
    removePoints = ~isnan(COM(:,1));
    rangeV = range(vPelvis(removePoints(1:end-1)));
    
    for i = 1:size(allStrides,1)
        if ~isnan(allStrides(i,1))
            if (allStrides(i,1)-10)>0 && (allStrides(i,1)+10)<length(vPelvis)
                v_now = mean(vPelvis(allStrides(i,1)-10:allStrides(i,1)+10));
            elseif (allStrides(i,1)-10)>0
                v_now = mean(vPelvis(allStrides(i,1)-10:end));
            elseif (allStrides(i,1)+10)<length(vPelvis)
                v_now = mean(vPelvis(1:allStrides(i,1)+10));
            else
                v_now = mean(vPelvis);
            end
            
            if v_now < trh*rangeV
                allStrides(i,9) = 1;
            else
                allStrides(i,9) = 0;
            end
        else
            allStrides(i,9) = 1;
        end
    end
    
    %% GEs to CHECK
    if isempty(find(isnan(allStrides(:,1)), 1))
        mdeltaT = mean(diff(allStrides(:,1)));
    else
        dt = nan(size(allStrides,1),1);
        for i = 1:size(allStrides,1)-1
            if ~isnan(allStrides(i,1)) &&  ~isnan(allStrides(i+1,1))
                dt(i) = allStrides(i+1,1)-allStrides(i,1);
            end
        end
        mdeltaT = mean(dt(~isnan(dt)));
    end
    i = 1;
    check = find(allStrides(:,9));
    if ~isempty(check)
        if size(check,1)==size(allStrides,1)
            elementsToRemove = false;
        else
            elementsToRemove = true;
        end
    else
        elementsToRemove = false;
    end
    if elementsToRemove
        while i < size(allStrides,1)-1
            if allStrides(i+1,9)== 1 && allStrides(i,9) == 1
                % two consecutive strides that could be spurious movements of both feet
                % movement of the feet while seating
                if ~isnan(allStrides(i+1,1))&&~isnan(allStrides(i,1))
                    if allStrides(i+1,1)- allStrides(i,1)< trh*mdeltaT
                        % these strides have to be removed
                        ALL_IC1l = [StrNowL.ICEvents]; ALL_IC1l = ALL_IC1l(1:2:end);
                        ALL_IC1r = [StrNowR.ICEvents]; ALL_IC1r = ALL_IC1r(1:2:end);
                        if allStrides(i,7) == 0
                            pos_toRem = find(ALL_IC1l == allStrides(i,1),1);
                            StrNowL(pos_toRem) = [];
                            pos_toRem = find(ALL_IC1r == allStrides(i+1,1),1);
                            StrNowR(pos_toRem) = [];
                            allStrides(i:i+1,:)=[];
                        else
                            pos_toRem = find(ALL_IC1l == allStrides(i+1,1),1);
                            StrNowL(pos_toRem) = [];
                            pos_toRem = find(ALL_IC1r == allStrides(i,1),1);
                            StrNowR(pos_toRem) = [];
                            allStrides(i:i+1,:)=[];
                        end
                    end
                else
                    ALL_IC1l = [StrNowL.ICEvents]; ALL_IC1l = ALL_IC1l(1:2:end);
                    ALL_IC1r = [StrNowR.ICEvents]; ALL_IC1r = ALL_IC1r(1:2:end);
                    if allStrides(i,7) == 0
                        pos_toRem = find(ALL_IC1l == allStrides(i,1),1);
                        StrNowL(pos_toRem) = [];
                        pos_toRem = find(ALL_IC1r == allStrides(i+1,1),1);
                        StrNowR(pos_toRem) = [];
                        allStrides(i:i+1,:)=[];
                    else
                        pos_toRem = find(ALL_IC1l == allStrides(i+1,1),1);
                        StrNowL(pos_toRem) = [];
                        pos_toRem = find(ALL_IC1r == allStrides(i,1),1);
                        StrNowR(pos_toRem) = [];
                        allStrides(i:i+1,:)=[];
                    end
                end
            end
            i = i+1;
        end
        
        while isnan(allStrides(1,1)) && size(allStrides,1)>=1
            allStrides(1,:) = [];
        end
    end
    
    if ~isempty(find(diff(allStrides(:,1))<10, 1))
        % there are still spurious movements
        i = 1;
        while i<size(allStrides,1)-1
            if allStrides(i+1,1)-allStrides(i,1)<10
                ALL_IC1l = [StrNowL.ICEvents]; ALL_IC1l = ALL_IC1l(1:2:end);
                ALL_IC1r = [StrNowR.ICEvents]; ALL_IC1r = ALL_IC1r(1:2:end);
                if allStrides(i,7) == 0
                    pos_toRem = find(ALL_IC1l == allStrides(i,1),1);
                    StrNowL(pos_toRem) = [];
                    pos_toRem = find(ALL_IC1r == allStrides(i+1,1),1);
                    StrNowR(pos_toRem) = [];
                    allStrides(i:i+1,:)=[];
                else
                    pos_toRem = find(ALL_IC1l == allStrides(i+1,1),1);
                    StrNowL(pos_toRem) = [];
                    pos_toRem = find(ALL_IC1r == allStrides(i,1),1);
                    StrNowR(pos_toRem) = [];
                    allStrides(i:i+1,:)=[];
                end
            end
            i = i+1;
        end
        
        while isnan(allStrides(1,1)) && size(allStrides,1)>=1
            allStrides(1,:) = [];
        end
    end
    
    WB_b = min(allStrides(:,1));
    WB_e = max(allStrides(:,2));
    allStrides(:,9) = [];
    allStrides(:,6) = 1:size(allStrides,1);
else
    WB_b = [];
    WB_e = [];
end
end