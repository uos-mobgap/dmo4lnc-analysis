function InclineWalkingOutputs = InclineWalkingDet(allStridesStruct, fs)
    
    allH = [allStridesStruct.MaxHeightFlag]';
    allSides = [allStridesStruct.side]';
    allHValues = [allStridesStruct.VerticalDispl]';
    TrueStrides = [allStridesStruct.TrueStrideFlag]';
    ElegibleStrides = zeros(size(allH));
    for i = 1:size(allH,1)
        if ~isnan(allH(i))
            if allH(i)~=0 && allStridesStruct(i).TrueStrideFlag
                ElegibleStrides(i) = i;
            end
        else
            ElegibleStrides(i) = allH(i);
        end
    end
    check = ~isempty(find(ElegibleStrides, 1));
    if check
        Incline_Start = [];
        Incline_End = [];
        Incline_Duration = [];
        Incline_NumberStrides = [];
        Incline_PositiveElevation = [];
        Incline_NegativeElevation = [];
        
        i = 1; i_now = i;
        while i_now<=length(ElegibleStrides) 
            while(ElegibleStrides(i)~=0) && i<=size(ElegibleStrides,1)-1
                i = i+1;
            end
            if i - i_now+1 >= 4
                if i - i_now+1 == 4
                    check = string(allSides(i_now:i));
                    selStrides = TrueStrides(i_now:i);
                    check1 = find(contains(check,'R')); check1R = selStrides(check1);
                    check2 = find(contains(check,'L')); check1L = selStrides(check2);
                    if size(check1,1)==2 && size(check2,1)==2 && isequal(check1R, ones(2,1)) && isequal(check1L, ones(2,1))
                        Incline_Start = [Incline_Start; allStridesStruct(i_now).ICEvents(1,1)/fs];
                        Incline_End = [Incline_End; allStridesStruct(i-1).ICEvents(1,2)/fs];
                        Incline_Duration = [Incline_Duration; Incline_End(end)-Incline_Start(end)];
                        Incline_NumberStrides = [Incline_NumberStrides; i - i_now];
                        allHValues_now = allHValues(i_now:i);
                        Incline_PositiveElevation = [Incline_PositiveElevation; sum(allHValues_now(allHValues_now>0))];
                        Incline_NegativeElevation = [Incline_NegativeElevation; sum(allHValues_now(allHValues_now<0))];
                        
                    elseif size(check1,1)>2 && (isempty(find(isnan(check1R), 1))  && isempty(find(check1R == 0, 1)))
                        Incline_Start = [Incline_Start; allStridesStruct(i_now).ICEvents(1,1)/fs];
                        Incline_End = [Incline_End; allStridesStruct(i-1).ICEvents(1,2)/fs];
                        Incline_Duration = [Incline_Duration; Incline_End(end)-Incline_Start(end)];
                        Incline_NumberStrides = [Incline_NumberStrides; i - i_now];
                        allHValues_now = allHValues(i_now:i);
                        Incline_PositiveElevation = [Incline_PositiveElevation; sum(allHValues_now(allHValues_now>0))];
                        Incline_NegativeElevation = [Incline_NegativeElevation; sum(allHValues_now(allHValues_now<0))];
                        
                    elseif size(check2,1)>2 && (isempty(find(isnan(check1L), 1))  && isempty(find(check1L == 0, 1)))
                        Incline_Start = [Incline_Start; allStridesStruct(i_now).ICEvents(1,1)/fs];
                        Incline_End = [Incline_End; allStridesStruct(i-1).ICEvents(1,2)/fs];
                        Incline_Duration = [Incline_Duration; Incline_End(end)-Incline_Start(end)];
                        Incline_NumberStrides = [Incline_NumberStrides; i - i_now];
                        allHValues_now = allHValues(i_now:i);
                        Incline_PositiveElevation = [Incline_PositiveElevation; sum(allHValues_now(allHValues_now>0))];
                        Incline_NegativeElevation = [Incline_NegativeElevation; sum(allHValues_now(allHValues_now<0))];
                    end
                else
                    Incline_Start = [Incline_Start; allStridesStruct(i_now).ICEvents(1,1)/fs];
                    Incline_End = [Incline_End; allStridesStruct(i-1).ICEvents(1,2)/fs];
                    Incline_Duration = [Incline_Duration; Incline_End(end)-Incline_Start(end)];
                    Incline_NumberStrides = [Incline_NumberStrides; i - i_now];
                    allHValues_now = allHValues(i_now:i);
                    Incline_PositiveElevation = [Incline_PositiveElevation; sum(allHValues_now(allHValues_now>0))];
                    Incline_NegativeElevation = [Incline_NegativeElevation; sum(allHValues_now(allHValues_now<0))];
                end
            end
            i_now = i+1;
            i = i+1;
        end
        InclineWalkingOutputs.Incline_Start = Incline_Start;
        InclineWalkingOutputs.Incline_End = Incline_End;
        InclineWalkingOutputs.Incline_Duration = Incline_Duration;
        InclineWalkingOutputs.Incline_Number = size(Incline_Start,1);
        InclineWalkingOutputs.Incline_NumberStrides = Incline_NumberStrides;
        InclineWalkingOutputs.Incline_PositiveElevation = Incline_PositiveElevation;
        InclineWalkingOutputs.Incline_NegativeElevation = Incline_NegativeElevation;
        
    else
        InclineWalkingOutputs.Incline_Start = [];
        InclineWalkingOutputs.Incline_End = [];
        InclineWalkingOutputs.Incline_Duration = [];
        InclineWalkingOutputs.Incline_Number = [];
        InclineWalkingOutputs.Incline_NumberStrides = [];
        InclineWalkingOutputs.Incline_PositiveElevation = [];
        InclineWalkingOutputs.Incline_NegativeElevation = [];
    end
end