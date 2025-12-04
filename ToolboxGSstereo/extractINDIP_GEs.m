function GE = extractINDIP_GEs(INDIP)

check = false;
if size(INDIP,2) > 0
    side = {'Right','Left'};
    side_2 = {'right','left'};
    for wb = 1:size(INDIP,2)
        ICs = INDIP(wb).InitialContact_Event';
        ICs_LR = INDIP(wb).InitialContact_LeftRight';
        
        if ~isempty(find(isnan(ICs),1))
            check =  true;
        end
        
        FCs = INDIP(wb).FinalContact_Event';
        FCs_LR = INDIP(wb).FinalContact_LeftRight';
        if wb ==1
            for s = 1:length(side)
                GE.HS.(side_2{s}) = ICs(strcmp(ICs_LR,side{s}));
            end
            for s = 1:length(side)
                GE.TO.(side_2{s}) = FCs(strcmp(FCs_LR,side{s}));
            end
            
        else
            for s = 1:length(side)
                GE.HS.(side_2{s}) = [GE.HS.(side_2{s}); ICs(strcmp(ICs_LR,side{s}))];
            end
            for s = 1:length(side)
                GE.TO.(side_2{s}) = [GE.TO.(side_2{s}); FCs(strcmp(FCs_LR,side{s}))];
            end
        end        
    end
    if check
        GE = [];
    end
else
    GE = [];
end
end