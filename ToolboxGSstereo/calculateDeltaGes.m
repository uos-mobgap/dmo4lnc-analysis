function deltaGE = calculateDeltaGes(GEs, RefGEs, fs, twindow)
side = fieldnames(GEs.HS);
for s = 1:size(side,1)
    IC = GEs.HS.(side{s});
    IC_Ref = round(RefGEs.HS.(side{s})*fs);
    deltaEventHS = nan(size(IC,1),3);
    for e = 1:size(IC,1)
        deltaEv_now = IC(e,1)-IC_Ref;
        [~, posToConsider] = min(abs(deltaEv_now));
        if abs(deltaEv_now(posToConsider))<twindow
            deltaEventHS(e) = deltaEv_now(posToConsider);
            deltaEventHS(e,2) = IC(e,1);
            deltaEventHS(e,3) = IC_Ref(posToConsider);
        end
    end
    if length(unique(deltaEventHS(~isnan(deltaEventHS(:,3)),3))) ~= length(deltaEventHS(~isnan(deltaEventHS(:,3)),3))
        % an event has been associated twice!
        for t = 1:size(deltaEventHS,1)-1
            if (deltaEventHS(t+1,3)-deltaEventHS(t,3))== 0
                if abs(deltaEventHS(t+1,1))< abs(deltaEventHS(t,1))
                    deltaEventHS(t,1) = nan;
                    deltaEventHS(t,3) = nan;
                else
                    deltaEventHS(t+1,1) = nan;
                    deltaEventHS(t+1,3) = nan;
                end
                
            else
            end
        end
    end
    deltaGE.HS.(side{s}) = deltaEventHS(:,1);
    
    FC = GEs.TO.(side{s});
    FC_Ref = round(RefGEs.TO.(side{s})*fs);
    deltaEventTO = nan(size(FC,1),3);
    for e = 1:size(FC,1)
        deltaEv_now = FC(e,1)-FC_Ref;
        [~, posToConsider] = min(abs(deltaEv_now));
        if abs(deltaEv_now(posToConsider))<twindow
            deltaEventTO(e) = deltaEv_now(posToConsider);
            deltaEventTO(e,2) = FC(e,1);
            deltaEventTO(e,3) = FC_Ref(posToConsider);
        end
    end
    if length(unique(deltaEventTO(~isnan(deltaEventTO(:,3)),3))) ~= length(deltaEventTO(~isnan(deltaEventTO(:,3)),3))
        % an event has been associated twice!
        for t = 1:size(deltaEventTO,1)-1
            if (deltaEventTO(t+1,3)-deltaEventTO(t,3))== 0
                if abs(deltaEventTO(t+1,1))< abs(deltaEventTO(t,1))
                    deltaEventTO(t,1) = nan;
                    deltaEventTO(t,3) = nan;
                else
                    deltaEventTO(t+1,1) = nan;
                    deltaEventTO(t+1,3) = nan;
                end
                
            else
            end
        end
    end
    deltaGE.TO.(side{s}) = deltaEventTO(:,1);
end
end