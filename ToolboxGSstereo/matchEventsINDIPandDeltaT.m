function [IC, FC, deltaEventHS, deltaEventTO] = matchEventsINDIPandDeltaT(IC, FC, GE_INDIP, fs, s, twindow)

%% ICs
% Keep only those events that have been found by the INDIP system
IC_r = round(GE_INDIP.HS.(s)*fs);
if ~isempty(IC)
    if ~isempty(find(IC(:,1)>IC_r(1)-twindow,1))
        IC = IC(IC(:,1)>IC_r(1)-twindow,1);
    else
        IC = [];
    end
else
    IC = [];
end
if ~isempty(IC)
    if ~isempty(find(IC(:,1)<IC_r(end)+twindow,1))
        IC = IC(IC(:,1)<IC_r(end)+twindow,1);
    else
        IC = [];
    end
else
    IC = [];
end
deltaEventHS = nan(size(IC,1),3);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for e = 1:size(IC,1)
    deltaEv_now = IC(e,1)-IC_r;
    [~, posFCConsider] = min(abs(deltaEv_now));
    if abs(deltaEv_now(posFCConsider))<twindow
        deltaEventHS(e,1) = deltaEv_now(posFCConsider);
        deltaEventHS(e,2) = IC(e,1);
        deltaEventHS(e,3) = IC_r(posFCConsider);
    end
end
if length(unique(deltaEventHS(~isnan(deltaEventHS(:,3)),3))) ~= length(deltaEventHS(~isnan(deltaEventHS(:,3)),3))
    % an event has been associated twice!
    for t = 1:size(deltaEventHS,1)-1
        if (deltaEventHS(t+1,3)-deltaEventHS(t,3))== 0
            if abs(deltaEventHS(t+1,1))< abs(deltaEventHS(t,1))
                deltaEventHS(t,1) = nan;
                deltaEventHS(t,3) = nan;
            elseif abs(deltaEventHS(t+1,1)) == abs(deltaEventHS(t,1))
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
deltaEventHS = deltaEventHS(:,1);

%% FCs
% Keep only those events that have been found by the INDIP system
FC_r = round(GE_INDIP.TO.(s)*fs);
if ~isempty(FC)
    if ~isempty(find(FC(:,1)>FC_r(1)-twindow,1))
        FC = FC(FC(:,1)>FC_r(1)-twindow,1);
    else
        FC = [];
    end
else
    FC = [];
end
if ~isempty(FC)
    if ~isempty(find(FC(:,1)<FC_r(end)+twindow,1))
        FC = FC(FC(:,1)<FC_r(end)+twindow,1);
    else
        FC = [];
    end
else
    FC = [];
end
deltaEventTO = nan(size(FC,1),3);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for e = 1:size(FC,1)
    deltaEv_now = FC(e,1)-FC_r;
    [~, posFCConsider] = min(abs(deltaEv_now));
    if abs(deltaEv_now(posFCConsider))<twindow
        deltaEventTO(e,1) = deltaEv_now(posFCConsider);
        deltaEventTO(e,2) = FC(e,1);
        deltaEventTO(e,3) = FC_r(posFCConsider);
    end
end
if length(unique(deltaEventTO(~isnan(deltaEventTO(:,3)),3))) ~= length(deltaEventTO(~isnan(deltaEventTO(:,3)),3))
    % an event has been associated twice!
    for t = 1:size(deltaEventTO,1)-1
        if (deltaEventTO(t+1,3)-deltaEventTO(t,3))== 0
            if abs(deltaEventTO(t+1,1))< abs(deltaEventTO(t,1))
                deltaEventTO(t,1) = nan;
                deltaEventTO(t,3) = nan;
            elseif abs(deltaEventTO(t+1,1)) == abs(deltaEventTO(t,1))
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
deltaEventTO = deltaEventTO(:,1);
end