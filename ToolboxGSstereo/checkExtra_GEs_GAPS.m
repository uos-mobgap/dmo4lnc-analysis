function check = checkExtra_GEs_GAPS(GE_now, Data, MrkGaps, GEname, side, versGE)

if strcmp(GEname,'HS')
    mrkNow = 'HEEL';
else
    mrkNow = 'TOE';
end
if ~isempty(find(strcmp(MrkGaps, strcat(side,mrkNow)), 1))    
    %% check where there are gaps
    mrksOfInt = strcat(side,mrkNow);
    mNow = Data.Mrks.(mrksOfInt)(:,1);
    
    if ~isempty(find(isnan(mNow), 1))
        gapsNow = find(isnan(mNow));
        startStopGaps = [];
        s = find(diff(gapsNow)>1);
        if ~isempty(s)
            for j = 1:size(s,1)
                if j == 1
                    startStopGaps = [startStopGaps; gapsNow(1),gapsNow(s(j))];
                else
                    startStopGaps = [startStopGaps; gapsNow(s(j-1)+1),gapsNow(s(j))];
                end
            end
            startStopGaps = [startStopGaps; gapsNow(s(j)+1),gapsNow(end)];
        else
            startStopGaps = [startStopGaps; gapsNow(1),gapsNow(end)];
        end
        
        PotentialAffectedEvents = GE_now;
        checkStartStop = zeros(size(startStopGaps,1),1);
        for k = 1:size(startStopGaps,1)
            if versGE == 10 % 10. Hrejac and Marshall - GEs
                if ~isempty(find(PotentialAffectedEvents>startStopGaps(k,1),1)) && ~isempty(find(PotentialAffectedEvents<startStopGaps(k,2),1))
                    checkStartStop(k,1) = 1;
                end
            else
                if ~isempty(intersect(startStopGaps(k,1):startStopGaps(k,2),PotentialAffectedEvents))
                    checkStartStop(k,1) = 1;
                end
            end
        end
        
        if ~isempty(find(checkStartStop,1))
            % there is a gap
            check = false;
        else
            % no gap
            check = true;
        end        
        
    else
        % no gap
        check = true;
    end
    
else
     % no gap
    check = true;
end

end