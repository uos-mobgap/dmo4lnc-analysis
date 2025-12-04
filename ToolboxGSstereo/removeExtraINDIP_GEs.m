function GE_INDIP = removeExtraINDIP_GEs(GE_INDIP, Data, MrkGaps)

Gaps = [];
fs = Data.Fs;
if ~isempty(find(strcmp(MrkGaps, 'LTOE'), 1))|| ~isempty(find(strcmp(MrkGaps, 'RTOE'), 1))||...
        ~isempty(find(strcmp(MrkGaps, 'LHEEL'), 1))|| ~isempty(find(strcmp(MrkGaps, 'RHEEL'), 1))
    
    %% check where there are gaps
    mrksOfInt = {'RTOE','LTOE','RHEEL','LHEEL'};
    for m = 1:size(mrksOfInt,2)
        mNow = Data.Mrks.(mrksOfInt{m})(:,1);
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
            Gaps.(mrksOfInt{m}) = startStopGaps;
        end
    end
end
if ~isempty(Gaps)
    toCheck = fieldnames(Gaps);
    for t = 1:size(toCheck,1)
        if strcmp(toCheck(t),'RHEEL')
            PotentialAffectedEvents = round(GE_INDIP.HS.right*fs);
            for k = 1:size(Gaps.RHEEL,1)
                if ~isempty(intersect(Gaps.RHEEL(k,1):Gaps.RHEEL(k,2),PotentialAffectedEvents))
                    [~, posToRem] = intersect(PotentialAffectedEvents,Gaps.RHEEL(k,1):Gaps.RHEEL(k,2));
                    GE_INDIP.HS.right(posToRem) = [];
                end
            end
            
        elseif strcmp(toCheck(t),'LHEEL')
            PotentialAffectedEvents = round(GE_INDIP.HS.left*fs);
            for k = 1:size(Gaps.LHEEL,1)
                if ~isempty(intersect(Gaps.LHEEL(k,1):Gaps.LHEEL(k,2),PotentialAffectedEvents))
                    posToRem = intersect(Gaps.LHEEL(k,1):Gaps.LHEEL(k,2),PotentialAffectedEvents);
                end
            end
            
        elseif strcmp(toCheck(t),'RTOE')
            PotentialAffectedEvents = round(GE_INDIP.TO.right*fs);
            for k = 1:size(Gaps.RTOE,1)
                if ~isempty(intersect(Gaps.RTOE(k,1):Gaps.RTOE(k,2),PotentialAffectedEvents))
                    posToRem = intersect(Gaps.RTOE(k,1):Gaps.RTOE(k,2),PotentialAffectedEvents);
                end
            end
            
        else % strcmp(toCheck(t),'LTOE')
            PotentialAffectedEvents = round(GE_INDIP.TO.left*fs);
            for k = 1:size(Gaps.LTOE,1)
                if ~isempty(intersect(Gaps.LTOE(k,1):Gaps.LTOE(k,2),PotentialAffectedEvents))
                    posToRem = intersect(Gaps.LTOE(k,1):Gaps.LTOE(k,2),PotentialAffectedEvents);
                end
            end
        end
    end
    
end
end