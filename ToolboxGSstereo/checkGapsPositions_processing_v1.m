function CheckOutputs = checkGapsPositions_processing_v1(Stereo_raw,Stereo, INDIP,mode, CheckOutputs, fs,h)

%% gaps on heel markers
side = {'R','L'};
for s = 1:2
    valNow = eval(['Stereo_raw.' side{s} 'HEEL(:,3)']);
    GapsOfInt = find(isnan(valNow));
    startEndGaps = [];
    if ~isempty(GapsOfInt)
        m = find(diff(GapsOfInt)>1);
        if isempty(m)
            startEndGaps = [startEndGaps; GapsOfInt(1),GapsOfInt(end)];
        else
            for i = 1:size(m,1)
                if i == 1
                    startEndGaps = [startEndGaps; GapsOfInt(1),GapsOfInt(m(i))];
                else
                    startEndGaps = [startEndGaps; GapsOfInt(m(i-1)+1),GapsOfInt(m(i))];
                end
            end
            startEndGaps = [startEndGaps; GapsOfInt(m(i)+1),GapsOfInt(end)];
        end
    end
    GapsStartStop.(strcat(side{s},'HEEL')) = startEndGaps;
end

%% gaps on DYN markers
dynMrks = {'0','REF','X','Y'};
for s = 1:length(dynMrks)
    valNow = eval(['Stereo_raw.DYNA' dynMrks{s} '(:,3)']);
    Gaps.(strcat('DYNA',dynMrks{s}))= ~isnan(valNow);
end
% Overall pelvic mrk quality 
Quality.Pelvis =(Gaps.DYNA0 + Gaps.DYNAREF + Gaps.DYNAX + Gaps.DYNAY)/4;

%% gaps on Feet
footMrks = {'HEEL','REF','TOE','INDIP'};
for f = 1:length(dynMrks)
    for s = 1:2
        valNow = eval(['Stereo_raw.' side{s} footMrks{f} '(:,3)']);
        Gaps.(strcat(side{s},footMrks{f}))= ~isnan(valNow);
    end
end
% Overall feet mrk quality 
for s = 1:2
    Quality.(strcat(side{s},'foot'))= eval(['(Gaps.' side{s} 'HEEL + Gaps.' side{s} 'REF + Gaps.' side{s} 'TOE + Gaps.' side{s} 'INDIP)/4']);
end

segments = fieldnames(Quality);

if isempty(GapsStartStop.LHEEL) && isempty(GapsStartStop.RHEEL)
    % no gaps in the trial;
    CheckOutputs.HeelMrkFlag = false;
    
    % check overall data quality
    CheckOutputs = checkGapsWBmatching(Stereo, INDIP, mode, fs, h, side, GapsStartStop, Quality, segments, CheckOutputs);
    
else
    CheckOutputs.HeelMrkFlag = '';
    CheckOutputs = checkGapsWBmatching(Stereo, INDIP, mode, fs, h, side, GapsStartStop, Quality, segments, CheckOutputs);
end % check if there are GAPs