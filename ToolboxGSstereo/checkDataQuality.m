function [OverallQuality, OverallQualityTest] = checkDataQuality(Stereo_raw, fs)

%% gaps on DYN markers
dynMrks = {'0','REF','X','Y'};
for s = 1:length(dynMrks)
    valNow = eval(['Stereo_raw.BACK' dynMrks{s} '(:,3)']);
    Gaps.(strcat('BACK',dynMrks{s}))= ~isnan(valNow);
end
% Overall pelvic mrk quality
Quality.Pelvis =(Gaps.BACK0 + Gaps.BACKREF + Gaps.BACKX + Gaps.BACKY)/4;

%% gaps on Feet
footMrks = {'HEEL','REF','TOE','MET5'};
side = {'L','R'};
for f = 1:length(dynMrks)
    for s = 1:2
        valNow = eval(['Stereo_raw.' side{s} footMrks{f} '(:,3)']);
        Gaps.(strcat(side{s},footMrks{f}))= ~isnan(valNow);
    end
end
% Overall feet mrk quality
for s = 1:2
    Quality.(strcat(side{s},'foot'))= eval(['(Gaps.' side{s} 'HEEL + Gaps.' side{s} 'REF + Gaps.' side{s} 'TOE + Gaps.' side{s} 'REF)/4']);
end

segments = fieldnames(Quality);

%% check overall mrk quality within the WB
for seg = 1:size(segments,1)
    valNow = Quality.(segments{seg});
    OverallQuality.(strcat('Quality_',segments{seg})).Mean = mean(valNow)*100;
    OverallQuality.(strcat('Quality_',segments{seg})).STD = std(valNow)*100;
    OverallQuality.(strcat('Quality_',segments{seg})).perc25 = prctile(valNow,25)*100;
    OverallQuality.(strcat('Quality_',segments{seg})).Median = median(valNow)*100;
    OverallQuality.(strcat('Quality_',segments{seg})).perc75 = prctile(valNow,75)*100;
end

%% Quality of heel mrk traj
N = size(valNow,1);
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

for s = 1:2
    OverallQuality.(strcat(side{s},'foot'))=(sum(Gaps.(strcat(side{s},'HEEL')))/N)*100;
    if ~isempty(GapsStartStop.(strcat(side{s},'HEEL')))
        OverallQuality.(strcat(side{s},'Gaps'))=(GapsStartStop.(strcat(side{s},'HEEL'))(:,2)-...
        GapsStartStop.(strcat(side{s},'HEEL'))(:,1))/fs;
        for k = 1:size(GapsStartStop.(strcat(side{s},'HEEL')),1)
            fprintf('Gaps %.1f [s] %s HEEL \n',OverallQuality.(strcat(side{s},'Gaps'))(k),side{s})
        end
    else
        OverallQuality.(strcat(side{s},'Gaps'))=[];
    end
end

fprintf('Overall Pelvis data quality: %.1f (%.1f) %% \n',OverallQuality.Quality_Pelvis.Mean,OverallQuality.Quality_Pelvis.STD)
fprintf('Overall R foot data quality: %.1f (%.1f) %% \n',OverallQuality.Quality_Rfoot.Mean,OverallQuality.Quality_Rfoot.STD)
fprintf('Overall L foot data quality: %.1f (%.1f) %% \n\n',OverallQuality.Quality_Rfoot.Mean,OverallQuality.Quality_Lfoot.STD)

OverallQualityTest.AllSegm = mean([OverallQuality.Quality_Pelvis.Mean, mean([OverallQuality.Quality_Rfoot.Mean,OverallQuality.Quality_Lfoot.Mean])]);
OverallQualityTest.LargeHeelGapsFlag = ~isempty(OverallQuality.LGaps)||~isempty(OverallQuality.RGaps);
