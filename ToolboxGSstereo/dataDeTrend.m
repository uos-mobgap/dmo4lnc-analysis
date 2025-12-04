function Data = dataDeTrend(Data, headerTraj, markers, fs)
%% check if there is a trend in the signal
mrkDyn = markers(contains(markers, 'DYN'));
fig = 1;
[Act, ~] = actRec_v3(Data, headerTraj, fs, mrkDyn, fig);
pos1 = find(Act,1,'first');
pos2 = find(Act,1,'last');
h1 = (Data.(headerTraj).RHEEL(pos1,3)+Data.(headerTraj).LHEEL(pos1,3))/2;
h2 = (Data.(headerTraj).RHEEL(pos2,3)+Data.(headerTraj).LHEEL(pos2,3))/2;
if abs(h2 - h1)>0.01 % delta > 1cm
    delta = h1-h2;
    N = size(Data.(headerTraj).RHEEL,1);
    lineToadd = nan(N,1);    
    lineToadd(1:pos1-1)=zeros(pos1-1,1);
    lineToadd(pos1:pos2) = (0:1/fs:(pos2-pos1)/fs)*(delta/(pos2-pos1))*(fs);
    lineToadd(pos2:end) = ones(size(lineToadd(pos2:end)))*lineToadd(pos2);
    figure
    for i = 1:size(markers,1)
        if ~strcmp(markers{i},'WRIST')
            subplot(size(markers,1)-1,1,i)
            %         figure
            hold on
            plot(Data.(headerTraj).(markers{i})(:,3),'b')
%             sigNow = detrend(Data.(headerTraj).(markers{i})(:,3), 'omitnan');
%             pos = find(~isnan(Data.(headerTraj).(markers{i})(:,3)),1);
%             if ~isempty(pos)
%                 deltaNow = -sigNow(1);
%                 sigNow = sigNow+deltaNow+Data.(headerTraj).(markers{i})(pos,3);
%             end
            sigNow = Data.(headerTraj).(markers{i})(:,3)+lineToadd;
            plot(sigNow,'r')
            Data.(headerTraj).(markers{i})(:,3) = sigNow;
        end
    end
end

end