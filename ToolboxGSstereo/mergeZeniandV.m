function [GEs, foreFootICs] = mergeZeniandV(GEs, GEs_v, Event, Side, TurnsNotDec, foreFootICs, fs)
    if ~isempty(TurnsNotDec)
        for i = 1:length(Event)
            for j = 1:length(Side)
                if ~isempty(GEs_v.(Event{i}).(Side{j}))
                    for k = 1:size(TurnsNotDec,1)
                        [~,d] = intersect(GEs_v.(Event{i}).(Side{j})(:,1),TurnsNotDec(k,1)-(fs/2):1:TurnsNotDec(k,2)+(fs/2));
                        for m = 1:size(d,1)
                            if ~isempty(GEs.(Event{i}).(Side{j}))
                                possibleGEs = GEs.(Event{i}).(Side{j})(GEs.(Event{i}).(Side{j})(:,1)<GEs_v.(Event{i}).(Side{j})(d(m),1));
                                if ~isempty(possibleGEs)
                                    if abs(possibleGEs(end)- GEs_v.(Event{i}).(Side{j})(d(m),1))/fs <.8*mean(diff(possibleGEs)/fs)
                                        d(m) = 0;
                                    end
                                    if size(possibleGEs,1)<size(GEs.(Event{i}).(Side{j}),1) && d(m)>0
                                        if abs(GEs.(Event{i}).(Side{j})(size(possibleGEs,1)+1,1)-GEs_v.(Event{i}).(Side{j})(d(m),1))/fs<.8*mean(diff(possibleGEs)/fs)
                                            d(m) = 0;
                                        end
                                    end
                                end
                            end
                        end
                        d = d(d>0);
                        GEs.(Event{i}).(Side{j})=[GEs.(Event{i}).(Side{j}); GEs_v.(Event{i}).(Side{j})(d,:)];
                    end
%                     if ~isempty(d)
                        [~,ia] = unique(GEs.(Event{i}).(Side{j})(:,1));
                        GEs.(Event{i}).(Side{j}) = GEs.(Event{i}).(Side{j})(ia,:);
%                     end
                    %% check inserted events!
                    checkNow = diff(GEs.(Event{i}).(Side{j})(:,1));
                    if ~isempty(find(checkNow<50, 1))
                        GEs.(Event{i}).(Side{j})(:,end+1) = ones(size(GEs.(Event{i}).(Side{j}),1),1);
                        posToCheck = find(checkNow<50);
                        for h = 1:length(posToCheck)
                            if GEs.(Event{i}).(Side{j})(posToCheck(h),2) == 0
                                GEs.(Event{i}).(Side{j})(posToCheck(h),end) = 0;
                            else
                                GEs.(Event{i}).(Side{j})(posToCheck(h)+1,end) = 0;
                            end
                        end
                        if strcmp(Event{i},'HS')
                            if size(find((GEs.(Event{i}).(Side{j})(:,3)==1)&(GEs.(Event{i}).(Side{j})(:,2)==1)),1)~=size(foreFootICs.(Side{j}),1)
                                keepZeniEvents = GEs.(Event{i}).(Side{j})((GEs.(Event{i}).(Side{j})(:,2)==1),:);
                                foreFootICs.(Side{j}) = foreFootICs.(Side{j})((keepZeniEvents(:,3)==1),1);
                            end
                        end
                        GEs.(Event{i}).(Side{j}) = GEs.(Event{i}).(Side{j})((GEs.(Event{i}).(Side{j})(:,3)==1),1:2);                        
                    end
                    %%
                end
            end
        end
    end
end