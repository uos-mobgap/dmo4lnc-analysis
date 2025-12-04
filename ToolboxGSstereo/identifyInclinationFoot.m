function InclinationFootStatic = identifyInclinationFoot(data,TM, tm, Test,testNow)
Trials = fieldnames(data.(TM{tm}).(Test{testNow}));
for tl = 1:1%size(Trials,1)
    standardsAll = fieldnames(data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards);
    pos = find(strcmp(standardsAll,'Stereophoto_raw'));
    if ~isempty(pos)
        close all;
        %% Skin-markers data Headers
        Data = data.(TM{tm}).(Test{testNow}).(Trials{tl}).Standards.(standardsAll{pos});
        headerTraj = 'Mrks'; %'mrks'
        markers = fieldnames(Data.(headerTraj));
        fs = Data.Fs; %Data.(headerTraj).fs
        N = length(Data.(headerTraj).(markers{1}));
        time = 1/fs:1/fs:N/fs;
        
        %% Filter data and TimeFLAG
        [Data, MrkGaps, ~, ~, ~, ~] = dataFilteringAndFlags(Data, headerTraj, markers, fs, time, N);
        if isempty(MrkGaps)
            sides = {'R','L'};
            for s = 1:length(sides)
                PosHeel_3d = eval(['Data.(headerTraj).' sides{s} 'HEEL(1,:)']);
                PosToe_3d = eval(['Data.(headerTraj).' sides{s} 'TOE(1,:)']);
                VectorFootNow = PosToe_3d-PosHeel_3d;
                AngleNow = dot(VectorFootNow,[0 0 1])/(norm(VectorFootNow));
                InclinationFootStatic.(sides{s}) = 90 - rad2deg(acos(AngleNow));
            end
        else
            warning('check for gaps in the static!')
            sides = {'R','L'};
            for s = 1:length(sides)
                PosHeel_3d = eval(['Data.(headerTraj).' sides{s} 'HEEL(1,:)']);
                PosToe_3d = eval(['Data.(headerTraj).' sides{s} 'TOE(1,:)']);
                if isempty(isnan(PosHeel_3d))&&isempty(isnan(PosToe_3d))
                    VectorFootNow = PosToe_3d-PosHeel_3d;
                    AngleNow = dot(VectorFootNow,[0 0 1])/(norm(VectorFootNow));
                    InclinationFootStatic.(sides{s}) = 90 - rad2deg(acos(AngleNow));
                else
                    PosHeel_3d_all = eval(['Data.(headerTraj).' sides{s} 'HEEL']);
                    PosToe_3d_all = eval(['Data.(headerTraj).' sides{s} 'TOE']);
                    posTocheckH = find(~isnan(PosHeel_3d_all(:,1)));
                    posTocheckT = find(~isnan(PosToe_3d_all(:,1)));
                    commonFrames = intersect(posTocheckH,posTocheckT);
                    VectorFootNow = PosHeel_3d_all(min(commonFrames),:)-PosToe_3d_all(min(commonFrames),:);
                    AngleNow = dot(VectorFootNow,[0 0 1])/(norm(VectorFootNow));
                    InclinationFootStatic.(sides{s}) = 90 - rad2deg(acos(AngleNow));
                end
            end
        end   
    else
        InclinationFootStatic = [];
    end
end
end