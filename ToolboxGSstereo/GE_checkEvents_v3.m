function [IC, FC, missingIC] = GE_checkEvents_v3(Data, headerTraj, time, FootMrk, side, fs, GEs, GEs_v, max_break, fig)
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% new version - info about GAPS is also included
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% vertical marker trajectories
H = Data.(headerTraj).(strcat(side,FootMrk{1}))(:,3);
T = Data.(headerTraj).(strcat(side,FootMrk{2}))(:,3);
StartStopGaps = [];

if range(H)+0.25*range(H)>range(fillgaps(H)) && range(H)-0.25<range(fillgaps(H))
    GapsNow = find(isnan(H));   
    if ~isempty(GapsNow)
        s = find(diff(GapsNow)>1);
        if ~isempty(s)
            for i = 1:size(s,1)
                if i == 1
                    StartStopGaps = [StartStopGaps; GapsNow(1),GapsNow(s(i))];
                else
                    StartStopGaps = [StartStopGaps; GapsNow(s(i-1)+1),GapsNow(s(i))];
                end
            end
            StartStopGaps = [StartStopGaps; GapsNow(s(i)+1),GapsNow(end)];
        else
            StartStopGaps = [StartStopGaps; GapsNow(1),GapsNow(end)];
        end
    end
    
    H = fillgaps(H);
end
if range(T)+0.25*range(T)>range(fillgaps(T)) && range(T)-0.25<range(fillgaps(T))
    T = fillgaps(T);
end
%

if side == 'R'
    s = 'right';
else
    s = 'left';
end

minDist = 0.5;
nMinDist = fs*minDist;
maxDist = 2;
nMaxDist = fs*maxDist;

IC_zeni = GEs.HS.(s);
FC_zeni = GEs.TO.(s);

IC_v = GEs_v.HS.(s);
FC_v = GEs_v.TO.(s);

IC = IC_zeni;
FC = FC_zeni;

%% find how many strides
[~,nStr] = findpeaks((H*1000),'MinPeakDistance',nMinDist,'MinPeakProminence',0.15*range(H)*1000);
[~,nStr2] = findpeaks((T*1000),'MinPeakDistance',nMinDist,'MinPeakProminence',0.25*range(T)*1000);

for i = 1:size(nStr,1)-1
    cond1 = find(nStr2<nStr(i+1));
    cond2 = find(nStr2>nStr(i),1);
    if ~isempty(cond1)
        if cond1(end) == cond2
            if (nStr(i+1)-nStr2(cond2)>nMinDist) && (nStr2(cond2)-nStr(i)>nMinDist)
                nStr = [nStr; nStr2(cond2)];
            end
        end
    end
end
nStr = sort(nStr);

%% check if there are GEs identified with the vertical algo that have been missed by the Zeni's one
for i = 1:size(IC_zeni,1)-1
    cond1 = find(IC_v(:,1)<IC_zeni(i+1,1));
    cond2 = find(IC_v(:,1)>IC_zeni(i,1),1);
    if ~isempty(cond1)
        if cond1(end) == cond2
            if (IC_zeni(i+1,1)-IC_v(cond2,1)>nMinDist) && (IC_v(cond2,1)-IC_zeni(i,1)>nMinDist)
                condA = nStr(nStr>IC_v(cond2,1));
                condB = nStr(nStr<IC_v(cond2,1));
                if ~isempty(condA)&& ~isempty(condB)
                    if condA(1)-condB(end)<nMaxDist
                        IC = [IC; IC_v(cond2,:)];
                    end
                end
            end
        elseif size(cond1,1)>1
            pos = intersect(cond1,cond2);
            if ~isempty(pos)
                if (IC_zeni(i+1,1)-IC_v(cond1(end),1)>nMinDist) && (IC_v(cond1(end),1)-IC_zeni(i,1)>nMinDist)
                    condA = nStr(nStr>IC_v(cond1(end),1));
                    condB = nStr(nStr<IC_v(cond1(end),1));
                    if ~isempty(condA)&& ~isempty(condB)
                        if condA(1)-condB(end)<nMaxDist
                            IC = [IC; IC_v(cond1(end),:)];
                        end
                    end
                end
            end
        end
    end
end

for i = 1:size(FC_zeni,1)-1
    cond1 = find(FC_v(:,1)<FC_zeni(i+1,1));
    cond2 = find(FC_v(:,1)>FC_zeni(i,1),1);
    if ~isempty(cond1)
        if cond1(end) == cond2
            if (FC_zeni(i+1,1)-FC_v(cond2,1)>nMinDist) && (FC_v(cond2,1)-FC_zeni(i,1)>nMinDist)
                condA = nStr(nStr>FC_v(cond2,1));
                condB = nStr(nStr<FC_v(cond2,1));
                if ~isempty(condA)&& ~isempty(condB)
                    if condA(1)-condB(end)<nMaxDist
                        FC = [FC; FC_v(cond2,:)];
                    end
                end
            end
        end
    end
end

%% check if there are additional GEs identified with the vertical algo
if size(IC_v,1) > size(IC,1) && size(IC,1)< size(nStr,1)
    for i = 1:size(IC_v,1)
        if i == 1
            condA = nStr(nStr>IC_v(i,1));
            condB = nStr(1) - IC_v(i,1);
            if ~isempty(condA)&& (condB>0)
                IC = [IC; IC_v(i,:)];
            end
            
        else
            condA = nStr(nStr>IC_v(i,1));
            condB = nStr(nStr<IC_v(i,1));
            if ~isempty(condA)&& ~isempty(condB)
                if condA(1)-condB(end)<nMaxDist
                    IC = [IC; IC_v(i,:)];
                end
            end
        end
    end
end

if size(FC_v,1) > size(FC,1) && size(FC,1)< size(nStr,1)
    for i = 1:size(FC_v,1)
        if i == 1
            condA = nStr(nStr>FC_v(i,1));
            condB = nStr(1) - FC_v(i,1);
            if ~isempty(condA)&& (condB>0)
                FC = [FC; FC_v(i,:)];
            end
            
        elseif  i == size(FC_v,1)
            condA = nStr(nStr>FC_v(i,1));
            condB = nStr(1) - FC_v(i,1);
            if ~isempty(condA)&& (condB>0)
                FC = [FC; FC_v(i,:)];
            end
            
        else
            condA = nStr(nStr>FC_v(i,1));
            condB = nStr(nStr<FC_v(i,1));
            if ~isempty(condA)&& ~isempty(condB)
                if condA(1)-condB(end)<nMaxDist
                    FC = [FC; FC_v(i,:)];
                end
            end
        end
    end
end

%%
if ~isempty(IC)
    [~,ia] = unique(IC(:,1));
    IC = IC(ia,:);
end
if ~isempty(FC)
    [~,ia] = unique(FC(:,1));
    FC = FC(ia,:);
end

%% check all data
missingIC = []; missingFC = [];
durMean = mean(diff(nStr)); durSTD = std(diff(nStr));

for i = 1:size(nStr,1)
    % initial contacts
    if i == 1
        durN = nStr(i)-1;
        [~,dIC] = intersect(IC(:,1),1:1:nStr(i));
    else
        durN = nStr(i)-nStr(i-1);
        [~,dIC] = intersect(IC(:,1),nStr(i-1):1:nStr(i));
    end
    if ~isempty(dIC)
        % at least an event has been found
        if length(dIC)>1 % an additional event has been found
            pos = find(IC(dIC,2)==1);
            if length(pos) == 1 && durN < durMean+2*durSTD % there is a single GE identified with the ZENI's algo
                posOK = dIC(pos);
                posN = find(dIC~=posOK);
                IC(dIC(posN),:)=[];
                
            elseif length(pos) > 1 && durN < durMean+2*durSTD % there more GEs identified with the ZENI's algo
                valuesN = H(IC(dIC(pos)));
                [~, valueToKeep] = min(valuesN);
                posOK = dIC(valueToKeep);
                posN = find(dIC~=posOK);
                IC(dIC(posN),:)=[];
                
            elseif length(pos) == 1  && durN >= durMean+2*durSTD % there is a single GE identified with the ZENI's algo
                posOK = dIC(pos);
                % all the values are kept!
                
            elseif length(pos) > 1 && durN >= durMean+2*durSTD % there more GEs identified with the ZENI's algo
                valuesN = H(IC(dIC(pos)));
                % all the values are kept!
                
            else % there more GEs identified with the vertical mrk traj
                valuesN = H(IC(dIC));
                [~, valueToKeep] = min(valuesN);
                posOK = dIC(valueToKeep);
                posN = find(dIC~=posOK);
                IC(dIC(posN),:)=[];
            end
        end
        
    else% find IC found with the vTry algo
        if i == 1
            [~,dIC] = intersect(IC_v(:,1),1:1:nStr(i));
        else
            [~,dIC] = intersect(IC_v(:,1),nStr(i-1):1:nStr(i));
        end        
        if ~isempty(dIC)
            if size(dIC,1) == 1
                posToAdd = find(IC(:,1)>IC_v(dIC,1),1);
                p = (1);
            elseif size(dIC,1)> 1 && i ==1
                [~, p] = min(abs(IC_v(dIC,1)-nStr(i)));
                posToAdd = find(IC(:,1)>IC_v(dIC(p),1),1); 
            else
                [~, p] = min(abs(IC_v(dIC,1)-nStr(i-1)));
                posToAdd = find(IC(:,1)>IC_v(dIC(p),1),1); 
            end
            
            % check value to add
            deltaT = abs(IC_v(dIC(p),1)-IC(posToAdd,1));
            if ~isempty(deltaT)
                if deltaT<durMean+2*std(diff(nStr)) && deltaT>durMean-2*std(diff(nStr))
                    IC(end+1,:) = nan(1,2);
                    IC(posToAdd+1:end,:)=IC(posToAdd:end-1,:);
                    IC(posToAdd,:) = IC_v(dIC(p),:); 
                end   
            end
        else
            warning('a correct IC cannot be assessed');
        end
    end
    
    % final contacts
    if i == 1
        [~,dFC] = intersect(FC(:,1),1:1:nStr(i));
    else
        [~,dFC] = intersect(FC(:,1),nStr(i-1):1:nStr(i));
    end
    
    if ~isempty(dFC)
        if length(dFC)>1
            pos = find(FC(dFC,2)==1);
            if length(pos) == 1 % there is a single GE identified with the ZENI's algo
                posOK = dFC(pos);
                posN = find(dFC~=posOK);
                FC(dFC(posN),:)=[];
                
            elseif length(pos) > 1 % there more GEs identified with the ZENI's algo
                valuesN = T(FC(dFC(pos)));
                [~, valueToKeep] = min(valuesN);
                posOK = dFC(valueToKeep);
                posN = find(dFC~=posOK);
                FC(dFC(posN),:)=[];
                
            else % there more GEs identified with the vertical mrk traj
                valuesN = T(FC(dFC));
                [~, valueToKeep] = min(valuesN);
                posOK = dFC(valueToKeep);
                posN = find(dFC~=posOK);
                FC(dFC(posN),:)=[];
            end
        end
    else
        if i == 1
            missingFC = [missingFC;1,nStr(i)];
        else
            missingFC = [missingFC;nStr(i-1),nStr(i)];
        end
    end
end
 
%%
nMinObsDist = 0.6*min(diff(nStr));
for i = 2:size(nStr,1)
    [~,dIC] = intersect(IC(:,1),nStr(i-1):1:nStr(i));
    if ~isempty(dIC)
        if length(dIC)>1
            selIC = IC(dIC,:);
            durIC = diff(selIC(:,1));
            pos = find(durIC<nMinObsDist);
            ICtoRem = [];
            if ~isempty(pos) % there are extra events
                for j = 1:size(pos,1)
                    if H(selIC(pos(j)+1,1))<H(selIC(pos(j),1))
                        ICtoRem = [ICtoRem; selIC(pos(j),:)];
                    else
                        ICtoRem = [ICtoRem; selIC(pos(j)+1,:)];
                    end
                end
            end
            for k = 1:size(ICtoRem,1)
                posN = find(IC(:,1) == ICtoRem(k,1));
                IC(posN,:) = [];
            end
        end
    else
        missingIC = [missingIC;nStr(i-1),nStr(i)];
    end
end

if ~isempty(IC)
    %% Check if an additional IC has been inserted
    check = true;
    while check
        if IC(1,1)-max_break*fs <= 50
            SigNow = H(1:IC(1,1));
        else
            SigNow = H(IC(1,1)-max_break*fs:IC(1,1));
        end
        if max(SigNow)< 0.45*max(H(IC(1,1)+1:end)) % && (max(SigNow)-H(1))< 0.1*(max(H(IC(1,1)+1:end))-H(1))
            if (IC(1,2)==0)
                if length(find(isnan(SigNow)))<= 0.7*length(SigNow)
                    if (max(SigNow)-H(1))<0.03
                        IC(1,:) = [];
                    else
                        check = false;
                    end
                else
                    check = false;
                end
            elseif (max(SigNow)-H(1))< 0.02 % max(SigNow)< 0.10 *max(H(IC(1,1)+1:end))
                IC(1,:) = []; % no movement!
            else
                check = false;
            end
        else
            check = false;
        end
    end
    
    %% Check if a final IC has been missed
    if IC(end,1)+ max_break*fs > length (H)
        SigNow = H(IC(end,1):end);
    else
        SigNow = H(IC(end,1):IC(end,1)+max_break*fs);
    end
    if max(SigNow)> 0.45*max(H(1:IC(end,1)))
        % look for the relevant IC
        possibleICs = IC_v(IC_v(:,1) >IC(end,1));
        if ~isempty(possibleICs)
            candidateIC_duration = min(possibleICs)-IC(end,1);
            if candidateIC_duration < durMean+4*durSTD && candidateIC_duration > durMean-4*durSTD
                IC(end+1,:) = [min(possibleICs), 0];
            elseif length(possibleICs)>1
                k = 2;
                while k <= length(possibleICs)
                    candidateIC_duration = possibleICs(k)-IC(end,1);
                    if candidateIC_duration < durMean+4*durSTD && candidateIC_duration > durMean-4*durSTD
                        IC(end+1,:) = [possibleICs(k), 0];
                    end
                    k = k+1;
                end
            end
        end
    end
end

%% CHECK SPURIOUS STRIDES
thr = 0.15;
thr2 = 0.2;
IC_temp = IC;
maxDeltaH_all = max(H)- H(IC(1),1);
maxDeltaH_intitialCond = max(H)- mean(H(fs:2*fs));
IC_temp(:,end+1) = ones(size(IC,1),1);
pos = zeros(1,size(IC,1)-1);
for k = 1:size(IC,1)-1
    [maxHNow, pos(k)] = max(H(IC(k):IC(k+1),1));
    delta1 = maxHNow - H(IC(k),1);
    delta2 = maxHNow - H(IC(k+1),1);
    if delta1<thr*maxDeltaH_all && delta2<thr*maxDeltaH_all
        %% condition that has been added on 13.11.2020
        if (IC_temp(k+1,2)==1 && IC_temp(k,2)==1) && (IC_temp(k+1,1) - IC_temp(k,1)) > 0.5*mean(diff(nStr))
            IC_temp(k+1,end)=1;
        else
            IC_temp(k+1,end)=0;
        end
        %%
    end
end
%% condition that has been added on 16.11.2020
stdPOS = std(pos);
if max(pos) < mean(pos)+2*stdPOS % only if data are normally distribuited - not true for SDA!
    extraIC = find(pos<thr2*mean(pos));
    window = 20;
    if ~isempty(extraIC)
        for k = 1:size(extraIC,2)
            signToCheck = H(IC(extraIC(k))-window:IC(extraIC(k)+1)+window);
            [~,pos1Min] = min(signToCheck); pos1Min=pos1Min +(IC(extraIC(k))-window);
            d1 = abs(pos1Min - IC(extraIC(k)));
            d2 = abs(pos1Min - IC(extraIC(k)+1));
            if IC_temp(extraIC(k)+1,end)==0 || IC_temp(extraIC(k),end)==0
                % this extraIC has been already considered
            elseif d1<d2
                IC_temp(extraIC(k)+1,end)=0;
            else
                IC_temp(extraIC(k),end)=0;
            end
        end
    end
end
%%
IC = IC(IC_temp(:,end)==1,:);

%% check strides before and after GAPS - 0-->1
% part added on 23.11.2020
stepD = diff(IC(:,1));
for i = 1:size(IC,1)-1
    if IC(i,2)~=IC(i+1,2)
        if i >1
            Hnow = H(IC(i-1,1)+.5*fs:IC(i+1,1)-.5*fs);
            delta = IC(i-1,1)+.5*fs;
        else
            Hnow = H(.5*fs:IC(i+1,1)-.5*fs);
            delta = .5*fs;
        end
        [~,posM]=max(Hnow);
        if abs(posM -length(Hnow))<10
            [~,posM]=max(Hnow(1:round(0.8*length(Hnow))));
        end
        if (H(IC(i),1)-min(Hnow(posM:end)))>0.02 % 2 cm higher than the min found!
            pos= find(islocalmin(Hnow(posM:end)),1); % min(Hnow(posM:end));
            if pos(1) ==1
                [~,pos]= min(Hnow(posM:end));
            end
            pos = pos(1)+delta+posM;
            IC(i,1) = pos; IC(i,2)=0;
        end
    end
    if i >1
        Hnow = H(IC(i-1,1)+.5*fs:IC(i+1,1)-.5*fs);
        delta = IC(i-1,1)+.5*fs;
    else
        Hnow = H(.5*fs:IC(i+1,1)-.5*fs);
        delta = .5*fs;
    end
    [~,posM]=max(Hnow);
    if (H(IC(i),1)-min(Hnow(posM:end)))>0.02 % 2 cm higher than the min found!
        pos= find(islocalmin(Hnow(posM:end)),1); % min(Hnow(posM:end));
        if pos(1) ==1
            [~,pos]= min(Hnow(posM:end));
        end
        pos = pos(1)+delta+posM;
        if IC(i,2) ==1 && (IC(i+1,1)-IC(i,1))>mean(stepD)-2*std(stepD) && (IC(i+1,1)-IC(i,1))<mean(stepD)+2*std(stepD)
            IC(i,1) = pos; IC(i,2)=0;
        elseif IC(i,2) ==1 && (IC(i+1,1)-pos)>mean(stepD)-2*std(stepD) && (IC(i+1,1)-pos)<mean(stepD)+2*std(stepD)
            IC(i,1) = pos; IC(i,2)=0;
        end
    end
end

%% check IF there are ICs that have been inserted where they should have NOT - there are GAPS!!!
IC(:,end+1) = ones(size(IC,1),1);
for i = 1:size(StartStopGaps,1)
    if ~isempty(intersect(find(IC(:,1)>StartStopGaps(i,1)), find(IC(:,1)<StartStopGaps(i,2))))
        posToRem = intersect(find(IC(:,1)>StartStopGaps(i,1)), find(IC(:,1)<StartStopGaps(i,2)));
        for k = 1:length(posToRem)
            IC(posToRem(k),end)=0;
        end
    end
end

IC = IC(IC(:,end)==1,:); IC = IC(:,1:2);

%% FIGURE
if fig == 1
    figure('Name',strcat(side,'_Corrected_GEs'),'NumberTitle','off')
    plot(time(1:end), T)
    hold on
    plot(time(1:end), H, 'r')
    
    scatter((FC(:,1)/fs),T(FC(:,1)),'v')
    if ~isempty(FC)
        text((FC(1,1)/fs),T(FC(1,1)),'TO')
    end
    scatter((IC(:,1)/fs),H(IC(:,1)),'^')
    if ~isempty(IC)
        text((IC(1,1)/fs),H(IC(1,1)),'HS')
    end
    ylabel(strcat(side, ' Vertical Mark Traj [m]'))
    xlabel('Time [s]')
end
end