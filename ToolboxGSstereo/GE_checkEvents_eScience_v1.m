function [IC, FC, missingIC] = GE_checkEvents_eScience_v1(Data, headerTraj, time, FootMrk, side, fs, GEs, GEs_v, fig)
    
    %% vertical marker trajectories
    H = Data.(headerTraj).(strcat(side,'_',FootMrk{1}))(:,3);
    T = Data.(headerTraj).(strcat(side,'_',FootMrk{2}))(:,3);
    if range(H)+0.25*range(H)>range(fillgaps(H)) && range(H)-0.25<range(fillgaps(H))      
        H = fillgaps(H);
    end
    if range(T)+0.25*range(T)>range(fillgaps(T)) && range(T)-0.25<range(fillgaps(T))      
        T = fillgaps(T);
    end
    %
    
    if side == 'r'
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
    [~,nStr] = findpeaks((H*1000),'MinPeakDistance',nMinDist,'MinPeakProminence',0.2*range(H)*1000);
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
            if length(dIC)>1
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
        else
            if i == 1
                missingIC = [missingIC;1,nStr(i)];
            else
                missingIC = [missingIC;nStr(i-1),nStr(i)];
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
    missingIC = [];    
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
    
    %% FIGURE
    if fig == 1
        figure('Name',strcat(side,'_Corrected_GEs'),'NumberTitle','off')
        plot(time(1:end), T)
        hold on
        plot(time(1:end), H, 'r')

        scatter((FC(:,1)/fs),T(FC(:,1)),'v')
        text((FC(1,1)/fs),T(FC(1,1)),'TO')
        scatter((IC(:,1)/fs),H(IC(:,1)),'^')
        text((IC(1,1)/fs),H(IC(1,1)),'HS')
        ylabel(strcat(side, ' Vertical Mark Traj [m]'))
        xlabel('Time [s]')   
    end   
end