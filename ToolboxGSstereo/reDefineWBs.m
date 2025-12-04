function CandidateWBs = reDefineWBs(WB_L, WB_R)
FlagL = zeros(1,size(WB_L,1));
if ~isempty(WB_R)
    for i = 1:size(WB_R,1)
        CandidateWBs(i).IndexR = i;
        CandidateWBs(i).IndexL = [];
    end
    for i = 1:size(WB_R,1)
        for j = 1:size(WB_L,1)
            if (~isempty(intersect(WB_R(i,1):1:WB_R(i,2),WB_L(j,1):1:WB_L(j,2))))
                CandidateWBs(i).IndexL = [CandidateWBs(i).IndexL j];
                FlagL(j)=1;
            else
            end
        end
    end
else
    CandidateWBs=[];
end

%% Save also WB_L with no intersections referring to FlagL
if ~isempty(WB_L)
    k=length(CandidateWBs)+1;
    for i=1:length(FlagL)
        if FlagL(i)==0
            %this is a WB_L with zero intersections
            CandidateWBs(k).IndexR= [];
            CandidateWBs(k).IndexL= i;
            k = k+1;
        else
        end
    end
end
end