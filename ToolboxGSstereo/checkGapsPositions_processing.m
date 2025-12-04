function CheckOutputs = checkGapsPositions_processing(Stereo_raw,Stereo, INDIP,mode, CheckOutputs,Test,testNow,Trials,tl,fs,h)

N = size(Stereo_raw.RHEEL,1);
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
    GapsStartStop.(strcat('DYNA',dynMrks{s})) = startEndGaps;
end

if isempty(GapsStartStop.LHEEL) && isempty(GapsStartStop.RHEEL)
    % no gaps in the trial;
    CheckOutputs.HeelMrkFlag = false;
else
    %% check where these gaps are
    if size(Stereo.(mode),1)>0
        % identify same WBs/CWPs
        for wb_SP = 1:size(Stereo.(mode),2)
            % find relevant INDIP WB
            for wb_INDIP = 1:size(INDIP.(mode),2)
                if ~isempty(intersect(Stereo.(mode)(wb_SP).Start*fs:Stereo.(mode)(wb_SP).End*fs,...
                        INDIP.(mode)(wb_INDIP).Start*fs:INDIP.(mode)(wb_INDIP).End*fs))
                    WB(wb_SP).IndexSP = wb_SP;
                    WB(wb_SP).IndexINDIP = wb_INDIP;
                else
                    % no match!
                end
            end
        end
        
        for wb = 1:size(WB,2)
            wb_s = WB(wb).IndexSP;
            wb_i = WB(wb).IndexINDIP;
            if ~isempty(wb_i)
                for s = 1:2
                    g = size(GapsStartStop.(strcat(side{s},'HEEL')),1);
                    Check_gaps = zeros(g,2);
                    for g = 1:size(GapsStartStop.(strcat(side{s},'HEEL')),1)
                        % start of the WB
                        if ~isempty(intersect(round(Stereo.(mode)(wb_s).Start-h)*fs:round(Stereo.(mode)(wb_s).Start+h)*fs,...
                                GapsStartStop.(strcat(side{s},'HEEL'))(g,1):GapsStartStop.(strcat(side{s},'HEEL'))(g,2)))
                            % check if there are differences with the INDIP
                            % system
                            
                            if Stereo.(mode)(wb_s).Start > INDIP.(mode)(wb_i).Start
                                deltaT = abs(Stereo.(mode)(wb_s).Start - INDIP.(mode)(wb_i).Start);
                                if deltaT>0.15
                                    Check_gaps(g,1)=true;
                                else
                                    Check_gaps(g,1)=false;
                                end
                            else
                                Check_gaps(g,1)=false;
                            end
                        else
                            Check_gaps(g,1) = false;
                        end % intersect @ initial edge
                        
                        % end of the WB
                        if ~isempty(intersect((round(Stereo.(mode)(wb_s).End-h))*fs:(round(Stereo.(mode)(wb_s).End+h))*fs,...
                                GapsStartStop.(strcat(side{s},'HEEL'))(g,1):GapsStartStop.(strcat(side{s},'HEEL'))(g,2)))
                            % check if there are differences with the INDIP
                            % system
                            
                            if Stereo.(mode)(wb_s).End < INDIP.(mode)(wb_i).End
                                deltaT = abs(Stereo.(mode)(wb_s).End - INDIP.(mode)(wb_i).End);
                                if deltaT>0.15
                                    Check_gaps(g,2)=true;
                                else
                                    Check_gaps(g,2)=false;
                                end
                            else
                                Check_gaps(g,2)=false;
                            end
                        else
                            Check_gaps(g,2)=false;
                        end % intersect @ final edge
                    end % for each gap start/stop
                    
                    % overall - start
                    if ~isempty(find(Check_gaps(:,1),1))
                        Check_gaps_all.(side{s})(1)= true;
                    else
                        Check_gaps_all.(side{s})(1)= false;
                    end
                    % overall - stop
                    if ~isempty(find(Check_gaps(:,2),1))
                        Check_gaps_all.(side{s})(2)= true;
                    else
                        Check_gaps_all.(side{s})(2)= false;
                    end
                end %sides
                
                % overall both sides - start
                if Check_gaps_all.R(1) || Check_gaps_all.L(1)
                    IssueEdges(wb).Start = true;
                else
                    IssueEdges(wb).Start = false;
                end
                % overall both sides - stop
                if Check_gaps_all.R(1) || Check_gaps_all.L(1)
                    IssueEdges(wb).End = true;
                else
                    IssueEdges(wb).End = false;
                end
            end % same wb
        end % different WBs
        
        if isempty(find([IssueEdges.Start], 1)) && isempty(find([IssueEdges.End], 1))
            CheckOutputs.HeelMrkFlag = false;
        else
            CheckOutputs.HeelMrkFlag = IssueEdges;
        end
        
    else
        %% ONLY INDIP WBs/CWPs
        % check if a WB was not found because there were mrks
        for wb = 1:size(INDIP.(mode),2)
            for s = 1:2
                g = size(GapsStartStop.(strcat(side{s},'HEEL')),1);                
                Check_gaps = zeros(g,2);
                for g = 1:size(GapsStartStop.(strcat(side{s},'HEEL')),1)
                    % start of the WB
                    if ~isempty(intersect(round(INDIP.(mode)(wb).Start-h)*fs:round(INDIP.(mode)(wb).Start+h)*fs,...
                            GapsStartStop.(strcat(side{s},'HEEL'))(g,1):GapsStartStop.(strcat(side{s},'HEEL'))(g,2)))
                        
                        Check_gaps(g,1)=true;
                    else
                        Check_gaps(g,1) = false;
                    end % intersect @ initial edge
                    
                    % end of the WB
                    if ~isempty(intersect((round(INDIP.(mode)(wb).End-h))*fs:(round(INDIP.(mode)(wb).End+h))*fs,...
                            GapsStartStop.(strcat(side{s},'HEEL'))(g,1):GapsStartStop.(strcat(side{s},'HEEL'))(g,2)))

                        Check_gaps(g,2)=true;
                    else
                        Check_gaps(g,2)=false;
                    end % intersect @ final edge
                end % for each gap start/stop
                
                % overall - start
                if ~isempty(find(Check_gaps(:,1),1))
                    Check_gaps_all.(side{s})(1)= true;
                else
                    Check_gaps_all.(side{s})(1)= false;
                end
                % overall - stop
                if ~isempty(find(Check_gaps(:,2),1))
                    Check_gaps_all.(side{s})(2)= true;
                else
                    Check_gaps_all.(side{s})(2)= false;
                end
            end %sides
            
            % overall both sides - start
            if Check_gaps_all.R(1) || Check_gaps_all.L(1)
                IssueEdges(wb).Start = true;
            else
                IssueEdges(wb).Start = false;
            end
            % overall both sides - stop
            if Check_gaps_all.R(1) || Check_gaps_all.L(1)
                IssueEdges(wb).End = true;
            else
                IssueEdges(wb).End = false;
            end
        end % different WBs
        
        if isempty(find([IssueEdges.Start], 1)) && isempty(find([IssueEdges.End], 1))
            CheckOutputs.HeelMrkFlag = false;
        else
            CheckOutputs.HeelMrkFlag = IssueEdges;
        end
    end
end % check if there are GAPs