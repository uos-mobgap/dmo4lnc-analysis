function CheckOutputs = checkGapsWBmatching(Stereo, INDIP, mode, fs, h, side, GapsStartStop, Quality, segments, CheckOutputs)

if ~isempty(INDIP)&& ~isempty(Stereo.(mode))
    %% check where these gaps are
    if size(Stereo.(mode),1)>0
        % identify same WBs/CWPs
        for wb_SP = 1:size(Stereo.(mode),2)
            WB(wb_SP).IndexSP = wb_SP;
            % find relevant INDIP WB
            for wb_INDIP = 1:size(INDIP.(mode),2)
                if ~isempty(intersect(round(Stereo.(mode)(wb_SP).Start*fs):round(Stereo.(mode)(wb_SP).End*fs),...
                        round(INDIP.(mode)(wb_INDIP).Start*fs):round(INDIP.(mode)(wb_INDIP).End*fs)))
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
                                if deltaT>0.25
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
                                if deltaT>0.25
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
                if Check_gaps_all.R(2) || Check_gaps_all.L(2)
                    IssueEdges(wb).End = true;
                else
                    IssueEdges(wb).End = false;
                end
                
                %% check overall mrk quality within the WB
                for seg = 1:size(segments,1)
                    valNow = Quality.(segments{seg})(round(Stereo.(mode)(wb_s).Start)*fs:round(Stereo.(mode)(wb_s).End)*fs);
                    CheckOutputs.(strcat('Quality_',segments{seg}))(wb).Mean = mean(valNow)*100;
                    CheckOutputs.(strcat('Quality_',segments{seg}))(wb).STD = std(valNow)*100;
                    CheckOutputs.(strcat('Quality_',segments{seg}))(wb).Min = min(valNow)*100;
                    CheckOutputs.(strcat('Quality_',segments{seg}))(wb).perc25 = prctile(valNow,25)*100;
                    CheckOutputs.(strcat('Quality_',segments{seg}))(wb).Median = median(valNow)*100;
                    CheckOutputs.(strcat('Quality_',segments{seg}))(wb).perc75 = prctile(valNow,75)*100;
                    CheckOutputs.(strcat('Quality_',segments{seg}))(wb).Max = max(valNow)*100;
                end
                % same wb
                
            else
                % no MATCH for the INDIP
                IssueEdges(wb).Start = false;
                IssueEdges(wb).End = false;                
            end 
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
            
            %% check overall mrk quality within the WB
            for seg = 1:size(segments,1)
                valNow = Quality.(segments{seg})(round(Stereo.(mode)(wb).Start)*fs:round(Stereo.(mode)(wb).End)*fs);
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb).Mean = mean(valNow)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb).STD = std(valNow)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb).Min = min(valNow)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb).perc25 = prctile(valNow,25)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb).Median = median(valNow)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb).perc75 = prctile(valNow,75)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb).Max = max(valNow)*100;
            end
            
        end % different WBs
        
        if isempty(find([IssueEdges.Start], 1)) && isempty(find([IssueEdges.End], 1))
            CheckOutputs.HeelMrkFlag = false;
        else
            CheckOutputs.HeelMrkFlag = IssueEdges;
        end
    end
    
elseif ~isempty(Stereo.(mode)) && isempty(INDIP)
    % no INDIP output
    %% check where these gaps are
    if size(Stereo.(mode),1)>0
        for wb_s = 1:size(Stereo.(mode),2)
            for s = 1:2
                g = size(GapsStartStop.(strcat(side{s},'HEEL')),1);
                Check_gaps = zeros(g,2);
                for g = 1:size(GapsStartStop.(strcat(side{s},'HEEL')),1)
                    % start of the WB
                    if ~isempty(intersect(round(Stereo.(mode)(wb_s).Start-h)*fs:round(Stereo.(mode)(wb_s).Start+h)*fs,...
                            GapsStartStop.(strcat(side{s},'HEEL'))(g,1):GapsStartStop.(strcat(side{s},'HEEL'))(g,2)))
                        
                        Check_gaps(g,1)=true;
                    else
                        Check_gaps(g,1) = false;
                    end % intersect @ initial edge
                    
                    % end of the WB
                    if ~isempty(intersect((round(Stereo.(mode)(wb_s).End-h))*fs:(round(Stereo.(mode)(wb_s).End+h))*fs,...
                            GapsStartStop.(strcat(side{s},'HEEL'))(g,1):GapsStartStop.(strcat(side{s},'HEEL'))(g,2)))
                        
                        Check_gaps(g,2)=false;
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
                IssueEdges(wb_s).Start = true;
            else
                IssueEdges(wb_s).Start = false;
            end
            % overall both sides - stop
            if Check_gaps_all.R(1) || Check_gaps_all.L(1)
                IssueEdges(wb_s).End = true;
            else
                IssueEdges(wb_s).End = false;
            end
            
            %% check overall mrk quality within the WB
            for seg = 1:size(segments,1)
                valNow = Quality.(segments{seg})(round(Stereo.(mode)(wb_s).Start)*fs:round(Stereo.(mode)(wb_s).End)*fs);
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb_s).Mean = mean(valNow)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb_s).STD = std(valNow)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb_s).Min = min(valNow)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb_s).perc25 = prctile(valNow,25)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb_s).Median = median(valNow)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb_s).perc75 = prctile(valNow,75)*100;
                CheckOutputs.(strcat('Quality_',segments{seg}))(wb_s).Max = max(valNow)*100;
            end
            
        end % different WBs
        
        if isempty(find([IssueEdges.Start], 1)) && isempty(find([IssueEdges.End], 1))
            CheckOutputs.HeelMrkFlag = false;
        else
            CheckOutputs.HeelMrkFlag = IssueEdges;
        end
    else
        %% check overall mrk quality within the acquisition
        for seg = 1:size(segments,1)
            valNow = Quality.(segments{seg});
            CheckOutputs.(strcat('Quality_',segments{seg}))(1).Mean = mean(valNow)*100;
            CheckOutputs.(strcat('Quality_',segments{seg}))(1).STD = std(valNow)*100;
            CheckOutputs.(strcat('Quality_',segments{seg}))(1).Min = min(valNow)*100;
            CheckOutputs.(strcat('Quality_',segments{seg}))(1).perc25 = prctile(valNow,25)*100;
            CheckOutputs.(strcat('Quality_',segments{seg}))(1).Median = median(valNow)*100;
            CheckOutputs.(strcat('Quality_',segments{seg}))(1).perc75 = prctile(valNow,75)*100;
            CheckOutputs.(strcat('Quality_',segments{seg}))(1).Max = max(valNow)*100;
        end
    end
    
elseif isempty(Stereo.(mode))
    % NO SP
    %% check overall mrk quality within the acquisition
    for seg = 1:size(segments,1)
        valNow = Quality.(segments{seg});
        CheckOutputs.(strcat('Quality_',segments{seg}))(1).Mean = mean(valNow)*100;
        CheckOutputs.(strcat('Quality_',segments{seg}))(1).STD = std(valNow)*100;
        CheckOutputs.(strcat('Quality_',segments{seg}))(1).Min = min(valNow)*100;
        CheckOutputs.(strcat('Quality_',segments{seg}))(1).perc25 = prctile(valNow,25)*100;
        CheckOutputs.(strcat('Quality_',segments{seg}))(1).Median = median(valNow)*100;
        CheckOutputs.(strcat('Quality_',segments{seg}))(1).perc75 = prctile(valNow,75)*100;
        CheckOutputs.(strcat('Quality_',segments{seg}))(1).Max = max(valNow)*100;
    end
end