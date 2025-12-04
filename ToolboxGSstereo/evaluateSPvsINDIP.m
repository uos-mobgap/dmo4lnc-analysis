function [checkErrors, flag] = evaluateSPvsINDIP(Stereo, INDIP, mode, fs, thr)

flag = false;
p_to_check = {'Start','End','AverageStepCadence','AverageStrideSpeed'};

if ~isempty(INDIP)
    if ~isempty(INDIP.(mode))
        % identify same WBs/CWPs
        for wb_SP = 1:size(Stereo.(mode),2)
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
        
        if size(Stereo.(mode),2)>=1
            for wb = 1:size(WB,2)
                wb_s = WB(wb).IndexSP;
                wb_i = WB(wb).IndexINDIP;
                for p = 1:length(p_to_check)
                    if ~isempty(wb_s)&&~isempty(wb_i)
                        if p<3
                            checkErrors(wb).(p_to_check{p}) = (abs(Stereo.(mode)(wb_s).(p_to_check{p}) - INDIP.(mode)(wb_i).(p_to_check{p}))/...
                                Stereo.(mode)(wb_s).Duration)*100;
                        else
                            checkErrors(wb).(p_to_check{p}) = (abs(Stereo.(mode)(wb_s).(p_to_check{p}) - INDIP.(mode)(wb_i).(p_to_check{p}))/...
                                Stereo.(mode)(wb_s).(p_to_check{p}))*100;
                        end
                        if checkErrors(wb).(p_to_check{p})> thr
                            flag = true;
                            fprintf('Check %d WB %s - percentage error of %.2f %% \n',wb, p_to_check{p}, checkErrors(wb).(p_to_check{p}))
                        end
                    end
                end
                
                % Match GEs
                if ~isempty(wb_s) && ~isempty(wb_i)
                    IC_SP = round(Stereo.(mode)(wb_s).InitialContact_Event*fs);
                    IC_SP_LR = Stereo.(mode)(wb_s).InitialContact_LeftRight;
                    IC_INDIP = round(INDIP.(mode)(wb_i).InitialContact_Event*fs);
                    IC_INDIP_LR = INDIP.(mode)(wb_i).InitialContact_LeftRight;
                    mEvents(:,1:4) = zeros(size(IC_SP,2),4);
                    h = 25;
                    for k = 1:size(IC_SP,2)
                        % find the relevant matching IC found with the INDIP
                        if ~isnan(IC_SP(k))
                            posToMatch = find((IC_INDIP>IC_SP(k)-h) & (IC_INDIP<IC_SP(k)+h));
                            if isempty(posToMatch)
                                %                 warning('No IC found!')
                            elseif length(posToMatch) == 1
                                % only one IC found
                                mEvents(k,1)=1;
                                mEvents(k,2)=IC_SP(k);
                                mEvents(k,3)=IC_INDIP(posToMatch);
                                side_Now = IC_SP_LR(k);
                                side_Now_IN = IC_INDIP_LR(posToMatch);
                                if strcmp(side_Now,side_Now_IN)
                                    if strcmp(side_Now,'Left')
                                        mEvents(k,4) = 0;
                                    else
                                        mEvents(k,4) = 1;
                                    end
                                else
                                    % events from different sides
                                    mEvents(k,2)=nan;
                                    mEvents(k,3)=nan;
                                    mEvents(k,4) = nan;
                                end
                                mPosEvents(k,1) = k;
                                mPosEvents(k,2) = posToMatch;
                            else % more ICs found
                                posToConsider = IC_INDIP(posToMatch);
                                [~,p] = min(abs(posToConsider-IC_SP(k)));
                                mEvents(k,1)=1;
                                mEvents(k,2)=IC_SP(k);
                                side_Now = IC_SP_LR(k);
                                side_Now_IN = IC_INDIP_LR(posToMatch(p));
                                if strcmp(side_Now,side_Now_IN)
                                    mEvents(k,3)=IC_INDIP(posToMatch(p));
                                    mPosEvents(k,1) = k;
                                    mPosEvents(k,2) = posToMatch(p);
                                    if strcmp(side_Now,'Left')
                                        mEvents(k,4) = 0;
                                    else
                                        mEvents(k,4) = 1;
                                    end
                                else
                                    % events from different sides
                                    mEvents(k,2)=nan;
                                    mEvents(k,3)=nan;
                                    mEvents(k,4) = nan;
                                end
                            end
                        else
                            mEvents(k,2)=nan(1);
                            mEvents(k,3)=nan(1);
                            mEvents(k,3)=nan(1);
                        end
                    end
                    mEvents = mEvents(mEvents(:,1)==1,2:4);
                    fprintf('All: Bias %.2f Precision %.2f Accuracy %.2f ICs \n',nanmean(mEvents(:,1)-mEvents(:,2)), ...
                        nanstd(mEvents(:,1)-mEvents(:,2)),nanmean(abs(mEvents(:,1)-mEvents(:,2))))
                    fprintf('R Bias %.2f Precision %.2f Accuracy %.2f ICs \n',nanmean(mEvents(mEvents(:,3)==1,1)-...
                        mEvents(mEvents(:,3)==1,2)), nanstd(mEvents(mEvents(:,3)==1,1)-mEvents(mEvents(:,3)==1,2)),...
                        nanmean(abs(mEvents(mEvents(:,3)==1,1)-mEvents(mEvents(:,3)==1,2))))
                    fprintf('L Bias %.2f Precision %.2f Accuracy %.2f ICs \n\n',nanmean(mEvents(mEvents(:,3)==0,1)-...
                        mEvents(mEvents(:,3)==0,2)), nanstd(mEvents(mEvents(:,3)==0,1)-mEvents(mEvents(:,3)==0,2)),...
                        nanmean(abs(mEvents(mEvents(:,3)==0,1)- mEvents(mEvents(:,3)==0,2))))
                    
                    clear mEvents;
                end
            end
        else
            % no WB identified with the SP
            for p = 1:length(p_to_check)
                checkErrors.(p_to_check{p}) = [];
            end
            
        end
        
    else
        for p = 1:length(p_to_check)
            checkErrors.(p_to_check{p}) = [];
        end
    end
else
    for p = 1:length(p_to_check)
        checkErrors.(p_to_check{p}) = [];
    end    
end
end
