function CheckOutputs = WB_GAPS_edges_check(Stereo, Stereo_raw, INDIP, mode,fs)

CheckOutputs.HeelMrkFlag = false;
h = 2; % [s] window to use for assessing whether a gap is at the edges of a WB/CWP

if ~isempty(INDIP)
    if size(Stereo.(mode),1) ~=0 && size(INDIP.(mode),1) ~=0
        %% both SP and INDIP AVAILABLE
        CheckOutputs.RefSystem = 'Both';
        CheckOutputs = checkGapsPositions_processing_v1(Stereo_raw,Stereo, INDIP,mode, CheckOutputs, fs,h);
        
    else
        if size(INDIP.(mode),1) ==0 && size(Stereo.(mode),1)==0
            CheckOutputs.RefSystem = 'No WB';
        elseif size(INDIP.(mode),1) ==0
            CheckOutputs.RefSystem = 'No WB - INDIP';
        else
            CheckOutputs.RefSystem = 'No WB - SP';
            CheckOutputs = checkGapsPositions_processing_v1(Stereo_raw,Stereo, INDIP,mode, CheckOutputs, fs,h);
        end
    end
else
    % no indip output
    if size(Stereo.(mode),1) ~=0
        %% SP AVAILABLE
        CheckOutputs.RefSystem = 'SP - WB';
        
    else
        %% SP not AVAILABLE
        CheckOutputs.RefSystem = 'No WB - SP';        
    end
    CheckOutputs = checkGapsPositions_processing_v1(Stereo_raw,Stereo, INDIP,mode, CheckOutputs, fs,h);
end
end