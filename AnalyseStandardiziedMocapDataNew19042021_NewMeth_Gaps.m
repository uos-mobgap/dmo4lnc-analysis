%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Extract Gold Standard gait parameters from mocap Data
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Toolbox v.5 - updated on 02.12.2020 based on CAU data and NEW OUTPUT standardization
% New data output structure as agreed on 16.09.2020 (T2.2 group)
% Correct/check GEs obtained either during steps/turns
% Toolbox developed for the MOBILISE-D Project (https://www.mobilise-d.eu/)
%
% Author code:  Dr Tecla Bonci
%               email: t.bonci@sheffield.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
clc

addpath('ToolboxGSstereo');
folderNow = pwd;

% version results DryRun
versResult = 10;
saveFileDiffFolder = false; % True/False for saving figures and files in a different folder
                        
originalZeni = false;       % True/False in using the original or modified Zeni's Algo
TVS = 1; % this flag is used in several places and should be 1 for our purposes
%~~~~~~~~~~~~~~~~~~~~~~~

%% Change this path as relevant
folderProcessed = false; % True/False if the data should be saved in a "Processed" subfolder

%% Some functions accept the input "fig" selecting if relevant figures are shown
%~~~~~~~
Fig_Settings.fig = 0;        % 1 == figures are shown; otherwise insert any other values
Fig_Settings.figAngle = 0;   % 1 == figures are shown; otherwise insert any other values
Fig_Settings.saveFigure = 0; % 1 == save matlab figure in the relevant folder
Fig_Settings.figV = 0;       % 1 == show figure about GEs found with the vertical mrk traj
%~~~~~~~

%% PARAMETERS TO SET
% STRIDE DEFINITION
strideParam.max_st = 3;         % Maximum stride duration accepted for each stride to be considered a correct stride [s]
strideParam.min_st = 0.2;       % Minimum stride duration accepted for each stride to be considered a correct stride [s]
strideParam.min_sl = 0.15;      % Minimum stride length accepted for each stride to be considered a correct stride [m]
strideParam.max_h = 0.15;       % Max height difference between initial and final stride instants [m]

% WB - MicroWB
WB_Param.n_min = 2;          % n_min: minimum number of consecutive strides which define a WB (for the same side - i.e., 2L + 2R)
WB_Param.max_break = 3;      % max_break: Maximum duration between two consecutive eligible strides [s]
WB_Param.n_strides = 4;      % n_strides: minimum number of consecutive strides in a WB 

% TURNS
TurnParam.turnThres = 45;     % angle used to identify the start/stop of turns [deg]

% parameters to identify "sharp turns"
TurnParam.maxTurn = 360;      % angle used to identify "sharp" turns that might be removed [deg]
TurnParam.averVel = 200;      % [deg/s] - mean angular velocity
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% PARAMETERS TO SET - CWP
% STRIDE DEFINITION
CWP_Param.max_h_CWP = 2;      % Max height difference between initial and final stride instants [m]

% TURNS
% parameters to identify "sharp turns"
CWP_Param.maxTurn_CWP = 360;      % angle used to identify "sharp" turns that might be removed [deg]
CWP_Param.averVel_CWP = 200;       % [deg/s] - mean angular velocity
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%
% High values have been inserted to avoid misclassification of turns as
% being "sharp" - porposed values were "maxTurn_CWP = 150" and "averVel_CWP = 50"
%%

%% Threshold for assessing whether there is a swing after an IC of not
thrs_ALL = 0.10:0.02:0.2;

% because we have colon slicing using fs, we get lots of warnings. However,
% we can't fix this because using ints (100 or 101) for fs causes errors 
% with the butterworth filters. so just suppress the warnings.
warning('off', 'MATLAB:colon:nonIntegerIndex');
warning('off', 'MATLAB:colon:nonScalarOperands');

% Base dataset folder
baseFolder = 'C:\Users\ac4jmi\Desktop\DMO4LNC\dmo4lnc-analysis\Dataset';

% Get all cohort folders (e.g., 'HA', 'HB', etc.)
cohortDirs = dir(baseFolder);
isCohort = [cohortDirs.isdir] & ~ismember({cohortDirs.name}, {'.', '..'});
cohortDirs = cohortDirs(isCohort);

for c = 1:length(cohortDirs)
    cohort = cohortDirs(c).name;
    cohortPath = fullfile(baseFolder, cohort);

    % Get all subject folders inside the cohort folder
    sbjDirs = dir(cohortPath);
    isSubject = [sbjDirs.isdir] & ~ismember({sbjDirs.name}, {'.', '..'});
    sbjDirs = sbjDirs(isSubject);

    for s = 1:length(sbjDirs)
        sbjNow = sbjDirs(s).name
        folder = fullfile(cohortPath, sbjNow, 'Lab\Laboratory');

        % Get all .mat files in the Laboratory folder
        filePattern = fullfile(folder, '*.mat');
        fileList = dir(filePattern);

        InclinationFootGait = [];

        for tr = 1:length(fileList)
            if strcmp(fileList(tr).name, 'data.mat')
                fileName = fullfile(folder, fileList(tr).name);
                data = load(fileName);

                name_now = fieldnames(data);
                data = data.(name_now{1});
                TM = fieldnames(data);

                for tm = 1:size(TM, 1)
                    Test = fieldnames(data.(TM{tm}));
                    % Standing still task
                    for testNow = 1:size(Test, 1)
                        if strcmp(Test{testNow}, 'Test1') % test 1 is always standing still
                            InclinationFootStatic = identifyInclinationFoot(data, TM, tm, Test, testNow);
                        end
                    end
                    % Walking tasks
                    for testNow = 1:size(Test, 1)
                        if ~(strcmp(Test{testNow}, 'Test1'))
                            Pers = 1;
                            [data, ~, InclinationFootGait, InclinationFootStatic] = ...
                                extractDMOsFromData(data, TM, tm, Test, testNow, ...
                                strideParam, WB_Param, TurnParam, CWP_Param, ...
                                Fig_Settings, thrs_ALL, InclinationFootStatic, ...
                                InclinationFootGait, Pers, versResult, folder, ...
                                sbjNow, originalZeni, saveFileDiffFolder, TVS);
                        end
                    end
                end
            end
        end

        % Save results
        cd(folder)
        if folderProcessed
            folderNow_processed = fullfile(folder, 'Processed');
            if ~exist(folderNow_processed, 'dir')
                mkdir(folderNow_processed);
            end
            cd(folderNow_processed);
        end
        save data data;

        % Go back to subject folder
        cd(folder)
    end
end