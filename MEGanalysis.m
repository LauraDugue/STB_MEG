% MEGanalysis.m
%
%      usage: MEGanalysis
%         by: laura
%       date: 21/05/15
%    purpose: this master program is used to run all sub-analysis program
%             for the MEG analysis

%% Add the path for fieldtrip and eeglab if necessary
addpath(genpath('/users2/purpadmin/Laura/MRI/GitRepo/fieldtrip'));
addpath(genpath('/users2/purpadmin/Laura/MRI/GitRepo/eeglab13_4_4b'));

%% Run the pre-processing of the data

ld_runMEGPreproc

%% Segment the .sqd file in half to run it through the manual inspection

exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];

dataDir = sprintf('%s/%s/preproc', exptDir, sessionDir);
preprocDir = sprintf('%s/%s/preprocmanual', exptDir, sessionDir);

%%% make the preproc dir if it doesn't exist
preprocmanualDir = sprintf('%s/preprocmanual', [exptDir sessionDir]);
if ~exist(preprocmanualDir,'dir')
    mkdir(preprocmanualDir)
end

runFiles = dir(sprintf('%s/%s*.sqd', preprocDir, fileBase));
if isempty(runFiles)
    %% segment original sqd into runs
    dataFile = sprintf('%s/%s_ebi.sqd', dataDir, fileBase);
    
    % check settings in ld_segmentSqd before running!
    nRuns = ld_segmentSqdHalf(dataFile);
    runs = 1:2;
    
    %% move run files into preproc directory
    runFiles = dir(sprintf('%s/*half*.sqd', dataDir));
    for iRun = runs
        movefile(sprintf('%s/%s', dataDir, runFiles(iRun).name), preprocDir)
    end
else
    % we have done preprocessing before, so find the number of runs
    runTag = getTag(runFiles(end).name,'run');
    nRuns = str2num(runTag(3:4));
    runs = 1:nRuns;
end

%% Manual Inspection through Fieldtrip
%RUN: ld_manualInspection

%% Epoching (cue, for just cue-only trials, or display, for just display-present trials)
%%% Make sure to recombine the two half of the data

%% Separating per condition - Saving as a mat file
 
%% ERP analysis

%% Time-frequency analysis

%% Phase-locking analysis
