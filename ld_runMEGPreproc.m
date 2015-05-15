% ld_runMEGPreproc
% addpath(genpath('/users2/purpadmin/Laura/MRI/GitRepo/fieldtrip'));
% addpath(genpath('/users2/purpadmin/Laura/MRI/GitRepo/eeglab13_4_4b'));
%% setup
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];

dataDir = sprintf('%s/%s', exptDir, sessionDir);
preprocDir = sprintf('%s/preproc', dataDir);
figDir = sprintf('%s/%s/%s', preprocDir, 'figures');

%% make the preproc dir if it doesn't exist
if ~exist(preprocDir,'dir')
    mkdir(preprocDir)
end
%% make the figDir dir if it doesn't exist
if ~exist(figDir,'dir')
    mkdir(figDir)
end

%% segment only if needed
runFiles = dir(sprintf('%s/%s*.sqd', preprocDir, fileBase));
if isempty(runFiles)
    %% segment original sqd into runs
    dataFile = sprintf('%s/%s.sqd', dataDir, fileBase);
    
    % check settings in ld_segmentSqd before running!
    nRuns = ld_segmentSqd(dataFile);
    runs = 1:nRuns
    
    %% move run files into preproc directory
    runFiles = dir(sprintf('%s/*run*.sqd', dataDir));
    for iRun = 1:nRuns
        movefile(sprintf('%s/%s', dataDir, runFiles(iRun).name), preprocDir)
    end
else
    % we have done preprocessing before, so find the number of runs
    runTag = getTag(runFiles(end).name,'run');
    nRuns = str2num(runTag(3:4));
    runs = 1:nRuns
end

%% view data
% % just from run 1
run1Data = sqdread(sprintf('%s/%s', preprocDir, runFiles(1).name));
run1Data  = run1Data(:,1:157)';

srate = 1000;
windowSize = [1 5 2560 1392];
eegplot(run1Data,'srate',srate,'winlength',20,'dispchans',80,'position',windowSize);

%% manually set bad channels
badChannels = []; % in matlab 1-indexing

%% run preproc for each run
for iRun = 1:nRuns
    run = runs(iRun);
    runFile = sprintf('%s/%s_run%02d.sqd', preprocDir, fileBase, run);
    preprocFileName = ld_MEGPreproc(runFile, figDir, badChannels);
end
keyboard
%% combine run files into preprocessed sqd
analStr = getTag(preprocFileName);
outFileName = sprintf('%s_%s.sqd', fileBase, analStr);
outfile = rd_combineSqd(preprocDir, outFileName, analStr)

%% view triggers for the combined file
rd_checkTriggers(outfile);