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
addpath(genpath('/users2/purpadmin/Laura/MRI/GitRepo/rufinLab'));

%% Run the pre-processing of the data

ld_runMEGPreproc

%% Segment the .sqd file in half to run it through the manual inspection
%%% directory
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];

dataDir = sprintf('%s/%s/preproc', exptDir, sessionDir);
preprocDir = sprintf('%s/%s/preprocmanual', exptDir, sessionDir);

%%% make the preproc dir if it doesn't exist
preprocmanualDir = sprintf('%s/preprocmanual', [exptDir '/' sessionDir]);
if ~exist(preprocmanualDir,'dir')
    mkdir(preprocmanualDir)
end

%%% segment original sqd into runs
dataFile = sprintf('%s/%s_ebi.sqd', dataDir, fileBase);

% check settings in ld_segmentSqd before running!
nRuns = ld_segmentSqdHalf(dataFile);
runs = 1:2;

%%% move run files into preproc directory
runFiles = dir(sprintf('%s/*half*.sqd', dataDir));
for iRun = runs
    movefile(sprintf('%s/%s', dataDir, runFiles(iRun).name), preprocDir)
end

%% Manual Inspection through Fieldtrip
% RUN: ld_manualInspection

%% Epoching (cueOnset, for just cue-only trials, or displayOnset, for just display-present trials)

%%% conditions
% 160 = Fixation (beginning of the trial)
% 161 = Cue onset
% 162 = target onset (target present trials)
% 166 = target onset (cue-only trials)
% 163 = Response cue onset
trigChan = 162; % epoched centered on cue onset
trigName = 'displayOnset';

%%% directory
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
preproc = 'preprocmanual';
matDir = 'matEpoch';

epochedData = [];
for half = 1:2
    data = ['ebi_half0' num2str(half)];
    
    ld_epoch(data,trigChan,trigName,exptDir,sessionDir,preproc,matDir,fileBase);
    
    load([exptDir '/' sessionDir '/' matDir '/epochedData_' data '_' trigName '.mat'])
    epochedData = cat(3,epochedData,trigData);
end
savename = sprintf('%s/%s/%s_epoch_workspace.mat', saveDir);
save([exptDir '/' sessionDir '/' matDir '/epochedData_' trigName '.mat'],'epochedData','nSamples','nChannels','nTrigs','trigName','-v7.3');

%% Separating per condition - Saving as a mat file
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
matDir = 'matEpoch';
trigName = 'cueOnset';

load([exptDir '/' sessionDir '/' matDir '/epochedData_' trigName '.mat'])
load([exptDir '/' sessionDir '/Log Files/Conditions/indexBehavior.mat'])

validData = epochedData(:,:,sessionValid);
save([exptDir '/' sessionDir '/' matDir '/validData_' trigName '.mat'],'validData','nSamples','nChannels','nTrigs','trigName','-v7.3');
invalidData = epochedData(:,:,sessionInvalid);
save([exptDir '/' sessionDir '/' matDir '/invalidData_' trigName '.mat'],'invalidData','nSamples','nChannels','nTrigs','trigName','-v7.3');
cueOnlyData = epochedData(:,:,sessioncueOnly);
save([exptDir '/' sessionDir '/' matDir '/cueOnlyData_' trigName '.mat'],'cueOnlyData','nSamples','nChannels','nTrigs','trigName','-v7.3');

%% ERP analysis - Load data
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
matDir = 'matEpoch';
trigName = 'cueOnset';

load([exptDir '/' sessionDir '/' matDir '/validData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/invalidData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/cueOnlyData_' trigName '.mat'])

load /Volumes/DRIVE1/DATA/laura/MEG/data_hdr.mat

%% Baseline correction - Average - plot data
tstart = -500;
tstop = 0;
baselinePeriod = tstart:tstop;
t = baselinePeriod;
inBaseline = ismember(t,baselinePeriod);

% Valid trials
baselineDC = squeeze(mean(mean(validData(inBaseline,:,:),1),3));
baselineTSeries = repmat(baselineDC,[size(validData,1),1,size(validData,3)]);
validData0 = validData;
validData1 = validData-baselineTSeries;
avgValid0 = mean(validData0,3);
avgValid1 = mean(validData1,3);

% Invalid trials
baselineDC = squeeze(mean(mean(invalidData(inBaseline,:,:),1),3));
baselineTSeries = repmat(baselineDC,[size(invalidData,1),1,size(invalidData,3)]);
invalidData0 = invalidData;
invalidData1 = invalidData-baselineTSeries;
avgInvalid0 = mean(invalidData0,3);
avgInvalid1 = mean(invalidData1,3);

% cueOnly trials
baselineDC = squeeze(mean(mean(cueOnlyData(inBaseline,:,:),1),3));
baselineTSeries = repmat(baselineDC,[size(cueOnlyData,1),1,size(cueOnlyData,3)]);
cueOnlyData0 = cueOnlyData;
cueOnlyData1 = cueOnlyData-baselineTSeries;
avgcueOnly0 = mean(cueOnlyData0,3);
avgcueOnly1 = mean(cueOnlyData1,3);

timePeriod = -500:1500;
figure;hold on;
plot(timePeriod,avgValid1(:,:));title('Valid trials')
figure;hold on;
plot(timePeriod,avgInvalid1(:,:));title('Invalid trials')
figure;hold on;
plot(timePeriod,avgcueOnly1(:,:));title('Cue-only trials')

badChannels = [];
inds = setdiff(0:156,badChannels)+1;;

figure;
fH = ssm_plotOnMesh(avgcueOnly1_157(1,:)', 'Cue-only', [], [], '2d');
% set(gca,'CLim',[0 4])
%% Time-frequency analysis

%% Phase-locking analysis
