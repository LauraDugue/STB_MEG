function nRuns = ld_segmentSqdHalf(fileName)
%
% function nRuns = ld_segmentSqd(fileName)
%
% Chunk data into smaller segments. Segments data between runs, based on 
% triggers.
%
% Rachel Denison (September 2014) and modified by Laura Dugue (2015)

%% setup
% file out
segmentLabel = 'half';

% trigger
trigChan = 163;

% 160 = Fixation (beginning of the trial)
% 161 = Cue onset
% 162 = target onset (target present trials)
% 166 = target onset (cue-only trials)
% 163 = Response cue onset

% segmentation options
segmentationOption = 'selectData'; % 'splitData' or 'selectData'
segmentCushion = 5; % used only if 'selectData'

% experiment info
nTrigsPerRun = 80; 
nRuns = 12;
nRunsPerSegment = 6;
trialDur = 3.4; %%%% TO ADJUST FOR EXO AND ENDO TRIAL TYPES

%% display sqd info
info = sqdread(fileName,'info')

%% read in all triggers
triggers = all_trigger(fileName, trigChan);
nTrigs = size(triggers,1);

%% check and inspect triggers
nTrigsExpected = nTrigsPerRun*nRuns;

if nTrigs~=nTrigsExpected
    fprintf('Found %d triggers, but expected %d!\n', nTrigs, nTrigsExpected)
    
    figure
    plot(triggers(:,1),triggers(:,2),'.')
    title('original triggers')
    
    fprintf('\nExiting to allow manual adjustments ...\n')
    return
end

%% exclude stray triggers
% manual (different for each subject)
if nTrigs>nTrigsExpected
    excludedTrigIdxs = nTrigsPerRun*6 + 1;
    triggers(excludedTrigIdxs,:) = [];
    nTrigs = size(triggers,1);

    fprintf('\nNow there are %d triggers (%d expected)\n\n', nTrigs, nTrigsExpected)
    figure
    plot(triggers(:,1),triggers(:,2),'.')
    title('after trigger exclusion')
end

%% specify data segments
Fs = info.SampleRate;

nTrigsPerSegment = nRunsPerSegment*nTrigsPerRun;
segmentFirstTrialIdxs = 1:nTrigsPerSegment:nTrigs;
segmentLastTrialIdxs = [segmentFirstTrialIdxs(2:end)-1 nTrigs];
segmentFirstTrialStartTimes = triggers(segmentFirstTrialIdxs,1);
segmentLastTrialEndTimes = triggers(segmentLastTrialIdxs,1) + trialDur*Fs;

switch segmentationOption
    case 'splitData'
        % % divide segments halfway through last trial end and first trial start of
        % % subsequent segments
        segmentDividePoints = mean([segmentLastTrialEndTimes(1:end-1) segmentFirstTrialStartTimes(2:end)],2)';
        segmentStartTimes = [int64(1) segmentDividePoints];
        segmentEndTimes = [segmentDividePoints+1 int64(info.SamplesAvailable)];
    case 'selectData'
        segmentStartTimes = segmentFirstTrialStartTimes - segmentCushion*Fs;
        segmentEndTimes = segmentLastTrialEndTimes + segmentCushion*Fs;
    otherwise
        error('segmentationOption not found')
end
nSegments = numel(segmentStartTimes);

%% read data segments and write new sqd file
for iSegment = 1:nSegments
    % read segment data from original file
	segment = sqdread(fileName, 'Samples', [segmentStartTimes(iSegment) segmentEndTimes(iSegment)]);
    
    % new segment file name
    segmentFileName = sprintf('%s_%s%02d.sqd', fileName(1:end-4), segmentLabel, iSegment);
    
    % write segment file 
    fprintf('Writing sqd: %s %d\n', segmentLabel, iSegment) 
    sqdwrite(fileName, segmentFileName, 'Data', segment);
    
    % check triggers in new file
    rd_checkTriggers(segmentFileName,[],0);
    title(sprintf('%s %d', segmentLabel, iSegment)) 
end

