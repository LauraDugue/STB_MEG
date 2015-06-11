% MEGanalysis.m
%
%      usage: MEGanalysis
%         by: laura
%       date: 21/05/15
%    purpose: this master program is used to run all sub-analysis program
%             for the MEG analysis

%% Obtain the ROIs in a nifti format for that observer to run MNE

% 1. Go in the right directory where the ROI is.
% 2. Load the T1
% 3. RUN: roi2nifti.m This will save 2 files (.hdr and .img) in the T1
% space that Noah is then converting in .mgz for MNE purposes

%% Add the path for fieldtrip and eeglab if necessary
%%% ATTENTION: eeglab and fieldtrip have redondable variables. DO NOT add both to the path at the same time!!!!
addpath(genpath('/users2/purpadmin/Laura/MRI/GitRepo/fieldtrip'));
% addpath(genpath('/users2/purpadmin/Laura/MRI/GitRepo/eeglab13_4_4b'));
% addpath(genpath('/users2/purpadmin/Laura/MRI/GitRepo/rufinLab'));

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
dataFile = sprintf('%s/%s_elbi.sqd', dataDir, fileBase);

% check settings in ld_segmentSqd before running!
nRuns = ld_segmentSqdHalf(dataFile);
runs = 1:2;

%%% move run files into preproc directory
runFiles = dir(sprintf('%s/*half*.sqd', dataDir));
for iRun = runs
    movefile(sprintf('%s/%s', dataDir, runFiles(iRun).name), preprocDir)
end

%% Manual Inspection through Fieldtrip
% RUN: ld_manualInspection.m

%% Epoching (cueOnset, for just cue-only trials, or displayOnset, for just display-present trials)

%%% conditions
% 160 = Fixation (beginning of the trial)
% 161 = Cue onset
% 162 = target onset (target present trials)
% 166 = target onset (cue-only trials)
% 163 = Response cue onset
trigChan = 161; % epoched centered on cue onset
trigName = 'cueOnset';

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
    data = ['elbi_half0' num2str(half)];
    
    %     ld_epoch(data,trigChan,trigName,exptDir,sessionDir,preproc,matDir,fileBase);
    
    load([exptDir '/' sessionDir '/' matDir '/epochedData_' data '_' trigName '.mat'])
    epochedData = cat(3,epochedData,trigData);
end
save([exptDir '/' sessionDir '/' matDir '/epochedData_elbi_' trigName '.mat'],'epochedData','nSamples','nChannels','nTrigs','trigName','-v7.3');

%% Obtain the indeces for each trial type

% RUN: behavior_analysis_meg.m

%% Separating per condition - Saving as a mat file
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
matDir = 'matEpoch';
trigName = 'cueOnset';

load([exptDir '/' sessionDir '/' matDir '/epochedData_elbi_' trigName '.mat'])
load([exptDir '/' sessionDir '/Log Files/Conditions/indexBehavior.mat'])

validCorrectLeftData = epochedData(:,:,sessionValidCorrectLeft);
save([exptDir '/' sessionDir '/' matDir '/validCorrectLeftData_' trigName '.mat'],'validCorrectLeftData','nSamples','nChannels','nTrigs','trigName','-v7.3');
validIncorrectLeftData = epochedData(:,:,sessionValidIncorrectLeft);
save([exptDir '/' sessionDir '/' matDir '/validIncorrectLeftData_' trigName '.mat'],'validIncorrectLeftData','nSamples','nChannels','nTrigs','trigName','-v7.3');
invalidCorrectLeftData = epochedData(:,:,sessionInvalidCorrectLeft);
save([exptDir '/' sessionDir '/' matDir '/invalidCorrectLeftData_' trigName '.mat'],'invalidCorrectLeftData','nSamples','nChannels','nTrigs','trigName','-v7.3');
invalidIncorrectLeftData = epochedData(:,:,sessionInvalidIncorrectLeft);
save([exptDir '/' sessionDir '/' matDir '/invalidIncorrectLeftData_' trigName '.mat'],'invalidIncorrectLeftData','nSamples','nChannels','nTrigs','trigName','-v7.3');

validCorrectRightData = epochedData(:,:,sessionValidCorrectRight);
save([exptDir '/' sessionDir '/' matDir '/validCorrectRightData_' trigName '.mat'],'validCorrectRightData','nSamples','nChannels','nTrigs','trigName','-v7.3');
validIncorrectRightData = epochedData(:,:,sessionValidIncorrectRight);
save([exptDir '/' sessionDir '/' matDir '/validIncorrectRightData_' trigName '.mat'],'validIncorrectRightData','nSamples','nChannels','nTrigs','trigName','-v7.3');
invalidCorrectRightData = epochedData(:,:,sessionInvalidCorrectRight);
save([exptDir '/' sessionDir '/' matDir '/invalidCorrectRightData_' trigName '.mat'],'invalidCorrectRightData','nSamples','nChannels','nTrigs','trigName','-v7.3');
invalidIncorrectRightData = epochedData(:,:,sessionInvalidIncorrectRight);
save([exptDir '/' sessionDir '/' matDir '/invalidIncorrectRightData_' trigName '.mat'],'invalidIncorrectRightData','nSamples','nChannels','nTrigs','trigName','-v7.3');

cueOnlyLeftData = epochedData(:,:,sessioncueOnlyLeft);
save([exptDir '/' sessionDir '/' matDir '/cueOnlyLeftData_' trigName '.mat'],'cueOnlyLeftData','nSamples','nChannels','nTrigs','trigName','-v7.3');
cueOnlyRightData = epochedData(:,:,sessioncueOnlyRight);
save([exptDir '/' sessionDir '/' matDir '/cueOnlyRightData_' trigName '.mat'],'cueOnlyRightData','nSamples','nChannels','nTrigs','trigName','-v7.3');

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
inds = setdiff(0:156,badChannels)+1;
cueOnly = squeeze(avgcueOnly1(500:1500,:,:));
avgcueOnly1_157 = to157chan(cueOnly,inds,'zeros');

figure;
fH = ssm_plotOnMesh(avgcueOnly1_157(1,:), 'Cue-only', [], data_hdr, '2d');
% set(gca,'CLim',[0 4])

%% Time-frequency analysis
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
matDir = 'matEpoch';
trigName = 'cueOnset';

load([exptDir '/' sessionDir '/' matDir '/validCorrectLeftData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/validIncorrectLeftData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/invalidCorrectLeftData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/invalidIncorrectLeftData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/validCorrectRightData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/validIncorrectRightData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/invalidCorrectRightData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/invalidIncorrectRightData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/cueOnlyLeftData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/cueOnlyRightData_' trigName '.mat'])

%%% make the timeFreq dir if it doesn't exist
timeFreqDir = sprintf('%s/timeFreq', [exptDir '/' sessionDir]);
if ~exist(timeFreqDir,'dir')
    mkdir(timeFreqDir)
end

for elec = 1:157
    % Valid-Correct-Left trials
    [tf, freqs, times] = timefreq(squeeze(validCorrectLeftData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_validCorrectLeft.mat'],'tf','freqs','times');
    % Valid-Incorrect-Left trials
    [tf, freqs, times] = timefreq(squeeze(validIncorrectLeftData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_validIncorrectLeft.mat'],'tf','freqs','times');
    % Invalid-Correct-Left trials
    [tf, freqs, times] = timefreq(squeeze(invalidCorrectLeftData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_invalidCorrectLeft.mat'],'tf','freqs','times');
    % Invalid-Incorrect-Left trials
    [tf, freqs, times] = timefreq(squeeze(invalidIncorrectLeftData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_invalidIncorrectLeft.mat'],'tf','freqs','times');
    % Valid-Correct-Right trials
    [tf, freqs, times] = timefreq(squeeze(validCorrectRightData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_validCorrectRight.mat'],'tf','freqs','times');
    % Valid-Incorrect-Right trials
    [tf, freqs, times] = timefreq(squeeze(validIncorrectRightData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_validIncorrectRight.mat'],'tf','freqs','times');
    % Invalid-Correct-Right trials
    [tf, freqs, times] = timefreq(squeeze(invalidCorrectRightData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_invalidCorrectRight.mat'],'tf','freqs','times');
    % Invalid-Incorrect-Right trials
    [tf, freqs, times] = timefreq(squeeze(invalidIncorrectRightData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_invalidIncorrectRight.mat'],'tf','freqs','times');
    % CueOnlyLeft trials
    [tf, freqs, times] = timefreq(squeeze(cueOnlyLeftData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_cueOnlyLeft.mat'],'tf','freqs','times');
    % CueOnlyRight trials
    [tf, freqs, times] = timefreq(squeeze(cueOnlyRightData(:,elec,:)), 1000, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_cueOnlyRight.mat'],'tf','freqs','times');
end

%% Amplitude analysis

exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];

%%% make the amplitude dir if it doesn't exist
amplitudeDir = sprintf('%s/amplitude', [exptDir '/' sessionDir]);
if ~exist(amplitudeDir,'dir')
    mkdir(amplitudeDir)
end

%%% run amplitude analysis
dirmat = [exptDir '/' sessionDir];
electrods = 1:157;
listCond = {'cueOnlyLeft','cueOnlyRight',...
    'validCorrectLeft','validIncorrectLeft','invalidCorrectLeft','invalidIncorrectLeft',...
    'validCorrectRight','validIncorrectRight','invalidCorrectRight','invalidIncorrectRight'};%

for cond = listCond
    for elec = electrods
        % load data
        dirElec = [dirmat '/timeFreq/' obs '_elec' num2str(elec) '_' cond{:}];
        load(dirElec)
        % compute amplitude
        amplitude = mean(abs(tf),3);
        save([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_' cond{:} '.mat'],'amplitude','freqs','times');
        % compute baseline corrected amplitude
        baseline_correct = mean(amplitude(:,1:61),2);
        baseline_correct = repmat(baseline_correct,1,size(times,2));
        amplitude_correct = amplitude - baseline_correct;
        save([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_correct_' cond{:} '.mat'],'amplitude_correct','freqs','times');
        % compute standard deviation
        amplitude_std_correct = std(abs(tf),[],3).^2/size(tf,3);
        save([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_std_correct_' cond{:} '.mat'],'amplitude_std_correct','freqs','times');
        
    end
end

%% Plot amplitudes - condition-by-condition
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];

%%% Condition to plot
listCond = {'cueOnlyLeft','cueOnlyRight',...
    'validCorrectLeft','validIncorrectLeft','invalidCorrectLeft','invalidIncorrectLeft',...
    'validCorrectRight','validIncorrectRight','invalidCorrectRight','invalidIncorrectRight'};

for cond = listCond
    dirmat = [exptDir '/' sessionDir];
    electrods = 1:157;
    amplitude_mean = zeros(50,512);
    amplitude_mean_correct = zeros(50,512);
    amplitude_mean_std_correct = zeros(50,512);
    
    for elec = electrods
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_' cond{:} '.mat'])
        amplitude_mean = amplitude_mean + amplitude;
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_correct_' cond{:} '.mat'])
        amplitude_mean_correct = amplitude_mean_correct  + amplitude_correct;
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_std_correct_' cond{:} '.mat'])
        amplitude_mean_std_correct = amplitude_mean_std_correct  + amplitude_std_correct;
    end
    amplitude_mean = amplitude_mean./length(electrods);
    amplitude_mean_correct = amplitude_mean_correct./length(electrods);
    amplitude_mean_std_correct = (sqrt(amplitude_mean_std_correct)./length(electrods))./sqrt(length(electrods));
    
    zscore = amplitude_mean_correct./amplitude_mean_std_correct;
    pvals = 1 - normcdf(double(zscore));
    pvals(pvals==0) = 0.0000000000000001;
    
    %% Plot amplitude
    figure;
    surf(times-500,freqs(1:50),double(amplitude_mean_correct));%
    hold on;
    view(0,90);
    set(gcf,'Renderer','Zbuffer');
    shading interp;
    colorbar
    % load('Ruf_colormap3.mat');colormap(mymap);
    set(gca,'clim',[-1000 1000])
    set(gca,'xlim',[-200 1200])
    set(gca,'ylim',[2 100])
    line([0 0], [0 100], [max(max(abs(amplitude_mean_correct))) max(max(abs(amplitude_mean_correct)))],'Color', 'w','LineStyle','--','LineWidth',2)
    if (strcmp(cond{:},{'invalidCorrectRight'}) || strcmp(cond{:},{'invalidIncorrectRight'})...
            || strcmp(cond{:},{'validCorrectRight'}) || strcmp(cond{:},{'validIncorrectRight'})...
            || strcmp(cond{:},{'invalidCorrectLeft'}) || strcmp(cond{:},{'invalidIncorrectLeft'})...
            || strcmp(cond{:},{'validCorrectLeft'}) || strcmp(cond{:},{'validIncorrectLeft'})) && strcmp(attCond,'exo')
        line([160 160], [0 100], [max(max(abs(amplitude_mean_correct))) max(max(abs(amplitude_mean_correct)))],'Color', 'b','LineStyle','--','LineWidth',2)
    end
    line([1010 1010], [0 100], [max(max(abs(amplitude_mean_correct))) max(max(abs(amplitude_mean_correct)))],'Color', 'b','LineStyle','--','LineWidth',2)
    title([cond{:} ' - all sensors'],'FontSize',16)
    hold off;
    
    namefig = sprintf(['/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/tfmap-' cond{:}]);
    %     print('-djpeg','-r600',namefig);
end

%% Plot amplitudes - difference between correct and incorrect trials - condition-by-condition
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];

%%% Condition to plot
listCond = {'validCorrectLeft','validIncorrectLeft','invalidCorrectLeft','invalidIncorrectLeft',...
    'validCorrectRight','validIncorrectRight','invalidCorrectRight','invalidIncorrectRight'};

for cond = listCond
    dirmat = [exptDir '/' sessionDir];
    electrods = 1:157;
    amplitude_mean = zeros(50,512);
    amplitude_mean_correct = zeros(50,512);
    amplitude_mean_std_correct = zeros(50,512);
    
    for elec = electrods
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_correct_' cond{:} '.mat'])
        amplitude_mean_correct = amplitude_mean_correct  + amplitude_correct;
    end
    
    if strcmp(cond{:},{'validCorrectLeft'})
        amplitude_mean_correct_validCorrectLeft = amplitude_mean_correct./length(electrods);
    elseif strcmp(cond{:},{'validIncorrectLeft'})
        amplitude_mean_correct_validIncorrectLeft = amplitude_mean_correct./length(electrods);
    elseif strcmp(cond{:},{'invalidCorrectLeft'})
        amplitude_mean_correct_invalidCorrectLeft = amplitude_mean_correct./length(electrods);
    elseif strcmp(cond{:},{'invalidIncorrectLeft'})
        amplitude_mean_correct_invalidIncorrectLeft = amplitude_mean_correct./length(electrods);
    elseif strcmp(cond{:},{'validCorrectRight'})
        amplitude_mean_correct_validCorrectRight = amplitude_mean_correct./length(electrods);
    elseif strcmp(cond{:},{'validIncorrectRight'})
        amplitude_mean_correct_validIncorrectRight = amplitude_mean_correct./length(electrods);
    elseif strcmp(cond{:},{'invalidCorrectRight'})
        amplitude_mean_correct_invalidCorrectRight = amplitude_mean_correct./length(electrods);
    elseif strcmp(cond{:},{'invalidIncorrectRight'})
        amplitude_mean_correct_invalidIncorrectRight = amplitude_mean_correct./length(electrods);
    end
end

%% Difference to plot
data = ((amplitude_mean_correct_invalidCorrectLeft+amplitude_mean_correct_invalidCorrectRight)/2) - ...
    ((amplitude_mean_correct_invalidIncorrectLeft+amplitude_mean_correct_invalidIncorrectRight)/2);
titleName = 'invalidCorrect - invalidIncorrect';
saveName = 'iC-iI';

figure;
surf(times-500,freqs(1:50),double(data));
hold on;
view(0,90);
set(gcf,'Renderer','Zbuffer');
shading interp;
colorbar
% load('Ruf_colormap3.mat');colormap(mymap);
set(gca,'clim',[-3000 3000])
set(gca,'xlim',[-200 1200])
set(gca,'ylim',[2 100])
line([0 0], [0 100], [max(max(abs(amplitude_mean_correct))) max(max(abs(amplitude_mean_correct)))],'Color', 'w','LineStyle','--','LineWidth',2)
line([160 160], [0 100], [max(max(abs(amplitude_mean_correct))) max(max(abs(amplitude_mean_correct)))],'Color', 'b','LineStyle','--','LineWidth',2)
line([1010 1010], [0 100], [max(max(abs(amplitude_mean_correct))) max(max(abs(amplitude_mean_correct)))],'Color', 'b','LineStyle','--','LineWidth',2)
title(titleName,'FontSize',16)
hold off;

namefig = sprintf(['/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/tfmap-dif' saveName]);
print('-djpeg','-r600',namefig);

%% Plot topography - condition-by-condition
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
load /Volumes/DRIVE1/DATA/laura/MEG/data_hdr.mat

%%% Times and Freqs to plot
indfreq = [8 12];%
indtime = [400 800];

%%% Condition to plot
listCond = {'cueOnlyLeft','cueOnlyRight',...
    'validCorrectLeft','validIncorrectLeft','invalidCorrectLeft','invalidIncorrectLeft',...
    'validCorrectRight','validIncorrectRight','invalidCorrectRight','invalidIncorrectRight'};

for cond = listCond
    dirmat = [exptDir '/' sessionDir];
    electrods = 1:157;
    amplitude_mean = zeros(1,length(electrods));
    amplitude_mean_correct = zeros(1,length(electrods));
    amplitude_mean_std_correct = zeros(1,length(electrods));
    
    for elec = electrods
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_' cond{:} '.mat'])
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_correct_' cond{:} '.mat'])
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_std_correct_' cond{:} '.mat'])
        
        time = find(times >= (indtime(1)+ 500) & times <= (indtime(2)+ 500));
        freq = find(freqs >= indfreq(1) & freqs <= indfreq(2));
        
        data = amplitude(freq,time); amplitude_mean(elec) = sum(data(:));
        data = amplitude_correct(freq,time); amplitude_mean_correct(elec) = sum(data(:));
        data = amplitude_std_correct(freq,time); amplitude_mean_std_correct(elec) = sum(data(:));
    end
    
    badChannels = [];
    inds = setdiff(0:156,badChannels)+1;
    data1_157 = to157chan(amplitude_mean_correct,inds,'zeros');
    
    figure;
    fH = ssm_plotOnMesh(data1_157, cond{:}, [], data_hdr, '2d');
    %     set(gca,'CLim',[-4000000 4000000])
    
    namefig = sprintf(['/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/topo-alpha-' cond{:}]);
    print('-djpeg','-r600',namefig);
end

%% Plot topography difference between cue ony left and right
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
load /Volumes/DRIVE1/DATA/laura/MEG/data_hdr.mat

%%% Times and Freqs to plot
indfreq = [8 12];%
indtime = [400 800];

%%% Condition to plot
% listCond = {'cueOnlyLeft','cueOnlyRight'};
% listCond = {'validCorrectLeft','validCorrectRight'};
% listCond = {'validIncorrectLeft','validIncorrectRight'};
% listCond = {'invalidCorrectLeft','invalidCorrectRight'};
listCond = {'invalidIncorrectLeft','invalidIncorrectRight'};

for cond = listCond
    dirmat = [exptDir '/' sessionDir];
    electrods = 1:157;
    amplitude_mean = zeros(1,length(electrods));
    amplitude_mean_correct = zeros(1,length(electrods));
    amplitude_mean_std_correct = zeros(1,length(electrods));
    
    for elec = electrods
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_' cond{:} '.mat'])
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_correct_' cond{:} '.mat'])
        load([dirmat '/amplitude/' obs '_elec' num2str(elec) '_amp_std_correct_' cond{:} '.mat'])
        
        time = find(times >= (indtime(1)+ 500) & times <= (indtime(2)+ 500));
        freq = find(freqs >= indfreq(1) & freqs <= indfreq(2));
        
        data = amplitude(freq,time); amplitude_mean(elec) = sum(data(:));
        data = amplitude_correct(freq,time); amplitude_mean_correct(elec) = sum(data(:));
        data = amplitude_std_correct(freq,time); amplitude_mean_std_correct(elec) = sum(data(:));
    end
    
    badChannels = [];
    inds = setdiff(0:156,badChannels)+1;
    data1_157 = to157chan(amplitude_mean_correct,inds,'zeros');
    
    if strcmp(cond{:},{'cueOnlyLeft'}) || strcmp(cond{:},{'validCorrectLeft'})...
            || strcmp(cond{:},{'validIncorrectLeft'}) || strcmp(cond{:},{'invalidCorrectLeft'}) || strcmp(cond{:},{'invalidIncorrectLeft'})
        cueLeft = data1_157;
    elseif strcmp(cond{:},{'cueOnlyRight'}) || strcmp(cond{:},{'validCorrectRight'})...
            || strcmp(cond{:},{'validIncorrectRight'}) || strcmp(cond{:},{'invalidCorrectRight'}) || strcmp(cond{:},{'invalidIncorrectRight'})
        cueRight = data1_157;
    end
end

data = cueLeft - cueRight;

figure;
fH = ssm_plotOnMesh(data, [listCond{1} ' - ' listCond{2}], [], data_hdr, '2d');
set(gca,'CLim',[-4000000 4000000])

namefig = sprintf('/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/topo-alpha-dif-iI');
print('-djpeg','-r600',namefig);


%% Phase-locking analysis
%%% set-up path 
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];

%%% make the amplitude dir if it doesn't exist
phaseLockingDir = sprintf('%s/phaseLocking', [exptDir '/' sessionDir]);
if ~exist(phaseLockingDir,'dir')
    mkdir(phaseLockingDir)
end

%%% run phase-locking analysis
dirmat = [exptDir '/' sessionDir];
electrods = 1:157;
% listCond = {'cueOnlyLeft','cueOnlyRight',...
%     'validCorrectLeft','validIncorrectLeft','invalidCorrectLeft','invalidIncorrectLeft',...
%     'validCorrectRight','validIncorrectRight','invalidCorrectRight','invalidIncorrectRight'};%
listCond = {'validCorrectLeft','invalidCorrectLeft','validCorrectRight','invalidCorrectRight'};%

for cond = listCond
    for elec = electrods
        % load data
        load([dirmat '/timeFreq/' obs '_elec' num2str(elec) '_' cond{:}]);
        tfTest = tf;
        phaseLocking = zeros(50,512);
        % compute phase-locking
        if strcmp(cond{:},{'validCorrectLeft'})% sub-sample to equate for the number of trials in correct and incorrect conditions
            load([dirmat '/timeFreq/' obs '_elec' num2str(elec) '_validIncorrectLeft']);
            tfIncorrect = tf;
            for a = 1:100 
                numTrials = randsample(size(tfTest,3),size(tfIncorrect,3));
                phaseLocking = phaseLocking + (abs(mean(tfTest(:,:,numTrials)./abs(tfTest(:,:,numTrials)),3)));
            end;
            phaseLocking = phaseLocking./100;
        elseif strcmp(cond{:},{'invalidCorrectLeft'})
            load([dirmat '/timeFreq/' obs '_elec' num2str(elec) '_invalidIncorrectLeft']);
            tfIncorrect = tf;
            for a = 1:100 
                numTrials = randsample(size(tfTest,3),size(tfIncorrect,3));
                phaseLocking = phaseLocking + (abs(mean(tfTest(:,:,numTrials)./abs(tfTest(:,:,numTrials)),3)));
            end;
            phaseLocking = phaseLocking./100;
        elseif strcmp(cond{:},{'validCorrectRight'})
            load([dirmat '/timeFreq/' obs '_elec' num2str(elec) '_validIncorrectRight']);
            tfIncorrect = tf;
            for a = 1:100 
                numTrials = randsample(size(tfTest,3),size(tfIncorrect,3));
                phaseLocking = phaseLocking + (abs(mean(tfTest(:,:,numTrials)./abs(tfTest(:,:,numTrials)),3)));
            end;
            phaseLocking = phaseLocking./100;
        elseif strcmp(cond{:},{'invalidCorrectRight'})
            load([dirmat '/timeFreq/' obs '_elec' num2str(elec) '_invalidIncorrectRight']);
            tfIncorrect = tf;
            for a = 1:100 
                numTrials = randsample(size(tfTest,3),size(tfIncorrect,3));
                phaseLocking = phaseLocking + (abs(mean(tfTest(:,:,numTrials)./abs(tfTest(:,:,numTrials)),3)));
            end;
            phaseLocking = phaseLocking./100;
        elseif strcmp(cond{:},{'validIncorrectLeft'}) || strcmp(cond{:},{'invalidIncorrectLeft'})...
                || strcmp(cond{:},{'validIncorrectRight'}) || strcmp(cond{:},{'invalidIncorrectRight'})...
                || strcmp(cond{:},{'cueOnlyLeft'}) || strcmp(cond{:},{'cueOnlyRight'})
            phaseLocking = abs(mean(tfTest./abs(tfTest),3));
        end
        save([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_' cond{:} 'sub.mat'],'phaseLocking','freqs','times')
    end
end

%% Plot Phase-locking
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
matDir = 'matEpoch';
trigName = 'cueOnset';

%%% Condition to plot
listCond = {'cueOnlyLeft','cueOnlyRight',...
    'validCorrectLeft','validIncorrectLeft','invalidCorrectLeft','invalidIncorrectLeft',...
    'validCorrectRight','validIncorrectRight','invalidCorrectRight','invalidIncorrectRight'};%

for cond = listCond
    dirmat = [exptDir '/' sessionDir];
    electrods = 1:157;
    phaseLocking_mean = zeros(50,512);
    
    for elec = electrods
        load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_' cond{:} '.mat'])
        phaseLocking_mean = phaseLocking_mean + phaseLocking;
    end
    phaseLocking_mean = phaseLocking_mean./length(electrods);
    
    figure;
    surf(times-500,freqs(1:50),double(phaseLocking_mean));
    hold on;
    view(0,90);
    set(gcf,'Renderer','Zbuffer');
    shading interp;
    colorbar
    % load('Ruf_colormap3.mat');colormap(mymap);
    set(gca,'clim',[0 .2])
    set(gca,'xlim',[-200 1200])
    set(gca,'ylim',[2 100])
    line([0 0], [0 100], [max(max(abs(phaseLocking_mean))) max(max(abs(phaseLocking_mean)))],'Color', 'w','LineStyle','--','LineWidth',2)
    if (strcmp(cond{:},{'invalidCorrect'}) || strcmp(cond{:},{'invalidIncorrect'})...
            || strcmp(cond{:},{'validCorrect'}) || strcmp(cond{:},{'validIncorrect'})) && strcmp(attCond,'exo')
        line([160 160], [0 100], [max(max(abs(phaseLocking_mean))) max(max(abs(phaseLocking_mean)))],'Color', 'b','LineStyle','--','LineWidth',2)
    end
    line([1010 1010], [0 100], [max(max(abs(phaseLocking_mean))) max(max(abs(phaseLocking_mean)))],'Color', 'b','LineStyle','--','LineWidth',2)
    title([cond{:} ' - all sensors'],'FontSize',16)
    hold off;
    
    namefig = sprintf(['/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/tfmap-PL-' cond{:}]);
    print('-djpeg','-r600',namefig);
    
end

%% Bootstrap of Phase-locking - Correct vs. Incorrect

exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];

%%% make the amplitude dir if it doesn't exist
bootPhaseLockingDir = sprintf('%s/bootPhaseLocking', [exptDir '/' sessionDir]);
if ~exist(bootPhaseLockingDir,'dir')
    mkdir(bootPhaseLockingDir)
end

dirmat = [exptDir '/' sessionDir];
electrods = 1:157;
repetition = 100;
listCond = {'validCorrectLeft','validIncorrectLeft','invalidCorrectLeft','invalidIncorrectLeft',...
    'validCorrectRight','validIncorrectRight','invalidCorrectRight','invalidIncorrectRight'};%

for elec = electrods
    for cond = listCond
        % load data
        dirElec = [dirmat '/timeFreq/' obs '_elec' num2str(elec) '_' cond{:}];
        load(dirElec)
        
        if strcmp(cond{:},{'validCorrectLeft'})
            tf_validCorrectLeft = tf;
        elseif strcmp(cond{:},{'validIncorrectLeft'})
            tf_validIncorrectLeft = tf;
        elseif strcmp(cond{:},{'invalidCorrectLeft'})
            tf_invalidCorrectLeft = tf;
        elseif strcmp(cond{:},{'invalidIncorrectLeft'})
            tf_invalidIncorrectLeft = tf;
        elseif strcmp(cond{:},{'validCorrectRight'})
            tf_validCorrectRight = tf;
        elseif strcmp(cond{:},{'validIncorrectRight'})
            tf_validIncorrectRight = tf;
        elseif strcmp(cond{:},{'invalidCorrectRight'})
            tf_invalidCorrectRight = tf;
        elseif strcmp(cond{:},{'invalidIncorrectRight'})
            tf_invalidIncorrectRight = tf;
        end
    end
        
    PL_hit_validLeft = zeros(50,512);PL_hit_std_validLeft = zeros(50,512);PL_mis_validLeft = zeros(50,512);PL_mis_std_validLeft = zeros(50,512);
    PL_hit_invalidLeft = zeros(50,512);PL_hit_std_invalidLeft = zeros(50,512);PL_mis_invalidLeft = zeros(50,512);PL_mis_std_invalidLeft = zeros(50,512);
    PL_hit_validRight = zeros(50,512);PL_hit_std_validRight = zeros(50,512);PL_mis_validRight = zeros(50,512);PL_mis_std_validRight = zeros(50,512);
    PL_hit_invalidRight = zeros(50,512);PL_hit_std_invalidRight = zeros(50,512);PL_mis_invalidRight = zeros(50,512);PL_mis_std_invalidRight = zeros(50,512);
    
    for rep = 1:repetition
        disp(['Sub_sampling = '  num2str(rep)])
        
        numTrials = randsample(size(tf_validCorrectLeft,3),size(tf_validIncorrectLeft,3));
        tf_validLeft = cat(3,tf_validCorrectLeft(:,:,numTrials),tf_validIncorrectLeft);
        numTrials = randsample(size(tf_invalidCorrectLeft,3),size(tf_invalidIncorrectLeft,3));
        tf_invalidLeft = cat(3,tf_invalidCorrectLeft(:,:,numTrials),tf_invalidIncorrectLeft);
        numTrials = randsample(size(tf_validCorrectRight,3),size(tf_validIncorrectRight,3));
        tf_validRight = cat(3,tf_validCorrectRight(:,:,numTrials),tf_validIncorrectRight);
        numTrials = randsample(size(tf_invalidCorrectRight,3),size(tf_invalidIncorrectRight,3));
        tf_invalidRight = cat(3,tf_invalidCorrectRight(:,:,numTrials),tf_invalidIncorrectRight);
        
        parfor rep2 = 1:repetition
            disp(['Surrogate = '  num2str(rep2)])
            
            %%% Valid Left
            ind = randperm(size(tf_validLeft,3));
            tfhit_validLeft = tf_validLeft(:,:,ind(1:size(tf_validIncorrectLeft,3)));
            tfmis_validLeft = tf_validLeft(:,:,ind((size(tf_validIncorrectLeft,3)+1):end));
            PL_hit_validLeft = PL_hit_validLeft+abs(mean(tfhit_validLeft./abs(tfhit_validLeft),3));
            PL_mis_validLeft = PL_mis_validLeft+abs(mean(tfmis_validLeft./abs(tfmis_validLeft),3));
            PL_hit_std_validLeft = PL_hit_std_validLeft+(abs(mean(tfhit_validLeft./abs(tfhit_validLeft),3))).^2;
            PL_mis_std_validLeft = PL_mis_std_validLeft+(abs(mean(tfmis_validLeft./abs(tfmis_validLeft),3))).^2;
            %%% Invalid Left
            ind = randperm(size(tf_invalidLeft,3));
            tfhit_invalidLeft = tf_invalidLeft(:,:,ind(1:size(tf_invalidIncorrectLeft,3)));
            tfmis_invalidLeft = tf_invalidLeft(:,:,ind((size(tf_invalidIncorrectLeft,3)+1):end));
            PL_hit_invalidLeft = PL_hit_invalidLeft+abs(mean(tfhit_invalidLeft./abs(tfhit_invalidLeft),3));
            PL_mis_invalidLeft = PL_mis_invalidLeft+abs(mean(tfmis_invalidLeft./abs(tfmis_invalidLeft),3));
            PL_hit_std_invalidLeft = PL_hit_std_invalidLeft+(abs(mean(tfhit_invalidLeft./abs(tfhit_invalidLeft),3))).^2;
            PL_mis_std_invalidLeft = PL_mis_std_invalidLeft+(abs(mean(tfmis_invalidLeft./abs(tfmis_invalidLeft),3))).^2;
            %%% Valid Right
            ind = randperm(size(tf_validRight,3));
            tfhit_validRight = tf_validRight(:,:,ind(1:size(tf_validIncorrectRight,3)));
            tfmis_validRight = tf_validRight(:,:,ind((size(tf_validIncorrectRight,3)+1):end));
            PL_hit_validRight = PL_hit_validRight+abs(mean(tfhit_validRight./abs(tfhit_validRight),3));
            PL_mis_validRight = PL_mis_validRight+abs(mean(tfmis_validRight./abs(tfmis_validRight),3));
            PL_hit_std_validRight = PL_hit_std_validRight+(abs(mean(tfhit_validRight./abs(tfhit_validRight),3))).^2;
            PL_mis_std_validRight = PL_mis_std_validRight+(abs(mean(tfmis_validRight./abs(tfmis_validRight),3))).^2;
            %%% Invalid Right
            ind = randperm(size(tf_invalidRight,3));
            tfhit_invalidRight = tf_invalidRight(:,:,ind(1:size(tf_invalidIncorrectRight,3)));
            tfmis_invalidRight = tf_invalidRight(:,:,ind((size(tf_invalidIncorrectRight,3)+1):end));
            PL_hit_invalidRight = PL_hit_invalidRight+abs(mean(tfhit_invalidRight./abs(tfhit_invalidRight),3));
            PL_mis_invalidRight = PL_mis_invalidRight+abs(mean(tfmis_invalidRight./abs(tfmis_invalidRight),3));
            PL_hit_std_invalidRight = PL_hit_std_invalidRight+(abs(mean(tfhit_invalidRight./abs(tfhit_invalidRight),3))).^2;
            PL_mis_std_invalidRight = PL_mis_std_invalidRight+(abs(mean(tfmis_invalidRight./abs(tfmis_invalidRight),3))).^2;
        end;
    end
    %%% Valid Left
    PL_hit_validLeft = PL_hit_validLeft./(repetition*repetition);
    PL_hit_std_validLeft = sqrt(PL_hit_std_validLeft./(repetition*repetition)-PL_hit_validLeft.^2);
    save([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_hitbootPL_validLeft.mat'],'PL_hit_validLeft','PL_hit_std_validLeft','freqs','times');
    PL_mis_validLeft = PL_mis_validLeft./(repetition*repetition);
    PL_mis_std_validLeft = sqrt(PL_mis_std_validLeft./(repetition*repetition)-PL_mis_validLeft.^2);
    save([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_misbootPL_validLeft.mat'],'PL_mis_validLeft','PL_mis_std_validLeft','freqs','times');
    %%% Invalid Left
    PL_hit_invalidLeft = PL_hit_invalidLeft./(repetition*repetition);
    PL_hit_std_invalidLeft = sqrt(PL_hit_std_invalidLeft./(repetition*repetition)-PL_hit_invalidLeft.^2);
    save([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_hitbootPL_invalidLeft.mat'],'PL_hit_invalidLeft','PL_hit_std_invalidLeft','freqs','times');
    PL_mis_invalidLeft = PL_mis_invalidLeft./(repetition*repetition);
    PL_mis_std_invalidLeft = sqrt(PL_mis_std_invalidLeft./(repetition*repetition)-PL_mis_invalidLeft.^2);
    save([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_misbootPL_invalidLeft.mat'],'PL_mis_invalidLeft','PL_mis_std_invalidLeft','freqs','times');
    %%% Valid Right
    PL_hit_validRight = PL_hit_validRight./(repetition*repetition);
    PL_hit_std_validRight = sqrt(PL_hit_std_validRight./(repetition*repetition)-PL_hit_validRight.^2);
    save([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_misbootPL_validRight.mat'],'PL_mis_validRight','PL_mis_std_validRight','freqs','times');
    PL_mis_validRight = PL_mis_validRight./(repetition*repetition);
    PL_mis_std_validRight = sqrt(PL_mis_std_validRight./(repetition*repetition)-PL_mis_validRight.^2);
    save([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_hitbootPL_validRight.mat'],'PL_hit_validRight','PL_hit_std_validRight','freqs','times');
    %%% Invalid Right
    PL_hit_invalidRight = PL_hit_invalidRight./(repetition*repetition);
    PL_hit_std_invalidRight = sqrt(PL_hit_std_invalidRight./(repetition*repetition)-PL_hit_invalidRight.^2);
    save([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_misbootPL_invalidRight.mat'],'PL_mis_invalidRight','PL_mis_std_invalidRight','freqs','times');
    PL_mis_invalidRight = PL_mis_invalidRight./(repetition*repetition);
    PL_mis_std_invalidRight = sqrt(PL_mis_std_invalidRight./(repetition*repetition)-PL_mis_invalidRight.^2);
    save([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_hitbootPL_invalidRight.mat'],'PL_hit_invalidRight','PL_hit_std_invalidRight','freqs','times');
end


%% Plot the Z-scores

exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
dirmat = [exptDir '/' sessionDir];

electrods = 1:157;
poolsumPL = zeros(50,512);
correct =zeros(50,512);
incorrect=zeros(50,512);

%%% Either
% list_cond = {'validLeft','validRight'};
%%% Or
list_cond = {'invalidLeft','invalidRight'};

PL_correct_cond1 = zeros(50,512);PL_incorrect_cond1 = zeros(50,512);PL_hit_cond1 = zeros(50,512);PL_mis_cond1 = zeros(50,512);PL_hit_cond1_std = zeros(50,512);PL_mis_cond1_std = zeros(50,512);
PL_correct_cond2 = zeros(50,512);PL_incorrect_cond2 = zeros(50,512);PL_hit_cond2 = zeros(50,512);PL_mis_cond2 = zeros(50,512);PL_hit_cond2_std = zeros(50,512);PL_mis_cond2_std = zeros(50,512);

for elec = electrods
    for cond = list_cond
        if strcmp(cond{:},{'validLeft'})
            % Load phase-locking data
            load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_validCorrectLeftsub.mat']);
            PL_correct_cond1 = PL_correct_cond1 + abs(phaseLocking);
            load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_validIncorrectLeft.mat']);
            PL_incorrect_cond1 = PL_incorrect_cond1 + abs(phaseLocking);
            % Load bootstrap data
            load([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_hitbootPL_validLeft.mat']);
            PL_hit_cond1 = PL_hit_cond1 + abs(PL_hit_validLeft);
            PL_hit_cond1_std = PL_hit_cond1_std + PL_hit_std_validLeft.^2;
            load([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_misbootPL_validLeft.mat']);
            PL_mis_cond1 = PL_mis_cond1 + abs(PL_mis_validLeft);
            PL_mis_cond1_std = PL_mis_cond1_std + PL_mis_std_validLeft.^2;
        elseif strcmp(cond{:},{'validRight'})
            % Load phase-locking data
            load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_validCorrectRightsub.mat']);
            PL_correct_cond2 = PL_correct_cond1 + abs(phaseLocking);
            load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_validIncorrectRight.mat']);
            PL_incorrect_cond2 = PL_incorrect_cond1 + abs(phaseLocking);
            % Load bootstrap data
            load([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_hitbootPL_validRight.mat']);
            PL_hit_cond2 = PL_hit_cond2 + abs(PL_hit_validRight);
            PL_hit_cond2_std = PL_hit_cond2_std + PL_hit_std_validRight.^2;
            load([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_misbootPL_validRight.mat']);
            PL_mis_cond2 = PL_mis_cond2 + abs(PL_mis_validRight);
            PL_mis_cond2_std = PL_mis_cond2_std + PL_mis_std_validRight.^2;
        elseif strcmp(cond{:},{'invalidLeft'})
            % Load phase-locking data
            load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_invalidCorrectLeftsub.mat']);
            PL_correct_cond1 = PL_correct_cond1 + abs(phaseLocking);
            load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_invalidIncorrectLeft.mat']);
            PL_incorrect_cond1 = PL_incorrect_cond1 + abs(phaseLocking);
            % Load bootstrap data
            load([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_hitbootPL_invalidLeft.mat']);
            PL_hit_cond1 = PL_hit_cond1 + abs(PL_hit_invalidLeft);
            PL_hit_cond1_std = PL_hit_cond1_std + PL_hit_std_invalidLeft.^2;
            load([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_misbootPL_invalidLeft.mat']);
            PL_mis_cond1 = PL_mis_cond1 + abs(PL_mis_invalidLeft);
            PL_mis_cond1_std = PL_mis_cond1_std + PL_mis_std_invalidLeft.^2;
        elseif strcmp(cond{:},{'invalidRight'})
            % Load phase-locking data
            load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_invalidCorrectRightsub.mat']);
            PL_correct_cond2 = PL_correct_cond1 + abs(phaseLocking);
            load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_invalidIncorrectRight.mat']);
            PL_incorrect_cond2 = PL_incorrect_cond1 + abs(phaseLocking);
            % Load bootstrap data
            load([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_hitbootPL_invalidRight.mat']);
            PL_hit_cond2 = PL_hit_cond2 + abs(PL_hit_invalidRight);
            PL_hit_cond2_std = PL_hit_cond2_std + PL_hit_std_invalidRight.^2;
            load([dirmat '/bootPhaseLocking/' obs '_elec' num2str(elec) '_misbootPL_invalidRight.mat']);
            PL_mis_cond2 = PL_mis_cond2 + abs(PL_mis_invalidRight);
            PL_mis_cond2_std = PL_mis_cond2_std + PL_mis_std_invalidRight.^2;
        end;
    end;
end

PL_correct_cond1 = PL_correct_cond1./length(electrods);
PL_incorrect_cond1 = PL_incorrect_cond1./length(electrods);
PL_correct_cond2 = PL_correct_cond2./length(electrods);
PL_incorrect_cond2 = PL_incorrect_cond2./length(electrods);

PL_hit_cond1 = PL_hit_cond1./length(electrods);
PL_mis_cond1 = PL_mis_cond1./length(electrods);
PL_hit_cond2 = PL_hit_cond2./length(electrods);
PL_mis_cond2 = PL_mis_cond2./length(electrods);

PL_hit_cond1_std = sqrt(PL_hit_cond1_std)./length(electrods);
PL_mis_cond1_std = sqrt(PL_mis_cond1_std)./length(electrods);
PL_hit_cond2_std = sqrt(PL_hit_cond2_std)./length(electrods);
PL_mis_cond2_std = sqrt(PL_mis_cond2_std)./length(electrods);

%%% Calcul of the z-score for PHASE OPPOSITION
sumPL = PL_correct_cond1 + PL_correct_cond2 + PL_incorrect_cond1 + PL_incorrect_cond2;
sumPLbootmean = PL_hit_cond1_std + PL_hit_cond2_std + PL_mis_cond1_std + PL_mis_cond2_std;
sumPLbootstd = sqrt(PL_hit_cond1_std.^2 + PL_hit_cond2_std.^2 + PL_mis_cond1_std.^2 + PL_mis_cond2_std.^2);
sumPLzOPP = (sumPL-sumPLbootmean)./sumPLbootstd;

%%% Calcul of the z-score for PHASE DIFFERENCE
sumPL = (PL_correct_cond1 + PL_correct_cond2) - (PL_incorrect_cond1 + PL_incorrect_cond2);
% sumPLzDIFF = mean(sumPL,3)./(std(sumPL,[],3)/sqrt(1));

pvalsOPP = 1 - tcdf(min(double(sumPLzOPP),8),(length(electrods)*length(obs))-1);
% pvalsDIFF = 1 - tcdf(min(double(sumPLzDIFF),8),(length(electrods)*length(obs))-1);

%% Phase-locking difference

figure;
surf(times-500,freqs(1:50),double(sumPL));
hold on;
view(0,90);
set(gcf,'Renderer','Zbuffer');
shading interp;
axis([times(1)-500 times(end)-500 freqs(1) freqs(end)]);
colorbar
set(gca,'ylim',[2 40])
% set(gca,'xlim',[-200 200])
set(gca,'clim',[-.2 .2])
title('InvalidCorrect - InvalidIncorrect','FontSize',18)
xlabel('Time from cue onset (ms)','FontSize',14)
ylabel('Frequency (Hz)','FontSize',14)
line([0 0], [0 100], [max(max(abs(double(sumPL)))) max(max(abs(double(sumPL))))], 'Color', 'k','LineWidth',2)
hold off;

namefig = sprintf('/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/tfmap-PLdiff-invalid');
print('-djpeg','-r600',namefig);

%% Phase opposition

thepvals = pvalsOPP;%(:,323:end);
[p_fdr p_masked] = fdr(thepvals,0.1);
fdr_map=zeros(size(thepvals));
if (~isempty(p_masked)); fdr_map((thepvals<p_fdr))=1; end
p = zeros(size(pvalsOPP));
% p(:,323:end) = fdr_map;

figure;
surf(times-500,freqs(1:50),-log10(pvalsOPP));%)
hold on;
view(0,90);
set(gcf,'Renderer','Zbuffer');
shading interp;
axis([times(1)-500 times(end)-500 freqs(1) freqs(end)]);
[a,b]=contour(times-500,freqs,1-sign(p*11),1,'Color',[0 0 0]);
colorbar
set(gca,'ylim',[2 40])
set(gca,'xlim',[-200 200])
% set(gca,'clim',[0.075 .08])
title('Phase opposition: InvalidCorrect/InvalidIncorrect')
line([0 0], [0 100], [max(max(abs(double(-log10(pvalsOPP))))) max(max(abs(double(-log10(pvalsOPP)))))], 'Color', 'k','LineWidth',2)
hold off;

children=get(b,'Children');
for k=1:length(children)
    set(children(k),'ZData',ones(1,length(get(children(k),'YData')))*14);
    set(children(k),'LineWidth',1.5);
end

namefig = sprintf('/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/tfmap-PLopp-invalid');
% print('-djpeg','-r600',namefig);
