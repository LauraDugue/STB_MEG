% MEGanalysis.m
%
%      usage: MEGanalysis
%         by: laura
%       date: 21/05/15
%    purpose: this master program is used to run all sub-analysis program
%             for the MEG analysis

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

% validData = epochedData(:,:,sessionValid);
% save([exptDir '/' sessionDir '/' matDir '/validData_' trigName '.mat'],'validData','nSamples','nChannels','nTrigs','trigName','-v7.3');
% invalidData = epochedData(:,:,sessionInvalid);
% save([exptDir '/' sessionDir '/' matDir '/invalidData_' trigName '.mat'],'invalidData','nSamples','nChannels','nTrigs','trigName','-v7.3');
% cueOnlyData = epochedData(:,:,sessioncueOnly);
% save([exptDir '/' sessionDir '/' matDir '/cueOnlyData_' trigName '.mat'],'cueOnlyData','nSamples','nChannels','nTrigs','trigName','-v7.3');

% validCorrectData = epochedData(:,:,sessionValidCorrect);
% save([exptDir '/' sessionDir '/' matDir '/validCorrectData_' trigName '.mat'],'validCorrectData','nSamples','nChannels','nTrigs','trigName','-v7.3');
% validIncorrectData = epochedData(:,:,sessionValidIncorrect);
% save([exptDir '/' sessionDir '/' matDir '/validIncorrectData_' trigName '.mat'],'validIncorrectData','nSamples','nChannels','nTrigs','trigName','-v7.3');
% invalidCorrectData = epochedData(:,:,sessionInvalidCorrect);
% save([exptDir '/' sessionDir '/' matDir '/invalidCorrectData_' trigName '.mat'],'invalidCorrectData','nSamples','nChannels','nTrigs','trigName','-v7.3');
% invalidIncorrectData = epochedData(:,:,sessionInvalidIncorrect);
% save([exptDir '/' sessionDir '/' matDir '/invalidIncorrectData_' trigName '.mat'],'invalidIncorrectData','nSamples','nChannels','nTrigs','trigName','-v7.3');

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

load([exptDir '/' sessionDir '/' matDir '/validCorrectData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/validIncorrectData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/invalidCorrectData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/invalidIncorrectData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/cueOnlyLeftData_' trigName '.mat'])
load([exptDir '/' sessionDir '/' matDir '/cueOnlyRightData_' trigName '.mat'])

%%% make the timeFreq dir if it doesn't exist
timeFreqDir = sprintf('%s/timeFreq', [exptDir '/' sessionDir]);
if ~exist(timeFreqDir,'dir')
    mkdir(timeFreqDir)
end

for elec = 1:157
    % Valid-Correct trials
    [tf, freqs, times] = timefreq(squeeze(validCorrectData(:,elec,:)), 512, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_validCorrect.mat'],'tf','freqs','times');
    % Valid-Incorrect trials
    [tf, freqs, times] = timefreq(squeeze(validIncorrectData(:,elec,:)), 512, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_validIncorrect.mat'],'tf','freqs','times');
    % Invalid-Correct trials
    [tf, freqs, times] = timefreq(squeeze(invalidCorrectData(:,elec,:)), 512, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_invalidCorrect.mat'],'tf','freqs','times');
    % Invalid-Incorrect trials
    [tf, freqs, times] = timefreq(squeeze(invalidIncorrectData(:,elec,:)), 512, 'cycles', [1 15], 'freqs', [2 100], 'ntimesout', 512, 'freqscale', 'log', 'nfreqs', 50);
    save([timeFreqDir,'/', obs,'_elec' num2str(elec), '_invalidIncorrect.mat'],'tf','freqs','times');
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
matDir = 'matEpoch';
trigName = 'cueOnset';

%%% make the amplitude dir if it doesn't exist
amplitudeDir = sprintf('%s/amplitude', [exptDir '/' sessionDir]);
if ~exist(amplitudeDir,'dir')
    mkdir(amplitudeDir)
end

%%% run amplitude analysis
dirmat = [exptDir '/' sessionDir];
electrods = 1:157;
listCond = {'cueOnlyLeft','cueOnlyRight','validCorrect','validIncorrect','invalidCorrect','invalidIncorrect'};%

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

%% Plot amplitudes
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
matDir = 'matEpoch';
trigName = 'cueOnset';

%%% Condition to plot
cond = {'validCorrect'};
% listCond = {'cueOnlyLeft','cueOnlyRight','validCorrect','validIncorrect','invalidCorrect','invalidIncorrect'};%
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
surf(times-500,freqs(1:50),double(amplitude_mean_correct));
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
if (strcmp(cond{:},{'invalidCorrect'}) || strcmp(cond{:},{'invalidIncorrect'})...
        || strcmp(cond{:},{'validCorrect'}) || strcmp(cond{:},{'validIncorrect'})) && strcmp(attCond,'exo')
    line([160 160], [0 100], [max(max(abs(amplitude_mean_correct))) max(max(abs(amplitude_mean_correct)))],'Color', 'b','LineStyle','--','LineWidth',2)
end
line([1010 1010], [0 100], [max(max(abs(amplitude_mean_correct))) max(max(abs(amplitude_mean_correct)))],'Color', 'b','LineStyle','--','LineWidth',2)
title([cond{:} ' - all sensors'],'FontSize',16)
hold off;

namefig = sprintf(['/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/tfmap-' cond{:}]);
print('-djpeg','-r600',namefig);

%% Plot p-values
thepvals = pvals;%(:,323:end);  %%%  Selection of the p_values of interest
[p_fdr p_masked] = fdr(thepvals,0.000000000000001); %%% Alpha = 10%
fdr_map=zeros(size(thepvals));
if (~isempty(p_masked)); fdr_map(thepvals<p_fdr)=1; end
p = zeros(size(pvals));%
p(:,:) = fdr_map;
z = zeros(size((1:100)));
z(:) = 20;

figure;
surf(times-500,freqs(1:50),-log10(pvals));
hold on;
view(0,90);
set(gcf,'Renderer','Zbuffer');
shading interp;
axis([times(1)-500 times(end)-500 freqs(1) freqs(end)]);
[a,b]=contour(times-500,freqs,1-sign(p*11),1,'Color',[0 0 0]);
colorbar
% load('Ruf_colormap3.mat'); colormap(mymap);
set(gca,'clim',[0 20])
set(gca,'xlim',[-200 1200])
set(gca,'ylim',[2 100])
line([0 0], [0 100], [max(max(abs(double(zscore)))) max(max(abs(double(zscore))))], 'Color', 'k','LineWidth',2)
hold off;

children=get(b,'Children');
for k=1:length(children)
    set(children(k),'ZData',ones(1,length(get(children(k),'YData')))*14);
    set(children(k),'LineWidth',1.5);
end

%% Plot topography
exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
matDir = 'matEpoch';
trigName = 'cueOnset';

%%% Times and Freqs to plot
indfreq = [8 12];%
indtime = [400 800];

%%% Condition to plot
cond = {'cueOnlyLeft'};
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

load /Volumes/DRIVE1/DATA/laura/MEG/data_hdr.mat
%%
% cueLeft = data1_157;
cueRight = data1_157;
data = cueLeft - cueRight;

%%
figure;
fH = ssm_plotOnMesh(data, 'cueLeft - cueRight', [], data_hdr, '2d');
set(gca,'CLim',[-4000000 4000000])

namefig = sprintf('/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/topo-cuelLeft-alpha');
print('-djpeg','-r600',namefig);

%%
figure;
fH = ssm_plotOnMesh(data, cond{:}, [], data_hdr, '2d');

%% Phase-locking analysis

exptDir = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot';
obs = 'id';
attCond = 'exo';
fileBase = 'R0947_STB_4.28.15';
sessionDir = [obs '/meg/' attCond '/' fileBase];
matDir = 'matEpoch';
trigName = 'cueOnset';

%%% make the amplitude dir if it doesn't exist
phaseLockingDir = sprintf('%s/phaseLocking', [exptDir '/' sessionDir]);
if ~exist(phaseLockingDir,'dir')
    mkdir(phaseLockingDir)
end

%%% run amplitude analysis
dirmat = [exptDir '/' sessionDir];
electrods = 1:157;
listCond = {'cueOnlyLeft','cueOnlyRight','validCorrect','validIncorrect','invalidCorrect','invalidIncorrect'};%

for cond = listCond
    for elec = electrods
        % load data
        dirElec = [dirmat '/timeFreq/' obs '_elec' num2str(elec) '_' cond{:}];
        load(dirElec)
        % compute phase-locking
        phaseLocking = abs(mean(tf./abs(tf),3));
        save([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_' cond{:} '.mat'],'phaseLocking','freqs','times')
        
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
cond = {'validIncorrect'};
dirmat = [exptDir '/' sessionDir];
electrods = 1:157;
phaseLocking_mean = zeros(50,512);

for elec = electrods
    load([dirmat '/phaseLocking/' obs '_elec' num2str(elec) '_PL_' cond{:} '.mat'])
    phaseLocking_mean = phaseLocking_mean + phaseLocking;
end
phaseLocking_mean = phaseLocking_mean./length(electrods);

%% Plot amplitude
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

namefig = sprintf('/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Figures/tfmap-PL-validIncorrect');
print('-djpeg','-r600',namefig);



