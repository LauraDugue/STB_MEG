function ld_epoch(analStr,trigChan,trigName,exptDir,sessionDir,preproc,matDir,fileBase)

%%% Epoch the data and sort them per trial condition

%% Setup Directory
dataDir = sprintf('%s/%s/%s', exptDir, sessionDir,preproc);
saveDir = sprintf('%s/%s/%s', exptDir, sessionDir,matDir);
filename = sprintf('%s/%s_%s.sqd', dataDir, fileBase, analStr);
savename = sprintf('%s/%s/%s_epoch_workspace.mat', saveDir);
% figDir = sprintf('%s/figures/%s', saveDir);

if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
% if ~exist(figDir,'dir')
%     mkdir(figDir)
% end

%% Setup channels
channelSets = {0:39,40:79,80:119,120:156};

% channels
badChannels = [];%[10 11 115 49 152];
highSNRChannelsL = [];%[26 60 14 92];
highSNRChannelsR = [];%[1 50 7 8];

tstart = -500;
tstop = 1500;

%% Get the data
trigMean = [];
trigData = [];
for iChSet = 1:numel(channelSets)
    allChannels = channelSets{iChSet};
    channels = setdiff(allChannels,badChannels);
    
    [trigM, triggers, Fs, trigD, trigEvents] =  rd_getData(filename, trigChan, channels, tstart, tstop);
    trigMean = cat(2,trigMean,trigM);
    trigData = cat(2,trigData,trigD);
end

nSamples = size(trigMean,1);
nChannels = size(trigMean,2);
nTrigs = size(trigMean,3);

%% Save the data

save([savename 'epochedData_' analStr '_' trigName '.mat'],'trigData','nSamples','nChannels','nTrigs','trigName');

end

