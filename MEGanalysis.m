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

%% Manual Inspection through Fieldtrip

ld_manualInspection