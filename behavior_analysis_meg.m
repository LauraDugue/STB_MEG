%% Analysis MEG Program

% This program analyse the behavioral data from the MEG experiment for both
% the exoegnous and the endogenous conditions

% Input: - initials of the observer
%        - attention condition
%        - directory of the behavioral data
clear all;clc
%% INPUT: set-up parameters of the runs to analyze
obs = 'id';       % initials of the observer
attCond = 'exo'; % 'exo' or 'endo' condition
dirData = '/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Log Files'; % define the directory where the data are for the given observer and condition

%% Find all the runs
data = struct2cell(dir([dirData '/*.mat'])); % take all the stim files that finish with .mat
data = data(1,:);                            % save the name of the stim files
numRuns = 1:length(data);                    % determine the number of files to analyze

%% Determine performance for all the runs
sessionValid = [];sessionInvalid = [];
sessionValidCorrect = [];sessionValidIncorrect = [];sessionInvalidCorrect = [];sessionInvalidIncorrect = [];
sessionValidCorrectLeft = [];sessionValidIncorrectLeft = [];sessionInvalidCorrectLeft = [];sessionInvalidIncorrectLeft = [];
sessionValidCorrectRight = [];sessionValidIncorrectRight = [];sessionInvalidCorrectRight = [];sessionInvalidIncorrectRight = [];
sessioncueOnly = [];sessioncueOnlyLeft = [];sessioncueOnlyRight = [];
for iRun = numRuns
    
    load([dirData '/' data{iRun}])
    
    idx = diff(response.correct)==1;
    roughAcc = response.correct(idx) + 1;
    perf(1,iRun) = round(mean(roughAcc)*100);
    
    timeSeq = (ceil(stimulus.seqtiming(:)*100))./100;
    
    for iTrial = 1:size(stimulus.trialSeq,1)
        if iTrial < size(stimulus.trialSeq,1)
            idx = find(timeSeq<stimulus.trialsStartTimes(iTrial+1) & timeSeq>=stimulus.trialsStartTimes(iTrial));
        elseif iTrial == size(stimulus.trialSeq,1)
            idx = find(timeSeq>=stimulus.trialsStartTimes(iTrial));
        end
        
        rt = response.secs(idx);
        respWindow = find(stimulus.fixSeq(idx)==1 | stimulus.fixSeq(idx)==2);
        rt = rt(respWindow);
        rt = min(find(rt>0))*16.5;
        
        resp = response.correct(idx);
        resp2 = diff(resp)==1;
        resp = resp(resp2) + 1;
        if isempty(max(find(resp~=0)))
            answer(iTrial) = 0;
            reactTime(iTrial) = 4;
        else
            answer(iTrial) = resp(max(find(resp~=0)));
            reactTime(iTrial) = rt;
        end
    end
    
    if strcmp(attCond,'endo')
        valid = find(stimulus.trialSeq(:,1) == 1 | stimulus.trialSeq(:,1) == 2 |...
            stimulus.trialSeq(:,1) == 3 | stimulus.trialSeq(:,1) == 4 |...
            stimulus.trialSeq(:,1) == 5 | stimulus.trialSeq(:,1) == 6);
        valid_left = find((stimulus.trialSeq(:,1) == 1 | stimulus.trialSeq(:,1) == 2 |...
            stimulus.trialSeq(:,1) == 3 | stimulus.trialSeq(:,1) == 4 |...
            stimulus.trialSeq(:,1) == 5 | stimulus.trialSeq(:,1) == 6)&stimulus.trialSeq(:,2) == 1);
        valid_right = find((stimulus.trialSeq(:,1) == 1 | stimulus.trialSeq(:,1) == 2 |...
            stimulus.trialSeq(:,1) == 3 | stimulus.trialSeq(:,1) == 4 |...
            stimulus.trialSeq(:,1) == 5 | stimulus.trialSeq(:,1) == 6)&stimulus.trialSeq(:,2) == 2);
        invalid = find(stimulus.trialSeq(:,1) == 7 | stimulus.trialSeq(:,1) == 8);
        invalid_left = find((stimulus.trialSeq(:,1) == 7 | stimulus.trialSeq(:,1) == 8)&stimulus.trialSeq(:,2) == 1);
        invalid_right = find((stimulus.trialSeq(:,1) == 7 | stimulus.trialSeq(:,1) == 8)&stimulus.trialSeq(:,2) == 2);
        cue_only = find(stimulus.trialSeq(:,1) == 9 | stimulus.trialSeq(:,1) == 10);
        cue_only_left = find((stimulus.trialSeq(:,1) == 9 | stimulus.trialSeq(:,1) == 10)&stimulus.trialSeq(:,2) == 1);
        cue_only_right = find((stimulus.trialSeq(:,1) == 9 | stimulus.trialSeq(:,1) == 10)&stimulus.trialSeq(:,2) == 2);
    elseif strcmp(attCond,'exo')
        valid = find(stimulus.trialSeq(:,1) == 1 | stimulus.trialSeq(:,1) == 2 |...
            stimulus.trialSeq(:,1) == 3 | stimulus.trialSeq(:,1) == 4);
        valid_left = find((stimulus.trialSeq(:,1) == 1 | stimulus.trialSeq(:,1) == 2 |...
            stimulus.trialSeq(:,1) == 3 | stimulus.trialSeq(:,1) == 4)&stimulus.trialSeq(:,2) == 1);
        valid_right = find((stimulus.trialSeq(:,1) == 1 | stimulus.trialSeq(:,1) == 2 |...
            stimulus.trialSeq(:,1) == 3 | stimulus.trialSeq(:,1) == 4)&stimulus.trialSeq(:,2) == 2);
        invalid = find(stimulus.trialSeq(:,1) == 5 | stimulus.trialSeq(:,1) == 6 | ...
            stimulus.trialSeq(:,1) == 7 | stimulus.trialSeq(:,1) == 8);
        invalid_left = find((stimulus.trialSeq(:,1) == 5 | stimulus.trialSeq(:,1) == 6 | ...
            stimulus.trialSeq(:,1) == 7 | stimulus.trialSeq(:,1) == 8)&stimulus.trialSeq(:,2) == 1);
        invalid_right = find((stimulus.trialSeq(:,1) == 5 | stimulus.trialSeq(:,1) == 6 | ...
            stimulus.trialSeq(:,1) == 7 | stimulus.trialSeq(:,1) == 8)&stimulus.trialSeq(:,2) == 2);
        cue_only = find(stimulus.trialSeq(:,1) == 9 | stimulus.trialSeq(:,1) == 10);
        cue_only_left = find((stimulus.trialSeq(:,1) == 9 | stimulus.trialSeq(:,1) == 10)&stimulus.trialSeq(:,2) == 1);
        cue_only_right = find((stimulus.trialSeq(:,1) == 9 | stimulus.trialSeq(:,1) == 10)&stimulus.trialSeq(:,2) == 2);
    end
    
    ansValid(iRun,:) = answer(valid);
    rtValid(iRun,:) = reactTime(valid);
    ansInvalid(iRun,:) = answer(invalid);
    rtInvalid(iRun,:) = reactTime(invalid);
    
    ansValidLeft(iRun,:) = answer(valid_left);
    rtValidLeft(iRun,:) = reactTime(valid_left);
    ansInvalidLeft(iRun,:) = answer(invalid_left);
    rtInvalidLeft(iRun,:) = reactTime(invalid_left);
    
    ansValidRight(iRun,:) = answer(valid_right);
    rtValidRight(iRun,:) = reactTime(valid_right);
    ansInvalidRight(iRun,:) = answer(invalid_right);
    rtInvalidRight(iRun,:) = reactTime(invalid_right);
    
    sessionValid = [sessionValid;valid+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionInvalid = [sessionInvalid;invalid+(size(stimulus.trialSeq,1)*(iRun-1))];
    
    sessionValidCorrect = [sessionValidCorrect;valid(ansValid(iRun,:)==1)+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionValidIncorrect = [sessionValidIncorrect;valid(ansValid(iRun,:)==0)+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionInvalidCorrect = [sessionInvalidCorrect;invalid(ansInvalid(iRun,:)==1)+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionInvalidIncorrect = [sessionInvalidIncorrect;invalid(ansInvalid(iRun,:)==0)+(size(stimulus.trialSeq,1)*(iRun-1))];
    
    sessionValidCorrectLeft = [sessionValidCorrectLeft;valid_left(ansValidLeft(iRun,:)==1)+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionValidIncorrectLeft = [sessionValidIncorrectLeft;valid_left(ansValidLeft(iRun,:)==0)+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionInvalidCorrectLeft = [sessionInvalidCorrectLeft;invalid_left(ansInvalidLeft(iRun,:)==1)+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionInvalidIncorrectLeft = [sessionInvalidIncorrectLeft;invalid_left(ansInvalidLeft(iRun,:)==0)+(size(stimulus.trialSeq,1)*(iRun-1))];
    
    sessionValidCorrectRight = [sessionValidCorrectRight;valid_right(ansValidRight(iRun,:)==1)+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionValidIncorrectRight = [sessionValidIncorrectRight;valid_right(ansValidRight(iRun,:)==0)+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionInvalidCorrectRight = [sessionInvalidCorrectRight;invalid_right(ansInvalidRight(iRun,:)==1)+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessionInvalidIncorrectRight = [sessionInvalidIncorrectRight;invalid_right(ansInvalidRight(iRun,:)==0)+(size(stimulus.trialSeq,1)*(iRun-1))];
    
    sessioncueOnly = [sessioncueOnly;cue_only+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessioncueOnlyLeft = [sessioncueOnlyLeft;cue_only_left+(size(stimulus.trialSeq,1)*(iRun-1))];
    sessioncueOnlyRight = [sessioncueOnlyRight;cue_only_right+(size(stimulus.trialSeq,1)*(iRun-1))];
end    
save('/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Log Files/Conditions/indexBehavior.mat',...
    'sessionValid','sessionInvalid','sessioncueOnly','sessionValidCorrect','sessionValidIncorrect','sessionInvalidCorrect','sessionInvalidIncorrect',...
    'sessioncueOnlyLeft','sessioncueOnlyRight','sessionValidCorrectLeft','sessionValidIncorrectLeft','sessionInvalidCorrectLeft','sessionInvalidIncorrectLeft',...
    'sessionValidCorrectRight','sessionValidIncorrectRight','sessionInvalidCorrectRight','sessionInvalidIncorrectRight')
%% Create averaged performance for Valid and invalid conditions for each run
avgPerfValid = (sum(ansValid,2) ./ size(ansValid, 2));
avgPerfInvalid = (sum(ansInvalid,2) ./ size(ansInvalid, 2));

%% Create averaged performance for Valid and invalid conditions averaged over run
avgValid = mean(avgPerfValid);
avgInvalid = mean(avgPerfInvalid);

%% Plot the performance

figure;
bar([avgValid avgInvalid]*100)
xlim([0 3])
title([obs ' (n = ' num2str(size(numRuns,2)) ' runs)'],'FontSize', 16, 'FontWeight','bold')
ylabel('Percent correct','FontSize', 14)
xlabel('Cueing condition','FontSize', 14)
ylim([40 100])
set(gca,'yTick',40:20:100)
set(gca,'xTick',1:2)
set(gca,'xTickLabel',{'Valid' 'Invalid'},'FontSize', 14)

namefig=sprintf(['/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Log Files/Figures/performance_' obs '_' attCond]);
print ('-djpeg', '-r500',namefig);

%% Average reaction time for valid and invalid (no 4)
greaterThanFourValid = find(rtValid > 4);
avgRtValid = mean(rtValid(greaterThanFourValid));

greaterThanFourInvalid = find(rtInvalid > 4);
avgRtInvalid  = mean(rtInvalid(greaterThanFourInvalid));

%% Plot the reaction time
figure;
bar([avgRtValid avgRtInvalid])
xlim([0 3])
ylim([0 500])
title([obs ' (n = ' num2str(size(numRuns,2)) ' runs)'],'FontSize', 16, 'FontWeight','bold')
ylabel('Reaction Time (ms)','FontSize', 14)
xlabel('Cueing condition','FontSize', 14)
set(gca,'yTick',0:100:400)
set(gca,'xTick',1:2)
set(gca,'xTickLabel',{'Valid' 'Invalid'},'FontSize', 14)

namefig=sprintf(['/Volumes/DRIVE1/DATA/laura/MEG/Pilot/id/meg/exo/R0947_STB_4.28.15/Log Files/Figures/reactTime_' obs '_' attCond]);
print ('-djpeg', '-r500',namefig);



 