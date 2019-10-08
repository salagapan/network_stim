% Script to preprocess iEEG file during stimulation session
% ----------------------------------------------------------------------- %
% 1. Epoch the data according to WM task
% 2. Extract the stimulation segments within the data and store in a new
%    EEG dataset
% 3. Run artifact removal
% 4. Stitch the cleaned segments back in the original dataset
% 5. Remove stimulation and seizure electrodes
% 6. Re-reference to common average and save the cleaned data in a separate
%    EEG dataset
% ----------------------------------------------------------------------- %
% Dependencies:
%   EEGLab toolbox (v 14.0.0)
%   StimArtifactRemoval_ICA - custom function
% ----------------------------------------------------------------------- %
% Author: Sankar Alagapan, Frohlich Lab, UNC Chapel Hill (2018)
% ----------------------------------------------------------------------- %

clearvars
clc

setFile = '/Data/Stimulation/P1_SMS_Stimulation.set';
outFile = '/Data/Stimulation/P1_SMS_Stimulation_Epoched.set';
load('/Data/Electrode Locations/P1_seizureLectrodes.mat');

% Working memory task timing parameters (to epoch the data appropiately)
fixationDuration = 1;
memSetSize = 5;
itemDuration = 0.5;
encodingDuration = memSetSize * itemDuration;
retentionDuration = 2;

taskDuration = cueDuration + encodingDuration + retentionDuration;
stimDuration = encodingDuration;

stimElectrodes = {'E33','E34','E51','E52'}; % For P1
% P2 - {'E30','E31','E46','E47'}
% P3 - {'E47','E48','E104','E105'};

EEG = pop_loadset(setFile);
channelLabelsSet = {iEEG.chanlocs.labels}; % Channel label in SET file

% Extract trials of WM task
epochedEEG = pop_epoch(EEG,{'DIN1'},[0, taskDuration]);
epochedEEG = eeg_checkset(epochedEEG);

% Extract stimulation epochs
stimEpochEEG = pop_epoch(epochedEEG,{'DIN8'},[0, stimDuration+.05]); % adding 0.05 to capture last stim completely
stimEpochEEG = eeg_checkset(stimEpochEEG);

% Get the data from struct
data = stimEpochEEG.data;
dataClean = StimArtifactRemoval_ICA(data, stimEpochEEG.srate, subName, outFilePath);

% Store in a different EEG struct
cleanedStimEEG = stimEpochEEG;
cleanedStimEEG.data = dataClean;
cleanedStimEEG = eeg_checkset(cleanedStimEEG);

% Plug in the cleaned EEG section in the whole EEG 
Fs = epochedEEG.srate;
events = {epochedEEG.event.type};
trialStartId = find(strcmp(events,'DIN1'));
nextEvents = events(trialStartId+1);
stimEventId = find(strcmp(nextEvents,'DIN8'));

epochData = epochedEEG.data;
stimCleanData = cleanedStimEEG.data;

% Compare the electrode labels in the original dataset and align them
% accordingly
epochChannels = {epochedEEG.chanlocs.labels};
cleanStimChannels = {cleanedStimEEG.chanlocs.labels};
chMapping = zeros(length(cleanStimChannels),1);
for iChannel = 1:length(cleanStimChannels)
    chMapping(iChannel) = find(strcmp(epochChannels,cleanStimChannels{iChannel}));
end
stimChannels = setdiff(1:length(epochChannels),chMapping);
stimDurationIndex = fixationDuration * Fs + 1 : (fixationDuration + stimDuration(iSub) + 0.05)*Fs;
assert(length(stimDurationIndex) == size(stimCleanData,2));

% Store the cleaned segments in the data
cleanData = epochData;
cleanData(chMapping, stimDurationIndex, stimEventId) = stimCleanData;
cleanData(stimChannels,:,:) = 0;

cleanEEG = epochedEEG;
cleanEEG.data = cleanData;
iEEG = eeg_checkset(cleanEEG);

% Get stimulation electrode index
if ~isempty(stimElectrodes)
    stimElectrodeIndex = zeros(length(stimElectrodes),1);
    for iChannel = 1:length(subStimElectrodes)
        stimElectrodeIndex(iChannel) = find(strcmp(channelLabelsSet, stimElectrodes{iChannel}));
    end
end

rejChannel = cat(1, stimElectrodeIndex, seizureElectrodes);
% Reject the seizure electrodes
iEEG = pop_select(iEEG, 'nochannel',rejChannel);
iEEG = eeg_checkset(iEEG);

% Re-reference to common average
iEEG = pop_reref(iEEG, []);
iEEG = eeg_checkset(iEEG);

pop_saveset(iEEG, outFile);


