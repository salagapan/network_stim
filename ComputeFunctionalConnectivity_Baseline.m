% Script to compute functional connectivity - coherence and debiased
% weighted phase lag index - from iEEG data collected during the baseline
% session of the Sternberg working memory task
% ----------------------------------------------------------------------- %
% -> The script is run on data that has been preprocessed (removal of
% stimulation channels, channels over seizure focus, non-eeg channels,
% epoched according to WM task trials and re-referenced
% to common average
% -> The data is cut into sections according to epochs within WM task
% trials and functional connectivity is computed using helper functions
% provided at the end
% ----------------------------------------------------------------------- %
% Dependencies: 
%   EEGLab toolbox (v 14.0.0)
%   Fieldtrip toolbox
% ----------------------------------------------------------------------- %
% Author: Sankar Alagapan, Frohlich Lab, UNC Chapel Hill (2018)
% ----------------------------------------------------------------------- %

clearvars
clc

% Parameters for estimation of spectra and functional connectivity
cfgSpectra.method      = 'mtmfft';
cfgSpectra.output      = 'powandcsd';
cfgSpectra.pad         = 'nextpow2' ;
cfgSpectra.taper       = 'hanning';
cfgSpectra.foilim      = [1 50];
cfgSpectra.channel     = 'all';
cfgSpectra.channelcmb  = {'all','all'};
cfgSpectra.includeauto = 1;
cfgSpectra.keeptrials  = 'yes';

cfgWPLI            = [];
cfgWPLI.method     = 'wpli_debiased';

cfgCoh = [];
cfgCoh.method = 'coh';

% Working memory task timing parameters (to epoch the data appropiately)
cueDuration = 1;
memSetSize = 5;
memItemDuration = 0.5;
encodingDuration = memSetSize * memItemDuration;
retentionDuration = 2;

setFileName = '/Data/Baseline/P1_SMS_Baseline_Epoched.set';

iEEG = pop_loadset(setFileName);

% Select the times corresponding to WM epochs
iEEG_Cue = pop_select(iEEG, 'time', [0 cueDuration]);
iEEG_Mem = pop_select(iEEG, 'time', [cueDuration cueDuration + encodingDuration]);
iEEG_Ret = pop_select(iEEG, 'time', [cueDuration + encodingDuration cueDuration + encodingDuration + retentionDuration]);

% Get functional connectivity for each epochs
cue = computeFC_Shuffle(iEEG_Cue, cfgSpectra, cfgWPLI, cfgCoh, nShuffle);
encoding = computeFC_Shuffle(iEEG_Mem, cfgSpectra, cfgWPLI, cfgCoh, nShuffle);
retention = computeFC_Shuffle(iEEG_Ret, cfgSpectra, cfgWPLI, cfgCoh, nShuffle);

% use both encoding and maintenance epoch to find frequencies of interest
wpliMat = cat(1, encoding.wpli, retention.wpli); 
f = encoding.f;
ind_thetaAlpha = f>=3 & f<=13; % Frequencies within which peaks need to be found
[pk,pkLoc] = findpeaks( nanmedian(wpliMat(:,ind_thetaAlpha),1),f(ind_thetaAlpha));
f_peak = pkLoc(pk == max(pk)); % Find the highest peak

% Get the functional connectivity values corresponding to f_peak
cue = computeFC_freq(cue, f_peak);
encoding = computeFC_freq(encoding, f_peak);
retention = computeFC_freq(retention, f_peak);

saveOutput(outFileName,cue,encoding,retention);


% ----------------------------------------------------------------------- %
% HELPER FUNCTIONS
% ----------------------------------------------------------------------- %
function funcConnStruct = computeFC_Shuffle(EEG, cfgSpectra, cfgWPLI, cfgCoh, nShuffle)
% Function to compute spectra, wpli and coherence; shuffle the data, get
% shuffled estimates and get p value from them

% Setting feedback to 'no' so that there are less outputs
cfgSpectra.feedback = 'no';
cfgWPLI.feedback = 'no';
cfgCoh.feedback = 'no';

data = convertEEGtoFT(EEG); % Convert EEG to Fieldtrip format
Spectra      = ft_freqanalysis(cfgSpectra, data); % Calculate Spectra

WPLI         = ft_connectivityanalysis(cfgWPLI, Spectra); % Calculate WPLI
wpliSpectrum = WPLI.wpli_debiasedspctrm;

Coh = ft_connectivityanalysis(cfgCoh, Spectra); % Calculate Coherence
cohSpectrum = Coh.cohspctrm;

% Initialize
wpliSpectrum_shuffle = zeros(size(wpliSpectrum,1),size(wpliSpectrum,2),nShuffle);
cohSpectrum_shuffle = zeros(size(cohSpectrum,1),size(cohSpectrum,2),nShuffle);

% Shuffle the data and calculate spectra/wpli/coherence on the shuffled
% data
for iShuffle = 1:nShuffle
    data = shuffleDataStruct(data);
    Spectra = ft_freqanalysis(cfgSpectra, data);
    
    WPLI = ft_connectivityanalysis(cfgWPLI, Spectra);
    wpliSpectrum_shuffle(:,:,iShuffle) = WPLI.wpli_debiasedspctrm;
    
    Coh = ft_connectivityanalysis(cfgCoh, Spectra);
    cohSpectrum_shuffle(:,:,iShuffle) = Coh.cohspctrm;
end

% Initialize
f = Spectra.freq;
p_wpli = zeros(size(wpliSpectrum));
p_coh = zeros(size(cohSpectrum));

% Get p-values
for iF = 1:f
    for iC = 1:size(wpliSpectrum,1)
        p_wpli(iC,iF) = computeEmpiricalPValue(wpliSpectrum(iC,iF),squeeze(wpliSpectrum_shuffle(iC,iF,:)));
        p_coh(iC,iF) = computeEmpiricalPValue(cohSpectrum(iC,iF),squeeze(cohSpectrum_shuffle(iC,iF,:)));
    end
end

% Get phaseLag
phaseLag = angle(squeeze(nanmean(Spectra.crsspctrm,1)));

% Store in struct
funcConnStruct.f = f;
funcConnStruct.spectra = squeeze(nanmean(Spectra.powspctrm,1));
funcConnStruct.coherence = cohSpectrum;
funcConnStruct.pCoherence = p_coh;
funcConnStruct.coherenceShuffle = cohSpectrum_shuffle;
funcConnStruct.phaseLag = phaseLag;
funcConnStruct.wpli = wpliSpectrum;
funcConnStruct.pWpli = p_wpli;
funcConnStruct.wpliShuffle = wpliSpectrum_shuffle;
end

function funcConnStruct = computeFC_freq(funcConnStruct, fInterest)
% Function to get adjacency matrices (coherence, wpli) in frequency band of
% interest (fInterest - 1.5, fInterest + 1.5)
f = funcConnStruct.f;
ind_fInterest = f > (fInterest - 1.5) & f < (fInterest + 1.5);

% Get estimates, shuffled estimates in freq band
coherence = nanmean(funcConnStruct.coherence(:,ind_fInterest),2);
coherenceShuffle = squeeze(nanmean(funcConnStruct.coherenceShuffle(:,ind_fInterest,:),2));

wpli = nanmean(funcConnStruct.wpli(:,ind_fInterest),2);
wpliShuffle = squeeze(nanmean(funcConnStruct.wpliShuffle(:,ind_fInterest,:),2));

pCoherence = nan(size(coherence,1),1);
pWpli = nan(size(coherence,1),1);

% compute pvalue
for i = 1:size(coherence,1)
    pCoherence(i) = computeEmpiricalPValue(coherence(i), coherenceShuffle(i,:));
    pWpli(i) = computeEmpiricalPValue(wpli(i), wpliShuffle(i,:));
end

% Store in same struct
funcConnStruct.fPeak = fInterest;
funcConnStruct.coherencefPeak = squareform(coherence);
funcConnStruct.pCoherencefPeak = squareform(pCoherence);
funcConnStruct.coherencefPeakShuffle = coherenceShuffle;
funcConnStruct.wplifPeak = squareform(wpli);
funcConnStruct.pWplifPeak = squareform(pWpli);
funcConnStruct.wplifPeakSuffle = wpliShuffle;
end

function data = convertEEGtoFT(EEG)
% Function to convert eeglab EEG structure to the data format fieldtrip
% likes. Includes code from EEGLAB_to_FieldTrip.m file from Sangtae Ahn

% FieldTrip Structure
data.fsample = EEG.srate;
data.trial = [];
data.time = [];
data.label = [];
timepoints = EEG.times/1000; % to make in Sec

for iTrial = 1:EEG.trials
    data.trial{1,iTrial} = EEG.data(:,:,iTrial);
    data.time{1,iTrial} = timepoints;
end

for iCh = 1:EEG.nbchan
    data.label{iCh} = EEG.chanlocs(iCh).labels;
end

data.label=data.label';
data.trialinfo = ones(EEG.trials,1);

% Re-define Trials
cfg=[];
cfg.begsample = 1;
cfg.endsample = length(timepoints);
data = ft_redefinetrial(cfg,data);
end

function pVal = computeEmpiricalPValue( estimate, shuffledEstimates)
% Function to get p value from shuffled data
pEstimate = length(find(shuffledEstimates > estimate)) / length(shuffledEstimates) ;
pVal = min(pEstimate, 1 - pEstimate);
end

function dataStruct = shuffleDataStruct(dataStruct)
% Function to shuffle the fieldtrip data structure

data = dataStruct.trial; % Get the data
nTrials = length(data);
nChannels = size(data{1},1);

% Reshape the data to get it ready for cell2mat cell array of {nCh, nT}
% into cell array of {1,nCh,nT}
tempData = cellfun(@(x)reshape(x,1,size(x,1),size(x,2)),data,'un',0);
tempData = cell2mat(tempData');

% Shuffle the data
dataShuffle = zeros(nTrials,nChannels,size(data{1},2));
for iChannel = 1:nChannels
    tempShuffle = randperm(nTrials);
    dataShuffle(:,iChannel,:) = tempData(tempShuffle,iChannel,:);
end
% Put it back into the dataStruct
for iTrial = 1:nTrials
    dataStruct.trial{iTrial} = squeeze(dataShuffle(iTrial,:,:));
end
end

function saveOutput(fileName,cue,encoding,retention)
% Function to save the output because saving inside parfor needs a separate
% function
fcStruct.cue = cue;
fcStruct.encoding = encoding;
fcStruct.retention = retention;
save(fileName,'-struct','fcStruct')
end
