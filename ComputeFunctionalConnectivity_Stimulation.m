% Script to compute functional connectivity differences between conditions
% weighted phase lag index - from iEEG data collected during the
% stimulation session of the Sternberg working memory task
% ----------------------------------------------------------------------- %
% -> The script is run on data that has been preprocessed (removal of
% stimulation channels, channels over seizure focus, non-eeg channels,
% epoched according to WM task trials, followed by stimulation artifact
% removal and re-referenced to common average
% -> The data is cut into sections according to epochs within WM task
% trials and functional connectivity is computed using helper functions
% -> The difference between the functional connectivity for the different
% stimulation conditions are computed followed by a shuffle procedure to
% determine statistically significant difference in functional connectivity
% between conditions
% ----------------------------------------------------------------------- %
% Dependencies: 
%   EEGLab toolbox (v 14.0.0)
%   Fieldtrip toolbox
% ----------------------------------------------------------------------- %
% Author: Sankar Alagapan, Frohlich Lab, UNC Chapel Hill (2018)
% ----------------------------------------------------------------------- %

clearvars
clc
subSetFile = '/Data/Stimulation/P1_SMS_Stimulation_Epoched.set';

behaviorFile = '/Data/Task Performance/P1_SternbergStimulation_Summary.csv';
behaviorSummary = readtable(behaviorFile);

stimFrequency = 4;

iEEG = pop_loadset(subSetFile);
% Find the trial type from the behavior information
inphaseTrials   = find(strcmp(behaviorSummary.Condition,'In Phase'));
antiphaseTrials = find(strcmp(behaviorSummary.Condition,'Anti Phase'));
shamTrials      = find(strcmp(behaviorSummary.Condition,'Sham'));

% Compute functional connectivity measures on Stimulated segment
inPhaseEEG = pop_select(iEEG,'trial',inphaseTrials,...
    'time',timeSegment(iSegment,:));
[inPhaseWPLI, inPhaseCoherence, inPhaseSpectra] = computeFC(inPhaseEEG, cfgSpectra, cfgWPLI, cfgCoh);

antiPhaseEEG = pop_select(iEEG,'trial',antiphaseTrials,...
    'time',timeSegment(iSegment,:));
[antiPhaseWPLI, antiPhaseCoherence, antiPhaseSpectra] = computeFC(antiPhaseEEG, cfgSpectra, cfgWPLI, cfgCoh);

shamEEG  = pop_select(iEEG,'trial',shamTrials,...
    'time',timeSegment(iSegment,:));
[shamWPLI, shamCoherence, shamSpectra] = computeFC(shamEEG, cfgSpectra, cfgWPLI, cfgCoh);


f = inPhaseSpectra.freq;
freqInd = f > stimFrequency - 1.5 & f < stimFrequency + 1.5 ;

% Compute measures in stim frequency band
inPhaseSpctrm   = squeeze(nanmean(inPhaseSpectra.powspctrm,1));
antiPhaseSpctrm = squeeze(nanmean(antiPhaseSpectra.powspctrm,1));
shamSpctrm      = squeeze(nanmean(shamSpectra.powspctrm,1));

inPhasePhaseLag = squareform(angle(squeeze(nanmean(nanmean(inPhaseSpectra.crsspctrm(:,:,freqInd),3),1))));
antiPhasePhaseLag = squareform(angle(squeeze(nanmean(nanmean(antiPhaseSpectra.crsspctrm(:,:,freqInd),3),1))));
shamPhaseLag = squareform(angle(squeeze(nanmean(nanmean(shamSpectra.crsspctrm(:,:,freqInd),3),1))));

% Compute difference in spectra (difference is measured in db)
Diff_Spectra_inPhSham = 10*log10(nanmean(inPhaseSpctrm(:,freqInd),2)./nanmean(shamSpctrm(:,freqInd),2));
Diff_Spectra_inPhAntiPh = 10*log10(nanmean(inPhaseSpctrm(:,freqInd),2)./nanmean(antiPhaseSpctrm(:,freqInd),2));
Diff_Spectra_antiPhSham = 10*log10(nanmean(antiPhaseSpctrm(:,freqInd),2)./nanmean(shamSpctrm(:,freqInd),2));

% Compute difference in WPLI
inPhaseWPLI_StimFreq   = squareform(mean(inPhaseWPLI.wpli_debiasedspctrm(:,freqInd),2));
antiPhaseWPLI_StimFreq = squareform(mean(antiPhaseWPLI.wpli_debiasedspctrm(:,freqInd),2));
shamWPLI_StimFreq      = squareform(mean(shamWPLI.wpli_debiasedspctrm(:,freqInd),2));

Diff_WPLI_inPhSham   = inPhaseWPLI_StimFreq - shamWPLI_StimFreq;
Diff_WPLI_inPhAntiPh = inPhaseWPLI_StimFreq - antiPhaseWPLI_StimFreq;
Diff_WPLI_antiPhSham = antiPhaseWPLI_StimFreq - shamWPLI_StimFreq;

% Compute difference in Coherence
inPhaseCoh_StimFreq   = squareform(mean(inPhaseCoherence.cohspctrm(:,freqInd),2));
antiPhaseCoh_StimFreq = squareform(mean(antiPhaseCoherence.cohspctrm(:,freqInd),2));
shamCoh_StimFreq      = squareform(mean(shamCoherence.cohspctrm(:,freqInd),2));

Diff_Coh_inPhSham   = inPhaseCoh_StimFreq - shamCoh_StimFreq;
Diff_Coh_inPhAntiPh = inPhaseCoh_StimFreq - antiPhaseCoh_StimFreq;
Diff_Coh_antiPhSham = antiPhaseCoh_StimFreq - shamCoh_StimFreq;

% Store in Structs
% Spectra
subSpectra.ChannelLabels = channelLabels_Final;

subSpectra.inPhaseSpectra = inPhaseSpctrm;
subSpectra.antiPhaseSpectra = antiPhaseSpctrm;
subSpectra.shamSpectra = shamSpctrm;

subSpectra.DiffSpectra_InSham = Diff_Spectra_inPhSham;
subSpectra.DiffSpectra_InAnti = Diff_Spectra_inPhAntiPh;
subSpectra.DiffSpectra_AntiSham = Diff_Spectra_antiPhSham;

% WPLI
subWPLI.ChannelLabels = channelLabels_Final;

subWPLI.inPhaseWPLI   = inPhaseWPLI_StimFreq;
subWPLI.antiPhaseWPLI = antiPhaseWPLI_StimFreq;
subWPLI.shamWPLI      = shamWPLI_StimFreq;

subWPLI.DiffWPLI_InSham   = Diff_WPLI_inPhSham;
subWPLI.DiffWPLI_InAnti   = Diff_WPLI_inPhAntiPh;
subWPLI.DiffWPLI_AntiSham = Diff_WPLI_antiPhSham;

% Coherence
subCoh.ChannelLabels = channelLabels_Final;

subCoh.inPhaseCoh   = inPhaseCoh_StimFreq;
subCoh.antiPhaseCoh = antiPhaseCoh_StimFreq;
subCoh.shamCoh      = shamCoh_StimFreq;

subCoh.DiffCoh_InSham   = Diff_Coh_inPhSham;
subCoh.DiffCoh_InAnti   = Diff_Coh_inPhAntiPh;
subCoh.DiffCoh_AntiSham = Diff_Coh_antiPhSham;

nChannels = size(iEEG.data,1);
nTrials = size(iEEG.data,3);

% Shuffle the data, get all the estimates and store in a shuffled
% matrix.

inPhaseWPLI_shuffle = zeros(nShuffle, nChannels, nChannels);
antiPhaseWPLI_shuffle = zeros(nShuffle, nChannels, nChannels);
shamWPLI_shuffle = zeros(nShuffle, nChannels, nChannels);

inPhaseCoh_shuffle = zeros(nShuffle, nChannels, nChannels);
antiPhaseCoh_shuffle = zeros(nShuffle, nChannels, nChannels);
shamCoh_shuffle = zeros(nShuffle, nChannels, nChannels);


Diff_Spectra_inPhSham_shuffle   = zeros(nShuffle,nChannels);
Diff_Spectra_inPhAntiPh_shuffle = zeros(nShuffle,nChannels);
Diff_Spectra_antiPhSham_shuffle = zeros(nShuffle,nChannels);

Diff_WPLI_inPhSham_shuffle   = zeros(nShuffle, nChannels, nChannels);
Diff_WPLI_inPhAntiPh_shuffle = zeros(nShuffle, nChannels, nChannels);
Diff_WPLI_antiPhSham_shuffle = zeros(nShuffle, nChannels, nChannels);

Diff_Coh_inPhSham_shuffle   = zeros(nShuffle, nChannels, nChannels);
Diff_Coh_inPhAntiPh_shuffle = zeros(nShuffle, nChannels, nChannels);
Diff_Coh_antiPhSham_shuffle = zeros(nShuffle, nChannels, nChannels);


for iter = 1:nShuffle
    sprintf('Running %s permutations. Currently permutation %d ...\n',subName,iter)
    iEEG.data = iEEG.data(:,:,randperm(nTrials));
    iEEG = eeg_checkset(iEEG);
    
    % Compute functional connectivity measures on the stimulated data part
    inPhaseEEG = pop_select(iEEG,'trial',inphaseTrials,...
        'time',timeSegment(iSegment,:));
    [inPhaseWPLI, inPhaseCoherence, inPhaseSpectra] = computeFC(inPhaseEEG, cfgSpectra, cfgWPLI, cfgCoh);
    
    antiPhaseEEG = pop_select(iEEG,'trial',antiphaseTrials,...
        'time',timeSegment(iSegment,:));
    [antiPhaseWPLI, antiPhaseCoherence, antiPhaseSpectra] = computeFC(antiPhaseEEG, cfgSpectra, cfgWPLI, cfgCoh);
    
    shamEEG  = pop_select(iEEG,'trial',shamTrials,...
        'time',timeSegment(iSegment,:));
    [shamWPLI, shamCoherence, shamSpectra] = computeFC(shamEEG, cfgSpectra, cfgWPLI, cfgCoh);
    
    inPhaseSpctrm   = squeeze(nanmean(inPhaseSpectra.powspctrm,1));
    antiPhaseSpctrm = squeeze(nanmean(antiPhaseSpectra.powspctrm,1));
    shamSpctrm      = squeeze(nanmean(shamSpectra.powspctrm,1));
    
    Diff_Spectra_inPhSham_shuffle(iter,:,:) = 10*log10(nanmean(inPhaseSpctrm(:,freqInd),2)./nanmean(shamSpctrm(:,freqInd),2));
    Diff_Spectra_inPhAntiPh_shuffle(iter,:,:) = 10*log10(nanmean(inPhaseSpctrm(:,freqInd),2)./nanmean(antiPhaseSpctrm(:,freqInd),2));
    Diff_Spectra_antiPhSham_shuffle(iter,:,:) = 10*log10(nanmean(antiPhaseSpctrm(:,freqInd),2)./nanmean(shamSpctrm(:,freqInd),2));
    
    inPhaseWPLI_StimFreq      = squareform(mean(inPhaseWPLI.wpli_debiasedspctrm(:,freqInd),2));
    antiPhaseWPLI_StimFreq    = squareform(mean(antiPhaseWPLI.wpli_debiasedspctrm(:,freqInd),2));
    shamWPLI_StimFreq         = squareform(mean(shamWPLI.wpli_debiasedspctrm(:,freqInd),2));
    
    inPhaseWPLI_shuffle(iter,:,:) = inPhaseWPLI_StimFreq;
    antiPhaseWPLI_shuffle(iter,:,:) = antiPhaseWPLI_StimFreq;
    shamWPLI_shuffle(iter,:,:) = shamWPLI_StimFreq;
    
    Diff_WPLI_inPhSham_shuffle(iter,:,:)   = inPhaseWPLI_StimFreq - shamWPLI_StimFreq;
    Diff_WPLI_inPhAntiPh_shuffle(iter,:,:) = inPhaseWPLI_StimFreq - antiPhaseWPLI_StimFreq;
    Diff_WPLI_antiPhSham_shuffle(iter,:,:) = antiPhaseWPLI_StimFreq - shamWPLI_StimFreq;
    
    inPhaseCoh_StimFreq   = squareform(mean(inPhaseCoherence.cohspctrm(:,freqInd),2));
    antiPhaseCoh_StimFreq = squareform(mean(antiPhaseCoherence.cohspctrm(:,freqInd),2));
    shamCoh_StimFreq      = squareform(mean(shamCoherence.cohspctrm(:,freqInd),2));
    
    inPhaseCoh_shuffle(iter,:,:) = inPhaseCoh_StimFreq;
    antiPhaseCoh_shuffle(iter,:,:) = antiPhaseCoh_StimFreq;
    shamCoh_shuffle(iter,:,:) = shamCoh_StimFreq;
    
    Diff_Coh_inPhSham_shuffle(iter,:,:)   = inPhaseCoh_StimFreq - shamCoh_StimFreq;
    Diff_Coh_inPhAntiPh_shuffle(iter,:,:) = inPhaseCoh_StimFreq - antiPhaseCoh_StimFreq;
    Diff_Coh_antiPhSham_shuffle(iter,:,:) = antiPhaseCoh_StimFreq - shamCoh_StimFreq;
end

pWPLI_inPhase = nan(nChannels,nChannels);
pWPLI_antiPhase = nan(nChannels,nChannels);
pWPLI_sham = nan(nChannels,nChannels);

pCoh_inPhase = nan(nChannels,nChannels);
pCoh_antiPhase = nan(nChannels,nChannels);
pCoh_sham = nan(nChannels,nChannels);

pWPLI_inSham = nan(nChannels,nChannels);
pWPLI_antiSham = nan(nChannels,nChannels);
pWPLI_inAnti = nan(nChannels,nChannels);

pCoh_inSham = nan(nChannels,nChannels);
pCoh_antiSham = nan(nChannels,nChannels);
pCoh_inAnti = nan(nChannels,nChannels);

pSpectra_inSham = nan(nChannels,1);
pSpectra_antiSham = nan(nChannels,1);
pSpectra_inAnti = nan(nChannels,1);

for iChannel = 1:nChannels
    
    pSpectra_inSham(iChannel) = ...
        getEmpiricalPValue(Diff_Spectra_inPhSham(iChannel) , Diff_Spectra_inPhSham_shuffle(iChannel));
    pSpectra_inAnti(iChannel) = ...
        getEmpiricalPValue(Diff_Spectra_inPhAntiPh(iChannel) , Diff_Spectra_inPhAntiPh_shuffle(iChannel));
    pSpectra_antiSham(iChannel) = ...
        getEmpiricalPValue(Diff_Spectra_antiPhSham(iChannel) , Diff_Spectra_antiPhSham_shuffle(iChannel));
    
    for jChannel = iChannel+1:nChannels
        % WPLI
        pWPLI_inSham(iChannel,jChannel) = getEmpiricalPValue(Diff_WPLI_inPhSham(iChannel,jChannel),...
            squeeze(Diff_WPLI_inPhSham_shuffle(:,iChannel,jChannel)));
        pWPLI_inSham(jChannel,iChannel) = getEmpiricalPValue(Diff_WPLI_inPhSham(iChannel,jChannel),...
            squeeze(Diff_WPLI_inPhSham_shuffle(:,iChannel,jChannel)));
        
        pWPLI_inAnti(iChannel,jChannel) = getEmpiricalPValue(Diff_WPLI_inPhAntiPh(iChannel,jChannel),...
            squeeze(Diff_WPLI_inPhAntiPh_shuffle(:,iChannel,jChannel)));
        pWPLI_inAnti(jChannel,iChannel) = getEmpiricalPValue(Diff_WPLI_inPhAntiPh(iChannel,jChannel),...
            squeeze(Diff_WPLI_inPhAntiPh_shuffle(:,iChannel,jChannel)));
        
        pWPLI_antiSham(iChannel,jChannel) = getEmpiricalPValue(Diff_WPLI_antiPhSham(iChannel,jChannel),...
            squeeze(Diff_WPLI_antiPhSham_shuffle(:,iChannel,jChannel)));
        pWPLI_antiSham(jChannel,iChannel) = getEmpiricalPValue(Diff_WPLI_antiPhSham(iChannel,jChannel),...
            squeeze(Diff_WPLI_antiPhSham_shuffle(:,iChannel,jChannel)));
        
        pWPLI_inPhase(iChannel,jChannel) = getEmpiricalPValue(inPhaseWPLI_StimFreq(iChannel,jChannel),...
            squeeze(inPhaseWPLI_shuffle(:,iChannel,jChannel)));
        pWPLI_inPhase(jChannel,iChannel) = getEmpiricalPValue(inPhaseWPLI_StimFreq(iChannel,jChannel),...
            squeeze(inPhaseWPLI_shuffle(:,iChannel,jChannel)));
        
        pWPLI_antiPhase(iChannel,jChannel) = getEmpiricalPValue(antiPhaseWPLI_StimFreq(iChannel,jChannel),...
            squeeze(antiPhaseWPLI_shuffle(:,iChannel,jChannel)));
        pWPLI_antiPhase(jChannel,iChannel) = getEmpiricalPValue(antiPhaseWPLI_StimFreq(iChannel,jChannel),...
            squeeze(antiPhaseWPLI_shuffle(:,iChannel,jChannel)));
        
        pWPLI_sham(iChannel,jChannel) = getEmpiricalPValue(shamWPLI_StimFreq(iChannel,jChannel),...
            squeeze(shamWPLI_shuffle(:,iChannel,jChannel)));
        pWPLI_sham(jChannel,iChannel) = getEmpiricalPValue(shamWPLI_StimFreq(iChannel,jChannel),...
            squeeze(shamWPLI_shuffle(:,iChannel,jChannel)));
        
        % Coherence
        pCoh_inSham(iChannel,jChannel) = getEmpiricalPValue(Diff_Coh_inPhSham(iChannel,jChannel),...
            squeeze(Diff_Coh_inPhSham_shuffle(:,iChannel,jChannel)));
        pCoh_inSham(jChannel,iChannel) = getEmpiricalPValue(Diff_Coh_inPhSham(iChannel,jChannel),...
            squeeze(Diff_Coh_inPhSham_shuffle(:,iChannel,jChannel)));
        
        pCoh_inAnti(iChannel,jChannel) = getEmpiricalPValue(Diff_Coh_inPhAntiPh(iChannel,jChannel),...
            squeeze(Diff_Coh_inPhAntiPh_shuffle(:,iChannel,jChannel)));
        pCoh_inAnti(jChannel,iChannel) = getEmpiricalPValue(Diff_Coh_inPhAntiPh(iChannel,jChannel),...
            squeeze(Diff_Coh_inPhAntiPh_shuffle(:,iChannel,jChannel)));
        
        pCoh_antiSham(iChannel,jChannel) = getEmpiricalPValue(Diff_Coh_antiPhSham(iChannel,jChannel),...
            squeeze(Diff_Coh_antiPhSham_shuffle(:,iChannel,jChannel)));
        pCoh_antiSham(jChannel,iChannel) = getEmpiricalPValue(Diff_Coh_antiPhSham(iChannel,jChannel),...
            squeeze(Diff_Coh_antiPhSham_shuffle(:,iChannel,jChannel)));
        
        pCoh_inPhase(iChannel,jChannel) = getEmpiricalPValue(inPhaseCoh_StimFreq(iChannel,jChannel),...
            squeeze(inPhaseCoh_shuffle(:,iChannel,jChannel)));
        pCoh_inPhase(jChannel,iChannel) = getEmpiricalPValue(inPhaseCoh_StimFreq(iChannel,jChannel),...
            squeeze(inPhaseCoh_shuffle(:,iChannel,jChannel)));
        
        pCoh_antiPhase(iChannel,jChannel) = getEmpiricalPValue(antiPhaseCoh_StimFreq(iChannel,jChannel),...
            squeeze(antiPhaseCoh_shuffle(:,iChannel,jChannel)));
        pCoh_antiPhase(jChannel,iChannel) = getEmpiricalPValue(antiPhaseCoh_StimFreq(iChannel,jChannel),...
            squeeze(antiPhaseCoh_shuffle(:,iChannel,jChannel)));
        
        pCoh_sham(iChannel,jChannel) = getEmpiricalPValue(shamCoh_StimFreq(iChannel,jChannel),...
            squeeze(shamCoh_shuffle(:,iChannel,jChannel)));
        pCoh_sham(jChannel,iChannel) = getEmpiricalPValue(shamCoh_StimFreq(iChannel,jChannel),...
            squeeze(shamCoh_shuffle(:,iChannel,jChannel)));
        
    end
end

subWPLI.p_DiffWPLI_InSham = pWPLI_inSham;
subWPLI.p_DiffWPLI_InAnti = pWPLI_inAnti;
subWPLI.p_DiffWPLI_AntiSham = pWPLI_antiSham;

subWPLI.p_inPhase = pWPLI_inPhase;
subWPLI.p_antiPhase = pWPLI_antiPhase;
subWPLI.p_sham = pWPLI_sham;

subCoh.p_DiffCoh_InSham = pCoh_inSham;
subCoh.p_DiffCoh_InAnti = pCoh_inAnti;
subCoh.p_DiffCoh_AntiSham = pCoh_antiSham;

subCoh.p_inPhase = pCoh_inPhase;
subCoh.p_antiPhase = pCoh_antiPhase;
subCoh.p_sham = pCoh_sham;

subSpectra.p_DiffSpectra_InSham = pSpectra_inSham;
subSpectra.p_DiffSpectra_InAnti = pSpectra_inAnti;
subSpectra.p_DiffSpectra_AntiSham = pSpectra_antiSham;

subSpectra.phaseLag.inPhase = inPhasePhaseLag;
subSpectra.phaseLag.antiPhase = antiPhasePhaseLag;
subSpectra.phaseLag.sham = shamPhaseLag;

save(wpliFileCell{iSegment},'-struct','subWPLI');
save(cohFileCell{iSegment},'-struct','subCoh');
save(spectraFileCell{iSegment},'-struct','subSpectra');

% ----------------------------------------------------------------------- %
% HELPER FUNCTIONS
% ----------------------------------------------------------------------- %

function [WPLI, Coh, Spectra] = computeFC(EEG, cfgSpectra, cfgWPLI, cfgCoh)

data = convertEEGtoFT(EEG);
Spectra = ft_freqanalysis(cfgSpectra, data);
WPLI = ft_connectivityanalysis(cfgWPLI, Spectra);
Coh = ft_connectivityanalysis(cfgCoh, Spectra);
end

function pVal = getEmpiricalPValue( estimate, shuffledEstimates)
pEstimate = length(find(shuffledEstimates > estimate)) / length(shuffledEstimates) ;
pVal = min(pEstimate, 1 - pEstimate);
end
