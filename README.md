# Network Targeted Stimulation to Enhance Working Memory
Codes for analyses used in Alagapan et al. 2019 Cell Reports, Network-targeted, multi-site direct cortical stimulation 
enhances working memory by modulating phase lag of low frequency oscillations

## Scripts:
- **Preprocessing.m** - Script to run preprocessing on the stimulation iEEG data
- **ComputeFunctionalConnectivity_Baseline** - Script to compute functional connectivity measures Coherence and DeBiased Weighted
      Phase Lag Index in the baseline dataset
- **ComputeFunctionalConnectivity_Stimulation** - Script to compute functional connectivity measures Coherence and DeBiased             Weighted Phase Lag Index in the stimulation dataset

## Functions:
- **StimArtifactRemoval_ICA.m** - Function to remove electrical stimulation artifacts - Used in Preprocessing.m
- **plot_timeseries.m** - Function to create stacked line plots of multichannel data - Used in StimArtifactRemoval_ICA.m

## Dependencies:
- EEGLab Toolbox (v 14.0.0) - https://github.com/sccn/eeglab
- Chronux Toolbox - https://github.com/jsiegle/chronux
- Fieldtrip Toolbox - https://github.com/fieldtrip/fieldtrip
