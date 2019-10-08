% Function to create a stacked line plot. Used in StimArtifactRemoval_ICA.m
% ----------------------------------------------------------------------- %
% Inputs:
%   signal      - nChannel x nTime matrix
%   t           - nTime x 1 (time corresponding to the data)
%   elec_labels - nChannel X 1 (Electrode labels)
%   path        - for saving output
%   stimlabel   - index of stimulation electrodes (to blank the channels)
%   siglevel    - peak to peak amplitude of signal (to set the separation
%                 between channels)
% ----------------------------------------------------------------------- %
% Author: Sankar Alagapan, Frohlich Lab, UNC Chapel Hill (2018)
% ----------------------------------------------------------------------- %
function plot_timeseries(signal_mat,t,elec_labels,stimindex,sigLevel)

if nargin<5
    sigLevel = 300;
end
[nc,ns]=size(signal_mat);

ind_1 = stimindex(1);
ind_2 = stimindex(2);
signal_mat(ind_1(1),:)=0;
signal_mat(ind_2(1),:)=0;

temp=repmat([1:nc]',1,ns);

signal_offset=signal_mat+sigLevel*temp;

figure(1)
plot(t,signal_offset)
set(gca,'YTick',sigLevel*(1:nc),'YTickLabel',elec_labels)
axis tight
if nargin>3
    title(stimindex)
end