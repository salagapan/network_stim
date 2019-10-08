% Function to remove stimulation artifacts using Independent Component
% Analysis
% ----------------------------------------------------------------------- %
% Inputs: 
%   data        - nChannel x nTime x nTrial matrix (typically from EEG struct)
%   Fs          - sampling rate
%   fileName    - for saving output in mat file
%   path        - for saving output
% Outputs:
%   dataClean   - nChannel x nTime x nTrial matrix
% Dependencies:
%   runica      - EEGLab toolbox (v 14.0.0)
%   mtspectrumc - Chronux toolbox
%   plot_timeseries - Custom function
% ----------------------------------------------------------------------- %
% Author: Sankar Alagapan, Frohlich Lab, UNC Chapel Hill (2018)
% ----------------------------------------------------------------------- %
function dataClean = StimArtifactRemoval_ICA(data, Fs, fileName, path)

% Parameters for computing spectra
params.Fs = Fs;
params.tapers = [3,5];

% Initialzing variables
[nChannels, nT, nTrials] = size(data);
dataClean = zeros(nChannels, nT, nTrials);
components = cell(nTrials,1);
componentsRejected = cell(nTrials,1);
icawinv = cell(nTrials,1);

% Iterate over trials
for iTrial = 1:nTrials
    
    tempData = squeeze(data(:,:,iTrial));    
    t = linspace(0,nT/Fs,nT);
    % Run ICA
    [weights,sphere,~,~,~,~,activations] = runica(tempData,'PCA',rank(tempData),'maxsteps',512);
    winv = pinv(weights*sphere);
    
    activationsFilt = activations;
    [S,f] = mtspectrumc(activations',params);
    
    % Plot data, component timeseries, component spectra
    h1 = figure(1);
    set(h1,'units','normalized','outerposition',[0 0 1 1]);
    plot_timeseries(tempData,t,1:length(activations),[],300)
    xlabel('Time [s]')
    ylabel('Channel #')
    title('Before ICA')
    
    h2 = figure(2);
    set(h2,'units','normalized','outerposition',[0 0 1 1]);
    plot_timeseries(activations,t,1:length(activations),[],20)
    xlabel('Time [s]')
    ylabel('Component #')
    title('Waveform of Components')
    
    h3 = figure(3);    
    set(h3,'units','normalized','outerposition',[0 0 1 1]);
    plot_timeseries(S',f,1:size(activations,1),[],0.5);
    xlabel('Time [s]')
    ylabel('Component #')
    title('Spectra of Components')    
    
    % Store the index of components in cell array for next step
    cellArray = cell(size(activations,1),1);
    for i = 1:size(activations,1)
        cellArray{i} = int2str(i);
    end
        
    pause;
    % Get manual input using a list dialog 
    [rejList,flag] = listdlg('PromptString','Select the components to reject:',...
                    'SelectionMode','multiple','ListString',cellArray);
    
    % If components were selected, set corresponding components to zero
    if flag == 1
        activationsFilt(rejList,:) = 0;
    end
    
    % Reconstruct data
    tempDataClean = winv*activationsFilt;
    
    h4 = figure(3);    
    set(h4,'units','normalized','outerposition',[0 0 1 1]);
    plot_timeseries(tempDataClean,t,1:size(tempDataClean,1),[],300)
    xlabel('Time [s]')
    ylabel('Channel #')
    title('After Artifact Suppression')    
    
    % Store in output mat
    dataClean(:,:,iTrial) = tempDataClean';
    
    components{iTrial} = activations;
    componentsRejected{iTrial} = rejList;
    icawinv{iTrial} = winv;
    
    % save output
    save(strcat(path,fileName),'components','componentsRejected','icawinv');
    
end
end