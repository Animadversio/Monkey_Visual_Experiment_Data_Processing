%% LFP background analysis 

% Find the trial starts time

% find pre-trial time window

% Analyze the band power and other features in that window! 

% Correlate the information in that window with the spiking activity in
% trial. 

%% 
meta.pathRAW = 'S:\Data-Ephys-Raw';
meta.ephysFN = 'Alfa-25062020-002';
meta.sdf = 'sdf'; % 'sdf' for convolved, 'raster' otherwise
%  [spikeChans,lfpChans,timeline,spikeID] = plxread_fullExperiment_v2(meta);
[spikeChans,lfpChans,timeline,spikeID] = plxread_fullExperiment_vcrp(meta);
%%
iCh = 3;
spectrogram(lfpChans(iCh,:), 200, 100, 50, 1000, 'Yaxis')