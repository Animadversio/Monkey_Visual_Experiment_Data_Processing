%% Movie Selectivity Analysis 
Animal="Beto";Set_Path;
ftr = contains(ExpRecord.ephysFN,"16102020");
[meta_new,rasters_new,~,Trials_new] = loadExperiments(find(ftr),Animal);
%%
ftr = contains(ExpRecord.ephysFN,"16102020-007");
[meta_new,rasters_new,~,Trials_new] = loadExperiments(find(ftr),Animal);
% 'rasterWindow',[-50 5000]
%% This one is at fovea. 
[meta_,rasters_,lfps_,Trials_] = loadData('Beto-16102020-007','expControlFN','201016_Beto_selectivity_movie(3)', 'rasterWindow',[-250 5000]); %'sdf', 'raster') ;
%% This one is at the receptive field of Chan 5 approximately.
[meta_,rasters_,lfps_,Trials_] = loadData('Beto-16102020-006','expControlFN','201016_Beto_selectivity_movie(2)', 'rasterWindow',[-250 5000]); %'sdf', 'raster') ;
%% Load Video and inspect its properties. 
moviedir = 'N:\Stimuli\2020-GANmovies\2020-10-16-Beto';
vid = VideoReader(fullfile(ExpRecord.stimuli{ftr}, Trials_.imageName{1}+".avi"));
vidLenth = vid.Duration * 1000;
frameN = vid.NumFrames; 
fps = vid.FrameRate; 
frameTick = (1000/fps)*(0:frameN); % Onset and offset for each frame. 
%%
figure(1);
for iCh = 1:numel(meta_.spikeID)
imagesc(squeeze(rasters_(iCh,:,:))')
vline(50,"r:","video onset")
title(num2str(meta_.spikeID(iCh)))
pause
end
%%
prefchan_id = find(meta_.spikeID==Trials_.TrialRecord.User.prefChan);
figure(2);
for iCh = prefchan_id'%1:numel(meta_.spikeID)
plot(movmean(squeeze(rasters_(iCh,1:4500,:)),100,1))
vline(250 + frameTick, repmat({"b:"},1,frameN+1), compose('F%d',0:frameN))
vline(250,"r:","VIDEO ON")
vline(250 + vidLenth, "r:", "VIDEO OFF")
vline(250 + (1000/fps)*(frameN+1)/2, "r-.", "MIDPOINT")
ylabel("Firing Rate")
title(num2str(meta_.spikeID(iCh))) 
pause
end
%% Visual and Eye Signal in LFP 
figure(3); 
for iCh = Trials_.TrialRecord.User.prefChan%prefchan_id'%1:numel(meta_.spikeID)
ax=subplot('Position', [0.05,0.05,0.9,0.9]);
plot(squeeze(lfps_(iCh,1:5000,:)))
vline(250 + frameTick, repmat({"b:"},1,frameN+1), compose('F%d',0:frameN))
vline(250,"r-.","VIDEO ON")
vline(250 + vidLenth,"r-.","VIDEO OFF")
vline(250 + (1000/fps)*(frameN+1)/2, "r-.", "MIDPOINT")
ylabel("LFP")
title(num2str(iCh))
pause
end
%% Trial by Trial analysis and visualization of LFPS
iCh = Trials_.TrialRecord.User.prefChan;
for triali = 1:size(lfps_,3)
ax=subplot('Position', [0.05,0.05,0.9,0.7]);
plot(squeeze(lfps_(iCh,1:5000,triali)))
vline(250,"r:","VIDEO ON")
vline(250+3300,"r:","VIDEO OFF")
ylabel("LFP")
title(num2str(iCh))

ax2=subplot('Position', [0.05,0.77,0.9,0.10]);
plot(Trials_.eyeXY{triali}(:,:))
vline(250,"r:","VIDEO ON")
vline(250+3300,"r:","VIDEO OFF")
ylabel("EyeXY")

ax3=subplot('Position', [0.05,0.89,0.9,0.10]);
plot(Trials_.eyePupil{triali})
vline(250,"r:","VIDEO ON")
vline(250+3300,"r:","VIDEO OFF")
ylabel("EyePupil")
pause
end
