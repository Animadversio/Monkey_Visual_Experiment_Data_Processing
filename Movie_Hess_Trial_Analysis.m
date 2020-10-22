
%% Movie Selectivity Analysis, Trial by Trial
Animal="Alfa";Set_Path;
ftr = find(contains(ExpRecord.ephysFN,"21102020"));
[meta_new,rasters_new,~,Trials_new] = loadExperiments(find(ftr),Animal);

%%
[meta_,rasters_,lfps_,Trials_] = loadData('Alfa-21102020-005','expControlFN','201021_Alfa_selectivity_movie(1)', 'rasterWindow',[-250 2500]); %'sdf', 'raster') ;
%%
[meta_,rasters_,lfps_,Trials_] = loadData('Alfa-21102020-006','expControlFN','201021_Alfa_selectivity_movie(2)', 'rasterWindow',[-250 2500]); %'sdf', 'raster') ;
%%
viewTime = Trials_.TrialRecord.User.viewTime;
cropWdw = [-250, 2500];
%% Sort the Trials
vidnms = string(Trials_.imageName);
uniqVid = unique(string(Trials_.imageName));
videoIdx = arrayfun(@(nm)find(vidnms == nm),uniqVid,'Uni',0);
%%
moviedir = 'N:\Stimuli\2020-GANmovies\2020-10-21-Alfa';
vid = VideoReader(fullfile(ExpRecord.stimuli(ftr(5)), Trials_.imageName{1}+".avi"));
vidLenth = vid.Duration * 1000;
fps = vid.FrameRate; 
frameN = int32(viewTime / (1000/fps) + 1);
frameTick = (1000/fps)*double(0:frameN); % Onset and offset for each frame. 
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
figure(2);clf;ax=subplot('Position', [0.05,0.05,0.9,0.85]);
for iCh = prefchan_id' % find(meta_.spikeID==6)' %1:numel(meta_.spikeID)
for iV = 1:numel(videoIdx)
% plot(movmean(squeeze(rasters_(iCh,1:2500,videoIdx{iV})),1,1))
cla(ax,'reset');
meanPSTH = mean(rasters_(iCh,:,videoIdx{iV}),3);
semPSTH = std(rasters_(iCh,:,videoIdx{iV}),1,3) / sqrt(numel(videoIdx{iV}));
% plot(cropWdw(1)+1:cropWdw(2), mean(rasters_(iCh,:,videoIdx{iV}),3))
shadedErrorBar(cropWdw(1)+1:cropWdw(2), meanPSTH, semPSTH, 'patchSaturation', 0.1)
vline(0 + frameTick, repmat({"b:"},1,frameN+1), compose('%d',0:frameN))
vline(0,"r:","VIDEO ON")
vline(0 + viewTime, "r:", "VIDEO OFF")
% vline(250 + (1000/fps)*(frameN+1)/2, "r-.", "MIDPOINT")
ylabel("Firing Rate");xlim(cropWdw);xlabel("Time from Onset")
title([num2str(meta_.spikeID(iCh)),strrep(uniqVid(iV),'_',' ')]) 
pause
end
end
%% Frame Triggered alignment

%% Visual and Eye Signal in LFP 
figure(3); 
for iCh = Trials_.TrialRecord.User.prefChan%prefchan_id'%1:numel(meta_.spikeID)
ax=subplot('Position', [0.05,0.05,0.9,0.85]);
plot(squeeze(lfps_(iCh,1:5000,:)))
vline(250 + frameTick, repmat({"b:"},1,frameN+1), compose('F%d',0:frameN))
vline(250,"r-.","VIDEO ON")
vline(250 + vidLenth,"r-.","VIDEO OFF")
vline(250 + (1000/fps)*(frameN+1)/2, "r-.", "MIDPOINT")
ylabel("LFP")
title(num2str(iCh))
pause
end
