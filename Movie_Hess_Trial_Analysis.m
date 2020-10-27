%% Basic analysis that could be applied to Movie -> neural response data.
%  
%  

%% Movie Selectivity Analysis, Trial by Trial
Animal="Alfa";Set_Path;setMatlabTitle("MovieStaticDynamics")
ftr = find(contains(ExpRecord.ephysFN,"21102020"));
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr(3:6),Animal);

%%
[meta_,rasters_,lfps_,Trials_] = loadData('Alfa-21102020-005','expControlFN','201021_Alfa_selectivity_movie(1)', 'rasterWindow',[-250 2500]); %'sdf', 'raster') ;
%%
[meta_,rasters_,lfps_,Trials_] = loadData('Alfa-21102020-006','expControlFN','201021_Alfa_selectivity_movie(2)', 'rasterWindow',[-250 2500]); %'sdf', 'raster') ;
%%
saveroot = "E:\OneDrive - Washington University in St. Louis\MovieDynamics";
expday = datetime(meta_.expControlFN(1:6),'InputFormat','yyMMdd');
fdrnm = compose("%s-%s-Chan%02d-2",datestr(expday,'yyyy-mm-dd'), Animal, Trials_.TrialRecord.User.prefChan);
figdir = fullfile(saveroot, fdrnm);
mkdir(figdir)
%%
viewTime = Trials_.TrialRecord.User.viewTime;
cropWdw = [-250, 2500];
unit_name_arr = generate_unit_labels(meta_.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta_.spikeID, rasters_);
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
prefchan_id = find(meta_.spikeID==Trials_.TrialRecord.User.prefChan);
figure(2);clf;ax=subplot('Position', [0.03,0.07,0.93,0.87]);
for iCh = find(meta_.spikeID==6)'%prefchan_id' %  %1:numel(meta_.spikeID)
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
ylabel("Firing Rate");xlabel("Time from Onset");xlim([-200,1600]);
title(compose("%s Chan %s %s",Animal,unit_name_arr(iCh),strrep(uniqVid(iV),'_',' '))) 
saveas(2, fullfile(figdir, compose("%s_%s_%s_psth.png",Animal,unit_name_arr(iCh),uniqVid(iV))))
pause
end
end
%%
color_seq = brewermap(numel(videoIdx), 'spectral');
%%
prefchan_id = find(meta_.spikeID==Trials_.TrialRecord.User.prefChan);
figure(2);clf;ax=subplot('Position', [0.03,0.07,0.93,0.87]);
for iCh = find(meta_.spikeID==6)'%prefchan_id' %  %1:numel(meta_.spikeID)
cla(ax,'reset');
for iV = 1:numel(videoIdx)
meanPSTH = mean(rasters_(iCh,:,videoIdx{iV}),3);
semPSTH = std(rasters_(iCh,:,videoIdx{iV}),1,3) / sqrt(numel(videoIdx{iV}));
% plot(cropWdw(1)+1:cropWdw(2), mean(rasters_(iCh,:,videoIdx{iV}),3))
shadedErrorBar(cropWdw(1)+1:cropWdw(2), meanPSTH, semPSTH, 'patchSaturation', 0.05, ...
                'lineProps', {'Color', color_seq(iV,:)})
end
vline(0 + frameTick, repmat({"b:"},1,frameN+1), compose('%d',0:frameN))
vline(0,"r:","VIDEO ON")
vline(0 + viewTime, "r:", "VIDEO OFF")
ylabel("Firing Rate");xlabel("Time from Onset");xlim([-200,1600]);
title(compose("%s Chan %s %s",Animal,unit_name_arr(iCh),strrep(uniqVid(iV),'_',' '))) 
saveas(2, fullfile(figdir, compose("%s_%s_psth_merge.png",Animal,unit_name_arr(iCh))))
end
%% Spike Density within frame duration
frameWdw = [-50,100];
frameInWdw = [-5:5]*(1000/fps);
frameInWdw = frameInWdw(frameInWdw<=frameWdw(2) & frameInWdw>=frameWdw(1));
figure(3);clf;ax=subplot('Position', [0.05,0.09,0.9,0.85]);
for iCh = prefchan_id' % find(meta_.spikeID==6)' %1:numel(meta_.spikeID)
for iV = 1:numel(videoIdx)
    frameRPSTH = cell2mat(arrayfun(@(fT)squeeze(rasters_(iCh,fT+frameWdw(1)+1:fT+frameWdw(2),videoIdx{iV}))', ...
            [int32(frameTick(1:end))-cropWdw(1)]', 'Uni', 0));
    frameRmean = mean(frameRPSTH,1);
    frameRsem = std(frameRPSTH,0, 1)/sqrt(size(frameRPSTH,1));
    cla(ax,'reset');
    plot(frameWdw(1)+1:frameWdw(2), mean(frameRPSTH,1))
    shadedErrorBar(frameWdw(1)+1:frameWdw(2), frameRmean, frameRsem, 'patchSaturation', 0.1)
    vline(frameInWdw, repmat({"b:"},1,numel(frameInWdw))) % , compose('%d',0:frameN)
    xlabel("Time to Frame Onset (ms)");ylabel("Firing Rate (mean over frames)")
    title(compose("%s Chan %s %s",Animal,unit_name_arr(iCh),strrep(uniqVid(iV),'_',' '))) 
    saveas(3, fullfile(figdir, compose("%s_%s_%s_frame_spk_density.png",Animal,unit_name_arr(iCh),uniqVid(iV))))
end
end
%% Spike Density within individual frame / Frame Triggered alignment individual frame
frameWdw = [-50,100];
frameInWdw = [-5:5]*(1000/fps);
frameInWdw = frameInWdw(frameInWdw<=frameWdw(2) & frameInWdw>=frameWdw(1));
figure;clf;ax=subplot('Position', [0.05,0.09,0.9,0.85]);
for iCh = prefchan_id' % find(meta_.spikeID==6)' %1:numel(meta_.spikeID)
for iV = 1:numel(videoIdx)
    frameRPSTH = cell2mat(arrayfun(@(fT)mean(rasters_(iCh,fT+frameWdw(1)+1:fT+frameWdw(2),videoIdx{iV}),3), ...
            [int32(frameTick(2:end))-cropWdw(1)]', 'Uni', 0));
    frameRPSTHsem = cell2mat(arrayfun(@(fT)std(rasters_(iCh,fT+frameWdw(1)+1:fT+frameWdw(2),videoIdx{iV}),0,3)/sqrt(numel(videoIdx{iV})), ...
            [int32(frameTick(2:end))-cropWdw(1)]', 'Uni', 0));
    cla(ax,'reset');hold on
    for fi = 1:size(frameRPSTH,1)
    plot(frameWdw(1)+1:frameWdw(2), frameRPSTH(fi,:))
%     shadedErrorBar(frameWdw(1)+1:frameWdw(2), frameRPSTH(fi,:), frameRPSTHsem(fi,:), 'patchSaturation', 0.1)
    end
    vline(frameInWdw, repmat({"b:"},1,numel(frameInWdw))) % , compose('%d',0:frameN)
    xlabel("Frame onset");ylabel("Firing Rate")
    pause
end
end
%% Visual and Eye Signal in LFP 
figure(4); clf; 
for iCh = Trials_.TrialRecord.User.prefChan % LFP channel and unit is equal %prefchan_id'%1:numel(meta_.spikeID)
for iV = 1:numel(videoIdx)
ax=subplot('Position', [0.05,0.12,0.9,0.81]);cla;hold on
% for iT = videoIdx{iV}'
plot(cropWdw(1)+1:cropWdw(2), squeeze(lfps_(iCh,:,videoIdx{iV})))
% end
vline(0 + frameTick, repmat({"b:"},1,frameN+1), compose('%d',0:frameN))
vline(0,"r:","VIDEO ON")
vline(0 + viewTime, "r:", "VIDEO OFF")
ylabel("LFP amplitude");xlabel("Time to Video Onset");xlim([-200,1700]);ylim([-3000,3000]);box off
title(compose("%s Chan %s %s",Animal,unit_name_arr(iCh),strrep(uniqVid(iV),'_',' '))) 
saveas(4, fullfile(figdir, compose("%s_%d_%s_LFP_trial.png",Animal,iCh,uniqVid(iV))))
end
% pause
end
%% Fourier Spectrum of PSTH. Check Periodicity
% [f, P1] = spect_ampl(meanPSTH(251:1450));
% figure;plot(f, P1);xlabel("frequency (Hz)")

for iCh = prefchan_id'% find(meta_.spikeID==6)'% %  %1:numel(meta_.spikeID)
psth_col = cell2mat(cellfun(@(idx)mean(rasters_(iCh,:,idx),3), videoIdx, 'Uni', 0));
[f, P1] = spect_ampl(psth_col(:,251:1450));
% for iV = 1:numel(videoIdx)
% meanPSTH = mean(rasters_(iCh,:,videoIdx{iV}),3);
% % semPSTH = std(rasters_(iCh,:,videoIdx{iV}),1,3) / sqrt(numel(videoIdx{iV}));
% end
end
figure;plot(f, mean(P1,1));xlabel("frequency (Hz)")
% figure;plot(f, P1');xlabel("frequency (Hz)")

function [f, P1] = spect_ampl(signal)
L = length(signal);
Fs = 1000; % sampling frequnecy
PSTHspect = fft(signal-mean(signal), L, 2);
P2 = abs(PSTHspect/L);
P1 = P2(:, 1:L/2+1);
P1(:, 2:end-1) = 2*P1(:, 2:end-1);
f = Fs*(0:(L/2))/L;
end