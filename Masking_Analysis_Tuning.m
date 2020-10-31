% Adapt the Masking Analysis Code 

Animal = "Alfa"; Set_Path;
ftr = find(contains(ExpRecord.ephysFN,"Alfa-29102020"));
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr(4),Animal);
%%
rasters = rasters_new{1};
meta = meta_new{1};
Trials = Trials_new{1};
wdw = meta.rasterWindow;
%%
% unit_name_arr = generate_unit_labels(meta.spikeID);
% [activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
% [activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);

imgnm_uniq = unique(Trials.imageName);
imgmsks = cellfun(@(nm)contains(Trials.imageName,nm),imgnm_uniq,'Uni',0);
%%
tOn_abs = nan(1,size(rasters,3));
tOn_abs = nan(1,size(rasters,3));
pre_offset = nan(1,size(rasters,3));
post_onset = nan(1,size(rasters,3));
duration = nan(1,size(rasters,3));
for iTr = 1:size(rasters,3)
curT0 = double(Trials.imageONtime(iTr));
iInTrial = Trials.imageInTrial(iTr);
iImgON = find(any(Trials.eventMarkers{iTr}(:,2) == (101:107),2));
iImgOFF = find(any(Trials.eventMarkers{iTr}(:,2) == (201:207),2)); % Note this may be trial start code, not image off code
% curT0 = Trials.eventMarkers{iTr}(iImgON(iInTrial),1);
tImgON = Trials.eventMarkers{iTr}(iImgON,1) - curT0;
tImgOFF = Trials.eventMarkers{iTr}(iImgOFF,1) - curT0;
trNum = Trials.trialNum(iTr);
tOn_abs(iTr) = curT0 + Trials.trialStart(trNum); % add the onset of the trial to the relative time of image to trial strat
tOff_abs(iTr) = min(tImgOFF(tImgOFF>0)) + curT0 + Trials.trialStart(trNum);
tLastOff = max(tImgOFF(tImgOFF<0));
tCurOff = min(tImgOFF(tImgOFF>0));
tNextOn = min(tImgON(tImgON>0.1));
if isempty(tNextOn), tNextOn = 600; end
if isempty(tLastOff), tLastOff = -250; end

pre_offset(iTr) = 0 - tLastOff;
post_onset(iTr) = tNextOn - tCurOff;
duration(iTr) = tCurOff;

if tNextOn - tCurOff<0, keyboard;end
end
%%
dur_msk = arrayfun(@(d) abs(duration - d)<2, [66.6], 'Uni', 0);
preOff_msk = cellfun(@(wdw) (pre_offset > wdw(1)) & (pre_offset < wdw(2)), {[0,30],[30,65],[65,100],[100,150],[150,250],[250,inf]}, 'Uni', 0);
postOn_msk = cellfun(@(wdw) (post_onset > wdw(1)) & (post_onset < wdw(2)), {[0,30],[30,65],[65,100],[100,150],[150,250],[250,inf]}, 'Uni', 0);
%%
clrseq = brewermap(numel(postOn_msk),'Spectral');
figdir = "E:\OneDrive - Washington University in St. Louis\VisualMask_PSTH\2020-10-29-Alfa-01";
mkdir(figdir)
smth = 9; dur = 66.66666; 
for iCh = 1:numel(meta.spikeID)
figure(9);clf;hold on;set(9,'pos',[1000         168         586         810]);
T = tiledlayout(numel(postOn_msk),1,'Padding','compact','TileSpacing','compact');
for subfi = 1:numel(postOn_msk)
nexttile(subfi)
msk = (preOff_msk{4} | preOff_msk{5} | preOff_msk{6}) & dur_msk{1} & postOn_msk{subfi}; 
plot(wdw(1)+1:wdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[clrseq(subfi,:),1],'LineWidth',2)
vline(0.0,'k-.');vline(dur,'k:'); vline(dur+median(post_onset(postOn_msk{subfi})),'r-.')
xlim([wdw(1),wdw(2)])
if subfi ~= numel(postOn_msk), xticklabels([]);end
end
xlabel("Time to Stimuli Onset")
title(T, compose("PSTH as a Function of Next Onset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(9,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.png",unit_name_arr(iCh))))
savefig(9,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.fig",unit_name_arr(iCh))))
end
%%
smth = 9; dur = 66.66666; 
for iCh = 1:numel(meta.spikeID)
figure(10);clf;hold on;set(11,'pos',[1000         168         586         810]);
T = tiledlayout(numel(preOff_msk),1,'Padding','compact','TileSpacing','compact');
for subfi = 1:numel(preOff_msk)
nexttile(subfi)
msk = (postOn_msk{4} | postOn_msk{5} | postOn_msk{6}) & dur_msk{1} & preOff_msk{subfi}; 
plot(wdw(1)+1:wdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[clrseq(subfi,:),1],'LineWidth',2)
vline(0.0,'k-.');vline(dur,'k:');vline( - median(pre_offset(preOff_msk{subfi})),'r-.')
xlim([wdw(1),wdw(2)])
if subfi ~= numel(preOff_msk), xticklabels([]);end
end
xlabel("Time to Stimuli Onset")
title(T, compose("PSTH as a Function of Next Onset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(10,fullfile(figdir,compose("Pre_offset_PSTH_allimg_Chan%s.png",unit_name_arr(iCh))))
savefig(10,fullfile(figdir,compose("Pre_offset_PSTH_allimg_Chan%s.fig",unit_name_arr(iCh))))
pause;
end
%% All Image plot separately pre offset
clrmap = brewermap(numel(imgmsks), 'Set3');
smth = 9;
for iCh = 1:numel(meta.spikeID)
figure(11);clf;hold on;set(11,'pos',[1000         168         586         810]);
T = tiledlayout(numel(preOff_msk),1,'Padding','compact','TileSpacing','compact');
for figi = 1:numel(preOff_msk)
nexttile(figi);hold on
for imgi = 1:numel(imgmsks)
msk = (preOff_msk{figi}) & dur_msk{1} & (postOn_msk{5} | postOn_msk{6}) & imgmsks{imgi}'; 
plot(wdw(1)+1:wdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[clrmap(imgi,:),0.5],'LineWidth',2)
end
vline(0.0,'k-.');vline(dur,'k:');vline(-median(pre_offset(preOff_msk{figi})),'r-.');
xlim([wdw(1),wdw(2)])
end
title(T, compose("PSTH as a Function of Previous Offset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(11, fullfile(figdir,compose("Pre_offset_PSTH_imgsep_Chan%s.png",unit_name_arr(iCh))))
savefig(11, fullfile(figdir,compose("Pre_offset_PSTH_imgsep_Chan%s.fig",unit_name_arr(iCh))))
% pause
end
%%
smth = 9;
for iCh = 1:numel(meta.spikeID)
figure(12);clf;hold on;set(12,'pos',[1000         168         586         810]);
T = tiledlayout(numel(postOn_msk),1,'Padding','compact','TileSpacing','compact');
for figi = 1:numel(postOn_msk)
nexttile(figi);hold on
for imgi = 1:numel(imgmsks)
msk = (preOff_msk{4} | preOff_msk{5} | preOff_msk{6}) & dur_msk{1} & (postOn_msk{figi}) & imgmsks{imgi}'; 
plot(wdw(1)+1:wdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[clrmap(imgi,:),0.5],'LineWidth',2)
end
vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{figi})),'r-.');
xlim([wdw(1),wdw(2)])
end
title(T, compose("PSTH as a Function of Previous Offset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(12, fullfile(figdir,compose("Post_onset_PSTH_imgsep_Chan%s.png",unit_name_arr(iCh))))
savefig(12, fullfile(figdir,compose("Post_onset_PSTH_imgsep_Chan%s.fig",unit_name_arr(iCh))))

end
