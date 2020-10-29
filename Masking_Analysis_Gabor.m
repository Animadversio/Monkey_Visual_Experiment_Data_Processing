Animal = "Alfa"; Set_Path;
ftr = find(contains(ExpRecord.ephysFN,"Alfa-27102020"));
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr(:),Animal);
%%
rasters = rasters_new{4};
meta = meta_new{4};
Trials = Trials_new{4};
%%
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);

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
dur_msk = arrayfun(@(d) abs(duration - d)<2, [50], 'Uni', 0);
preOff_msk = cellfun(@(wdw) (pre_offset > wdw(1)) & (pre_offset < wdw(2)), {[0,30],[30,100],[100,150],[150,250],[250,inf]}, 'Uni', 0);
postOn_msk = cellfun(@(wdw) (post_onset > wdw(1)) & (post_onset < wdw(2)), {[0,30],[30,100],[100,150],[150,250],[250,inf]}, 'Uni', 0);
%%
figdir = "E:\OneDrive - Washington University in St. Louis\VisualMask_PSTH\2020-10-27-Alfa-01";
mkdir(figdir)
smth = 9;dur = 50;
for iCh = 1:76
figure(9);clf;hold on;set(9,'pos',[1000         168         586         810]);
T = tiledlayout(5,1,'Padding','compact','TileSpacing','compact');
nexttile(1)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{1}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[1,0,0,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{1})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(2)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{2}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[1,0,.5,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{2})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(3)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{3}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0.5,0,.5,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{3})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(4)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{4}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0,0,.5,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{4})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(5)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{5}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0,0,1,0.4],'LineWidth',2)
vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{5})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
xlabel("Time to Stimuli Onset")
title(T, compose("PSTH as a Function of Next Onset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(9,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.png",unit_name_arr(iCh))))
savefig(9,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.fig",unit_name_arr(iCh))))
end
%%
figdir = "E:\OneDrive - Washington University in St. Louis\VisualMask_PSTH\2020-10-27-Alfa-01";
mkdir(figdir)
smth = 9;dur = 50;
for iCh = 1:76
figure(10);clf;hold on;set(11,'pos',[1000         168         586         810]);
T = tiledlayout(5,1,'Padding','compact','TileSpacing','compact');
nexttile(1)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{1}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[1,0,0,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{1})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(2)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{2}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[1,0,.5,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{2})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(3)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{3}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0.5,0,.5,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{3})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(4)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{4}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0,0,.5,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{4})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(5)
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & postOn_msk{5}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0,0,1,0.4],'LineWidth',2)
vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{5})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
xlabel("Time to Stimuli Onset")
title(T, compose("PSTH as a Function of Next Onset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(10,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.png",unit_name_arr(iCh))))
savefig(10,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.fig",unit_name_arr(iCh))))
end
%% All Image pre offset
clrmap = brewermap(5,'RdYlBu');
smth = 9;
for iCh = 1:76
figure(11);clf;hold on;set(11,'pos',[1000         168         586         810]);
T = tiledlayout(5,1,'Padding','compact','TileSpacing','compact');
for figi = 1:5
nexttile(figi);hold on
for imgi = 1:5
msk = (preOff_msk{figi}) & dur_msk{1} & (postOn_msk{4}|postOn_msk{5}) & imgmsks{imgi}'; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[clrmap(imgi,:),0.5],'LineWidth',2)
end
vline(0.0,'k-.');vline(dur,'k:');vline(-median(pre_offset(preOff_msk{figi})),'r-.');
xlim([rasterWdw(1),rasterWdw(2)])
end
title(T, compose("PSTH as a Function of Previous Offset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(11, fullfile(figdir,compose("Pre_offset_PSTH_imgsep_Chan%s.png",unit_name_arr(iCh))))
savefig(11, fullfile(figdir,compose("Pre_offset_PSTH_imgsep_Chan%s.fig",unit_name_arr(iCh))))
pause
end
%%
smth = 9;
for iCh = 1:76
figure(12);clf;hold on;set(12,'pos',[1000         168         586         810]);
T = tiledlayout(5,1,'Padding','compact','TileSpacing','compact');
for figi = 1:5
nexttile(figi);hold on
for imgi = 1:5
msk = (preOff_msk{4} | preOff_msk{5}) & dur_msk{1} & (postOn_msk{figi}) & imgmsks{imgi}'; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[clrmap(imgi,:),0.5],'LineWidth',2)
end
vline(0.0,'k-.');vline(dur,'k:');vline(dur+median(post_onset(postOn_msk{figi})),'r-.');
xlim([rasterWdw(1),rasterWdw(2)])
end
title(T, compose("PSTH as a Function of Previous Offset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(12, fullfile(figdir,compose("Post_onset_PSTH_imgsep_Chan%s.png",unit_name_arr(iCh))))
savefig(12, fullfile(figdir,compose("Post_onset_PSTH_imgsep_Chan%s.fig",unit_name_arr(iCh))))
pause
end
