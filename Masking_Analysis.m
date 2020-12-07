% Visual Masking 
searchstr = string([char(Animal),'-',datestr(datetime(),"ddmmyyyy"),'*.pl2']);
ephysFNs = string(ls("N:\Data-Ephys-Raw\"+searchstr));
% batchProcessPL2(ephysFNs)
%%
Animal = "Alfa"; Set_Path;
ftr = find(contains(ExpRecord.ephysFN,"Alfa-29102020"));
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr(:),Animal);
%%
rasters = rasters_new{9};
meta = meta_new{9};
Trials = Trials_new{9};
%%
figure(4);
Chans = 1:76;
for iTr = 1:size(rasters,3)
curT0 = double(Trials.imageONtime(iTr));
iInTrial = Trials.imageInTrial(iTr);
iImgON = find(any(Trials.eventMarkers{iTr}(:,2) == (101:107),2));
iImgOFF = find(any(Trials.eventMarkers{iTr}(:,2) == (201:207),2)); % Note this may be trial start code, not image off code
% curT0 = Trials.eventMarkers{iTr}(iImgON(iInTrial),1);
tImgON = Trials.eventMarkers{iTr}(iImgON,1) - curT0;
tImgOFF = Trials.eventMarkers{iTr}(iImgOFF,1) - curT0;
hold off
imagesc(-249:600,meta.spikeID(Chans),rasters(Chans,:,iTr));hold on
vline(tImgON, '-.r')
vline(tImgOFF, '-.yellow')
pause;
end
%%
pre_offset = nan(1,size(rasters,3));
post_onset = nan(1,size(rasters,3));
pre_offset(2:end) = Trials.imageONtime(2:end) - Trials.imageOFFtime(1:end-1);
post_onset(1:end-1) = Trials.imageONtime(2:end) - Trials.imageOFFtime(1:end-1);
%% Get the delay to next image onset and the distance from last exp..
tOn_abs = zeros(1,size(rasters,3));
tOff_abs = zeros(1,size(rasters,3));
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
tOn_abs(iTr) = curT0;
tOff_abs(iTr) = min(tImgOFF(tImgOFF>0)) + curT0;
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
rasterWdw = meta.rasterWindow;
%%
imgnm_uniq = unique(Trials.imageName);
imgmsk = contains(Trials.imageName, imgnm_uniq{1});
%%
dur_msk = arrayfun(@(d) abs(duration - d)<2, [50, 66.6], 'Uni', 0);
preOff_msk = cellfun(@(wdw) (pre_offset > wdw(1)) & (pre_offset < wdw(2)), {[0,30],[30,100],[100,150],[150,250],[250,inf]}, 'Uni', 0);
postOn_msk = cellfun(@(wdw) (post_onset > wdw(1)) & (post_onset < wdw(2)), {[0,30],[30,100],[100,150],[150,250],[250,inf]}, 'Uni', 0);
%%
iCh = 10;
smth = 9;
figure(6);clf;hold on
% msk = (preOff_msk{3} | preOff_msk{4}) & dur_msk{1}; 
% plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),7),'Color',[0,0,0,0.4],'LineWidth',2)
msk = (preOff_msk{3} | preOff_msk{4}) & dur_msk{1} & postOn_msk{1}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[1,0,0,0.4],'LineWidth',2)
msk = (preOff_msk{3} | preOff_msk{4}) & dur_msk{1} & postOn_msk{2}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[1,0,.5,0.4],'LineWidth',2)
msk = (preOff_msk{3} | preOff_msk{4}) & dur_msk{1} & postOn_msk{3}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0.5,0,.5,0.4],'LineWidth',2)
msk = (preOff_msk{3} | preOff_msk{4}) & dur_msk{1} & postOn_msk{4}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0,0,.5,0.4],'LineWidth',2)
msk = (preOff_msk{3} | preOff_msk{4}) & dur_msk{1} & postOn_msk{5}; 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0,0,1,0.4],'LineWidth',2)
%%
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
%%
figdir = "E:\OneDrive - Washington University in St. Louis\VisualMask_PSTH\2020-10-27-Alfa-02";

smth = 9;dur = 50;
for iCh = 1:76
figure(7);clf;hold on;set(7,'pos',[1000         168         586         810]);
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
saveas(7,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.png",unit_name_arr(iCh))))
savefig(7,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.fig",unit_name_arr(iCh))))
end


%%

smth = 9;
for iCh = 1:76
figure(8);clf;hold on;set(8,'pos',[1000         168         586         810]);
T = tiledlayout(5,1,'Padding','compact','TileSpacing','compact');
nexttile(1)
msk = (preOff_msk{1}) & dur_msk{1} & (postOn_msk{4}|postOn_msk{5}); 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[1,0,0,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(-median(pre_offset(preOff_msk{1})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(2)
msk = (preOff_msk{2}) & dur_msk{1} & (postOn_msk{4}|postOn_msk{5}); 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[1,0,.5,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(-median(pre_offset(preOff_msk{2})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(3)
msk = (preOff_msk{3}) & dur_msk{1} & (postOn_msk{4}|postOn_msk{5}); 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0.5,0,.5,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(-median(pre_offset(preOff_msk{3})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(4)
msk = (preOff_msk{4}) & dur_msk{1} & (postOn_msk{4}|postOn_msk{5}); 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0,0,.5,0.4],'LineWidth',2)
xticklabels([]);vline(0.0,'k-.');vline(-median(pre_offset(preOff_msk{4})),'r-.')
xlim([rasterWdw(1),rasterWdw(2)])
nexttile(5)
msk = (preOff_msk{5}) & dur_msk{1} & (postOn_msk{4}|postOn_msk{5}); 
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),smth),'Color',[0,0,1,0.4],'LineWidth',2)
vline(0.0,'k-.');vline(-median(pre_offset(preOff_msk{5})),'r-.');
xlim([rasterWdw(1),rasterWdw(2)])
title(T, compose("PSTH as a Function of Previous Offset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(8, fullfile(figdir,compose("Pre_offset_PSTH_all_Chan%s.png",unit_name_arr(iCh))))
savefig(8, fullfile(figdir,compose("Pre_offset_PSTH_all_Chan%s.fig",unit_name_arr(iCh))))
end
%% 



%%
iCh = 11;
figure(6);clf;hold on
for iImg = 1:3
imgmsk = contains(Trials.imageName', imgnm_uniq{iImg}) & abs(duration - 50)<1;
% msk = imgmsk' & (post_onset < 30) & (pre_offset<30);
% plot(rasterWdw(1)+1:rasterWdw(2), mean(rasters(iCh,:,msk),3),'r','LineWidth',2)
msk = imgmsk & (post_onset > 30) & (post_onset < 100) & (pre_offset > 30) & (pre_offset < 100);
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),9),'Color',[0,0,1,0.4],'LineWidth',2)
msk = imgmsk & (post_onset > 150) & (pre_offset > 150);
plot(rasterWdw(1)+1:rasterWdw(2), movmean(mean(rasters(iCh,:,msk),3),9),'Color',[0,0,0,0.4],'LineWidth',2)
end
%%
iCh = 11;
figure(7);clf;hold on
msk = imgmsk' & (post_onset < 30) & (pre_offset<30);
plot(rasterWdw(1)+1:rasterWdw(2), mean(rasters(iCh,:,msk),3))
msk = imgmsk' & (post_onset > 30) & (post_onset < 100) & (pre_offset > 30) & (pre_offset < 100);
plot(rasterWdw(1)+1:rasterWdw(2), mean(rasters(iCh,:,msk),3))
msk = imgmsk' & (post_onset > 150) & (pre_offset > 150);
plot(rasterWdw(1)+1:rasterWdw(2), mean(rasters(iCh,:,msk),3))
