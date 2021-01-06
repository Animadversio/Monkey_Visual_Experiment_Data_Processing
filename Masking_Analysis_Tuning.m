Animal = "Alfa"; Set_Path;
ftr = find(contains(ExpRecord.ephysFN,["Alfa-29122020","Alfa-30122020"]));  % "Alfa-06112020" "Alfa-04112020"
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr(:),Animal,false);
%%
clearvars -except Trials* meta* rasters* ExpRecord*
%%
% Adapt the Masking Analysis Code 
Animal = "Alfa"; Set_Path;
ftr = find(contains(ExpRecord.Exp_collection,'Masking'));  % "Alfa-06112020" "Alfa-04112020"
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr(:),Animal);
%%
% resultroot = "OneDrive - Washington University in St. Louis";
figroot = "E:\OneDrive - Washington University in St. Louis\VisualMask_PSTH";
flag.plot_chan = true;
flag.plot_imgsep = false;
flag.plot_tuning = false;
for Triali = 7:numel(meta_new)
rasters = rasters_new{Triali};
meta = meta_new{Triali};
Trials = Trials_new{Triali};
wdw = meta.rasterWindow;
fprintf("Experiment: %s, %s:\n%s\n",meta.ephysFN,meta.expControlFN,meta.comments)
% naming convention for the thing.
fdrnm = figdir_naming(meta, ExpRecord);
figdir = fullfile(figroot, fdrnm);
if exist(figdir,'dir')
    fprintf("Is %s OK?", fdrnm)
%     keyboard
end
mkdir(figdir)
%%
% unit_name_arr = generate_unit_labels(meta.spikeID);
% [activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
imgnm_uniq = unique(Trials.imageName);
imgmsks = cellfun(@(nm)contains(Trials.imageName,nm),imgnm_uniq,'Uni',0);
% Parse the time interval between images trial by trial
%% Parse the delays and ISIs. And print them out for better examination. 
[pre_ISI, post_ISI, duration, preISI_msk, postISI_msk, dur_msk, preISI_arr, postISI_arr, dur_arr] = parser_stimuli_timing(Trials);
fprintf("Experiment Configuration:\n")
fprintf("Pre-Stimuli ISIs (med of category): %s \n",join(compose("%.1f(%d) ",preISI_arr',cellfun(@sum,preISI_msk)')));
fprintf("Post-Stimuli ISIs (med of category): %s \n",join(compose("%.1f(%d) ",postISI_arr',cellfun(@sum,postISI_msk)')));
fprintf("Stimuli Duration (med of category): %s \n",join(compose("%.1f(%d) ",dur_arr',cellfun(@sum,dur_msk)')));
% warning(numel(dur_arr)==1, "Multiple stimuli duration detected, examine your code!\n")
%% Collect Statics for this 
set(groot,'defaultAxesTickLabelInterpreter','none'); 
if flag.plot_tuning
figure(2); set(2,'pos',[1092, 339, 960, 480]);
figure(3); set(3,'pos',[125, 339, 960, 480]);
figure(4); set(4,'pos',[125, 339, 960, 480]);
nonoverlap_bnd = 150;
preISI_far_msk = any(cat(1,preISI_msk{preISI_arr > nonoverlap_bnd}),1);
stats_all = [];
for chid = 1:numel(meta.spikeID)
    stats = []; stats_bsl = []; stats_del = [];
    respmat = [];respsem = [];
    respmat_bsl = [];respsem_bsl = [];
    respmat_del = [];respsem_del = [];    
    for postmsk = postISI_msk
    resp_dist_col = cellfun(@(msk) squeeze(mean(rasters(chid, 250+[100:200], msk' & postmsk{1} & preISI_far_msk),2)), imgmsks, 'uni',0);
    resp_bsl_dist_col = cellfun(@(msk) squeeze(mean(rasters(chid, 250+[-50:50], msk' & postmsk{1} & preISI_far_msk),2)), imgmsks, 'uni',0);
    resp_del_dist_col = cellfun(@(msk) squeeze(mean(rasters(chid, 250+[250:350], msk' & postmsk{1} & preISI_far_msk),2)), imgmsks, 'uni',0);
    resp_col = cellfun(@mean, resp_dist_col);
    resp_sem_col = cellfun(@(rsps) std(rsps)/sqrt(numel(rsps)), resp_dist_col);
    resp_bsl_col = cellfun(@mean, resp_bsl_dist_col);
    resp_bsl_sem_col = cellfun(@(rsps) std(rsps)/sqrt(numel(rsps)), resp_bsl_dist_col);
    resp_del_col = cellfun(@mean, resp_del_dist_col);
    resp_del_sem_col = cellfun(@(rsps) std(rsps)/sqrt(numel(rsps)), resp_del_dist_col);
    stats = [stats, anova_cells(resp_dist_col)];
    stats_bsl = [stats_bsl, anova_cells(resp_bsl_dist_col)];
    stats_del = [stats_del, anova_cells(resp_del_dist_col)];
    respmat = [respmat, resp_col];
    respsem = [respsem, resp_sem_col];
    respmat_bsl = [respmat_bsl, resp_bsl_col];
    respsem_bsl = [respsem_bsl, resp_bsl_sem_col];
    respmat_del = [respmat_del, resp_del_col];
    respsem_del = [respsem_del, resp_del_sem_col];
    end
    stats_all = [stats_all; stats];
    figure(2); clf;hold on
    for i = 1:numel(postISI_arr)
    errorbar(respmat(:,i), respsem(:,i),'LineWidth',1.5)
    end
    xticks(1:numel(imgnm_uniq));xlim([0.5,numel(imgnm_uniq)+0.5])
    xticklabels(imgnm_uniq);xtickangle(50)
    title([unit_name_arr(chid)+" response period 100,200 ms","ANOVA F="+join(compose("%.2f ",struct2table(stats).F'))])
    legend(compose("postISI %.1fms",postISI_arr))
    box off
    figure(3); clf;hold on
    for i = 1:numel(postISI_arr)
    errorbar(respmat_bsl(:,i), respsem_bsl(:,i),'LineWidth',1.5)
    end
    xticks(1:numel(imgnm_uniq));xlim([0.5,numel(imgnm_uniq)+0.5])
    xticklabels(imgnm_uniq);xtickangle(50);box off
    title([unit_name_arr(chid)+" baseline period -50,50 ms","ANOVA F="+join(compose("%.2f ",struct2table(stats_bsl).F'))])
    legend(compose("postISI %.1fms",postISI_arr))
    figure(4); clf;hold on
    for i = 1:numel(postISI_arr)
    errorbar(respmat_del(:,i), respsem_del(:,i),'LineWidth',1.5)
    end
    xticks(1:numel(imgnm_uniq));xlim([0.5,numel(imgnm_uniq)+0.5])
    xticklabels(imgnm_uniq);xtickangle(50); box off
    title([unit_name_arr(chid)+" baseline period 250,350 ms","ANOVA F="+join(compose("%.2f ",struct2table(stats_del).F'))])
    legend(compose("postISI %.1fms",postISI_arr))
    align_figures_ylim([2,3,4])
    saveas(2, fullfile(figdir,compose("tuningcurv_postISI_%s.png", unit_name_arr(chid))))
    saveas(3, fullfile(figdir,compose("tuningcurv_postISI_bsl_%s.png", unit_name_arr(chid))))
    saveas(4, fullfile(figdir,compose("tuningcurv_postISI_delay_%s.png", unit_name_arr(chid))))
end
figure(5);imagesc(arrayfun(@(S)S.F, stats_all)) % Summary Heatmap
yticks(1:numel(meta.spikeID));yticklabels(unit_name_arr)
%% 
nonoverlap_bnd = 150;
postISI_far_msk = any(cat(1,postISI_msk{postISI_arr > nonoverlap_bnd}),1);
stats_all = [];
for chid = 1:numel(meta.spikeID)
    stats = [];
    stats_bsl = [];
    respmat = [];
    respsem = [];
    respmat_bsl = [];
    respsem_bsl = [];
    for premsk = preISI_msk
    resp_dist_col = cellfun(@(msk) squeeze(mean(rasters(chid, 250+[100:200], msk' & premsk{1} & postISI_far_msk),2)), imgmsks, 'uni',0);
    resp_bsl_dist_col = cellfun(@(msk) squeeze(mean(rasters(chid, 250+[-50:50], msk' & premsk{1} & postISI_far_msk),2)), imgmsks, 'uni',0);
    resp_col = cellfun(@(msk) squeeze(mean(rasters(chid, 250+[100:200], msk' & premsk{1} & postISI_far_msk),[2,3])), imgmsks);
    resp_sem_col = cellfun(@(msk) squeeze(std(mean(rasters(chid, 250+[100:200], msk' & premsk{1} & postISI_far_msk),[2]),1,3)/sqrt(sum(msk' & premsk{1} & postISI_far_msk))), imgmsks);
    resp_bsl_col = cellfun(@(msk) squeeze(mean(rasters(chid, 250+[-50:50], msk' & premsk{1} & postISI_far_msk),[2,3])), imgmsks);
    resp_bsl_sem_col = cellfun(@(msk) squeeze(std(mean(rasters(chid, 250+[-50:50], msk' & premsk{1} & postISI_far_msk),[2]),1,3))/sqrt(sum(msk' & premsk{1} & postISI_far_msk)), imgmsks);
    stats = [stats, anova_cells(resp_dist_col)];
    stats_bsl = [stats_bsl, anova_cells(resp_bsl_dist_col)];
    respmat = [respmat, resp_col];
    respsem = [respsem, resp_sem_col];
    respmat_bsl = [respmat_bsl, resp_bsl_col];
    respsem_bsl = [respsem_bsl, resp_bsl_sem_col];
    end
    stats_all = [stats_all; stats]; 
    
    figure(7); set(7,'pos',[1092, 339, 960, 480]);clf;hold on
    for i = 1:numel(preISI_msk)
    errorbar(respmat(:,i), respsem(:,i),'LineWidth',1.5)
    end
    xticks(1:numel(imgnm_uniq));xlim([0.5,numel(imgnm_uniq)+0.5])
    xticklabels(imgnm_uniq);xtickangle(50)
    title([unit_name_arr(chid)+" response period 100, 200 ms","ANOVA F="+join(compose("%.2f ",struct2table(stats).F'))])
    legend(compose("preISI %.1fms",preISI_arr))
    box off
    figure(8); set(8,'pos',[1092, 339, 960, 480]);clf;hold on
    %plot(respmat,'LineWidth',1.5)
    for i = 1:numel(preISI_msk)
    errorbar(respmat_bsl(:,i), respsem_bsl(:,i),'LineWidth',1.5)
    end
    xticks(1:numel(imgnm_uniq));xlim([0.5,numel(imgnm_uniq)+0.5])
    xticklabels(imgnm_uniq);xtickangle(50)
    title([unit_name_arr(chid)+" baseline period -50,50 ms","ANOVA F="+join(compose("%.2f ",struct2table(stats_bsl).F'))])
    legend(compose("preISI %.1fms",preISI_arr))
    box off
    align_figures_ylim([7,8])
    saveas(7, fullfile(figdir,compose("tuningcurv_preISI_%s.png", unit_name_arr(chid))))
    saveas(8, fullfile(figdir,compose("tuningcurv_preISI_bsl_%s.png", unit_name_arr(chid))))
end
figure(6); imagesc(arrayfun(@(S)S.F, stats_all));colorbar() % Summary Heatmap
yticks(1:numel(meta.spikeID));yticklabels(unit_name_arr)
end
%%
if flag.plot_chan % Plot PSTH for each chan
%% Func of post onset
clrseq = brewermap(numel(postISI_msk),'Spectral');
smth = 9; nonoverlap_bnd = 150; % dur = dur_arr(1); 
for iCh = 1:numel(meta.spikeID)
figure(9);clf;hold on;set(9,'pos',[1000         168         586         810]);
T = tiledlayout(numel(postISI_msk),1,'Padding','compact','TileSpacing','compact');
preISI_far_msk = any(cat(1,preISI_msk{preISI_arr > nonoverlap_bnd}),1); % (preISI_msk{4} | preISI_msk{5} | preISI_msk{6})
for subfi = 1:numel(postISI_msk)
nexttile(subfi);
for duri = 1:numel(dur_msk)
dur = dur_arr(duri);
msk = preISI_far_msk & dur_msk{duri} & postISI_msk{subfi}; hold on
plot(wdw(1)+1:wdw(2), movmean(mean(rasters(iCh,:,msk),3),smth), 'Color',[clrseq(duri,:),1], 'LineWidth',2)
vline(0.0,'k-.');vline(dur,'k:'); vline(dur+postISI_arr(subfi),'r-.')
end
hold off
xlim([wdw(1),wdw(2)])
if subfi ~= numel(postISI_msk), xticklabels([]);end
end
align_tile_ylim(T);
xlabel("Time to Stimuli Onset")
title(T, compose("PSTH as a Function of Next Onset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(9,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.png",unit_name_arr(iCh))))
savefig(9,fullfile(figdir,compose("Post_onset_PSTH_allimg_Chan%s.fig",unit_name_arr(iCh))))
end
%% Func of Pre offset
smth = 9; dur = dur_arr; nonoverlap_bnd = 250;
for iCh = 1:numel(meta.spikeID)
figure(10);clf;hold on;set(10,'pos',[1000         168         586         810]);
T = tiledlayout(numel(preISI_msk),1,'Padding','compact','TileSpacing','compact');
postISI_far_msk = any(cat(1, postISI_msk{postISI_arr > nonoverlap_bnd}),1); % (postISI_msk{4} | postISI_msk{5} | postISI_msk{6})
for subfi = 1:numel(preISI_msk)
nexttile(subfi)
for duri = 1:numel(dur_msk)
dur = dur_arr(duri);
msk = postISI_far_msk & dur_msk{duri} & preISI_msk{subfi}; hold on
plot(wdw(1)+1:wdw(2), movmean(mean(rasters(iCh,:,msk),3),smth), 'Color',[clrseq(duri,:),1], 'LineWidth',2)
vline(0.0,'k-.');vline(dur,'k:');vline( - preISI_arr(subfi),'r-.')
end
hold off
xlim([wdw(1),wdw(2)])
if subfi ~= numel(preISI_msk), xticklabels([]);end
end
align_tile_ylim(T);
xlabel("Time to Stimuli Onset")
title(T, compose("PSTH as a Function of Next Onset\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(10,fullfile(figdir,compose("Pre_offset_PSTH_allimg_Chan%s.png",unit_name_arr(iCh))))
savefig(10,fullfile(figdir,compose("Pre_offset_PSTH_allimg_Chan%s.fig",unit_name_arr(iCh))))
% pause;
end
end
%% All Image plot separately pre offset
if flag.plot_imgsep
clrmap = brewermap(numel(imgmsks), 'Set3');
smth = 9;
preISI_far_msk = any(cat(1, preISI_msk{preISI_arr > nonoverlap_bnd}),1);
for iCh = 1:numel(meta.spikeID)
figure(12);clf;hold on;set(12,'pos',[1000         168         586         810]);
T = tiledlayout(numel(postISI_msk),1,'Padding','compact','TileSpacing','compact');
for figi = 1:numel(postISI_msk)
nexttile(figi);hold on
for duri = 1:numel(dur_msk)
dur = dur_arr(duri);
for imgi = 1:numel(imgmsks)
msk = preISI_far_msk & postISI_msk{figi} & dur_msk{duri} & imgmsks{imgi}'; 
plot(wdw(1)+1:wdw(2), movmean(mean(rasters(iCh,:,msk),3),smth), 'Color',[clrmap(imgi,:),0.5], 'LineWidth',2)
end
vline(0.0,'k-.');vline(dur,'k:');vline(dur + postISI_arr(figi),'r-.');
xlim([wdw(1),wdw(2)])
end
end
align_tile_ylim(T);
xlabel("Time to Stimuli Onset")
title(T, compose("PSTH as a Function of Post Offset ISI\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(12, fullfile(figdir,compose("Post_onset_PSTH_imgsep_Chan%s.png",unit_name_arr(iCh))))
savefig(12, fullfile(figdir,compose("Post_onset_PSTH_imgsep_Chan%s.fig",unit_name_arr(iCh))))
end
%%
smth = 9;
postISI_far_msk = any(cat(1, postISI_msk{postISI_arr > nonoverlap_bnd}),1); % (postISI_msk{5} | postISI_msk{6})
for iCh = 1:numel(meta.spikeID)
figure(11);clf;hold on;set(11,'pos',[1000         168         586         810]);
T = tiledlayout(numel(preISI_msk),1,'Padding','compact','TileSpacing','compact');
for figi = 1:numel(preISI_msk)
nexttile(figi);hold on
for duri = 1:numel(dur_msk)
dur = dur_arr(duri);
for imgi = 1:numel(imgmsks)
msk = preISI_msk{figi} & postISI_far_msk & dur_msk{duri} & imgmsks{imgi}'; 
plot(wdw(1)+1:wdw(2), movmean(mean(rasters(iCh,:,msk),3),smth), 'Color',[clrmap(imgi,:),0.5], 'LineWidth',2)
end
xlim([wdw(1),wdw(2)])
vline(0.0,'k-.');vline(dur,'k:');vline(- preISI_arr(figi),'r-.');
end
end
align_tile_ylim(T);
xlabel("Time to Stimuli Onset")
title(T, compose("PSTH as a Function of Pre Onset ISI\n%s Ch %s",Animal,unit_name_arr(iCh)))
saveas(11, fullfile(figdir,compose("Pre_offset_PSTH_imgsep_Chan%s.png",unit_name_arr(iCh))))
savefig(11, fullfile(figdir,compose("Pre_offset_PSTH_imgsep_Chan%s.fig",unit_name_arr(iCh))))
end
end
end




%%
[pre_ISI, post_ISI, duration, preISI_msk, postISI_msk, dur_msk, preISI_arr, postISI_arr, dur_arr] = parser_stimuli_timing(Trials);
function [pre_ISI, post_ISI, duration, preISI_msk, postISI_msk, dur_msk, preISI_arr, postISI_arr, dur_arr] = parser_stimuli_timing(Trials)
% Key function for the masking experiment analysis.
trialNum = numel(Trials.trialNum);
tOn_abs = nan(1,trialNum);
pre_ISI = nan(1,trialNum);
post_ISI = nan(1,trialNum);
duration = nan(1,trialNum);
for iTr = 1:trialNum
curT0 = double(Trials.imageONtime(iTr)); % image onset time 
iInTrial = Trials.imageInTrial(iTr); % image numbering in this RSVP trial

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
if isempty(tLastOff), tLastOff = -500; end

pre_ISI(iTr) = 0 - tLastOff;
post_ISI(iTr) = tNextOn - tCurOff;
duration(iTr) = tCurOff; % Trials.imageOFFtime(iTr) - Trials.imageONtime(iTr);

if tNextOn - tCurOff<0, keyboard;end
end

% Quantize the trials by the frame number (8.33333ms)
delay_limit = 250;
frameD = 100/12;
preISI_qu = round(pre_ISI/frameD);
postISI_qu = round(post_ISI/frameD);
% Sort the Quantized trials. 
preISI_categ = unique(round(pre_ISI(pre_ISI < delay_limit)/frameD)); 
postISI_categ = unique(round(post_ISI(post_ISI < delay_limit)/frameD)); 
% Create mask and array for it 
preISI_msk = arrayfun(@(frN)preISI_qu==frN, preISI_categ, 'Uni',0); 
preISI_msk{end+1} = pre_ISI >= delay_limit;
postISI_msk = arrayfun(@(frN)postISI_qu==frN, postISI_categ, 'Uni',0); 
postISI_msk{end+1} = post_ISI >= delay_limit;
preISI_arr = cellfun(@(idx) median(pre_ISI(idx)), preISI_msk);
postISI_arr = cellfun(@(idx) median(post_ISI(idx)), postISI_msk);

minor_categ = cellfun(@sum, preISI_msk) < trialNum / 20; 
preISI_msk(minor_categ) = [];
preISI_arr(minor_categ) = [];
minor_categ = cellfun(@sum, postISI_msk) < trialNum / 20; 
postISI_msk(minor_categ) = [];
postISI_arr(minor_categ) = [];

dur_qu = round(duration / frameD);
dur_categ = unique(dur_qu);
dur_msk = arrayfun(@(frN) dur_qu==frN, dur_categ, 'Uni', 0);
dur_arr = cellfun(@(idx) median(duration(idx)), dur_msk);
% get rid of some category with few trials (minority)
minor_categ = cellfun(@sum, dur_msk) < trialNum / 20; 
dur_msk(minor_categ) = [];
dur_arr(minor_categ) = [];
% hand coded categories
% preISI_msk = cellfun(@(wdw) (pre_ISI > wdw(1)) & (pre_ISI < wdw(2)), {[0,30],[30,65],[65,100],[100,150],[150,250],[250,inf]}, 'Uni', 0);
% postISI_msk = cellfun(@(wdw) (post_ISI > wdw(1)) & (post_ISI < wdw(2)), {[0,30],[30,65],[65,100],[100,150],[150,250],[250,inf]}, 'Uni', 0);
end
function [pre_ISI, post_ISI, duration, preISI_msk, postISI_msk, dur_msk] = postproc_timing(...
    Trials,pre_ISI, post_ISI, duration, preISI_msk, postISI_msk, dur_msk, preISI_arr, postISI_arr, dur_arr)

end
function fdrnm = figdir_naming(meta, ExpRecord)
Animal = meta.ephysFN(1:4);
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
chari = strfind(meta.expControlFN, 'Masking');
rowi_day = find(contains(ExpRecord.expControlFN, meta.expControlFN(1:chari+6)));
rowi = find(strcmp(ExpRecord.expControlFN, meta.expControlFN));
numday = find(rowi_day == rowi); % the numbering of experiment among exp in that day. 
fdrnm = compose("%s-%s-%02d",datestr(expday,'yyyy-mm-dd'), Animal, numday);
fprintf("Proposed folder name %s\n",fdrnm);
end

function align_tile_ylim(T)
ax_arr = get(T,'Children');
ylim_arr = [];
for i = 1:numel(ax_arr)
    ylim_arr(i,:) = ylim(ax_arr(i));
end
YLIM(1) = min(ylim_arr(:,1)); YLIM(2) = max(ylim_arr(:,2));
for i = 1:numel(ax_arr)
    ylim(ax_arr(i), YLIM);
    alllines = allchild(ax_arr(i));
    for ln = alllines' 
    if strcmp(ln.Tag,'vline')
        ln.YData = YLIM; % Change the y limit for the Vlines
    end
    end
end
end

function align_figures_ylim(figs)
ax_arr = [];
for fignum = figs
    ax = findobj(fignum,'Type','Axes');
    ax_arr = [ax_arr, ax];
end

ylim_arr = [];
for i = 1:numel(ax_arr)
    ylim_arr(i,:) = ylim(ax_arr(i));
end
YLIM(1) = min(ylim_arr(:,1)); YLIM(2) = max(ylim_arr(:,2));
for i = 1:numel(ax_arr)
    ylim(ax_arr(i), YLIM);
    alllines = allchild(ax_arr(i));
    for ln = alllines' 
    if strcmp(ln.Tag,'vline')
        ln.YData = YLIM; % Change the y limit for the Vlines
    end
    end
end
end

