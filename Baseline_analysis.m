
%% Get a elongated window version of formatted file
[meta_,rasters_,lfps_,Trials_] = loadData(meta.ephysFN,'expControlFN',meta.expControlFN,'rasterWindow',[-100 1000]) ;
%%
channel_j = 45;
figure, imagesc(squeeze(rasters(channel_j,:,:)))
yticks(0:100:1100);yticklabels(compose("%d",-100:100:1000))
%%
figure
plot(-99:1000, mean(rasters(45,:,:),3))
xticks(-100:100:1000);xticklabels(compose("%d",-100:100:1000))
fig_vertline(0,axis);
fig_vertline(50,axis);
fig_vertline(200,axis);
%%
% std(rasters_(45,:,:),0,3);
figure(10);set(10,'position',[827         274        1193         566])
channel_j = 45;
for channel_j = 1:size(rasters,1)
if ~ activ_msk(channel_j), continue, end
%plot(-99:1000, mean(rasters_(45,:,:),3))
set(0,"CurrentFigure",10);clf;
shadedErrorBar(-99:1000, mean(rasters(channel_j,:,:),3), std(rasters(channel_j,:,:),0,3))
%xticks(-100:100:1000);xticklabels(compose("%d",-100:100:1000))
fig_vertline(0,axis);
fig_vertline(50,axis);
fig_vertline(200,axis);
xlim([-100,1000]);xlabel("Time relative to Stimuli onset (ms)")
title([Exp_label_str,sprintf("Channel %s", unit_name_arr(channel_j))])
saveas(10,fullfile(exp_dir, sprintf("mean_rsp_curve_ch%s.jpg", unit_name_arr(channel_j))))
%shadedErrorBar
end
%%
figure, imagesc(1:size(rasters,3),-99:1000,squeeze(lfps(30,:,:)))
%yticks(0:100:1100);yticklabels(compose("%d",-100:100:1000))
%%
prebsl_bsl_col = nan(size(rasters,1),1);
prebsl_rsp_col = nan(size(rasters,1),1);
bsl_evk_col = nan(size(rasters,1),1);
bsl_score_col = nan(size(rasters,1),1);
bsl_preevk_col = nan(size(rasters,1),1);
evk_preevk_col = nan(size(rasters,1),1);
prebsl_preevk_col = nan(size(rasters,1),1);
exp_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Baseline_Seq_Effect\Optim_tuning_Exp4";
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
tmp = fopen(fullfile(exp_dir, "corr_stats_add.txt"),'w');
%channel_j = pref_chan_id(threadi);
for channel_j = 1:size(rasters,1)
prebsl_pref_ch = squeeze(mean(rasters(channel_j, 1:100, :), 2)); % pre image onset 
bsl_pref_ch = squeeze(mean(rasters(channel_j, 101:140, :), 2)); % baseline after image onset 
rsp_pref_ch = squeeze(mean(rasters(channel_j, 150:300, :), 2)); % response part 
scores_pref_ch = rsp_pref_ch - bsl_pref_ch;
prebsl_bsl = corr(prebsl_pref_ch, bsl_pref_ch); prebsl_bsl_col(channel_j) = prebsl_bsl;
prebsl_rsp = corr(prebsl_pref_ch, rsp_pref_ch); prebsl_rsp_col(channel_j) = prebsl_rsp;
bsl_evk   =  corr(bsl_pref_ch, rsp_pref_ch);  bsl_evk_col(channel_j) = bsl_evk;
bsl_score =  corr(bsl_pref_ch, scores_pref_ch); bsl_score_col(channel_j) = bsl_score;

bsl_preevk    = corr(bsl_pref_ch(2:end),rsp_pref_ch(1:end-1));  bsl_preevk_col(channel_j) = bsl_preevk;
evk_preevk    = corr(rsp_pref_ch(2:end),rsp_pref_ch(1:end-1));  evk_preevk_col(channel_j) = evk_preevk;
prebsl_preevk = corr(prebsl_pref_ch(2:end),rsp_pref_ch(1:end-1));  prebsl_preevk_col(channel_j) = prebsl_preevk;
if ~ activ_msk(channel_j), continue, end
fprintf("Chan%s\t bsl~prebsl:%.3f, evk~prebsl:%.3f, prebsl~evk(last):%.3f, bsl~evk(last):%.3f, evk~evk(last):%.3f, bsl~evk:%.3f, bsl~score:%.3f\n",...
    unit_name_arr(channel_j), prebsl_bsl, prebsl_rsp, prebsl_preevk, bsl_preevk, evk_preevk, bsl_evk, bsl_score)
fprintf(tmp, "Chan%s\t bsl~prebsl:%.3f, evk~prebsl:%.3f, prebsl~evk(last):%.3f, bsl~evk(last):%.3f, evk~evk(last):%.3f, bsl~evk:%.3f, bsl~score:%.3f\n",...
    unit_name_arr(channel_j), prebsl_bsl, prebsl_rsp, prebsl_preevk, bsl_preevk, evk_preevk, bsl_evk, bsl_score);
end
fclose(tmp);
%% Summary Figure of all statistics across channel (outliers can be spotted in the text files)
Exp_label_str = "Optim Tuning Exp4 pref chan 30";
figure(1);clf;set(1,'position',[ 214         122        1587         856]);
subplot(231)
histogram(prebsl_bsl_col,30)
xlabel("correlation coefficient");title(["pre Bsl (-100-0ms) ~"," Bsl (0-50ms) rate"]);
subplot(232)
histogram(prebsl_rsp_col,30)
xlabel("correlation coefficient");title(["pre Bsl (-100-0ms) ~"," reponse (50-200ms) rate"]);
subplot(233)
histogram(bsl_evk_col,30)
xlabel("correlation coefficient");title(["Bsl (0-50ms) ~"," reponse (50-200ms) rate"]);
subplot(234)
histogram(bsl_preevk_col,30)
xlabel("correlation coefficient");title(["bsl (0-50ms) ~"," reponse (50-200ms) rate (last trial)"]);
subplot(235)
histogram(evk_preevk_col,30)
xlabel("correlation coefficient");title(["response (50-200ms) ~"," reponse (50-200ms) rate (last trial)"]);
subplot(236)
histogram(prebsl_preevk_col,30)
xlabel("correlation coefficient");title(["pre Bsl (-100-0ms) ~"," reponse (50-200ms) rate (last trial)"]);
% subplot(231)
% histogram(bsl_score,bin=30)
% xlabel("correlation coefficient");title("bsl_score");
suptitle([Exp_label_str, "Additional Corr Coef and Effect of Last Trial Response"])
saveas(1,fullfile(exp_dir, sprintf("corr_coef_seq_summary.jpg")))
%%
figure(11);clf;set(11,'position',[725         143        1033         835])
for channel_j = 1:size(rasters,1)
if ~ activ_msk(channel_j), continue, end
prebsl_pref_ch = squeeze(mean(rasters(channel_j, 1:100, :), 2)); % pre image onset 
bsl_pref_ch = squeeze(mean(rasters(channel_j, 101:140, :), 2)); % baseline after image onset 
rsp_pref_ch = squeeze(mean(rasters(channel_j, 150:300, :), 2)); % response part 
scores_pref_ch = rsp_pref_ch - bsl_pref_ch;
prebsl_bsl = corr(prebsl_pref_ch, bsl_pref_ch); prebsl_bsl_slp = regress(bsl_pref_ch, [prebsl_pref_ch, ones(length(prebsl_pref_ch),1)]);
prebsl_rsp = corr(prebsl_pref_ch, rsp_pref_ch); prebsl_rsp_slp = regress(rsp_pref_ch, [prebsl_pref_ch, ones(length(prebsl_pref_ch),1)]);
bsl_evk   =  corr(bsl_pref_ch, rsp_pref_ch);  bsl_evk_slp = regress(rsp_pref_ch, [bsl_pref_ch, ones(length(bsl_pref_ch),1)]);
bsl_score =  corr(bsl_pref_ch, scores_pref_ch); bsl_score_slp = regress(scores_pref_ch, [bsl_pref_ch, ones(length(bsl_pref_ch),1)]);
bsl_preevk    = corr(bsl_pref_ch(2:end),rsp_pref_ch(1:end-1));  bsl_preevk_slp = regress(bsl_pref_ch(2:end), [rsp_pref_ch(1:end-1), ones(length(rsp_pref_ch)-1,1)]);
evk_preevk    = corr(rsp_pref_ch(2:end),rsp_pref_ch(1:end-1));  evk_preevk_slp = regress(rsp_pref_ch(2:end), [rsp_pref_ch(1:end-1), ones(length(rsp_pref_ch)-1,1)]);
prebsl_preevk = corr(prebsl_pref_ch(2:end),rsp_pref_ch(1:end-1));  prebsl_preevk_slp = regress(prebsl_pref_ch(2:end), [rsp_pref_ch(1:end-1), ones(length(rsp_pref_ch)-1,1)]);

set(0,"CurrentFigure",11);clf;
subplot(221)
scatter(prebsl_pref_ch, bsl_pref_ch),xlabel("pre bsl (-100-0)"),ylabel("baseline (0-40)")
title(sprintf("CorrCoef  %.3f, (slope %.3f, interc %.1f)", prebsl_bsl, prebsl_bsl_slp(1), prebsl_bsl_slp(2)))
subplot(223)
scatter(rsp_pref_ch(1:end-1), prebsl_pref_ch(2:end)),ylabel("pre bsl (-100-0)"),xlabel("evoked rsp (last trial)")
title(sprintf("CorrCoef  %.3f, (slope %.3f, interc %.1f)", prebsl_preevk, prebsl_preevk_slp(1), prebsl_preevk_slp(2)))
subplot(222)
scatter(rsp_pref_ch(1:end-1), bsl_pref_ch(2:end)),ylabel("baseline (0-40)"),xlabel("evoked rsp (last trial)")
title(sprintf("CorrCoef  %.3f, (slope %.3f, interc %.1f)",bsl_preevk, bsl_preevk_slp(1), bsl_preevk_slp(2)))
subplot(224)
scatter(rsp_pref_ch(1:end-1), rsp_pref_ch(2:end)),ylabel("evoked rsp"),xlabel("evoked rsp (last trial)")
title(sprintf("CorrCoef  %.3f, (slope %.3f, interc %.1f)",evk_preevk, evk_preevk_slp(1), evk_preevk_slp(2)))
suptitle([Exp_label_str, sprintf("Chan%s",unit_name_arr(channel_j))])
saveas(11,fullfile(exp_dir, sprintf("bsl_seq_corr_add_ch%s.jpg",unit_name_arr(channel_j))))
%shadedErrorBar
end