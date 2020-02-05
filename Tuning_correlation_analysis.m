% Correlation analysis between channels
Set_Exp_Specs;
%Expi = 8;
for Expi = 10
rasters = storedStruct.rasters{Expi};
meta = storedStruct.meta{Expi};
Trials = storedStruct.Trials{Expi};
ang_step = 18;
sphere_norm = norm_arr(Expi);
pref_chan = pref_chan_arr(Expi);
savepath = sprintf("C:\\Users\\ponce\\OneDrive\\Desktop\\OneDrive_Binxu\\OneDrive\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
global unit_label_arr
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
n_chan = length(unit_label_arr);
exp_str = sprintf("Exp%d preferred channel %d", Expi, pref_chan);
%% Calculate the baseline and response rate 
bsl_all = squeeze(mean(rasters(:, 1:50, :), 2));
score_all = squeeze(mean(rasters(:, 51:200, :), 2)) - bsl_all;
%% Calculate the mean subtracted residue rate 
res_all = zeros(size(score_all));
img_ids = Trials.imageName;
pattern_arr = {'norm_%d_PC2_%d_PC3_%d', 'norm_%d_PC49_%d_PC50_%d', 'norm_%d_RND1_%d_RND2_%d'};
for k=1:3
    name_pattern = pattern_arr{k};
    for i =-5:5
        cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, 90);
        img_idx = find(contains(Trials.imageName, cur_fn));
        for t =1:length(img_idx)
            img_ids{img_idx(t)} = sprintf(name_pattern, sphere_norm, 0, 90);
        end
        cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, -90);
        img_idx = find(contains(Trials.imageName, cur_fn));
        for t =1:length(img_idx)
            img_ids{img_idx(t)} = sprintf(name_pattern, sphere_norm, 0, -90);
        end
    end
end
uniq_img_idx = unique(img_ids);
for j = 1:length(uniq_img_idx)
    cur_fn = uniq_img_idx{j};
    img_idx = find(contains(img_ids, cur_fn));
    res_all(:, img_idx) = score_all(:, img_idx) - mean(score_all(:, img_idx), 2);
end
%%
figure(9);set(9,'position',[0          40        2560         963])
ax1 = subplot(131);set(ax1,'position',[0.05, 0.05, 0.28, 0.90])
bsl_corr = corr_mat_plot(ax1, bsl_all);
title([exp_str, "Baseline Firing Rate Correlation"])
ax2 = subplot(132);set(ax2,'position',[0.35, 0.05, 0.28, 0.90])
rsp_corr = corr_mat_plot(ax2, score_all);
title([exp_str, "Rsp - Baseline Rate Correlation"])
ax3 = subplot(133);set(ax3,'position',[0.65, 0.05, 0.28, 0.90])
res_corr = corr_mat_plot(ax3, res_all);
title([exp_str, "Stimuli Avg Subtracted (Residue) Rate Correlation"])
%%
diff_corr = (rsp_corr - res_corr) ;
figure
imagesc(diff_corr)%'AlphaData',double(eye(n_chan) == 0))
xticks(1:n_chan);xticklabels(unit_label_arr);
yticks(1:n_chan);yticklabels(unit_label_arr);
%%
figure(10);set(10,'position',[0          40        2560         963])
ax1 = subplot(131);set(ax1,'position',[0.05, 0.05, 0.28, 0.90])
img_in_space_idx = find(and(contains(Trials.imageName, "PC2"), contains(Trials.imageName, "PC3")));
rsp_corr1 = corr_mat_plot(ax1, score_all(:,img_in_space_idx));
title([exp_str, "Score Correlation in PC23 space"])
ax2 = subplot(132);set(ax2,'position',[0.35, 0.05, 0.28, 0.90])
img_in_space_idx = find(and(contains(Trials.imageName, "PC49"), contains(Trials.imageName, "PC50")));
rsp_corr2 = corr_mat_plot(ax2, score_all(:,img_in_space_idx));
title([exp_str, "Score Correlation in PC49 50 space"])
ax3 = subplot(133);set(ax3,'position',[0.65, 0.05, 0.28, 0.90])
img_in_space_idx = find(and(contains(Trials.imageName, "RND1"), contains(Trials.imageName, "RND2")));
rsp_corr3 = corr_mat_plot(ax3, score_all(:,img_in_space_idx));
title([exp_str, "Score Correlation in RND12 space"])
%%
savefig( 9, fullfile(savepath, 'corr_mat_bsl_rsp_res_cmp.fig'))
savefig(10, fullfile(savepath, 'corr_mat_subspace_cmp.fig'))
saveas( 9, fullfile(savepath, 'corr_mat_bsl_rsp_res_cmp.bmp'))
saveas(10, fullfile(savepath, 'corr_mat_subspace_cmp.bmp'))
end
% Baseline firing events

%% Demo code
figure(11)
scatter(score_all(37,:) + bsl_all(37,:),score_all(47,:) + bsl_all(47,:))
% 
[~,P,CI] = ttest(bsl_all(47,:) , score_all(47,:))
%%
rndidx = 1:size(rasters,3);%randsample(size(rasters,3), 100);
traces = squeeze(mean(rasters(:, 51:200, rndidx), 2))-squeeze(mean(rasters(:, 1:50, rndidx), 2));
figure(11);set(11,'position',[ 1380          42         829         954]);ax = subplot(1,1,1);
corr_mat_plot(ax, traces);
%
figure(12);
imagesc(traces)
yticks(1:n_chan);yticklabels(unit_label_arr);
%%
function rsp_corr = corr_mat_plot(ax, traces)
global unit_label_arr
n_chan = size(traces, 1);
rsp_corr = corr(traces');
cmax_rsp = max(rsp_corr(rsp_corr<1),[],"all");
cmin_rsp = min(rsp_corr(rsp_corr<1),[],"all");
imagesc(ax, rsp_corr,'AlphaData',double(eye(n_chan) == 0))
xticks(1:n_chan);xticklabels(unit_label_arr);
yticks(1:n_chan);yticklabels(unit_label_arr);
caxis([cmin_rsp, cmax_rsp]);colorbar;axis image
end