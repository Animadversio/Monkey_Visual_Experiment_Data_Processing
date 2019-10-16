% Correlation analysis between channels

Expi = 8;
rasters = storedStruct.rasters{Expi};
meta = storedStruct.meta{Expi};
Trials = storedStruct.Trials{Expi};
n_chan = length(unit_label_arr);
%%
unit_label_arr = generate_unit_labels(meta.spikeID);
bsl_all = squeeze(mean(rasters(:, 1:50, :), 2));
score_all = squeeze(mean(rasters(:, 51:200, :), 2)) - bsl_all;
%%
img_ids = Trials.imageName;

%%
bsl_corr = corr(bsl_all');
rsp_corr = corr(score_all');
cmax_bsl = max(bsl_corr(bsl_corr<1),[],"all");
cmin_bsl = min(bsl_corr(bsl_corr<1),[],"all");
cmax_rsp = max(rsp_corr(rsp_corr<1),[],"all");
cmin_rsp = min(rsp_corr(rsp_corr<1),[],"all");
%%
figure(9)
subplot(131)
imagesc(bsl_corr,'AlphaData',double(eye(n_chan) == 0))
caxis([cmin_bsl, cmax_bsl]); colorbar;axis image
xticks(1:n_chan);xticklabels(unit_label_arr);
yticks(1:n_chan);yticklabels(unit_label_arr);
title(["Baseline Firing Rate Correlation"])
subplot(132)
imagesc(rsp_corr,'AlphaData',double(eye(n_chan) == 0))
caxis([cmin_rsp, cmax_rsp]);colorbar;axis image
xticks(1:n_chan);xticklabels(unit_label_arr);
yticks(1:n_chan);yticklabels(unit_label_arr);
title(["Rsp - Baseline Rate Correlation"])
subplot(133)
imagesc(rsp_corr,'AlphaData',double(eye(n_chan) == 0))
caxis([cmin_rsp, cmax_rsp]);colorbar;axis image
xticks(1:n_chan);xticklabels(unit_label_arr);
yticks(1:n_chan);yticklabels(unit_label_arr);
title(["Rsp - Baseline Rate Correlation"])
%%
figure(10)
subplot(121)
imagesc(bsl_corr,'AlphaData',double(eye(n_chan) == 0))
caxis([cmin_bsl, cmax_bsl]); colorbar;axis image
xticks(1:n_chan);xticklabels(unit_label_arr);
yticks(1:n_chan);yticklabels(unit_label_arr);

subplot(122)
imagesc(rsp_corr,'AlphaData',double(eye(n_chan) == 0))
caxis([cmin_rsp, cmax_rsp]);colorbar;axis image
xticks(1:n_chan);xticklabels(unit_label_arr);
yticks(1:n_chan);yticklabels(unit_label_arr);