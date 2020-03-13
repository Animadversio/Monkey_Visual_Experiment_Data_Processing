%% Sqeuence Effect and Baseline Subtraction

bsl_score_col = zeros(size(rasters,1), 1);
bsl_evk_col = zeros(size(rasters,1), 1);
threadi=1;
channel_j = pref_chan_id(threadi);
tmp = fopen(fullfile(exp_dir, "corr_stats.txt"),'w');
for channel_j = 1:size(rasters,1)
% scores_pref_ch = scores_tsr(channel_j, :);
bsl_pref_ch = squeeze(mean(rasters(channel_j, 1:40, :), 2));
rsp_pref_ch = squeeze(mean(rasters(channel_j, 50:200, :), 2));
scores_pref_ch = rsp_pref_ch - bsl_pref_ch;
bsl_RL = corr(bsl_pref_ch,max(0,scores_pref_ch) );
bsl_score = corr(bsl_pref_ch,scores_pref_ch); bsl_score_col(channel_j) = bsl_score;
bsl_evk = corr(bsl_pref_ch,rsp_pref_ch); bsl_evk_col(channel_j) = bsl_evk;
if ~ activ_msk(channel_j), continue, end
fprintf("Chan%s\tbsl~score:%.3f, bsl~ReLU(score):%.3f, bsl~evoked:%.3f\n",unit_name_arr(channel_j), bsl_score,bsl_RL,bsl_evk)
fprintf(tmp, "Chan%s\tbsl~score:%.3f, bsl~ReLU(score):%.3f, bsl~evoked:%.3f\n",unit_name_arr(channel_j), bsl_score,bsl_RL,bsl_evk);
end
fclose(tmp);
%%
exp_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Baseline_Seq_Effect\Optim_tuning_Exp4";
figure(8);clf;set(8,"position",[ 524         477        1084         501]);
subplot(121)
scatter(scores_pref_ch', bsl_pref_ch),xlabel("score"),ylabel("baseline")
title(sprintf("CorrCoef  %.3f",bsl_score(1,2)))
subplot(122)
scatter(rsp_pref_ch, bsl_pref_ch),xlabel("evoked rsp"),ylabel("baseline")
title(sprintf("CorrCoef  %.3f",bsl_evk(1,2)))
suptitle([Exp_label_str, sprintf("Chan%s",unit_name_arr(channel_j))])
set(8,'Units','pixels');
saveas(8,fullfile(exp_dir, sprintf("bsl_score_corr_ch%s.jpg",unit_name_arr(channel_j))))
%save_to_pdf(8, fullfile(exp_dir, sprintf("bsl_score_corr_ch%s.pdf",unit_name_arr(channel_j))))
% saveas(8,fullfile(exp_dir, sprintf("bsl_score_corr_ch%s.eps",unit_name_arr(channel_j))))
%% Plot the summry histogram
figure(1);clf
subplot(121)
histogram(bsl_score_col,30)
xlabel("correlation coefficient")
title("corr: baseline rate ~ (evk - baseline) rate")
subplot(122)
histogram(bsl_evk_col,30)
xlabel("correlation coefficient")
title("corr: baseline rate ~ evoked rate")
suptitle([Exp_label_str, "baseline correlation summary"])
saveas(1,fullfile(exp_dir, sprintf("corr_coef_summary.jpg")))

%% Plot the individual scatter for each channel. 
figure(8);clf;set(8,"position",[ 524         477        1084         501]);
for channel_j = 1:size(rasters,1)
    if ~ activ_msk(channel_j)
        continue
    end
    % scores_pref_ch = scores_tsr(channel_j, :);
    bsl_pref_ch = squeeze(mean(rasters(channel_j, 1:40, :), 2));
    rsp_pref_ch = squeeze(mean(rasters(channel_j, 50:200, :), 2));
    scores_pref_ch = rsp_pref_ch - bsl_pref_ch;
    bsl_score = corr(bsl_pref_ch,scores_pref_ch); bsl_score_col(channel_j) = bsl_score;
    bsl_evk = corr(bsl_pref_ch,rsp_pref_ch); bsl_evk_col(channel_j) = bsl_evk;
    % regression coeficient
    bsl_score_slp = regress(scores_pref_ch, [bsl_pref_ch, ones(length(bsl_pref_ch),1)]);
    bsl_evk_slp = regress(rsp_pref_ch, [bsl_pref_ch, ones(length(bsl_pref_ch),1)]);
    
    set(0,"CurrentFigure",8)
    subplot(121)
    scatter(scores_pref_ch', bsl_pref_ch),xlabel("score"),ylabel("baseline")
    title(sprintf("CorrCoef  %.3f, (slope %.3f, interc %.1f)",bsl_score, bsl_score_slp(1), bsl_score_slp(2)))
    subplot(122)
    scatter(rsp_pref_ch, bsl_pref_ch),xlabel("evoked rsp"),ylabel("baseline")
    title(sprintf("CorrCoef  %.3f, (slope %.3f, interc %.1f)",bsl_evk, bsl_evk_slp(1), bsl_evk_slp(2)))
    suptitle([Exp_label_str, sprintf("Chan%s",unit_name_arr(channel_j))])
    set(8,'Units','pixels');
    saveas(8,fullfile(exp_dir, sprintf("bsl_score_corr_ch%s.jpg",unit_name_arr(channel_j))))
end