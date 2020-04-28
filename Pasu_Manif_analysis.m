clearvars -except PasuStats
%%
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'))
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
%%

N = samp_num;comp_n = 12;
FFTidx = [1:comp_n + 1, N-comp_n+1:N];
coef_col = zeros(length(FFT_col), length(FFTidx));
for imgi = 1:length(FFT_col)
    coef_col(imgi, :) = FFT_col{imgi}(FFTidx);
end
tsne_coord=tsne([real(coef_col),imag(coef_col)],"NumDimensions",2);
[~, umap]=run_umap([real(coef_col),imag(coef_col)]);
save(fullfile(result_dir,"PasuPatchStat.mat"),'umap','tsne_coord','coef_col','FFT_col')
%%
pref_chan_arr = arrayfun(@(c) c.units.pref_chan, Stats);
V1_msk = pref_chan_arr >=33 & pref_chan_arr <49;
V4_msk = pref_chan_arr >=49;
V4Expi = find(V4_msk);
Stats(19).ref.pasu_psths
%%
Expi = 20;
act_map = cellfun(@(psth) nanmean(psth(1,51:200,:),'all'), Stats(Expi).ref.pasu_psths);
act_vec = reshape(act_map',1,[]);
act_vec = act_vec(~isnan(act_vec));
figure;
scatter(umap.embedding(:,1),umap.embedding(:,2),act_vec/3,act_vec,'filled')
title(compose("Exp %d Channel %d unit 1",Expi, Stats(Expi).units.pref_chan))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir, Animal+"_Manif_PasuStats.mat"))
%%
save(fullfile(mat_dir, Animal+"_Manif_PasuStats.mat"), PasuStats)
%%
% img_dir = 'N:\Stimuli\2019-Manifold\pasupathy-wg-f-4-ori';
img_dir = 'E:/Github_Projects/Fourier_Curve_Motions/imgs/pasupathy-wg-f-4-ori/';
imlist = string(ls(fullfile(pasu_path,"*.jpg")));
pasuimg_list = cell(51, 4);%{};
for j = 1:4
    for i = 1:51
        cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
        img_idx = find(contains(imlist, cur_fn));
        if length(img_idx) == 0 % empty space
            pasuimg_list{i,j} = [];
        else
            pasuimg_list{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
        end
    end
end
%% Main loop
result_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Pasupathy_Space";
figure(10);set(10,'position',[826   241   732   732])
figure(11);set(11,'position',[826   241   732   732])
for Expi = 11:45
savedir = fullfile(result_dir,compose("%s_Exp%d_chan%02d", Animal, Expi, PasuStats(Expi).units.pref_chan));
mkdir(savedir)
PasuStats(Expi).ref.pasu_stats = repmat(summary,0,0);
PasuStats(Expi).ref.pasu_strs =  repmat("",0,0);
for chan_j = 1:length(PasuStats(Expi).units.spikeID)
    act_cell = cellfun(@(c)squeeze(c(chan_j,:,:)),PasuStats(Expi).ref.pasu_act,'UniformOutput',false);
    bsl_cell = cellfun(@(c)squeeze(c(chan_j,:,:)),PasuStats(Expi).ref.pasu_bsl,'UniformOutput',false);
    [summary, stat_str] = calc_pasu_tuning_stats(act_cell, bsl_cell);
    PasuStats(Expi).ref.pasu_stats(chan_j) = summary;
    PasuStats(Expi).ref.pasu_strs(chan_j) = stat_str;
    act_mat = cellfun(@(c)nanmean(c(chan_j,:,:),"all"),PasuStats(Expi).ref.pasu_act);
%     frame_img_list = score_frame_image_arr(pasuimg_list, hl_mat, ...
%         caxis(ax1), colormap(ax1), 50);
    act_vec = reshape(act_mat',1,[]);
    act_vec = act_vec(~isnan(act_vec));
    set(0,'CurrentFigure',10);clf
    scatter(umap.embedding(:,1),umap.embedding(:,2),40,act_vec,'filled');
    title(compose("%s Exp %d Channel %s\nUMAP embedding\n%s", Animal, Expi, PasuStats(Expi).units.unit_name_arr(chan_j),stat_str))
    CLIM = caxis(gca);
    set(0,'CurrentFigure',11);clf
    scatter(tsne_coord(:,1),tsne_coord(:,2),40,act_vec,'filled');
    title(compose("%s Exp %d Channel %s\ntSNE embedding\n%s", Animal, Expi, PasuStats(Expi).units.unit_name_arr(chan_j),stat_str))
    saveas(10, fullfile(savedir, sprintf("Pasu_UMAP_chan%s.png", PasuStats(Expi).units.unit_name_arr{chan_j})))
    saveas(11, fullfile(savedir, sprintf("Pasu_tSNE_chan%s.png", PasuStats(Expi).units.unit_name_arr{chan_j})))
    
end

end
%%
idx_mat = reshape(1:4*51,51,4);
idx_mat = arrayfun(@(idx){idx},idx_mat);

act_cell = cellfun(@(c)squeeze(c(chan_j,:,:)),PasuStats(Expi).ref.pasu_act,'UniformOutput',false);
bsl_cell = cellfun(@(c)squeeze(c(chan_j,:,:)),PasuStats(Expi).ref.pasu_bsl,'UniformOutput',false);
idx_cell = cellfun(@(idx,act)repmat(idx,length(act),1),idx_mat,act_cell,'UniformOutput',false); % create the idx cell array of same size
act_vec = cat(1, act_cell{:});
bsl_vec = cat(1, bsl_cell{:});
idx_vec = cat(1, idx_cell{:});
[H,P,CI,STATS] = ttest(act_vec, bsl_vec);
[P,ANOVATAB,STATS] = anova1(act_vec,idx_vec,'off');
%%
frame_img_list = score_frame_image_arr(pasuimg_list, act_mat, CLIM, colormap(gca), 50);
%%
[summary, stat_str] = calc_pasu_tuning_stats(act_cell, bsl_cell)
%%
idx_mat = reshape(1:4*51,51,4);
idx_mat = arrayfun(@(idx){idx},idx_mat);
function [summary, stat_str] = calc_pasu_tuning_stats(act_cell, bsl_cell) 
idx_mat = reshape(1:4*51,51,4);
idx_mat = arrayfun(@(idx){idx},idx_mat);   
    idx_cell = cellfun(@(idx,act)repmat(idx,length(act),1),idx_mat,act_cell,'UniformOutput',false); % create the idx cell array of same size
    act_vec = cat(1, act_cell{:});
    bsl_vec = cat(1, bsl_cell{:});
    idx_vec = cat(1, idx_cell{:});
    % Do statistics
    [p,tbl,stats] = anova1(act_vec,idx_vec,'off');
    stats.F = tbl{2,5};
    stats.p = p;
    summary.anova_F = stats.F;
    summary.anova_p = stats.p; 
    [p,tbl,stats_bsl] = anova1(bsl_vec,idx_vec,'off');
    stats_bsl.F = tbl{2,5};
    stats_bsl.p = p;
    summary.anova_F_bsl = tbl{2,5};
    summary.anova_p_bsl = p; 
    %
    [~,P,CI,STATS] = ttest(act_vec, bsl_vec);
    summary.t = STATS.tstat;
    summary.t_p = P;
    summary.t_CI = CI;
    % visualize
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.1f,%.1f] t=%.2f(%.2e)\n' ...
            'Evoked Modulation: All image, F=%.2f(%.2e)\n' ...
            'Baseline Modulation: All image, F=%.2f(%.2e)'],CI(1), CI(2), STATS.tstat, P, stats.F, stats.p,stats_bsl.F,stats_bsl.p);
end