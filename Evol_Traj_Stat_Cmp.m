%% Evolution Experiments compare statistics
%  For the paired experiments 1 deg evolution vs 3 deg evolution. 
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Evolution_Exp";
savepath = fullfile(result_dir, compose("Manifold_Evol_Resize_Cmp"));
mkdir(savepath);
global  Trials rasters 
final_scores_col = cell(1,16);
cnt = 1;
for groupi = 1:9
rel_i = (groupi - 1) * 2;
%     MAX_BLOCK_NUM = max(max(cell2mat(Trials_new{rel_i+1}.block)), max(cell2mat(Trials_new{rel_i+2}.block)));
%     color_seq = brewermap(MAX_BLOCK_NUM, 'spectral');
    % MAX_BLOCK_NUM = 50
for Triali = 1:2
meta = meta_new{Triali + rel_i};
rasters = rasters_new{Triali + rel_i};
Trials = Trials_new{Triali + rel_i};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi_tab = ExpSpecTable_Aug.Expi(exp_rowi);
%assert(Expi_tab == Expi, "Data Expi doesn't match that in exp record! Check your indexing or record.")
Expi = Expi_tab; 
fprintf("Processing Exp %d, %s\n", Expi, meta.comments)
%% Sort channel id
pref_chan = Trials.TrialRecord.User.prefChan;
pref_chan_id = find(meta.spikeID==pref_chan); % the id in the raster and lfps matrix 
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan);
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
%% Sort the block id and the identity of image. 
block_arr = cell2mat(Trials.block);
imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
gen_list = min(block_arr):max(block_arr);
%%
scores_tsr = squeeze(mean(rasters(:, 151:200, :), 2) - mean(rasters(:, 1:40, :), 2));
meanscore_syn = zeros(size(rasters, 1), length(gen_list)); % []
stdscore_syn = zeros(size(rasters, 1), length(gen_list));
meanscore_nat = zeros(size(rasters, 1), length(gen_list));
stdscore_nat = zeros(size(rasters, 1), length(gen_list));
for blocki = gen_list
    meanscore_syn(:,blocki) = mean(scores_tsr(:, row_gen & block_arr == blocki), 2);%cat(2, meanscore_syn, tmpscore_syn);
    meanscore_nat(:,blocki) = mean(scores_tsr(:, row_nat & block_arr == blocki), 2);%cat(2, meanscore_nat, tmpscore_nat);
    stdscore_syn(:,blocki)  = std(scores_tsr(:, row_gen & block_arr == blocki), 1, 2) / sqrt(sum(row_gen & block_arr == blocki));%cat(2, stdscore_syn, tmpstdscore_syn);
    stdscore_nat(:,blocki)  = std(scores_tsr(:, row_nat & block_arr == blocki), 1, 2) / sqrt(sum(row_gen & block_arr == blocki));%cat(2, stdscore_nat, tmpstdscore_nat);
end
% 
evol_stim_fr = zeros(size(rasters, 1), size(rasters, 2), length(gen_list));
evol_stim_sem = zeros(size(rasters, 1), size(rasters, 2), length(gen_list));
for blocki = 1:length(gen_list)
    evol_stim_fr(:, :, blocki) = mean(rasters(:,:, row_gen & block_arr == blocki),3);
    evol_stim_sem(:, :, blocki) = std(rasters(:,:, row_gen & block_arr == blocki),1,3)/sqrt(sum(row_gen & block_arr == blocki));
end
%%
channel_j = pref_chan_id;
%%
img_select = (block_arr >= max(block_arr)-5) & (block_arr <= max(block_arr)-1) & row_gen;
final_scores_col{cnt} = scores_tsr(channel_j, img_select);
fprintf("Mean score %.1f, ste %.1f\n", mean(final_scores_col{cnt}), std(final_scores_col{cnt})/sqrt(sum(img_select)))
cnt = cnt + 1;
end
end
%% Plot the statistics
h1 = figure('Visible','on');clf; h1.Position = [  19         235        1779         743];
hold on 
Exp_labels = [];
for groupi = 1:9
    rel_i = (groupi - 1) * 2;
    Expi_1 = ExpSpecTable_Aug.Expi(contains(ExpSpecTable_Aug.ephysFN, meta_new{rel_i + 1}.ephysFN));
    Expi_2 = ExpSpecTable_Aug.Expi(contains(ExpSpecTable_Aug.ephysFN, meta_new{rel_i + 2}.ephysFN));
    pair_label = compose("Exp%02d,%02d Ch%02d", Expi_1, Expi_2, Trials_new{rel_i + 1}.TrialRecord.User.prefChan);
    Exp_labels = [Exp_labels, pair_label];
    deg1scores = final_scores_col{rel_i + 1};
    deg3scores = final_scores_col{rel_i + 2};
    [~,P,CI] = ttest2(deg1scores, deg3scores);
    fprintf("P=%f,CI=[%.1f,%.1f]\n",P,CI(1),CI(2));
    text(rel_i + 1, 200, compose("P=%f\nCI=[%.1f,%.1f]",P,CI(1),CI(2)),'FontSize',12);
    mean_vec = [mean(deg1scores), mean(deg3scores)];
    ste_vec = [std(deg1scores)/sqrt(length(deg1scores)), ...
               std(deg3scores)/sqrt(length(deg3scores))];
    errorbar(rel_i + [1:2], mean_vec, ste_vec, 'LineWidth',2,'Color','k') % discard rel_i to get single comparison
end
xticks(1.5 + [0:2:16] )
xticklabels(Exp_labels)
xlabel("Experimental Pairs")
title("Comparing the Final Generations (end-5:end-1) Scores of 1 deg Evolution and 3 deg Evolutions")