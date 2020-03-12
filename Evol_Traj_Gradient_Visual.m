clearvars -except meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug ExpRecord
%% Visualizing Code Evolution Traj

Animal = "Beto"; Set_Path;
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Optimizer_Tuning";
% load Generator
G = FC6Generator("matlabGANfc6");
%% Loading the Exp data
Triali = 1;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);
fprintf("Processing  Exp %d:\n",Expi)
disp(ExpRecord.comments(exp_rowi))
% Fetch basic info of the experimebnt
pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
savepath = fullfile(result_dir, compose("%s_Evol%02d_chan%02d", Animal, Expi, pref_chan(1)));
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, strrep(string(Trials.TrialRecord.User.evoConfiguration{i,end}),"_"," ")];
end
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
for i = 1:thread_num
    pref_chan_id(i) = find(meta.spikeID==pref_chan(i) & ... % the id in the raster and lfps matrix 
                    unit_num_arr==unit_in_pref_chan(i)); % match for unit number
end

imgnm = Trials.imageName;
% seperate the thread natural images and generated images 
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm) & ...% first 2 characters are not digits
        cellfun(@(c) ~contains(c(end-4:end), "_nat"), imgnm) ; % last 4 characters are not `_nat`
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the thread 1 and 2 
row_thread0 = contains(imgnm, compose("thread%03d", 0));
row_thread1 = contains(imgnm, compose("thread%03d", 1));
assert(sum(row_thread0)+sum(row_thread1) == length(imgnm))
thread_msks = {row_thread0, row_thread1}; % store masks in a structure for the ease to iterate
% get the generation number 
block_arr = cell2mat(Trials.block);

% Sort image name and get scores in the given window
block_list = min(block_arr):max(block_arr);
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
meanscore_syn = nan(size(rasters, 1), length(block_list), 2); % [unit_num, gen_nums, threads]
stdscore_syn = nan(size(rasters, 1), length(block_list), 2); 
meanscore_nat = nan(size(rasters, 1), length(block_list), 2);
stdscore_nat = nan(size(rasters, 1), length(block_list), 2);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        meanscore_syn(:, blocki, threadi) = mean(scores_tsr(:, gen_msk), 2);
        meanscore_nat(:, blocki, threadi) = mean(scores_tsr(:, nat_msk), 2);
        stdscore_syn(:, blocki, threadi)  = std(scores_tsr(:, gen_msk), 1, 2) / sqrt(sum(gen_msk));
        stdscore_nat(:, blocki, threadi)  = std(scores_tsr(:, nat_msk), 1, 2) / sqrt(sum(nat_msk));
    end
end
%%
threadi=2;
scores_pref_ch = scores_tsr(pref_chan_id(threadi), :);
bsl_pref_ch = squeeze(mean(rasters(pref_chan_id(threadi), 1:40, :), 2));
rsp_pref_ch = squeeze(mean(rasters(pref_chan_id(threadi), 50:200, :), 2));
corrcoef(bsl_pref_ch',scores_pref_ch) 
corrcoef(bsl_pref_ch',rsp_pref_ch) 
% note the correlation between the baseline and score is significantly negative -0.8938! 
% Correlation bettwen rsponse and baseline rate is -0.1112, not
% significant. 
% The change of score is majorly (negatively) driven by fluctuation in baseline firing
% rate
%% 
% Find and load all codes 
[codes_all, img_ids, code_geni] = load_codes_all(meta.stimuli, threadi);
code_scores = nan(1, length(img_ids));
% img_ids_wsfx = cellfun(@(c) c(1:end-4), img_ids, 'UniformOutput', false);
%% For each generation in the experiment, sort the scores onto codes
for blocki = 2:max(block_arr)
    gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
    %nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
    gen_imgnms = imgnm(gen_msk);
    code_idx = zeros(1, length(gen_imgnms));
    for imgi = 1:numel(gen_imgnms)
        code_idx(imgi) = find(contains(img_ids, gen_imgnms{imgi}));
    end
    code_scores(code_idx) = scores_pref_ch(gen_msk);
end
%% get all codes in this gen, basis 
% do local PCA on codes? Or norm ? 
% Plot the image at the given location 
figure(1);clf;
for geni = 1:max(code_geni) - 1
    codes_cur = codes_all(code_geni <= geni+3 & code_geni >= geni-2, :);
    scores_cur = code_scores(code_geni <= geni+3 & code_geni >= geni-2); 
    [codeRD, RDmap] = pca(codes_cur, 3);
    sz_arr = 50*ones(1, length(scores_cur)); sz_arr(1:41:end) = 150;
    scatter3(codeRD(:,1),codeRD(:,2), codeRD(:,3), sz_arr, scores_cur,"filled")%,"MarkerFaceColor",0.5)
    xlabel("PC1");ylabel("PC2");zlabel("PC3")
    colorbar()
    pause
end

% Visualize the gradient / the next sample 
