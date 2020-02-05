Set_Path;
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Evolution_Exp";
ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx");
expftr = ExpSpecTable_Aug.Expi==20 & ...
     contains(ExpSpecTable_Aug.Exp_collection, "Manifold");
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr)); 
%%
h  = figure('Visible','on');h.Position = [828 42 1026 954];
for Triali = 2%1:length(meta_new)
% Fetch the trial info
%Triali = Expi - 26;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
Expi = ExpSpecTable_Aug.Expi(exp_rowi);
fprintf("Exp %d:\n",Expi)
disp(ExpSpecTable_Aug.comments(exp_rowi))
% Prepare the relevant folder and info
pref_chan = Trials.TrialRecord.User.prefChan;
pref_chan_id = find(meta.spikeID==pref_chan);
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan);
savepath = fullfile(result_dir, compose("Manifold_Evol%02d_chan%02d", Expi, pref_chan));
mkdir(savepath);
unit_name_arr = generate_unit_labels(meta.spikeID, savepath); % Generate readable labels for each channel
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, spikeID, rasters, savepath);
% Compute the block structure from imagenames
imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
block_arr = cell2mat(Trials.block);
%% Generate Gradual Changing Color Labels 
block_list = min(block_arr):max(block_arr);
% fun = @(m)srgb_to_Lab(m);
% color_seq = maxdistcolor(30,fun); % Use this to generate maximally distinguishable color sequence
color_seq = brewermap(length(block_list), 'spectral'); % Use this to generate gradual changing sequence
%% 
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
meanscore_syn = nan(size(rasters, 1), length(block_list)); % [unit_num, gen_nums, threads]
stdscore_syn = nan(size(rasters, 1), length(block_list)); 
meanscore_nat = nan(size(rasters, 1), length(block_list));
stdscore_nat = nan(size(rasters, 1), length(block_list));
for blocki = min(block_arr):max(block_arr)
    gen_msk = row_gen & block_arr == blocki;% & thread_msks{threadi}; 
    nat_msk = row_nat & block_arr == blocki;% & thread_msks{threadi};
    meanscore_syn(:, blocki) = mean(scores_tsr(:, gen_msk), 2);
    meanscore_nat(:, blocki) = mean(scores_tsr(:, nat_msk), 2);
    stdscore_syn(:, blocki)  = std(scores_tsr(:, gen_msk), 1, 2) / sqrt(sum(gen_msk));
    stdscore_nat(:, blocki)  = std(scores_tsr(:, nat_msk), 1, 2) / sqrt(sum(nat_msk));
end
%%
block_list = min(block_arr):max(block_arr);
evol_stim_fr = zeros(size(rasters, 1), size(rasters, 2), length(block_list));
evol_stim_sem = zeros(size(rasters, 1), size(rasters, 2), length(block_list));
for blocki = block_list
    gen_msk = row_gen & block_arr == blocki;
    evol_stim_fr(:, :, blocki) = mean(rasters(:,:, gen_msk),3);
    evol_stim_sem(:, :, blocki) = std(rasters(:,:, gen_msk),1,3)/sqrt(sum(gen_msk));
end
% nat_stim_fr = zeros(max(natural_stim_i), size(rasters, 2), size(rasters, 3));
% nat_stim_fr_std = zeros(max(natural_stim_i), size(rasters, 2), size(rasters, 3));
% nat_stim_fr_sem = zeros(max(natural_stim_i), size(rasters, 2), size(rasters, 3));
% for i = 1:max(natural_stim_i)
%     nat_stim_fr(i,:,:) =    mean(rasters(sort_idx(natural_stim_i==i), :, :),1);
%     nat_stim_fr_std(i,:,:) = std(rasters(sort_idx(natural_stim_i==i), :, :),1,1);
%     nat_stim_fr_sem(i,:,:) = nat_stim_fr_std(i,:,:) / sqrt(sum(natural_stim_i==i));
% end
%% Get the Image FileName Sequence
channel_j = pref_chan_id(1);
imgColl = repmat("", length(block_list),1);
scoreColl = zeros(length(block_list),1);
for blocki = block_list
    gen_msk = row_gen & block_arr == blocki; 
    [maxScore, maxIdx] = max(scores_tsr(channel_j, gen_msk));
    tmpimgs = imgnm(gen_msk);
    imgfullfn = ls(fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]))
    imgColl(blocki) = fullfile(meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % solved @ Jan.29
    scoreColl(blocki) = maxScore;
end
% Montage the images
set(0,'CurrentFigure',h); %clf; %
montage(imgColl(:))
title([Exp_label_str, compose('Best Image per Generation')])
saveas(h, fullfile(savepath, "EvolImageSeq.png"))
end