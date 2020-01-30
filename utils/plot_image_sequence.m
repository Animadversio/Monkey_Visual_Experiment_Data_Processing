
h = figure();set(h,'position',[1          41        2560         963]);
axs{1} = subplot(1,2,1);axs{2} = subplot(1,2,2);
%%
for Triali = 1
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpSpecTable_Aug.Expi(exp_rowi);
% 
pref_chan = Trials.TrialRecord.User.prefChan;
assert(pref_chan(1) == pref_chan(2))
pref_chan_id = find(meta.spikeID==pref_chan(1)); % the id in the raster and lfps matrix
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
% 
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, string(Trials.TrialRecord.User.evoConfiguration{i,end})];
end
% 
savepath = fullfile(result_dir, compose("Evol%02d_chan%02d", Expi, pref_chan(1)));
mkdir(savepath);
% Plot Image sequence
imgnm = Trials.imageName;
% seperate the thread natural images and generated images 
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the thread 1 and 2 
row_thread0 = contains(imgnm, compose("thread%03d", 0));
row_thread1 = contains(imgnm, compose("thread%03d", 1));
assert(sum(row_thread0)+sum(row_thread1) == length(imgnm))
thread_msks = {row_thread0, row_thread1}; % store masks in a structure for the ease to iterate
% get the generation number 
channel_j = pref_chan_id;
block_arr = cell2mat(Trials.block);
block_list = min(block_arr):max(block_arr);
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
%% Get the Image FileName Sequence
imgColl = repmat("", length(block_list), thread_num);
scoreColl = zeros(length(block_list), thread_num);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        [maxScore, maxIdx] = max(scores_tsr(channel_j, gen_msk));
        tmpimgs = imgnm(gen_msk);
        imgfullfn = ls(fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]));
        assert(~isempty(imgfullfn), "Image not found %s",fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]))
        imgColl(blocki, threadi) = fullfile(meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % solved @ Jan.29
        scoreColl(blocki, threadi) = maxScore;
    end
end
%% Montage the images
set(0,'CurrentFigure',h); %clf; %
set(gcf, "CurrentAxes", axs{1}); cla(axs{1},'reset'); 
montage(imgColl(:,1))
title([Exp_label_str, compose('Best Image per Generation'), compose("Optimizer %s", Optim_names(1))])
set(gcf, "CurrentAxes", axs{2}); cla(axs{2},'reset'); 
montage(imgColl(:,2))
title([Exp_label_str, compose('Best Image per Generation'), compose("Optimizer %s", Optim_names(2))])
saveas(h, fullfile(savepath, "EvolImageSeq_cmp.png"))

%% t test the last 5 generations 
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
scores_thread1 = scores_tsr(pref_chan_id, row_gen & (block_arr > max(block_arr)-6) & thread_msks{1});
scores_thread2 = scores_tsr(pref_chan_id, row_gen & (block_arr > max(block_arr)-6) & thread_msks{2});
[~,Pval,CI] = ttest2(scores_thread1, scores_thread2);
end