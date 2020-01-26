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
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
%%
imgColl = repmat("", length(block_list), thread_num);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        % nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        % meanscore_syn(:, blocki, threadi) = mean(scores_tsr(:, gen_msk), 2);
        [maxScore, maxIdx] = max(scores_tsr(channel_j, gen_msk));
        tmpimgs = imgnm(gen_msk);
        imgColl(blocki, threadi) = fullfile(meta.stimuli, [tmpimgs(maxIdx)+".bmp"]);
    end
end
%%
figure(5)
subplot(1,2,1)
montage(imgColl(:,1))
title([Exp_label_str, compose('Best Image per Generation'), compose("Optimizer %s", Optim_names(1))])
subplot(1,2,2)
montage(imgColl(:,2))
title([Exp_label_str, compose('Best Image per Generation'), compose("Optimizer %s", Optim_names(2))])
