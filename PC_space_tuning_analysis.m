%load(fullfile("\\storage1.ris.wustl.edu\crponce\Active\Data-Ephys-MAT","Beto64chan-02102019-003_formatted.mat"))
load(fullfile("\\storage1.ris.wustl.edu\crponce\Active\Data-Ephys-MAT","Beto64chan-03102019-002_formatted.mat"))
% meta,rasters,lfps,Trials
%% Loading code from "D:\Poncelab_Github\office-main\Project_Selectivity_Beto_loadRaw.m"

img_names = unique(Trials.imageName);

norm = 326;
ang_step = 18;
pref_chan = 6;
savepath = "C:\Users\ponce\OneDrive\Desktop\OneDrive_Binxu\OneDrive\PC_space_tuning\Exp2_chan06";
for pref_chan = 1:size(rasters,1)
figure(6);clf
set(gcf, 'position', [131         338        2109         616]);
%%
score_mat = zeros(11,11,5);
cnt_mat = zeros(11,11);
for i = -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters(pref_chan, :, img_idx);
        scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
        score_mat(i+6, j+6, 1:length(img_idx)) = scores;
    end
end
figure(3)
imagesc(-90:18:90, -90:18:90, sum(score_mat,3)./cnt_mat)
ylabel("PC 2 degree")
xlabel("PC 3 degree")
title("Tuning map on PC2 3 subspace")
shading flat
axis image
colorbar

figure(6)
subplot(131)
imagesc(-90:18:90, -90:18:90, sum(score_mat,3)./cnt_mat)
ylabel("PC 2 degree")
xlabel("PC 3 degree")
title("Tuning map on PC2 3 subspace")
shading flat
axis image
colorbar
%%
score_mat = zeros(11,11,5);
cnt_mat = zeros(11,11);
for i = -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_PC49_%d_PC50_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters(pref_chan, :, img_idx);
        scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
        score_mat(i+6, j+6, 1:length(img_idx)) = scores;
    end
end
figure(4)
imagesc(-90:18:90, -90:18:90, sum(score_mat,3)./cnt_mat)
ylabel("PC 49 degree")
xlabel("PC 50 degree")
title("Tuning map on PC49 50 subspace")
shading flat
axis image
colorbar

figure(6)
subplot(132)
imagesc(-90:18:90, -90:18:90, sum(score_mat,3)./cnt_mat)
ylabel("PC 49 degree")
xlabel("PC 50 degree")
title("Tuning map on PC49 50 subspace")
shading flat
axis image
colorbar
%%
score_mat = zeros(11,11,5);
cnt_mat = zeros(11,11);
for i = -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_RND1_%d_RND2_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters(pref_chan, :, img_idx);
        scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
        score_mat(i+6, j+6, 1:length(img_idx)) = scores;
    end
end
figure(5)
imagesc(-90:18:90, -90:18:90, mean(score_mat,3))
ylabel("Rand vector 1 degree")
xlabel("Rand vector 2 degree")
title("Tuning map on Random Vector (outside first 50 PCs) subspace")
shading flat
axis image
colorbar

figure(6)
subplot(133)
imagesc(-90:18:90, -90:18:90, mean(score_mat,3))
ylabel("Rand vector 1 degree")
xlabel("Rand vector 2 degree")
title("Tuning map on Random Vector (outside first 50 PCs) subspace")
shading flat
axis image
colorbar
%%

saveas(6, fullfile(savepath, sprintf("chan%02d_PC_tune_cmp.png", pref_chan)))
saveas(3, fullfile(savepath, sprintf("chan%02d_PC23_tune.png", pref_chan)))
saveas(4, fullfile(savepath, sprintf("chan%02d_PC4950_tune.png", pref_chan)))
saveas(5, fullfile(savepath, sprintf("chan%02d_RND_tune.png", pref_chan)))

end


% %%
% Unit_id = [meta.spikeID];
% % char(97) = 'a'
% cnt = 0;
% cur_chan = 1;
% for i = 1:length(Unit_id)
%     if cur_chan == Unit_id(i)
%         cnt = cnt + 1;
%     else
%         cur_chan = Unit_id(i);
%         cnt = 1;
%     end
%     sprintf()
% end
%% Modulate the contrast by the score of firing
img_folder = "\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-03a-beto";
cnt = 1;
img_list={};
norm_score_mat = sum(score_mat,3)./cnt_mat;
norm_score_mat = (norm_score_mat - min(norm_score_mat,[],'all')) / (max(norm_score_mat,[],'all') - min(norm_score_mat,[],'all'));
for i = -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        % cnt_mat(i+6, j+6) = length(img_idx);
        cur_img = single(imread(fullfile(img_folder, [cur_fn, '.jpg']))) / 255;
        %norm_score_mat(i+6, j+6)
%         psths = rasters(channel, :, img_idx);
%         scores = squeeze(mean(psths(1, 50:150, :))- mean(psths(1, 1:50, :)) );
        score = norm_score_mat(i+6, j+6);
        img_list{cnt} = score * cur_img + (1 - score) * 0.5 * ones(size(cur_img)); % imresize seems not working
        cnt = cnt + 1;
    end
end
figure(11)
montage(img_list, 'Size', [11 11])
%% Modulate the size by the score of firing