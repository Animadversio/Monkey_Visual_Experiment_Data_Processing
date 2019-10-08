%% Correlated PC space tuning
%%
load('D:\Beto64chan-03102019-001_formatted.mat')
%%
img_names = [Trials.imageName];
syn_image_idx = [];
nat_image_idx = [];
syn_img_cnt = 1;
nat_img_cnt = 1;
cur_gen = 0;
for i=1:numel(img_names)
    id_arr = regexp(img_names{i}, 'block(\d+)_thread\d+_gen_gen(\d+)_(\d+)','tokens');
    if isempty(id_arr)
        nat_image_idx = [nat_image_idx, i];
        % nat_block_id(nat_img_cnt) = block_id(syn_img_cnt);
        nat_gen_id(nat_img_cnt) = cur_gen;
        nat_img_cnt = nat_img_cnt + 1; 
    else
        syn_image_idx = [syn_image_idx, i];
        block_id(syn_img_cnt) = str2num(id_arr{1}{1});
        gen_id(syn_img_cnt) = str2num(id_arr{1}{2});
        img_id(syn_img_cnt) = str2num(id_arr{1}{3});
        cur_gen = gen_id(syn_img_cnt);
        syn_img_cnt = syn_img_cnt + 1;
    end
end
%%
channel = 6;
syn_scores = squeeze(mean(rasters(channel, 10:200, syn_image_idx), 2));
nat_scores = squeeze(mean(rasters(channel, 10:200, nat_image_idx), 2));
figure;hold on 
scatter(gen_id, syn_scores, 'filled', 'MarkerFaceAlpha', 0.4);
%ftmp.MarkerFaceAlpha = 0.2;
scatter(nat_gen_id, nat_scores,'filled', 'MarkerFaceAlpha', 0.4);
%ftmp.MarkerFaceAlpha = 0.2;
%%
figure
subplot(131)
channel = 5;
imagesc(squeeze(rasters(channel, :, nat_image_idx))', [0, 500])
title(num2str(channel))
subplot(132)
channel = 6;
imagesc(squeeze(rasters(channel, :, nat_image_idx))', [0, 300])
title(num2str(channel))
subplot(133)
channel = 7;
imagesc(squeeze(rasters(channel, :, nat_image_idx))', [0, 400])
title(num2str(channel))
%% Scatter plot the scores for nature stimuli and synthesized stimuli
figure;
subplot(131);hold on 
channel = 5;
syn_scores = squeeze(mean(rasters(channel, 10:200, syn_image_idx), 2));
nat_scores = squeeze(mean(rasters(channel, 10:200, nat_image_idx), 2));
scatter(gen_id, syn_scores, 'filled', 'MarkerFaceAlpha', 0.4)
scatter(nat_gen_id, nat_scores, 'filled', 'MarkerFaceAlpha', 0.4)
title(num2str(channel))
subplot(132);hold on 
channel = 6;
syn_scores = squeeze(mean(rasters(channel, 10:200, syn_image_idx), 2));
nat_scores = squeeze(mean(rasters(channel, 10:200, nat_image_idx), 2));
scatter(gen_id, syn_scores, 'filled', 'MarkerFaceAlpha', 0.4)
scatter(nat_gen_id, nat_scores, 'filled', 'MarkerFaceAlpha', 0.4)
title(num2str(channel))
subplot(133);hold on 
channel = 7;
syn_scores = squeeze(mean(rasters(channel, 10:200, syn_image_idx), 2));
nat_scores = squeeze(mean(rasters(channel, 10:200, nat_image_idx), 2));
scatter(gen_id, syn_scores, 'filled', 'MarkerFaceAlpha', 0.4)
scatter(nat_gen_id, nat_scores, 'filled', 'MarkerFaceAlpha', 0.4)
title(num2str(channel))
