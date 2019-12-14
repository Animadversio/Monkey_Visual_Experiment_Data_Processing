% Code to do basic analysis on PC + Pasupathy patch manifold experiment
% analysis
global  Trials rasters channel sphere_norm ang_step Reps
%storedStruct = load("D:\\Manifold_Exps.mat");
% Load code from "D:\Poncelab_Github\office-main\Project_Selectivity_Beto_loadRaw.m"
Reps = 15; % constant for maximum number of repetitions (as long as it's larger than the maximum, it's fine)
Set_Exp_Specs; 
%%
for Expi = 44:45 % universal manifold experiment identifier
% Load the dataset 
% meta = storedStruct.meta{Expi};
% rasters = storedStruct.rasters{Expi};
% Trials = storedStruct.Trials{Expi};
% Load new dataset
% meta = meta_new{2*(Expi-24)-1};
% rasters = rasters_new{2*(Expi-24)-1};
% Trials = Trials_new{2*(Expi-24)-1};
% sphere_norm = Pasu_norm_arr(Expi-10); % Load the specific information
Triali = 3*(Expi - 43)-2;%Expi - 34;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi_tab = ExpSpecTable_Aug.Expi(exp_rowi);
assert(Expi_tab == Expi, "Data Expi doesn't match that in exp record! Check your indexing or record.")
fprintf("Processing Exp %d, %s\n", Expi, meta.comments)

sphere_norm = norm_arr(Expi);
pref_chan = pref_chan_arr(Expi);

ang_step = 18;
% Save basic info
savepath = sprintf("C:\\Users\\ponce\\OneDrive - Washington University in St. Louis\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
mkdir(savepath);
unit_name_arr = generate_unit_labels(meta.spikeID, savepath); % Generate readable labels for each channel

% Load and visualize the Pasupathy shapes
img_dir = 'N:\Stimuli\2019-Manifold\pasupathy-wg-f-4-ori';
img_list = cell(51, 4);%{};
for j = 1:4
    for i = 1:51
        cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
        img_idx = find(contains(Trials.imageName, cur_fn));
        if length(img_idx) == 0 % empty space
            img_list{i,j} = [];
        else
            img_list{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
        end
    end
end
clear cur_fn i j img_idx
figure(1);set(1,'position', [ 326         629        2089         254])
montage(img_list, 'Size', [4, 51]);
saveas(1, fullfile(savepath, "pasupathy_images.jpg"))
% Load and visualize the Pasupathy shapes
img_dir = meta.stimuli; % image storage path online
evolv_img_list = cell(11, 11);%{};
for j = 1:11
    for i = 1:11
        cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', sphere_norm, ...
            (i - 6)*ang_step, (j - 6)*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        if length(img_idx) == 0 % empty space
            evolv_img_list{i,j} = [];
        else
            evolv_img_list{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
        end
    end
end
clear cur_fn i j img_idx
figure(2);set(2,'position', [ 326         629        2089         254])
montage(evolv_img_list, 'Size', [11, 11]);
title(sprintf("Exp%d pref chan%02d PC23 hemisphere", Expi, pref_chan))
ylabel("PC3 axis");xlabel("PC2 axis");
saveas(2, fullfile(savepath, sprintf("norm_%d_PC23.jpg", sphere_norm)))
% Record the basic statisitcs or load them from the saved file. 
Stat_summary = {};
% load(fullfile(savepath, "KentFit_Stats.mat"), 'Param_summary') % Load computed Kent Fit 
% load(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary') % Load Stats
for channel = 1:size(rasters, 1)
chan_label_str = sprintf("Exp%d Channel %s", Expi, unit_name_arr{channel});

figure(3);clf; set(3, 'position', [ 805         197        1559         781]);
figure(4);clf; set(4, 'position', [  73         181        2479         593]);
%%
[score_mat, bsl_mat, summary, stat_str] = get_stats_from_result('norm_%d_PC2_%d_PC3_%d');
Stat_summary{channel, 1} = summary;
figure(3);clf;
ax1 = subplot(1,2,1);
imagesc(-90:18:90, -90:18:90, nanmean(score_mat, 3)) % sum(score_mat,3)./cnt_mat
ylabel("PC 2 degree");xlabel("PC 3 degree")
title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
shading flat;axis image;colorbar
% Get the colormap from the heatmap in order to visualize the map below. 
frame_img_list = score_frame_image_arr(evolv_img_list, nanmean(score_mat, 3)...
    , caxis(ax1), colormap(ax1), 20);
ax2 = subplot(1,2,2);
montage(frame_img_list', 'Size', [11, 11]);
set(ax2, 'Position', [0.5303    0.0500    0.4347    0.9049]);
set(ax1, 'position', [0.0500    0.0700    0.4018    0.8050]);
%%
[score_mat, bsl_mat, summary, stat_str] = get_pasu_tuning_stats();
Stat_summary{channel, 2} = summary;
% visualize the pasupathy tuning 
hl_mat = nanmean(score_mat,3);
figure(4);clf % 2 panel plot for response and images
ax1 = subplot(2,1,1);
imagesc(hl_mat');
xlabel("Shape id");ylabel("orientation")
yticks(1:4);yticklabels(0:90:270)
title([chan_label_str, "Tuning map on Pasupathy patches", stat_str])
shading flat;axis image;colorbar
% Get the colormap from the heatmap in order to visualize the map below. 
% cmap = colormap(); 
% [Cmin, Cmax] = caxis();
frame_img_list = score_frame_image_arr(img_list, hl_mat, ...
    caxis(ax1), colormap(ax1), 50);
ax2 = subplot(2,1,2);
montage(frame_img_list, 'Size', [4, 51]);
set(ax1,'position',[0.0500    0.5038    0.9050    0.2977])
set(ax2,'position',[0.0500    0.0800    0.9050    0.2977])
%
saveas(3, fullfile(savepath, sprintf("chan%02d_PC23_tune.png", channel)))
saveas(4, fullfile(savepath, sprintf("chan%02d_Pasu_tune.png", channel)))
% save_to_pdf(3, fullfile(savepath, sprintf("chan%02d_PC23_tune.pdf", channel)))
% save_to_pdf(4, fullfile(savepath, sprintf("chan%02d_Pasu_tune.pdf", channel)))
end
save(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary')
end
% 
% %% Modulate the contrast by the score of firing
% img_folder = "\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-03a-beto";
% cnt = 1;
% img_list={};
% norm_score_mat = sum(score_mat,3)./cnt_mat;
% norm_score_mat = (norm_score_mat - min(norm_score_mat,[],'all')) / (max(norm_score_mat,[],'all') - min(norm_score_mat,[],'all'));
% for i = -5:5
%     for j = -5:5
%         cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', sphere_norm, i*ang_step, j*ang_step);
%         img_idx = find(contains(Trials.imageName, cur_fn));
%         % cnt_mat(i+6, j+6) = length(img_idx);
%         cur_img = single(imread(fullfile(img_folder, [cur_fn, '.jpg']))) / 255;
%         %norm_score_mat(i+6, j+6)
% %         psths = rasters(channel, :, img_idx);
% %         scores = squeeze(mean(psths(1, 50:150, :))- mean(psths(1, 1:50, :)) );
%         score = norm_score_mat(i+6, j+6);
%         img_list{cnt} = score * cur_img + (1 - score) * 0.5 * ones(size(cur_img)); % imresize seems not working
%         cnt = cnt + 1;
%     end
% end
% figure(11)
% montage(img_list, 'Size', [11 11])
%% Getting Stats
function [score_mat, bsl_mat, summary, stat_str] = get_stats_from_result(name_pattern)
    global  Trials rasters channel sphere_norm ang_step Reps
    score_mat = nan(11,11,Reps); 
    bsl_mat = nan(11,11,Reps);  % make sure the 3 dimension has larger size than repitition number! 
    cnt_mat = zeros(11,11); 
    id_mat = zeros(11,11); % record the id correspond to i,j
    for i =  -5:5
        for j =  -5:5
            cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
            img_idx = find(contains(Trials.imageName, cur_fn));
            cnt_mat(i+6, j+6) = length(img_idx);
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 50:150, :))- mean(psths(1, 1:50, :)) );
            baseline = squeeze(mean(psths(1, 50:150, :)));
            score_mat(i+6, j+6, 1:length(img_idx)) = scores;
            bsl_mat(i+6, j+6, 1:length(img_idx)) = baseline;
            if j ~= 5 && j ~= -5
                id = 11 * i + j;
            elseif j == 5 
                id = 5;
            elseif j == -5
                id = -5;
            end
            id_mat(i+6, j+6) = id;
        end
    end
    [theta_mat, phi_mat] = meshgrid(ang_step*(-5:5), ang_step*(-5:5));
    mean_fr_mat = bsl_mat + score_mat;
    id_vec_nan = reshape(repmat(id_mat, 1,1, Reps), 1, []);
    score_vec_nan = reshape(score_mat, 1, []);
    bsl_vec_nan = reshape(bsl_mat, 1, []);
    mean_fr_vec_nan = bsl_vec_nan + score_vec_nan;
    % Do statistics
    [p,tbl,stats] = anova1(score_vec_nan, id_vec_nan, 'off');
    stats.F = tbl{2,5};
    stats.p = p;
    summary.anova_F = stats.F;
    summary.anova_p = stats.p; 
    %
    [p2,tbl2,stats2] = anovan(score_vec_nan, {reshape(repmat(theta_mat, 1,1,Reps),1,[]), ...
                              reshape(repmat(phi_mat, 1,1,Reps),1,[])}, 'model', 'interaction','display' ,'off');
    stats2.p = p2;
    stats2.F = [tbl2{2:4,6}];
    summary.anova2_p = p2; % p for theta, phi and interaction 
    summary.anova2_F = [tbl2{2:4,6}]; % F for theta, phi and interaction
    %
    [~,P,CI] = ttest(mean_fr_vec_nan, bsl_vec_nan);
    summary.t_p = P;
    summary.t_CI = CI;
    % visualize
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],CI(1), CI(2), P, stats.F, stats.p, ...
            stats2.F(1),stats2.p(1), stats2.F(2),stats2.p(2),stats2.F(3),stats2.p(3));

end

function [score_mat, bsl_mat, summary, stat_str] = get_pasu_tuning_stats() %name_pattern is defaultly 'pasu_%02d_ori_%02d_wg_f'
    global  Trials rasters channel sphere_norm ang_step Reps
    score_mat = nan(51,4,Reps); 
    bsl_mat = nan(51,4,Reps);  % make sure the 3 dimension has larger size than repitition number! 
    cnt_mat = zeros(51,4); 
    id_mat = zeros(51,4); % record the id correspond to i,j
    for i = 1:51
        for j = 1:4
            cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
            img_idx = find(contains(Trials.imageName, cur_fn));
            cnt_mat(i, j) = length(img_idx);
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
            score_mat(i, j, 1:length(img_idx)) = scores;
            baseline = squeeze(mean(psths(1, 50:150, :)));
            bsl_mat(i, j, 1:length(img_idx)) = baseline;
            if i==1 || i==2 || i==3
                id = 4 * (i-1) + 1;
            else 
                id = 4 * (i-1) + j;
            end
            id_mat(i, j) = id;
        end
    end
    [shape_mat, rot_mat] = meshgrid(1:51, 1:4);
    mean_fr_mat = bsl_mat + score_mat;
    id_vec_nan = reshape(repmat(id_mat, 1,1, Reps), 1, []);
    score_vec_nan = reshape(score_mat, 1, []);
    bsl_vec_nan = reshape(bsl_mat, 1, []);
    mean_fr_vec_nan = bsl_vec_nan + score_vec_nan;
    % Do statistics
    [p,tbl,stats] = anova1(score_vec_nan, id_vec_nan, 'off');
    stats.F = tbl{2,5};
    stats.p = p;
    summary.anova_F = stats.F;
    summary.anova_p = stats.p; 
    %
    [p2,tbl2,stats2] = anovan(score_vec_nan, {reshape(repmat(shape_mat, 1,1,Reps),1,[]), ...
                              reshape(repmat(rot_mat, 1,1,Reps),1,[])}, 'model', 'interaction','display' ,'off');
    stats2.p = p2;
    stats2.F = [tbl2{2:4,6}];
    summary.anova2_p = p2; % p for theta, phi and interaction 
    summary.anova2_F = [tbl2{2:4,6}]; % F for theta, phi and interaction
    %
    [~,P,CI] = ttest(mean_fr_vec_nan, bsl_vec_nan);
    summary.t_p = P;
    summary.t_CI = CI;
    % visualize
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Shape, F=%.2f(%.3f), Rotation, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],CI(1), CI(2), P, stats.F, stats.p, ...
            stats2.F(1),stats2.p(1), stats2.F(2),stats2.p(2),stats2.F(3),stats2.p(3));

end

% Plot the padded image montage 
% Note this can be a seperate function. 
function frame_img_list = score_frame_image_arr(img_list, score_mat, clim, cmap, LineWidth)
% Use the cmap and clim to map values in hl_mat to color, and form color frame
% for corresponding image in the img_list. 
% Argument 
% img_list is a image array, same shape as hl_mat
% hl_mat is a score matrix, with nan is images with no observation.
% clim is the limit of value mapped to color
% cmap is a K-by-3 matrix coding the RGB values in 1:64. 
Cmin = clim(1); Cmax = clim(2);
assert(all(size(img_list) == size(score_mat)), ...
    "Score matrix and image cell array size doesn't matach")
% LineWidth = 50; % Key parameter controlling the width of margin, can be different for different
frame_img_list = cell(size(img_list)); % the list storing the padde image
for j = 1:size(img_list, 2)
    for i = 1:size(img_list, 1)
%         cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
%         img_idx = find(contains(Trials.imageName, cur_fn));
        if isnan(score_mat(i,j)) % check if there is any score
            frame_img_list{i,j} = [];
        else
            scale_val = (score_mat(i,j) - Cmin) / (Cmax - Cmin);
            c = interp1(cmap, scale_val * (size(cmap, 1) - 1) + 1); % Note, interpolation can be done from 1-64, not from 0
            
            pad_img = padarray(img_list{i,j}, [2*LineWidth, 2*LineWidth], 0);
            tmp_img = insertShape(pad_img, ...
                'Rectangle', [LineWidth,LineWidth,...
                            size(img_list{i,j},2)+2*LineWidth,...
                            size(img_list{i,j},1)+2*LineWidth], ...
                'LineWidth', 2 * LineWidth, 'Color', 256 * c);
            frame_img_list{i,j} = tmp_img;
        end
    end
end

end