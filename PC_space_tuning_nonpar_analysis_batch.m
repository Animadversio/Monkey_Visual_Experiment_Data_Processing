% Batch processing code  calculating the statistics (t, 1way ANOVA, 2way ANOVA) 
% and generate annotated figures for each and every experiment 
%% 
storedStruct = load("D:\\Manifold_Exps.mat");
%%
%load("D:\\Beto64chan-02102019-003_formatted");
% norm_arr = [328, 326, 269, 329, 401, 386, 300];
% pref_chan_arr = [29, 6, 5, 20, 19, 13, 28]; 
Set_Exp_Specs
Vis_Images = true;
global sphere_norm Trials channel rasters ang_step Reps meta
%%
for Expi=6:10
    rasters = storedStruct.rasters{Expi};
    Trials = storedStruct.Trials{Expi};
    meta = storedStruct.meta{Expi};
Stat_summary = {};
pref_chan = pref_chan_arr(Expi);
sphere_norm = norm_arr(Expi);
ang_step = 18;
Reps = 11; % can be any number LARGER than the largest repitition. or there may be problems caused by NAN and 0 filling
% savepath = "C:\Users\ponce\OneDrive\Desktop\OneDrive_Binxu\OneDrive\PC_space_tuning\Exp1_chan29";
savepath = sprintf("C:\\Users\\ponce\\OneDrive - Washington University in St. Louis\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
mkdir(savepath);
unit_name_arr = generate_unit_labels(meta.spikeID, savepath);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters, savepath);

if Vis_Images
    PC23_img_list = get_image_arr('norm_%d_PC2_%d_PC3_%d');
    PC4950_img_list = get_image_arr('norm_%d_PC49_%d_PC50_%d');
    RND12_img_list = get_image_arr('norm_%d_RND1_%d_RND2_%d');

    figure(5);set(5,'position', [ 326         629        2089         254])
    montage(PC23_img_list, 'Size', [11, 11]);
    title(sprintf("Exp%d pref chan%02d PC23 hemisphere", Expi, pref_chan))
    ylabel("PC3 axis");xlabel("PC2 axis");
    saveas(5, fullfile(savepath, sprintf("norm_%d_PC23.jpg", sphere_norm)))

    figure(6);set(6,'position', [ 326         629        2089         254])
    montage(PC23_img_list, 'Size', [11, 11]);
    title(sprintf("Exp%d pref chan%02d PC4950 hemisphere", Expi, pref_chan))
    ylabel("PC50 axis");xlabel("PC49 axis");
    saveas(6, fullfile(savepath, sprintf("norm_%d_PC4950.jpg", sphere_norm)))

    figure(7);set(7,'position', [ 326         629        2089         254])
    montage(PC23_img_list, 'Size', [11, 11]);
    title(sprintf("Exp%d pref chan%02d RND12 hemisphere", Expi, pref_chan))
    ylabel("RND2 axis");xlabel("RND1 axis");
    saveas(7, fullfile(savepath, sprintf("norm_%d_RND12.jpg", sphere_norm)))
end
for channel = 1:size(rasters,1)
chan_label_str = sprintf("Exp%d Channel %s", Expi, unit_name_arr{channel});
if Vis_Images 
    figure(1);clf;set(1, 'position', [ 805         197        1559         781]);
    figure(2);clf;set(2, 'position', [ 805         197        1559         781]);
    figure(3);clf;set(3, 'position', [ 805         197        1559         781]);
    figure(4);clf;set(4, 'position', [73  -40  2418  705]);
else
    figure(1);clf;set(1, 'position', [304    12   560   577]);
    figure(2);clf;set(2, 'position', [304    12   560   577]);
    figure(3);clf;set(3, 'position', [304    12   560   577]);
    figure(4);clf;set(4, 'position', [73  -40  2418  705]);
end
%% PC12
[score_mat, ~, summary, stat_str] = get_stats_from_result('norm_%d_PC2_%d_PC3_%d');
Stat_summary{channel, 1} = summary;
disp(summary)
% visualize
figure(1);
if Vis_Images, ax1 = subplot(1,2,1); end
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("PC 2 degree");xlabel("PC 3 degree")
title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
shading flat;axis image;colorbar
if Vis_Images
    frame_img_list = score_frame_image_arr(PC23_img_list, nanmean(score_mat, 3)...
        , caxis(ax1), colormap(ax1), 20);
    ax2 = subplot(1,2,2);
    montage(frame_img_list', 'Size', [11, 11]);
    set(ax2, 'Position', [0.5303    0.0500    0.4347    0.9049]);
    set(ax1, 'position', [0.0500    0.0700    0.4018    0.8050]);
end

figure(4)
if 0 ax1=subplot(2,3,1); else subplot(131), end 
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("PC 2 degree");xlabel("PC 3 degree")
title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
shading flat;axis image;colorbar
% if Vis_Images
%     ax2 = subplot(2,3,4);
%     montage(frame_img_list', 'Size', [11, 11]);
% end
%% PC4950
[score_mat, ~, summary, stat_str] = get_stats_from_result('norm_%d_PC49_%d_PC50_%d');
Stat_summary{channel, 2} = summary;
disp(summary)
    
figure(2);clf;
if Vis_Images, ax1 = subplot(1,2,1); end
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("PC 49 degree");xlabel("PC 50 degree")
title([chan_label_str, "Tuning map on PC49 50 subspace", stat_str])
shading flat;axis image;colorbar
if Vis_Images
    frame_img_list = score_frame_image_arr(PC4950_img_list, nanmean(score_mat, 3)...
        , caxis(ax1), colormap(ax1), 20);
    ax2 = subplot(1,2,2);
    montage(frame_img_list', 'Size', [11, 11]);
    set(ax2, 'Position', [0.5303    0.0500    0.4347    0.9049]);
    set(ax1, 'position', [0.0500    0.0700    0.4018    0.8050]);
end

figure(4)
if 0 ax1=subplot(2,3,2); else subplot(132), end 
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("PC 49 degree");xlabel("PC 50 degree")
title([chan_label_str, "Tuning map on PC49 50 subspace", stat_str])
shading flat;axis image;colorbar
% if Vis_Images
%     ax2 = subplot(2,3,5);
%     montage(frame_img_list', 'Size', [11, 11]);
% end
%% RND12
[score_mat, ~, summary, stat_str] = get_stats_from_result('norm_%d_RND1_%d_RND2_%d');
Stat_summary{channel, 3} = summary;
disp(summary)

figure(3);clf;
if Vis_Images, ax1 = subplot(1,2,1); end
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("Rand vector 1 degree");xlabel("Rand vector 2 degree")
title([chan_label_str, "Tuning map on Random Vector (outside first 50 PCs) subspace", stat_str])
shading flat;axis image;colorbar
if Vis_Images
    frame_img_list = score_frame_image_arr(RND12_img_list, nanmean(score_mat, 3)...
        , caxis(ax1), colormap(ax1), 20);
    ax2 = subplot(1,2,2);
    montage(frame_img_list', 'Size', [11, 11]);
    set(ax2, 'Position', [0.5303    0.0500    0.4347    0.9049]);
    set(ax1, 'position', [0.0500    0.0700    0.4018    0.8050]);
end

figure(4)
if 0 ax1=subplot(2,3,3); else subplot(133), end 
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("Rand vector 1 degree");xlabel("Rand vector 2 degree")
title([chan_label_str, "Tuning map on Random Vector (outside first 50 PCs) subspace", stat_str])
shading flat;axis image;colorbar
% if Vis_Images
%     ax2 = subplot(2,3,6);
%     montage(frame_img_list', 'Size', [11, 11]);
% end
%%
saveas(4, fullfile(savepath, sprintf("chan%02d_PC_tune_cmp_stat.png", channel)))
saveas(1, fullfile(savepath, sprintf("chan%02d_PC23_tune_stat.png", channel)))
saveas(2, fullfile(savepath, sprintf("chan%02d_PC4950_tune_stat.png", channel)))
saveas(3, fullfile(savepath, sprintf("chan%02d_RND_tune_stat.png", channel)))
end
save(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary')
end
%%
function img_list = get_image_arr(name_pattern)
global Trials meta ang_step sphere_norm
img_dir = meta.stimuli; % image storage path online
img_list = cell(11, 11);% 
for j = 1:11
    for i = 1:11
        cur_fn = sprintf(name_pattern, sphere_norm, ...
            (i - 6)*ang_step, (j - 6)*ang_step); % 'norm_%d_PC2_%d_PC3_%d'
        img_idx = find(contains(Trials.imageName, cur_fn));
        if length(img_idx) == 0 % empty space
            img_list{i,j} = [];
        else
            img_list{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
        end
    end
end
end

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

% Have moved to a separate function in its file. 
% function unit_name_arr = generate_unit_labels()
% % Generate the unit labels 17B from the spikeID variable
% global meta
% Unit_id = meta.spikeID;
% unit_name_arr = {}; % name tag for each unit 
% for i = 1:length(Unit_id)
%     cur_chan = Unit_id(i);
%     if sum(Unit_id == cur_chan) == 1
%         unit_name_arr{i} = num2str(cur_chan);
%     else
%         cur_chan = Unit_id(i);
%         rel_idx = find(find(Unit_id == cur_chan) == i);
%         unit_name_arr{i} = [num2str(cur_chan), char(64+rel_idx)];
%     end
% end
% end

