%load(fullfile("\\storage1.ris.wustl.edu\crponce\Active\Data-Ephys-MAT","Beto64chan-02102019-003_formatted.mat"))
%load(fullfile("\\storage1.ris.wustl.edu\crponce\Active\Data-Ephys-MAT","Beto64chan-03102019-002_formatted.mat"))
% meta,rasters,lfps,Trials
meta = meta_new{5};
rasters = rasters_new{5};
Trials = Trials_new{5};
img_names = unique(Trials.imageName);
%% Loading code from "D:\Poncelab_Github\office-main\Project_Selectivity_Beto_loadRaw.m"
global  Trials rasters channel sphere_norm ang_step Reps
Reps = 11;
Set_Exp_Specs;
%% Set 
img_dir = 'S:\Stimuli\2019-Manifold\pasupathy-wg-f-4-ori';
img_list = {};
ctr = 1;
for j = 1:4
    for i = 1:51
        cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
        img_idx = find(contains(Trials.imageName, cur_fn));
        if length(img_idx) == 0
            img_list{i,j} = [];
        else
            img_list{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
        end
        ctr = ctr+1;
    end
end
figure(1);set(1,'position', [ 326         629        2089         254])
montage(img_list, 'Size', [4, 51]);

for Expi = 11:13
Stat_summary = {};
ang_step = 18;
meta = meta_new{2*(Expi-10)-1};
rasters = rasters_new{2*(Expi-10)-1};
Trials = Trials_new{2*(Expi-10)-1};
sphere_norm = Pasu_norm_arr(Expi-10); % 269 Day3 % 326 Daye % 328 Day 1
pref_chan = Pasu_pref_chan_arr(Expi-10);

savepath = sprintf("C:\\Users\\ponce\\OneDrive\\Desktop\\OneDrive_Binxu\\OneDrive\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
mkdir(savepath);
saveas(1, fullfile(savepath, ["pasupathy_images.jpg"]))
unit_name_arr = generate_unit_labels(meta.spikeID);
% load(fullfile(savepath, "KentFit_Stats.mat"), 'Param_summary') % Load computed Kent Fit 
% load(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary') % Load Stats
for channel = 1:size(rasters, 1)
    chan_label_str = sprintf("Exp%d Channel %s", Expi, unit_name_arr{channel});
% figure(6);clf
% set(gcf, 'position', [131         338        2109         616]);
figure(3);clf; set(3, 'position', [1000         558         560         420]);
figure(4);clf; set(4, 'position', [ 73         181        2479         593]);
%%
score_mat = zeros(11,11,11);
cnt_mat = zeros(11,11);
for i = -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters(channel, :, img_idx);
        scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
        score_mat(i+6, j+6, 1:length(img_idx)) = scores;
    end
end
[score_mat, bsl_mat, summary, stat_str] = get_stats_from_result('norm_%d_PC2_%d_PC3_%d');
Stat_summary{channel, 1} = summary;
figure(3);clf
imagesc(-90:18:90, -90:18:90, nanmean(score_mat, 3)) % sum(score_mat,3)./cnt_mat
ylabel("PC 2 degree")
xlabel("PC 3 degree")
title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
shading flat
axis image
colorbar
%%
score_mat = zeros(51,4,11);
cnt_mat = zeros(51,4);
for i = 1:51
    for j = 1:4
        cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
        img_idx = find(contains(Trials.imageName, cur_fn));
        cnt_mat(i, j) = length(img_idx);
        psths = rasters(channel, :, img_idx);
        scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
        score_mat(i, j, 1:length(img_idx)) = scores;
    end
end
[score_mat, bsl_mat, summary, stat_str] = get_pasu_tuning_stats();
Stat_summary{channel, 2} = summary;
% figure(4);clf
% imagesc(1:51, 1:4, nanmean(score_mat, 3)');% (sum(score_mat,3)./cnt_mat)'
% xlabel("Shape id")
% ylabel("orientation")
% title([chan_label_str, "Tuning map on Pasupathy patches", stat_str])
% shading flat
% axis image
% colorbar
% hl_mat = nanmean(score_mat,3);
figure(4);clf
subplot(2,1,1)
h = imagesc(hl_mat');
xlabel("Shape id")
ylabel("orientation")
title([chan_label_str, "Tuning map on Pasupathy patches", stat_str])
shading flat
colorbar
axis image;
cmap = colormap();
[Cmin, Cmax] = caxis();
mod_img_list = cell(51,4);
for j = 1:4
    for i = 1:51
        cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
        img_idx = find(contains(Trials.imageName, cur_fn));
        if length(img_idx) == 0
            mod_img_list{i,j} = [];
        else
            scale_val = (hl_mat(i,j) - Cmin) / (Cmax - Cmin);
            c = interp1(cmap, scale_val * (size(cmap, 1) - 1) + 1);
            LineWidth = 50;
            pad_img = padarray(img_list{i,j}, [2*LineWidth, 2*LineWidth], 0);
            tmp_img = insertShape(pad_img, ...
                'Rectangle', [LineWidth,LineWidth,size(img_list{i,j})+2*LineWidth], ...
                'LineWidth', 2 * LineWidth, 'Color', 256 * c);
            mod_img_list{i,j} = tmp_img;
        end
    end
end
subplot(2,1,2)
montage(mod_img_list, 'Size', [4, 51]);

%%
saveas(3, fullfile(savepath, sprintf("chan%02d_PC23_tune.png", channel)))
saveas(4, fullfile(savepath, sprintf("chan%02d_Pasu_tune.png", channel)))
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