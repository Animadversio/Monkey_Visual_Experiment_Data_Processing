clearvars -except meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug  ExpSpecTable_Aug_alfa ExpRecord
%%
% Code to do basic analysis on PC + Pasupathy patch manifold experiment
% analysis
global  Trials rasters channel ang_step Reps
Reps = 15; % constant for maximum number of repetitions (as long as it's larger than the maximum, it's fine)
%%
%storedStruct = load("D:\\Manifold_Exps.mat");
% Load code from "D:\Poncelab_Github\office-main\Project_Selectivity_Beto_loadRaw.m"
% Set_Exp_Specs; 
Set_Path;
Animal = "Alfa";
expftr = ExpRecord.Expi>=34 & ...%ExpRecord.Expi<=40 & 
    contains(ExpRecord.expControlFN,"selectivity") & ...
     contains(ExpRecord.Exp_collection, "Manifold");
[meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr),Animal); 
%%
figure(1);clf; set(1,'position', [ 326         629        2089         254]); % all pasu images montage 
figure(2);clf; set(2,'position', [ 1070          48         910         820]); % all manifold images montaged
figure(5);clf; set(5,'position', [ 1070          48         910         820]); % all manifold images montaged
figure(6);clf; set(6,'position', [ 1070          48         910         820]); % all manifold images montaged
figure(3);clf; set(3, 'position', [ 805         197        1559         781]); % neural response to manifold images
figure(4);clf; set(4, 'position', [  73         181        2479         593]); % neural response to pasu images
figure(7);clf; set(7, 'position', [ 805         197        1559         781]); % neural response to manifold images PC 49 50 
figure(8);clf; set(8, 'position', [ 805         197        1559         781]); % neural response to manifold images RND 1 2 
%%
Reps = 15;
for Triali = 8:length(meta_new) % universal manifold experiment identifier
% Load the dataset 
% meta = storedStruct.meta{Expi};
% rasters = storedStruct.rasters{Expi};
% Trials = storedStruct.Trials{Expi};
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
Expi = ExpRecord.Expi(exp_rowi);% Check the Expi match number
fprintf("Processing %s Exp %d:\n", Animal, Expi)
disp(ExpRecord.comments(exp_rowi))

% assert(Expi_tab == Expi, "Data Expi doesn't match that in exp record! Check your indexing or record.")
sphere_norm = infer_norm_from_imgname(Trials.imageName, "PC2");
pref_chan = Trials.TrialRecord.User.prefChan;
if Expi == 16, pref_chan = 20; end % Change the mistake in input this variable.
unit_in_pref_chan = 1; % TODO find this number somewhere
% unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
% thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);

% sphere_norm = norm_arr(Expi);
% pref_chan = pref_chan_arr(Expi);

ang_step = 18;
% Save basic info
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
savepath = sprintf("C:\\Users\\ponce\\OneDrive - Washington University in St. Louis\\PC_space_tuning\\%s_Exp%d_chan%02d", Animal, Expi, pref_chan);
mkdir(savepath);
unit_name_arr = generate_unit_labels(meta.spikeID, savepath); % Generate readable labels for each channel
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters, savepath);

% Load and visualize the Pasupathy shapes
% Pre load the pasu images and montage them in figure(1)
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
set(0,'CurrentFigure',1); %clf; %
montage(img_list, 'Size', [4, 51]);
saveas(1, fullfile(savepath, "pasupathy_images.jpg"))
%% Load and visualize the Synthesized images
% if Expi ==18,  meta.stimuli ="N:\Stimuli\2019-Manifold\alfa-191209a-selectivity";  end
% if Expi ==19,  meta.stimuli ="N:\Stimuli\2019-Manifold\alfa-191210a-selectivity";  end
% if Expi ==20,  meta.stimuli ="N:\Stimuli\2019-Manifold\alfa-191211a-selectivity";  end
% if Expi ==21,  meta.stimuli ="N:\Stimuli\2019-Manifold\alfa-191212a-selectivity";  end
img_dir = meta.stimuli; % image storage path online
% if ~contains(img_dir, "N:\") && ~contains(img_dir, "\\storage1.ris.wustl.edu\crponce\Active\") && contains(img_dir(1:8), "Stimuli")
%     img_dir = fullfile("N:\", img_dir);
% end
sphere_norm = infer_norm_from_imgname(Trials.imageName, "PC2");
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
set(0,'CurrentFigure',2); %clf; %
montage(evolv_img_list, 'Size', [11, 11]);
title(sprintf("Exp%d pref chan%02d PC23 hemisphere", Expi, pref_chan))
ylabel("PC3 axis");xlabel("PC2 axis");
saveas(2, fullfile(savepath, sprintf("norm_%d_PC23.jpg", sphere_norm)))
%
sphere_norm = infer_norm_from_imgname(Trials.imageName, "PC49");
did_PC49 = ~isnan(sphere_norm);
if did_PC49
evolv_img_list1 = cell(11, 11);%{};
for j = 1:11
    for i = 1:11
        cur_fn = sprintf('norm_%d_PC49_%d_PC50_%d', sphere_norm, ...
            (i - 6)*ang_step, (j - 6)*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        if length(img_idx) == 0 % empty space
            evolv_img_list1{i,j} = [];
        else
            evolv_img_list1{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
        end
    end
end
clear cur_fn i j img_idx
set(0,'CurrentFigure',5); clf; %
montage(evolv_img_list1, 'Size', [11, 11]);
title(sprintf("Exp%d pref chan%02d PC4950 hemisphere", Expi, pref_chan))
ylabel("PC50 axis");xlabel("PC49 axis");
saveas(5, fullfile(savepath, sprintf("norm_%d_PC4950.jpg", sphere_norm)))
end
% 
sphere_norm = infer_norm_from_imgname(Trials.imageName, "RND1");
did_RND1 = ~isnan(sphere_norm);
if did_RND1
evolv_img_list2 = cell(11, 11);%{};
for j = 1:11
    for i = 1:11
        cur_fn = sprintf('norm_%d_RND1_%d_RND2_%d', sphere_norm, ...
            (i - 6)*ang_step, (j - 6)*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        if length(img_idx) == 0 % empty space
            evolv_img_list2{i,j} = [];
        else
            evolv_img_list2{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
        end
    end
end
clear cur_fn i j img_idx
set(0,'CurrentFigure',6); clf; %
montage(evolv_img_list2, 'Size', [11, 11]);
title(sprintf("Exp%d pref chan%02d RND12 hemisphere", Expi, pref_chan))
ylabel("RND2 axis");xlabel("RND1 axis");
saveas(6, fullfile(savepath, sprintf("norm_%d_RND12.jpg", sphere_norm)))
end


did_Gabor = false;
did_Pasu = false;
if sum(contains(Trials.imageName, "pasu")) > 186, did_Pasu = true; end
if sum(contains(Trials.imageName, "gab")) > 12, did_Gabor = true; end
%%
% Record the basic statisitcs or load them from the saved file. 
Stat_summary = {};
% load(fullfile(savepath, "KentFit_Stats.mat"), 'Param_summary') % Load computed Kent Fit 
% load(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary') % Load Stats
for channel = 1:size(rasters, 1)
chan_label_str = sprintf("%s Exp%d Channel %s", Animal, Expi, unit_name_arr{channel});
%% PC23 plot
sphere_norm = infer_norm_from_imgname(Trials.imageName, "PC2");
[score_mat, bsl_mat, summary, stat_str] = get_stats_from_result('norm_%d_PC2_%d_PC3_%d',sphere_norm);
Stat_summary{channel, 1} = summary; % statistics of PC23
set(0,'CurrentFigure',3); clf;
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
%% Pasupathy plot
if did_Pasu
[score_mat, bsl_mat, summary, stat_str] = get_pasu_tuning_stats();
Stat_summary{channel, 2} = summary; % statistics of Pasu images
% visualize the pasupathy tuning 
hl_mat = nanmean(score_mat,3);
set(0,'CurrentFigure',4); clf;% 2 panel plot for response and images
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
end

if did_PC49
sphere_norm = infer_norm_from_imgname(Trials.imageName, "PC49");
[score_mat, bsl_mat, summary, stat_str] = get_stats_from_result('norm_%d_PC49_%d_PC50_%d',sphere_norm);
Stat_summary{channel, 3} = summary; % statistics of PC49 50 images
set(0,'CurrentFigure',7); clf;
ax1 = subplot(1,2,1);
imagesc(-90:18:90, -90:18:90, nanmean(score_mat, 3)) % sum(score_mat,3)./cnt_mat
ylabel("PC 49 degree");xlabel("PC 50 degree")
title([chan_label_str, "Tuning map on PC49 50 subspace", stat_str])
shading flat;axis image;colorbar
% Get the colormap from the heatmap in order to visualize the map below. 
frame_img_list = score_frame_image_arr(evolv_img_list, nanmean(score_mat, 3)...
    , caxis(ax1), colormap(ax1), 20);
ax2 = subplot(1,2,2);
montage(frame_img_list', 'Size', [11, 11]);
set(ax2, 'Position', [0.5303    0.0500    0.4347    0.9049]);
set(ax1, 'position', [0.0500    0.0700    0.4018    0.8050]);
end

if did_RND1
sphere_norm = infer_norm_from_imgname(Trials.imageName, "RND1");
[score_mat, bsl_mat, summary, stat_str] = get_stats_from_result('norm_%d_RND1_%d_RND2_%d',sphere_norm);
Stat_summary{channel, 4} = summary; % statistics of PC49 50 images
set(0,'CurrentFigure',8); clf;
ax1 = subplot(1,2,1);
imagesc(-90:18:90, -90:18:90, nanmean(score_mat, 3)) % sum(score_mat,3)./cnt_mat
ylabel("RND 1 degree");xlabel("RND 2 degree")
title([chan_label_str, "Tuning map on RND1 2 subspace", stat_str])
shading flat;axis image;colorbar
% Get the colormap from the heatmap in order to visualize the map below. 
frame_img_list = score_frame_image_arr(evolv_img_list, nanmean(score_mat, 3)...
    , caxis(ax1), colormap(ax1), 20);
ax2 = subplot(1,2,2);
montage(frame_img_list', 'Size', [11, 11]);
set(ax2, 'Position', [0.5303    0.0500    0.4347    0.9049]);
set(ax1, 'position', [0.0500    0.0700    0.4018    0.8050]);
end

%
saveas(3, fullfile(savepath, sprintf("PC23_tune_chan%s.png", unit_name_arr{channel})))
if did_Pasu, saveas(4, fullfile(savepath, sprintf("Pasu_tune_chan%s.png", unit_name_arr{channel}))); end
if did_PC49, saveas(7, fullfile(savepath, sprintf("PC4950_tune_chan%s.png", unit_name_arr{channel}))); end
if did_RND1, saveas(8, fullfile(savepath, sprintf("RND12_tune_chan%s.png", unit_name_arr{channel}))); end
% Obsolete 
% save_to_pdf(3, fullfile(savepath, sprintf("chan%02d_PC23_tune.pdf", channel)))
% save_to_pdf(4, fullfile(savepath, sprintf("chan%02d_Pasu_tune.pdf", channel)))
end
save(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary')
PC23Tab = struct2table([Stat_summary{:,1}]);
PC23Tab.unit_name = unit_name_arr;
PC23Tab.channel = meta.spikeID;
PC23Tab.unit_num = unit_num_arr;
writetable(PC23Tab, fullfile(savepath, "PC23_Stats.xlsx"))
if did_Pasu
PasuTab = struct2table([Stat_summary{:,2}]);
PasuTab.unit_name = unit_name_arr;
PasuTab.channel = meta.spikeID;
PasuTab.unit_num = unit_num_arr;
writetable(PasuTab, fullfile(savepath, "Pasu_Stats.xlsx"))
end
if did_PC49
PC49Tab = struct2table([Stat_summary{:,3}]);
PC49Tab.unit_name = unit_name_arr;
PC49Tab.channel = meta.spikeID;
PC49Tab.unit_num = unit_num_arr;
writetable(PC49Tab, fullfile(savepath, "PC49_Stats.xlsx"))
end
if did_RND1
RNDTab = struct2table([Stat_summary{:,4}]);
RNDTab.unit_name = unit_name_arr;
RNDTab.channel = meta.spikeID;
RNDTab.unit_num = unit_num_arr;
writetable(RNDTab, fullfile(savepath, "RND1_Stats.xlsx"))
end
end
%%
PC23Tab = struct2table([Stat_summary{:,1}]);
PC23Tab.unit_name = unit_name_arr;
PC23Tab.channel = meta.spikeID;
PC23Tab.unit_num = unit_num_arr;
writetable(PC23Tab, fullfile(savepath, "PC23_Stats.xlsx"))
if did_Pasu
PasuTab = struct2table([Stat_summary{:,2}]);
PasuTab.unit_name = unit_name_arr;
PasuTab.channel = meta.spikeID;
PasuTab.unit_num = unit_num_arr;
writetable(PasuTab, fullfile(savepath, "Pasu_Stats.xlsx"))
end
if did_PC49
PC49Tab = struct2table([Stat_summary{:,3}]);
PC49Tab.unit_name = unit_name_arr;
PC49Tab.channel = meta.spikeID;
PC49Tab.unit_num = unit_num_arr;
writetable(PC49Tab, fullfile(savepath, "PC49_Stats.xlsx"))
end
if did_RND1
RNDTab = struct2table([Stat_summary{:,4}]);
RNDTab.unit_name = unit_name_arr;
RNDTab.channel = meta.spikeID;
RNDTab.unit_num = unit_num_arr;
writetable(RNDTab, fullfile(savepath, "RND1_Stats.xlsx"))
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
%%

%% Getting Stats
function sphere_norm = infer_norm_from_imgname(imgnames, pattern)
    if nargin == 1
        pattern = "PC2"; % by default it's the PC2 PC3 space
    end
    tmp = cellfun(@(c) regexp(c,strcat("norm_(?<norm>\d*)_", pattern),'names'), imgnames, 'UniformOutput', false);      
    extnorms = cellfun(@(c) str2num(c.norm), tmp(~cellfun('isempty',tmp)));
    sphere_norm = mode(extnorms);
end

function [score_mat, bsl_mat, summary, stat_str] = get_stats_from_result(name_pattern, sphere_norm)
    global  Trials rasters channel ang_step Reps
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
            baseline = squeeze(mean(psths(1, 1:50, :)));
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
    global  Trials rasters channel Reps
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
            baseline = squeeze(mean(psths(1, 1:50, :)));
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

function [score_mat, bsl_mat, summary, stat_str] = get_gabor_tuning_stats() %name_pattern is defaultly 'pasu_%02d_ori_%02d_wg_f'
    global  Trials rasters channel Reps
    ori = 0:30:150; sf = [0.5, 1];
    score_mat = nan(2,6,Reps); 
    bsl_mat = nan(2,6,Reps);  % make sure the 3 dimension has larger size than repitition number! 
    cnt_mat = zeros(2,6); 
    id_mat = zeros(2,6); % record the id correspond to i,j
    for i = 1:2
        for j = 1:6
            cur_fn = sprintf('gab_ori_%.1f_%.1f', ori(i), sf(j));
            img_idx = find(contains(Trials.imageName, cur_fn));
            cnt_mat(j, i) = length(img_idx);
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
            score_mat(i, j, 1:length(img_idx)) = scores;
            baseline = squeeze(mean(psths(1, 50:150, :)));
            bsl_mat(i, j, 1:length(img_idx)) = baseline;
            id_mat(i, j) = 6 * (i-1) + j;
        end
    end
    [ori_mat, sf_mat] = meshgrid(1:6, 1:2);
    mean_fr_mat = bsl_mat + score_mat;
    id_vec_nan = reshape(repmat(id_mat, 1, 1, Reps), 1, []);
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
    [p2,tbl2,stats2] = anovan(score_vec_nan, {reshape(repmat(ori_mat, 1,1,Reps),1,[]), ...
                              reshape(repmat(sf_mat, 1,1,Reps),1,[])}, 'model', 'interaction','display' ,'off');
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
            'Orient, F=%.2f(%.3f), SpatFreq, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],CI(1), CI(2), P, stats.F, stats.p, ...
            stats2.F(1),stats2.p(1), stats2.F(2),stats2.p(2),stats2.F(3),stats2.p(3));
end

% Plot the padded image montage 
% Note this is a seperate function in utils
% function frame_img_list = score_frame_image_arr(img_list, score_mat, clim, cmap, LineWidth)
% end