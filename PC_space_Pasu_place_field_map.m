%% PC space tuning summary 
global  Trials rasters channel sphere_norm ang_step Reps
storedStruct = load("D:\\Manifold_Exps.mat");
%%
Reps = 11; % constant for maximum number of repetitions (as long as it's larger than the maximum, it's fine)
Set_Exp_Specs;

for Expi = 24
meta = meta_new{2*(Expi-23)-1};
rasters = rasters_new{2*(Expi-23)-1};
Trials = Trials_new{2*(Expi-23)-1};
% meta = storedStruct.meta{Expi}; % loading data from long term store
% rasters = storedStruct.rasters{Expi};
% Trials = storedStruct.Trials{Expi};
sphere_norm = Pasu_norm_arr(Expi-10); % Load the specific information
pref_chan = Pasu_pref_chan_arr(Expi-10);
ang_step = 18;
% Save basic info
savepath = sprintf("C:\\Users\\ponce\\OneDrive - Washington University in St. Louis\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
mkdir(savepath);
unit_name_arr = generate_unit_labels(meta.spikeID); % Generate readable labels for each channel
load(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary') % Load Stats
%
signif_chan = [];
rsp_mat_arr = {};
signif_F = [];
signif_P = [];
fprintf("=======\nExp %02d Significantly modulated channels:\n", Expi)
for channel = 1:size(rasters,1)
    chan_label_str = sprintf("Channel %s  ", unit_name_arr{channel});
    if Stat_summary{channel,1}.anova_p < 0.01
        signif_chan(end+1) = channel;
        signif_F(end+1) = Stat_summary{channel,1}.anova_F;
        signif_P(end+1) = Stat_summary{channel,1}.anova_p;
        [score_mat, bsl_mat] = get_score_mat('norm_%d_PC2_%d_PC3_%d');
        rsp_mat_arr{end+1} = score_mat;
        fprintf("%s\t", unit_name_arr{channel})
    else
        
    end
end
fprintf('\n')
fprintf(num2str(signif_F,'%.2f\t'))
fprintf("\n%d units in total \n", length(signif_chan))
coln = ceil(length(signif_chan)/4);
figwidth =  1000 / 4 * coln +100;
%%
figure(8);clf;set(8,'position',[42,          42,        figwidth,       1000 ])
t = tiledlayout(4,coln,'TileSpacing','Compact');
Exp_label = sprintf("Exp %d Evolv channel %d", Expi, pref_chan);
for i = 1:numel(signif_chan)
    channel = signif_chan(i);
    chan_label_str = sprintf("Ch %s F%.2f(p=%.1e)", unit_name_arr{channel}, signif_F(i), signif_P(i));
    score_mat = rsp_mat_arr{i};
    maximum = max(nanmean(score_mat, 3), [], 'all');
    ax = nexttile(t);
    imagesc(-90:18:90, -90:18:90, nanmean(score_mat, 3)) % sum(score_mat,3)./cnt_mat
    colormap('parula')
    hold on 
    % Add contour to the plot for visibility
    contour(-90:18:90, -90:18:90, nanmean(score_mat, 3), [2 0.9] * maximum,...
        'LineColor',[1 0 0],'LineWidth',1) % Note single value will cause it to interpret that value as contour line number! 
    contour(-90:18:90, -90:18:90, nanmean(score_mat, 3), [2 0.75] * maximum,...
        'LineColor',[0.6350 0.0780 0.1840],'LineWidth',1) 
    contour(-90:18:90, -90:18:90, nanmean(score_mat, 3), [2 0.5] * maximum,...
        'LineColor','w','LineWidth',1) 
    %contourcmap('hot')
    %ylabel("PC 2 degree");xlabel("PC 3 degree")
    title([chan_label_str])
    axis off;shading flat;axis image;colorbar
end
% suptitle(Exp_label)
title(t,Exp_label)
%%
saveas(8, fullfile(savepath, sprintf("Exp%02d_PC23_placemap.jpg", Expi)))
saveas(8, fullfile(result_dir, sprintf("Exp%02d_PC23_placemap.jpg", Expi)))
end
%%
% % Load and visualize the Pasupathy shapes
% img_dir = 'S:\Stimuli\2019-Manifold\pasupathy-wg-f-4-ori';
% img_list = cell(51, 4);%{};
% for j = 1:4
%     for i = 1:51
%         cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
%         img_idx = find(contains(Trials.imageName, cur_fn));
%         if length(img_idx) == 0 % empty space
%             img_list{i,j} = [];
%         else
%             img_list{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
%         end
%     end
% end
% clear cur_fn i j img_idx
% figure(1);set(1,'position', [ 326         629        2089         254])
% montage(img_list, 'Size', [4, 51]);
% saveas(1, fullfile(savepath, "pasupathy_images.jpg"))
% % Load and visualize the Pasupathy shapes
% img_dir = meta.stimuli; % image storage path online
% evolv_img_list = cell(11, 11);%{};
% for j = 1:11
%     for i = 1:11
%         cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', sphere_norm, ...
%             (i - 6)*ang_step, (j - 6)*ang_step);
%         img_idx = find(contains(Trials.imageName, cur_fn));
%         if length(img_idx) == 0 % empty space
%             evolv_img_list{i,j} = [];
%         else
%             evolv_img_list{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
%         end
%     end
% end
% clear cur_fn i j img_idx
% figure(2);set(2,'position', [ 326         629        2089         254])
% montage(evolv_img_list, 'Size', [11, 11]);
% title(sprintf("Exp%d pref chan%02d PC23 hemisphere", Expi, pref_chan))
% ylabel("PC3 axis");xlabel("PC2 axis");
% saveas(2, fullfile(savepath, sprintf("norm_%d_PC23.jpg", sphere_norm)))
function [score_mat, bsl_mat] = get_score_mat(name_pattern)
    global  Trials rasters channel sphere_norm ang_step Reps
    score_mat = nan(11,11,Reps); 
    bsl_mat = nan(11,11,Reps);  % make sure the 3 dimension has larger size than repitition number! 
    for i =  -5:5
        for j =  -5:5
            cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
            img_idx = find(contains(Trials.imageName, cur_fn));
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 50:150, :))- mean(psths(1, 1:50, :)) );
            baseline = squeeze(mean(psths(1, 50:150, :)));
            score_mat(i+6, j+6, 1:length(img_idx)) = scores;
            bsl_mat(i+6, j+6, 1:length(img_idx)) = baseline;
        end
    end
end

