%% PC space tuning summary on cortex 
global  Trials rasters channel sphere_norm ang_step Reps
%storedStruct = load("D:\\Manifold_Exps.mat");
%%
Reps = 11; % constant for maximum number of repetitions (as long as it's larger than the maximum, it's fine)
Set_Exp_Specs;

for Expi = 1:10
% meta = meta_new{2*(Expi-23)-1};
% rasters = rasters_new{2*(Expi-23)-1};
% Trials = Trials_new{2*(Expi-23)-1};
meta = storedStruct.meta{Expi}; % loading data from long term store
rasters = storedStruct.rasters{Expi};
Trials = storedStruct.Trials{Expi};
sphere_norm = norm_arr(Expi); % Load the specific information
pref_chan = pref_chan_arr(Expi);
ang_step = 18;
% Save basic info
savepath = sprintf("C:\\Users\\ponce\\OneDrive - Washington University in St. Louis\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
mkdir(savepath);
unit_name_arr = generate_unit_labels(meta.spikeID); % Generate readable labels for each channel
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
load(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary') % Load Stats
sign_filter = true;

name_pattern_list = {'norm_%d_PC2_%d_PC3_%d', 'norm_%d_PC49_%d_PC50_%d', 'norm_%d_RND1_%d_RND2_%d'};
space_list = {'PC23', 'PC4950', 'RND12'};
for space_i = 1:3
% assert(sum(contains(unit_name_arr, "C"))==0, "Too many units in a channel! FIX THE CODE TO HANDLE ME! ")
% Select the significant channels
signif_chan = [];
rsp_mat_arr = {};
signif_F = [];
signif_P = [];
fprintf("=======\nExp %02d Significantly modulated channels in %s:\n", Expi, space_list{space_i})
for channel = 1:size(rasters,1)
    chan_label_str = sprintf("Channel %s  ", unit_name_arr{channel});
    if Stat_summary{channel,1}.anova_p < 0.01
        signif_chan(end+1) = channel;
        signif_F(end+1) = Stat_summary{channel,space_i}.anova_F;
        signif_P(end+1) = Stat_summary{channel,space_i}.anova_p;
        [score_mat, bsl_mat] = get_score_mat( name_pattern_list{space_i} );
        rsp_mat_arr{end+1} = score_mat;
        fprintf("%s\t", unit_name_arr{channel})
    else
        
    end
end
fprintf('\n')
fprintf(num2str(signif_F,'%.2f\t'))
fprintf("\n%d units in total \n", length(signif_chan))
% coln = ceil(length(signif_chan)/4);
% figwidth =  1000 / 4 * coln +100;
%% Plotting only the significant channels 
Exp_label = sprintf("Exp %d Evolv channel %d %s space", Expi, pref_chan, space_list{space_i});
figIT = figure(9);clf;hold on
figV1 = figure(10);clf;hold on
figV4 = figure(11);clf;hold on
%set(figA,'visible','off')
%[tA, ax_arrA] = Cortex_Channel_Tile_Layout("IT", figA); % index to array channel should be the absolute array channel index 
[ax_arrA,tIT,tV1,tV4] = Cortex_Channel_Tile_Layout_All(figIT, figV1, figV4);
[chan_idxA, chan_idxB] = unit_id2_chan_idx(1:64, meta.spikeID); 
if all(chan_idxA==chan_idxB) % each channel has only 1 unit exactly
    plot2unit = false;
else % some channels have more than 1 unit
    %figB = figure(10);clf;hold on;%set(figB,'visible','off')
    figITB = figure(12);clf;hold on
    figV1B = figure(13);clf;hold on
    figV4B = figure(14);clf;hold on
    %[tB, ax_arrB] = Cortex_Channel_Tile_Layout("IT", figB);
    [ax_arrB,tITB,tV1B,tV4B] = Cortex_Channel_Tile_Layout_All(figITB, figV1B, figV4B);
    plot2unit = true;
end
for arr_chan = 1:64 % array channel! not number in the resulting array
    channel = chan_idxA(arr_chan); % the 
    if isnan(channel)
        set(ax_arrA{arr_chan},'Visible','off')
        set(ax_arrB{arr_chan},'Visible','off')
        continue
    end
    chan_label_str = sprintf("Ch %s F%.2f(p=%.1e)", unit_name_arr{channel}, ...
        Stat_summary{channel,1}.anova_F, Stat_summary{channel,1}.anova_p);
    [score_mat, bsl_mat] = get_score_mat( name_pattern_list{space_i} );
    plot_contour_heatmap(nanmean(score_mat, 3), ax_arrA{arr_chan})
    title([chan_label_str])
    
    if plot2unit
        channel = chan_idxB(arr_chan); % the c
        chan_label_str = sprintf("Ch %s F%.2f(p=%.1e)", unit_name_arr{channel}, ...
            Stat_summary{channel,1}.anova_F, Stat_summary{channel,1}.anova_p);
        [score_mat, bsl_mat] = get_score_mat( name_pattern_list{space_i} );
        plot_contour_heatmap(nanmean(score_mat, 3), ax_arrB{arr_chan})
        title([chan_label_str])
    end
    
end
% suptitle(Exp_label)
title(tIT,strcat(Exp_label, " IT array"))
title(tV1,strcat(Exp_label, " V1V2 array"))
title(tV4,strcat(Exp_label, " V4 array"))
if plot2unit
    title(tITB,strcat(Exp_label, " IT array"))
    title(tV1B,strcat(Exp_label, " V1V2 array"))
    title(tV4B,strcat(Exp_label, " V4 array"))
end
%%
saveas(figIT, fullfile(savepath, sprintf("IT_array_Exp%02d_%s_placemapA.jpg", Expi, space_list{space_i})))
saveas(figIT, fullfile(result_dir, sprintf("IT_array_Exp%02d_%s_placemapA.jpg", Expi, space_list{space_i})))
saveas(figV1, fullfile(savepath, sprintf("V1_array_Exp%02d_%s_placemapA.jpg", Expi, space_list{space_i})))
saveas(figV1, fullfile(result_dir, sprintf("V1_array_Exp%02d_%s_placemapA.jpg", Expi, space_list{space_i})))
saveas(figV4, fullfile(savepath, sprintf("V4_array_Exp%02d_%s_placemapA.jpg", Expi, space_list{space_i})))
saveas(figV4, fullfile(result_dir, sprintf("V4_array_Exp%02d_%s_placemapA.jpg", Expi, space_list{space_i})))
if plot2unit
saveas(figITB, fullfile(savepath, sprintf("IT_array_Exp%02d_%s_placemapB.jpg", Expi, space_list{space_i})))
saveas(figITB, fullfile(result_dir, sprintf("IT_array_Exp%02d_%s_placemapB.jpg", Expi, space_list{space_i})))
saveas(figV1B, fullfile(savepath, sprintf("V1_array_Exp%02d_%s_placemapB.jpg", Expi, space_list{space_i})))
saveas(figV1B, fullfile(result_dir, sprintf("V1_array_Exp%02d_%s_placemapB.jpg", Expi, space_list{space_i})))
saveas(figV4B, fullfile(savepath, sprintf("V4_array_Exp%02d_%s_placemapB.jpg", Expi, space_list{space_i})))
saveas(figV4B, fullfile(result_dir, sprintf("V4_array_Exp%02d_%s_placemapB.jpg", Expi, space_list{space_i})))
end
%%
if sign_filter
    for arr_chan = 1:64 % array channel! not number in the resulting array
    channel = chan_idxA(arr_chan);
    if isnan(channel) % The channel is not in the data! 
        continue
    end
    if ~ (Stat_summary{channel,1}.anova_p < 0.01) % Note not < is different from > because there are NaN p values
        set(get(ax_arrA{arr_chan},'Children'),'Visible','off');
        colorbar(ax_arrA{arr_chan},'off')
    end
    if plot2unit
    channel = chan_idxB(arr_chan);
    if ~ (Stat_summary{channel,1}.anova_p < 0.01)
        set(get(ax_arrB{arr_chan},'Children'),'Visible','off');
        colorbar(ax_arrB{arr_chan},'off')
    end
    end
    end
    saveas(figIT, fullfile(savepath, sprintf("IT_array_Exp%02d_%s_placemap_signifA.jpg", Expi, space_list{space_i})))
    saveas(figIT, fullfile(result_dir, sprintf("IT_array_Exp%02d_%s_placemap_signifA.jpg", Expi, space_list{space_i})))
    saveas(figV1, fullfile(savepath, sprintf("V1_array_Exp%02d_%s_placemap_signifA.jpg", Expi, space_list{space_i})))
    saveas(figV1, fullfile(result_dir, sprintf("V1_array_Exp%02d_%s_placemap_signifA.jpg", Expi, space_list{space_i})))
    saveas(figV4, fullfile(savepath, sprintf("V4_array_Exp%02d_%s_placemap_signifA.jpg", Expi, space_list{space_i})))
    saveas(figV4, fullfile(result_dir, sprintf("V4_array_Exp%02d_%s_placemap_signifA.jpg", Expi, space_list{space_i})))
    if plot2unit
    saveas(figITB, fullfile(savepath, sprintf("IT_array_Exp%02d_%s_placemap_signifB.jpg", Expi, space_list{space_i})))
    saveas(figITB, fullfile(result_dir, sprintf("IT_array_Exp%02d_%s_placemap_signifB.jpg", Expi, space_list{space_i})))
    saveas(figV1B, fullfile(savepath, sprintf("V1_array_Exp%02d_%s_placemap_signifB.jpg", Expi, space_list{space_i})))
    saveas(figV1B, fullfile(result_dir, sprintf("V1_array_Exp%02d_%s_placemap_signifB.jpg", Expi, space_list{space_i})))
    saveas(figV4B, fullfile(savepath, sprintf("V4_array_Exp%02d_%s_placemap_signifB.jpg", Expi, space_list{space_i})))
    saveas(figV4B, fullfile(result_dir, sprintf("V4_array_Exp%02d_%s_placemap_signifB.jpg", Expi, space_list{space_i})))
    end
end
end
end
%%

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
function plot_contour_heatmap(score_mat_avg, ax)
    % Given a averaged score matrix, plot a contour heatmap to certain
    % axis! 
    maximum = max(score_mat_avg, [], 'all');
    axes(ax);
    imagesc(-90:18:90, -90:18:90, score_mat_avg) % sum(score_mat,3)./cnt_mat
    colormap('parula')
    hold on 
    % Add contour to the plot for visibility
    contour(-90:18:90, -90:18:90, score_mat_avg, [2 0.9] * maximum,...
        'LineColor',[1 0 0],'LineWidth',1) % Note single value will cause it to interpret that value as contour line number! 
    contour(-90:18:90, -90:18:90, score_mat_avg, [2 0.75] * maximum,...
        'LineColor',[0.6350 0.0780 0.1840],'LineWidth',1) 
    contour(-90:18:90, -90:18:90, score_mat_avg, [2 0.5] * maximum,...
        'LineColor','w','LineWidth',1) 
    %contourcmap('hot')
    %ylabel("PC 2 degree");xlabel("PC 3 degree")
    axis off;shading flat;axis image;colorbar
end
