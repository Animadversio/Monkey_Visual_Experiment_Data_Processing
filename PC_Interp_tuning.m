% General Code for analyze tuning in the interpolation space between 2
% evolutions
global sphere_norm Trials rasters ang_step Reps meta
%%
Expi = 1;
sphere_norm = 294; % 269 Day3 % 326 Daye % 328 Day 1
ang_step = 18;
pref_chan = 26;
Reps = 11;
rasters = rasters_new{1};
meta = meta_new{1};
Trials = Trials_new{1};
%savepath = "C:\Users\ponce\OneDrive\Desktop\OneDrive_Binxu\OneDrive\PC_space_tuning\Interp_Exp1_chan26";
savepath = sprintf("C:\\Users\\ponce\\OneDrive\\Desktop\\OneDrive_Binxu\\OneDrive\\PC_space_tuning\\Interp_Exp%d_chan%02d", Expi, pref_chan);
mkdir(savepath);
%%
Stat_summary = {};%cell(size(rasters,1), 2);
unit_name_arr = generate_unit_labels(meta.spikeID,savepath);
global channel
for channel = 1:size(rasters,1)
    
set(1, 'position', [ 971         335        1025         643])
set(2, 'position', [ 971         335        1025         643])
set(3, 'position', [1283          42         776         954])
chan_label_str = sprintf("Exp%d Channel %s", Expi, unit_name_arr{channel});

[score_mat, bl_mat, cnt_mat] = get_score_from_pattern('norm_%d_interp_%d_PC3_%d');
theta_arr = -90:ang_step:180;
phi_arr = -90:ang_step:90;
[summary, stat_str] = calc_tuning_stats(score_mat, bl_mat, theta_arr, phi_arr);
Stat_summary{channel, 1} = summary;
figure(1)
imagesc(theta_arr,phi_arr,nanmean(score_mat,3)',...
    'AlphaData',double(cnt_mat' ~= 0))
axis image;colorbar;box off
xlabel("Interpolation (thata)");ylabel("Deviation in PC3*(phi)");
xticks(theta_arr);yticks(phi_arr)
title([chan_label_str, "Tuning map on Interpolation PC3 subspace", stat_str]);%, param_str])

figure(3)
subplot(211)
imagesc(theta_arr,phi_arr,nanmean(score_mat,3)',...
    'AlphaData',double(cnt_mat' ~= 0))
axis image;colorbar;box off
xlabel("Interpolation (thata)");ylabel("Deviation in PC3*(phi)");
xticks(theta_arr);yticks(phi_arr)
title([chan_label_str, "Tuning map on Interpolation PC3 subspace", stat_str]);%, param_str])

[score_mat, bl_mat, cnt_mat] = get_score_from_pattern('norm_%d_interp_%d_PC4_%d');
[summary, stat_str] = calc_tuning_stats(score_mat, bl_mat, theta_arr, phi_arr);
Stat_summary{channel, 2} = summary;
theta_arr = -90:ang_step:180;
phi_arr = -90:ang_step:90;
figure(2)
imagesc(theta_arr,phi_arr,nanmean(score_mat,3)',...
    'AlphaData',double(cnt_mat' ~= 0))
axis image;colorbar;box off
xlabel("Interpolation (thata)");ylabel("Deviation in PC4*(phi)");
xticks(theta_arr);yticks(phi_arr)
title([chan_label_str, "Tuning map on Interpolation PC4 subspace", stat_str]);%, param_str])

figure(3)
subplot(212)
imagesc(theta_arr,phi_arr,nanmean(score_mat,3)',...
    'AlphaData',double(cnt_mat' ~= 0))
axis image;colorbar;box off
xlabel("Interpolation (thata)");ylabel("Deviation in PC4*(phi)");
xticks(theta_arr);yticks(phi_arr)
title([chan_label_str, "Tuning map on Interpolation PC4 subspace", stat_str]);%, param_str])

saveas(1, fullfile(savepath, sprintf("chan%02d_Interp_PC3_tune.png", channel)))
saveas(2, fullfile(savepath, sprintf("chan%02d_Interp_PC4_tune.png", channel)))
saveas(3, fullfile(savepath, sprintf("chan%02d_Interp_tune_cmp.png", channel)))
end
save(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary')
Full_Data_table = cell2table([cellfun(@(c) {c.t_CI}, Stat_summary),...
                              cellfun(@(c) {c.t_p}, Stat_summary), ...
                              cellfun(@(c) {c.anova_F}, Stat_summary), ...
                              cellfun(@(c) {c.anova_p}, Stat_summary)], ...
                              'VariableNames',{'PC3_CI' 'PC4_CI' 'PC3_p' 'PC4_p' 'PC3_ANOVA_F' 'PC4_ANOVA_F' 'PC3_ANOVA_p' 'PC4_ANOVA_p'}, ...
                              'RowNames', unit_name_arr); % Can be appended to more data! 
writetable(Full_Data_table, fullfile(savepath,'StatsTable.csv'),'QuoteStrings',true,'WriteRowNames',true);

figure(8);clf;set(8, 'position', [977    42   392   954]);
imagesc(cell2mat(cellfun(@(c) {c.anova_F}, Stat_summary)))
title(sprintf("Exp%d F Stats, ANOVA1", Expi))
colorbar;xticks(1:2)
yticks(0.5:1:length(unit_name_arr)-0.5);yticklabels(unit_name_arr)
saveas(8, fullfile(savepath, sprintf("Exp%d_ANOVA_stat.bmp", Expi)))


%%
function [score_mat, bl_mat, cnt_mat] = get_score_from_pattern(name_pattern)
    global Reps rasters Trials sphere_norm ang_step channel
    score_mat = nan(16,11,Reps);
    bl_mat = nan(16,11,Reps);
    cnt_mat = zeros(16,11);
    for i = -5:10
        for j = -5:5
            cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
            img_idx = find(contains(Trials.imageName, cur_fn));
            cnt_mat(i+6, j+6) = length(img_idx);
            if isempty(img_idx)
                continue
            end
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
            score_mat(i+6, j+6, 1:length(img_idx)) = scores;
            bl_mat(i+6, j+6, 1:length(img_idx)) = squeeze(mean(psths(1, 1:50, :)) );
        end
    end
end
