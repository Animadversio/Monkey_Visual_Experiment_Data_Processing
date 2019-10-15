storedStruct = load("D:\\Manifold_Exps.mat");
%%
%load("D:\\Beto64chan-02102019-003_formatted");
Set_Exp_Specs;% load 
global sphere_norm Trials channel rasters ang_step Reps meta
%%
for Expi=1
rasters = storedStruct.rasters{Expi};
Trials = storedStruct.Trials{Expi};
meta = storedStruct.meta{Expi};
unit_name_arr = generate_unit_labels();
% Param_summary = cell(size(rasters,1),3);
pref_chan = pref_chan_arr(Expi);
sphere_norm = norm_arr(Expi);
ang_step = 18;
Reps = 11; % can be any number LARGER than the largest repitition. or there may be problems caused by NAN and 0 filling
savepath = sprintf("C:\\Users\\ponce\\OneDrive\\Desktop\\OneDrive_Binxu\\OneDrive\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
mkdir(savepath);
load(fullfile(savepath, "KentFit_Stats.mat"), 'Param_summary') % Load computed Kent Fit 
load(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary') % Load Stats
%fid = fopen(fullfile(savepath, "Unit_Label.txt"),'w');
%fprintf(fid,'%s\n', unit_name_arr{:});
%fclose(fid);

for channel = 1%:size(rasters,1)
    figure(1);set(1, 'position', [573         294        1377         591]);
    figure(2);set(2, 'position', [573         294        1377         591]);
    figure(3);set(3, 'position', [573         294        1377         591]);
    chan_label_str = sprintf("Exp%d Channel %s", Expi, unit_name_arr{channel});
    
    theta_arr = -90:ang_step:90;
    phi_arr   = -90:ang_step:90;
    % PC1 2
    Parameter = Param_summary{channel, 1};
    summary = Stat_summary{channel, 1};
    [phi_grid, theta_grid] = meshgrid(theta_arr, phi_arr);
    [score_mat, bl_mat] = get_score_mat_from_result('norm_%d_PC2_%d_PC3_%d');
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],summary.t_CI(1), summary.t_CI(2), summary.t_p, summary.anova_F, summary.anova_p, ...
            summary.anova2_F(1),summary.anova2_p(1),summary.anova2_F(2),summary.anova2_p(2),summary.anova2_F(3),summary.anova2_p(3));
    V = coeffvalues(Parameter);
    CI = confint(Parameter);
    param_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
                    V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
    figure(1)
    ax1 = subplot(121);
    sphere_plot(ax1, theta_grid, phi_grid, nanmean(score_mat, 3));
    ylabel("PC2(theta)");zlabel("PC3(phi)");cl = caxis;
    title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
    ax2 = subplot(122);
    Kent_sphere_plot(ax2, coeffvalues(Parameter));
    ylabel("PC2(theta)");zlabel("PC3(phi)");caxis(cl);
    title([chan_label_str, "Kent Fit Tuning map on PC2 3 subspace", param_str])
    % PC49 50
    Parameter = Param_summary{channel, 2};
    summary = Stat_summary{channel, 2};
    [phi_grid, theta_grid] = meshgrid(theta_arr, phi_arr);
    [score_mat, bl_mat] = get_score_mat_from_result('norm_%d_PC49_%d_PC50_%d');
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],summary.t_CI(1), summary.t_CI(2), summary.t_p, summary.anova_F, summary.anova_p, ...
            summary.anova2_F(1),summary.anova2_p(1),summary.anova2_F(2),summary.anova2_p(2),summary.anova2_F(3),summary.anova2_p(3));
    V = coeffvalues(Parameter);
    CI = confint(Parameter);
    param_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
                    V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
    figure(2)
    ax1 = subplot(121);
    sphere_plot(ax1, theta_grid, phi_grid, nanmean(score_mat, 3));
    ylabel("PC49(theta)");zlabel("PC50(phi)");cl = caxis;
    title([chan_label_str, "Tuning map on PC49 50 subspace", stat_str])
    ax2 = subplot(122);
    Kent_sphere_plot(ax2, coeffvalues(Parameter));
    ylabel("PC49(theta)");zlabel("PC50(phi)");caxis(cl);
    title([chan_label_str, "Kent Fit Tuning map on PC49 50subspace", param_str])
    % RND1 2
    Parameter = Param_summary{channel, 3};
    summary = Stat_summary{channel, 3};
    [phi_grid, theta_grid] = meshgrid(theta_arr, phi_arr);
    [score_mat, bl_mat] = get_score_mat_from_result('norm_%d_RND1_%d_RND2_%d');
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],summary.t_CI(1), summary.t_CI(2), summary.t_p, summary.anova_F, summary.anova_p, ...
            summary.anova2_F(1),summary.anova2_p(1),summary.anova2_F(2),summary.anova2_p(2),summary.anova2_F(3),summary.anova2_p(3));
    V = coeffvalues(Parameter);
    CI = confint(Parameter);
    param_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
                    V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
    figure(3)
    ax1 = subplot(121);
    sphere_plot(ax1, theta_grid, phi_grid, nanmean(score_mat, 3));
    ylabel("RND1(theta)");zlabel("RND2(phi)");cl = caxis;
    title([chan_label_str, "Tuning map on RND 1 2 subspace", stat_str])
    ax2 = subplot(122);
    Kent_sphere_plot(ax2, coeffvalues(Parameter));
    ylabel("RND1(theta)");zlabel("RND2(phi)");caxis(cl);
    title([chan_label_str, "Kent Fit Tuning map on RND 1 2 subspace", param_str])
    
    saveas(1, fullfile(savepath, sprintf("chan%02d_PC23_tune_sfit_sph.png", channel)))
    saveas(2, fullfile(savepath, sprintf("chan%02d_PC4950_stune_sfit_sph.png", channel)))
    saveas(3, fullfile(savepath, sprintf("chan%02d_RND_tune_sfit_sph.png", channel)))
end
end
function [score_mat, bl_mat] = get_score_mat_from_result(name_pattern)
    global  Trials rasters channel sphere_norm ang_step Reps
    score_mat = nan(11,11,Reps);
    bl_mat = nan(11,11,Reps);
    cnt_mat = zeros(11,11);
    for i = -5:5
        for j = -5:5
            cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
            img_idx = find(contains(Trials.imageName, cur_fn));
            cnt_mat(i+6, j+6) = length(img_idx);
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
            score_mat(i+6, j+6, 1:length(img_idx)) = scores;
            bl_mat(i+6, j+6, 1:length(img_idx)) = squeeze(mean(psths(1, 1:50, :)));
        end
    end
end
function unit_name_arr = generate_unit_labels()
% Generate the unit labels 17B from the spikeID variable
global meta
Unit_id = meta.spikeID;
unit_name_arr = {}; % name tag for each unit 
for i = 1:length(Unit_id)
    cur_chan = Unit_id(i);
    if sum(Unit_id == cur_chan) == 1
        unit_name_arr{i} = num2str(cur_chan);
    else
        cur_chan = Unit_id(i);
        rel_idx = find(find(Unit_id == cur_chan) == i);
        unit_name_arr{i} = [num2str(cur_chan), char(64+rel_idx)];
    end
end
end

