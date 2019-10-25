% Batch processing code  calculating the statistics (t, 1way ANOVA, 2way ANOVA) 
% and generate annotated figures for each and every experiment 
%% 
storedStruct = load("D:\\Manifold_Exps.mat");
%%
%load("D:\\Beto64chan-02102019-003_formatted");
Set_Exp_Specs;
background_plot = 0;
global sphere_norm Trials channel rasters ang_step Reps meta
%%
background_plot = 0;
for Expi=10
rasters = storedStruct.rasters{Expi};
Trials = storedStruct.Trials{Expi};
meta = storedStruct.meta{Expi};
unit_name_arr = generate_unit_labels();
Param_summary = cell(size(rasters,1),3);
pref_chan = pref_chan_arr(Expi);
sphere_norm = norm_arr(Expi);
ang_step = 18;
Reps = 11; % can be any number LARGER than the largest repitition. or there may be problems caused by NAN and 0 filling
% savepath = "C:\Users\ponce\OneDrive\Desktop\OneDrive_Binxu\OneDrive\PC_space_tuning\Exp1_chan29";
savepath = sprintf("C:\\Users\\ponce\\OneDrive\\Desktop\\OneDrive_Binxu\\OneDrive\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
mkdir(savepath);
%fid = fopen(fullfile(savepath, "Unit_Label.txt"),'w');
%fprintf(fid,'%s\n', unit_name_arr{:});
%fclose(fid);
load(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary')
for channel = 1:size(rasters,1)
    figure(1);clf;set(1, 'position', [555         120        1507         671]);
    figure(2);clf;set(2, 'position', [555         120        1507         671]);
    figure(3);clf;set(3, 'position', [555         120        1507         671]);
    figure(4);clf;set(4, 'position', [131          42        1596         954]);
    figure(5);set(5, 'position', [573         294        1377         591]);
    figure(6);set(6, 'position', [573         294        1377         591]);
    figure(7);set(7, 'position', [573         294        1377         591]);
    chan_label_str = sprintf("Exp%d Channel %s", Expi, unit_name_arr{channel});
    theta_arr = -90:ang_step:90;
    phi_arr   = -90:ang_step:90;
    [phi_grid, theta_grid] = meshgrid(theta_arr, phi_arr);
    %% PC12
    [score_mat, Parameter, param_str, score_estim] = get_fitting_from_result('norm_%d_PC2_%d_PC3_%d');
    Param_summary{channel, 1} = Parameter;
    summary = Stat_summary{channel, 1};
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],summary.t_CI(1), summary.t_CI(2), summary.t_p, summary.anova_F, summary.anova_p, ...
            summary.anova2_F(1),summary.anova2_p(1),summary.anova2_F(2),summary.anova2_p(2),summary.anova2_F(3),summary.anova2_p(3));
    disp(Parameter)
    % visualize
    avg_score = nanmean(score_mat,3); scorelim = [min(avg_score,[],'all'), max(avg_score,[],'all')+0.1]; 
    % Use this to match color bar
    figure(1);if background_plot, set(1, 'Visible', 'off'); end
    subplot(121)
    imagesc(-90:18:90, -90:18:90, avg_score, scorelim)
    ylabel("PC 2 degree")
    xlabel("PC 3 degree")
    title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
    shading flat;axis image;colorbar
    subplot(122)
    imagesc(-90:18:90, -90:18:90, score_estim, scorelim)
    ylabel("PC 2 degree")
    xlabel("PC 3 degree")
    title([chan_label_str, "Tuning map on PC2 3 subspace", param_str])
    shading flat;axis image;colorbar
    
    % visualize fitting on sphere
    figure(5);if background_plot, set(5, 'Visible', 'off'); end
    ax1 = subplot(121);
    sphere_plot(ax1, theta_grid, phi_grid, nanmean(score_mat, 3));
    ylabel("PC2(theta)");zlabel("PC3(phi)");cl = caxis;
    title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
    ax2 = subplot(122);
    Kent_sphere_plot(ax2, coeffvalues(Parameter));
    ylabel("PC2(theta)");zlabel("PC3(phi)");caxis(cl);
    title([chan_label_str, "Kent Fit Tuning map on PC2 3 subspace", param_str])
    
    figure(4);if background_plot, set(4, 'Visible', 'off'); end
    subplot(231)
    imagesc(-90:18:90, -90:18:90, avg_score, scorelim)
    ylabel("PC 2 degree")
    xlabel("PC 3 degree")
    title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
    shading flat;axis image;colorbar
    subplot(234)
    imagesc(-90:18:90, -90:18:90, score_estim, scorelim)
    ylabel("PC 2 degree")
    xlabel("PC 3 degree")
    title([chan_label_str, "Tuning map on PC2 3 subspace", param_str])
    shading flat;axis image;colorbar
    %% PC4950
    [score_mat, Parameter, param_str, score_estim] = get_fitting_from_result('norm_%d_PC49_%d_PC50_%d');
    Param_summary{channel, 2} = Parameter;
    disp(Parameter)
    avg_score = nanmean(score_mat,3); scorelim = [min(avg_score,[],'all'), max(avg_score,[],'all')+0.1]; 
    summary = Stat_summary{channel, 2};
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],summary.t_CI(1), summary.t_CI(2), summary.t_p, summary.anova_F, summary.anova_p, ...
            summary.anova2_F(1),summary.anova2_p(1),summary.anova2_F(2),summary.anova2_p(2),summary.anova2_F(3),summary.anova2_p(3));
    
    % Use this to match color bar
    figure(2);if background_plot, set(2, 'Visible', 'off'); end
    subplot(121)
    imagesc(-90:18:90, -90:18:90, avg_score, scorelim)
    ylabel("PC 49 degree")
    xlabel("PC 50 degree")
    title([chan_label_str, "Tuning map on PC49 50 subspace", stat_str])
    shading flat;axis image;colorbar
    subplot(122)
    imagesc(-90:18:90, -90:18:90, score_estim, scorelim)
    ylabel("PC 49 degree")
    xlabel("PC 50 degree")
    title([chan_label_str, "Tuning map on PC2 3 subspace", param_str])
    shading flat;axis image;colorbar
    
    figure(6);if background_plot, set(6, 'Visible', 'off'); end
    ax1 = subplot(121);
    sphere_plot(ax1, theta_grid, phi_grid, nanmean(score_mat, 3));
    ylabel("PC49(theta)");zlabel("PC50(phi)");cl = caxis;
    title([chan_label_str, "Tuning map on PC49 50 subspace", stat_str])
    ax2 = subplot(122);
    Kent_sphere_plot(ax2, coeffvalues(Parameter));
    ylabel("PC49(theta)");zlabel("PC50(phi)");caxis(cl);
    title([chan_label_str, "Kent Fit Tuning map on PC49 50subspace", param_str])
    
    figure(4);if background_plot, set(6, 'Visible', 'off'); end
    subplot(232)
    imagesc(-90:18:90, -90:18:90, avg_score, scorelim)
    ylabel("PC 49 degree")
    xlabel("PC 50 degree")
    title([chan_label_str, "Tuning map on PC49 50 subspace", stat_str])
    shading flat;axis image;colorbar
    subplot(235)
    imagesc(-90:18:90, -90:18:90, score_estim, scorelim)
    ylabel("PC 49 degree")
    xlabel("PC 50 degree")
    title([chan_label_str, "Tuning map on PC2 3 subspace", param_str])
    shading flat;axis image;colorbar
    %% RND12
    [score_mat, Parameter, param_str, score_estim] = get_fitting_from_result('norm_%d_RND1_%d_RND2_%d');
    Param_summary{channel, 3} = Parameter;
    disp(Parameter)
    avg_score = nanmean(score_mat,3); scorelim = [min(avg_score,[],'all'), max(avg_score,[],'all')+0.1]; 
    summary = Stat_summary{channel, 3};
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],summary.t_CI(1), summary.t_CI(2), summary.t_p, summary.anova_F, summary.anova_p, ...
            summary.anova2_F(1),summary.anova2_p(1),summary.anova2_F(2),summary.anova2_p(2),summary.anova2_F(3),summary.anova2_p(3));
    
    % Use this to match color bar
    figure(3);if background_plot, set(3, 'Visible', 'off'); end
    subplot(121)
    imagesc(-90:18:90, -90:18:90, avg_score, scorelim)
    ylabel("Rand vector 1 degree")
    xlabel("Rand vector 2 degree")
    title([chan_label_str, "Tuning map on Random Vector (outside first 50 PCs) subspace", stat_str])
    shading flat;axis image;colorbar
    subplot(122)
    imagesc(-90:18:90, -90:18:90, score_estim, scorelim)
    ylabel("Rand vector 1 degree")
    xlabel("Rand vector 2 degree")
    title([chan_label_str, "Tuning map on Random Vector (outside first 50 PCs) subspace", param_str])
    shading flat;axis image;colorbar
    
    figure(7);if background_plot, set(7, 'Visible', 'off'); end
    ax1 = subplot(121);
    sphere_plot(ax1, theta_grid, phi_grid, nanmean(score_mat, 3));
    ylabel("RND1(theta)");zlabel("RND2(phi)");cl = caxis;
    title([chan_label_str, "Tuning map on RND 1 2 subspace", stat_str])
    ax2 = subplot(122);
    Kent_sphere_plot(ax2, coeffvalues(Parameter));
    ylabel("RND1(theta)");zlabel("RND2(phi)");caxis(cl);
    title([chan_label_str, "Kent Fit Tuning map on RND 1 2 subspace", param_str])
    
    figure(4);if background_plot, set(4, 'Visible', 'off'); end
    subplot(233)
    imagesc(-90:18:90, -90:18:90, avg_score, scorelim)
    ylabel("Rand vector 1 degree")
    xlabel("Rand vector 2 degree")
    title([chan_label_str, "Tuning map on Random Vector (outside first 50 PCs) subspace", stat_str])
    shading flat;axis image;colorbar
    subplot(236)
    imagesc(-90:18:90, -90:18:90, score_estim, scorelim)
    ylabel("Rand vector 1 degree")
    xlabel("Rand vector 2 degree")
    title([chan_label_str, "Tuning map on Random Vector (outside first 50 PCs) subspace", param_str])
    shading flat;axis image;colorbar
    %%
    saveas(4, fullfile(savepath, sprintf("chan%02d_PC_tune_cmp_sfit.png", channel)))
    saveas(1, fullfile(savepath, sprintf("chan%02d_PC23_tune_sfit.png", channel)))
    saveas(2, fullfile(savepath, sprintf("chan%02d_PC4950_stune_sfit.png", channel)))
    saveas(3, fullfile(savepath, sprintf("chan%02d_RND_tune_sfit.png", channel)))
    saveas(5, fullfile(savepath, sprintf("chan%02d_PC23_tune_sfit_sph.png", channel)))
    saveas(6, fullfile(savepath, sprintf("chan%02d_PC4950_stune_sfit_sph.png", channel)))
    saveas(7, fullfile(savepath, sprintf("chan%02d_RND_tune_sfit_sph.png", channel)))
end
save(fullfile(savepath, "KentFit_Stats.mat"), 'Param_summary')
end

%%
function [score_mat, Parameter, stat_str, festim] = get_fitting_from_result(name_pattern)
    global  Trials rasters channel sphere_norm ang_step Reps
    score_mat = nan(11,11,Reps);
    cnt_mat = zeros(11,11);
    for i = -5:5
        for j = -5:5
            cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
            img_idx = find(contains(Trials.imageName, cur_fn));
            cnt_mat(i+6, j+6) = length(img_idx);
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
            score_mat(i+6, j+6, 1:length(img_idx)) = scores;
        end
    end
    theta_arr = -90:ang_step:90;
    phi_arr   = -90:ang_step:90;
    [phi_grid, theta_grid] = meshgrid(theta_arr, phi_arr);
    theta_grid = theta_grid / 180 * pi;
    phi_grid = phi_grid / 180 * pi;
    theta_grid = repmat(theta_grid, 1,1,Reps);
    phi_grid = repmat(phi_grid, 1,1,Reps);
    theta_data = theta_grid(~isnan(score_mat(:)));
    phi_data = phi_grid(~isnan(score_mat(:)));
    score_data = score_mat(~isnan(score_mat(:)));
    %
    ft = fittype( @(theta, phi, psi, kappa, beta, A, x, y) KentFunc(theta, phi, psi, kappa, beta, A, x, y), ...
        'independent', {'x', 'y'},'dependent',{'z'});
    Parameter = fit([theta_data(:), phi_data(:)], score_data, ft, ...
                    'StartPoint', [0, 0, pi/2, 0.1, 0.1, 0.1], ...
                    'Lower', [-pi, -pi/2,  0, -Inf,   0,   0], ...
                    'Upper', [ pi,  pi/2,  pi,  Inf, Inf, Inf],...)%, ...
                    'Robust', 'LAR' );
    V = coeffvalues(Parameter);
    CI = confint(Parameter);
    stat_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
                    V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
    [phi_grid, theta_grid] = meshgrid(theta_arr, phi_arr);
    festim = Parameter(theta_grid(:)/180*pi, phi_grid(:)/180*pi);
    festim = reshape(festim, size(theta_grid));
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

