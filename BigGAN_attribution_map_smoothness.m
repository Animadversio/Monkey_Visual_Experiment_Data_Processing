load('N:\Stimuli\Project_BigGAN_paired_evolutions\Table_with_SPAM.mat');
Table_Evos_SPAM
%%
SPAM_stats = repmat(struct(), 1, length(Table_Evos_SPAM));
for iTr = 1:length(Table_Evos_SPAM)
    % assert(Table_Evos_SPAM{iTr,3} == Table_Evos_SPAM{iTr+1,3},"Experiment %d and %d Expid are not the same!",iTr,iTr+1)
    SPAM = Table_Evos_SPAM{iTr,6};
    norm_SPAM = SPAM./max(SPAM(:));
    SPAM_stats(iTr).SPAM = SPAM;
    SPAM_stats(iTr).max_R2 = max(SPAM(:));
    SPAM_stats(iTr).mean_R2 = mean(SPAM(:));
    % Smooth
    [total_variation, dirichlet_energy] = quantify_smoothness(norm_SPAM);
    SPAM_stats(iTr).total_variation_norm = total_variation;
    SPAM_stats(iTr).dirichlet_energy_norm = dirichlet_energy;
    [total_variation, dirichlet_energy] = quantify_smoothness(SPAM);
    SPAM_stats(iTr).total_variation_raw = total_variation;
    SPAM_stats(iTr).dirichlet_energy_raw = dirichlet_energy;
    [centrality_weighted_avg, centrality_variance] = quantify_centrality(SPAM, 0.0);
    % Floor at 0, this one seems NOT to be a good measure
    SPAM_stats(iTr).centrality_weighted_avg_zerofloor = centrality_weighted_avg;
    SPAM_stats(iTr).centrality_variance_zerofloor = centrality_variance;
    % Floor at min R2, this one seems to be a good measure
    [centrality_weighted_avg, centrality_variance] = quantify_centrality(SPAM, min(SPAM(:)));
    SPAM_stats(iTr).centrality_weighted_avg_minfloor = centrality_weighted_avg;
    SPAM_stats(iTr).centrality_variance_minfloor = centrality_variance;
end
%%
netfigdir = "N:\Stimuli\Project_BigGAN_paired_evolutions\Figures and Analyses\Evol_spatial_attribution_smoothness";
figdir = "E:\OneDrive - Harvard University\Manuscript_BigGAN\Figures\Evol_spatial_attribution_smoothness";
savefast(fullfile(netfigdir, 'SPAM_smoothness_stats.mat'), 'SPAM_stats')
savefast(fullfile(figdir, 'SPAM_smoothness_stats.mat'), 'SPAM_stats')
%%
figure('position',[1440         913         986         425])
subplot(1,2,1)
scatter([SPAM_stats.max_R2],[SPAM_stats.total_variation_norm])
[r, p] = corr([SPAM_stats.max_R2]', [SPAM_stats.total_variation_norm]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Total Variation (Normalized SPAM)")
xlabel("Max R2")
subplot(1,2,2)
scatter([SPAM_stats.max_R2],[SPAM_stats.total_variation_raw])
[r, p] = corr([SPAM_stats.max_R2]', [SPAM_stats.total_variation_raw]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Total Variation (Raw SPAM)")
xlabel("Max R2")
sgtitle("Spatial Attribution Map Smoothness vs Max R2 of Response Prediction")
saveallform([figdir, netfigdir], 'SPAM_TotalVariation_vs_max_R2', gcf)
%%
figure('position',[1440         913         986         425])
subplot(1,2,1)
scatter([SPAM_stats.mean_R2],[SPAM_stats.total_variation_norm])
[r, p] = corr([SPAM_stats.mean_R2]', [SPAM_stats.total_variation_norm]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Total Variation (Normalized SPAM)")
xlabel("Mean R2")
subplot(1,2,2)
scatter([SPAM_stats.mean_R2],[SPAM_stats.total_variation_raw])
[r, p] = corr([SPAM_stats.mean_R2]', [SPAM_stats.total_variation_raw]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Total Variation (Raw SPAM)")
xlabel("Mean R2")
sgtitle("Spatial Attribution Map Smoothness vs Mean R2 of Response Prediction")
saveallform([figdir, netfigdir], 'SPAM_TotalVariation_vs_mean_R2', gcf)
%%
figure('position',[1440         913         986         425])
subplot(1,2,1)
scatter([SPAM_stats.max_R2],[SPAM_stats.dirichlet_energy_norm])
[r, p] = corr([SPAM_stats.max_R2]', [SPAM_stats.dirichlet_energy_norm]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Dirichlet Energy (Normalized SPAM)")
xlabel("Max R2")
subplot(1,2,2)
scatter([SPAM_stats.max_R2],[SPAM_stats.dirichlet_energy_raw])
[r, p] = corr([SPAM_stats.max_R2]', [SPAM_stats.dirichlet_energy_raw]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Dirichlet Energy (Raw SPAM)")
xlabel("Max R2")
sgtitle("Spatial Attribution Map Smoothness vs Max R2 of Response Prediction")
saveallform([figdir, netfigdir], 'SPAM_DirichletEnergy_vs_max_R2', gcf)
%%
figure('position',[1440         913         986         425])
subplot(1,2,1)
scatter([SPAM_stats.mean_R2],[SPAM_stats.dirichlet_energy_norm])
[r, p] = corr([SPAM_stats.mean_R2]', [SPAM_stats.dirichlet_energy_norm]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Dirichlet Energy (Normalized SPAM)")
xlabel("Mean R2")
subplot(1,2,2)
scatter([SPAM_stats.mean_R2],[SPAM_stats.dirichlet_energy_raw])
[r, p] = corr([SPAM_stats.mean_R2]', [SPAM_stats.dirichlet_energy_raw]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Dirichlet Energy (Raw SPAM)")
xlabel("Mean R2")
sgtitle("Spatial Attribution Map Smoothness vs Mean R2 of Response Prediction")
saveallform([figdir, netfigdir], 'SPAM_DirichletEnergy_vs_mean_R2', gcf)
%%
figure('position',[1440         913         986         425])
subplot(1,2,1)
scatter([SPAM_stats.mean_R2],[SPAM_stats.centrality_weighted_avg_minfloor])
[r, p] = corr([SPAM_stats.mean_R2]', [SPAM_stats.centrality_weighted_avg_minfloor]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Weighted Avg Distance to Center of Mass (floor at min R2)")
xlabel("Mean R2")
subplot(1,2,2)
scatter([SPAM_stats.mean_R2],[SPAM_stats.centrality_weighted_avg_zerofloor])
[r, p] = corr([SPAM_stats.mean_R2]', [SPAM_stats.centrality_weighted_avg_zerofloor]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Weighted Avg Distance to Center of Mass (floor at 0)")
xlabel("Mean R2")
sgtitle("Spatial Attribution Map Centrality Weighted Avg vs Mean R2 of Response Prediction")
saveallform([figdir, netfigdir], 'SPAM_CentralityWeightedAvg_vs_mean_R2', gcf)

%%
figure('position',[1440         913         986         425])
subplot(1,2,1)
scatter([SPAM_stats.max_R2],[SPAM_stats.centrality_weighted_avg_minfloor])
[r, p] = corr([SPAM_stats.max_R2]', [SPAM_stats.centrality_weighted_avg_minfloor]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Weighted Avg Distance to Center of Mass (floor at min R2)")
xlabel("Max R2")
subplot(1,2,2)
scatter([SPAM_stats.mean_R2],[SPAM_stats.centrality_weighted_avg_minfloor])
[r, p] = corr([SPAM_stats.mean_R2]', [SPAM_stats.centrality_weighted_avg_minfloor]', 'Type', 'Pearson');
title(sprintf("Pearson Corr: %.2f, p val=%.1e", r, p))
ylabel("Weighted Avg Distance to Center of Mass (floor at min R2)")
xlabel("Mean R2")
sgtitle("Spatial Attribution Map Centrality Weighted Avg vs Mean R2 of Response Prediction")
saveallform([figdir, netfigdir], 'SPAM_CentralityWeightedAvg_vs_mean_max_R2', gcf)
%%
function [total_variation, dirichlet_energy] = quantify_smoothness(matrix)
    % This function quantifies the smoothness of a 13x13 matrix using total variation and Dirichlet energy.
    % Input:
    %   matrix - a 13x13 matrix
    % Output:
    %   total_variation - the total variation of the matrix
    %   dirichlet_energy - the Dirichlet energy of the matrix

    % Ensure the input is a 13x13 matrix
    % assert(all(size(matrix) == [13, 13]), 'Input matrix must be 13x13.');

    % Compute the gradients
    [dx, dy] = gradient(matrix);
    % Compute the total variation
    total_variation = sum(sum(sqrt(dx.^2 + dy.^2)));
    % Compute the Dirichlet energy
    dirichlet_energy = sum(sum(dx.^2 + dy.^2));
end

function [centrality_weighted_avg, centrality_variance] = quantify_centrality(matrix, floor_thresh)
    % This function quantifies the centrality of a matrix using weighted average distance to the center of mass.
    % Input:
    %   matrix - a matrix of any size
    % Output:
    %   centrality_weighted_avg - the weighted average distance to the center of mass
    %   centrality_variance - the variance of distances from the weighted average distance
    % Compute the center of mass
    if nargin < 2
        floor_thresh = min(matrix(:));
    end
    [rows, cols] = size(matrix);
    [X, Y] = meshgrid(1:cols, 1:rows);
    matrix_thresh = matrix - floor_thresh;
    matrix_thresh(matrix_thresh < 0) = 0;

    total = sum(matrix_thresh(:));
    center_x = sum(X(:) .* matrix_thresh(:)) / total;
    center_y = sum(Y(:) .* matrix_thresh(:)) / total;

    % Compute weighted average distance to the center of mass
    distances = sqrt((X - center_x).^2 + (Y - center_y).^2);
    weighted_avg_distance = sum(distances(:) .* matrix_thresh(:)) / total;

    % Compute variance of distances from the weighted average distance
    distance_variance = sum(((distances(:) - weighted_avg_distance).^2) .* matrix_thresh(:)) / total;

    % Assign the centrality measures to output variables
    centrality_weighted_avg = weighted_avg_distance;
    centrality_variance = distance_variance;

    % (Optional) Display the centrality measures
    fprintf('Weighted Average Distance to Center of Mass: %.4f\n', centrality_weighted_avg);
    fprintf('Variance of Distances: %.4f\n', centrality_variance);
end

