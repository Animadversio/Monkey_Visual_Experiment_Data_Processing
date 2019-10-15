
%meta,rasters,lfps,Trials

Expi=1;
pref_chan = 6;
channel = 29;
Reps = 10; 
norm = 328;
ang_step = 18;
score_mat = nan(11,11,Reps);
cnt_mat = zeros(11,11);
for i = -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials{Expi}.imageName, cur_fn));
        cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters{Expi}(channel, :, img_idx);
        scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
        score_mat(i+6, j+6, 1:length(img_idx)) = scores;
    end
end
% score_mat, row num go with PC2 value (theta), col num go with PC3 value
% (phi)
%%
theta_arr = -90:ang_step:90;
phi_arr   = -90:ang_step:90;
[theta_grid, phi_grid] = meshgrid(theta_arr, phi_arr);
theta_grid = theta_grid / 180 * pi;
phi_grid = phi_grid / 180 * pi;
theta_grid = repmat(theta_grid, 1,1,Reps);
phi_grid = repmat(phi_grid, 1,1,Reps);
theta_data = theta_grid(~isnan(score_mat(:)));
phi_data = phi_grid(~isnan(score_mat(:)));
score_data = score_mat(~isnan(score_mat(:)));
%%
ft = fittype( @(theta, phi, psi, kappa, beta, A, x, y) KentFunc(theta, phi, psi, kappa, beta, A, x, y), ...
    'independent', {'x', 'y'},'dependent',{'z'});
Parameter = fit([theta_data(:), phi_data(:)], score_data, ft, ...
                'StartPoint', [0, 0, pi/2, 0.1, 0.1, 0.1], ...
                'Lower', [-pi, -pi/2,  0, -Inf,   0,   0], ...
                'Upper', [ pi,  pi/2,  pi,  Inf, Inf, Inf],...)%, ...
                'Robust', 'LAR' )
V = coeffvalues(Parameter);
CI = confint(Parameter);
stat_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
                V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
[theta_grid, phi_grid] = meshgrid(theta_arr, phi_arr);
festim = Parameter(theta_grid(:)/180*pi, phi_grid(:)/180*pi);
%%
figure(1)
imagesc(theta_arr, phi_arr, reshape(festim, size(theta_grid)))
axis image
ylabel("PC3 phi")
xlabel("PC2 theta")
title(stat_str)