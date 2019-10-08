%% Statistics on the 2d tuning map 
pref_chan = 7;
score_mat = zeros(11,11,5);
cnt_mat = zeros(11,11);
for i = -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters(pref_chan, :, img_idx);
        scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
        score_mat(i+6, j+6, 1:length(img_idx)) = scores;
    end
end
% figure(3)
% imagesc(-90:18:90, -90:18:90, sum(score_mat,3)./cnt_mat)
% ylabel("PC 2 degree")
% xlabel("PC 3 degree")
% title("Tuning map on PC2 3 subspace")
% shading flat
% axis image
% colorbar
%
opts.tilted = true;
[PC2_ang, PC3_ang] = meshgrid(18*(-5:5), 18*(-5:5));
results = autoGaussianSurf(PC2_ang, PC3_ang, sum(score_mat,3)./cnt_mat, opts)
%
fprintf("x0=%.1f, y0=%.1f, theta=%.1f, sgmx=%.1f, sgmy=%.1f\n", results.x0, results.y0, results.angle, results.sigmax, results.sigmay)
