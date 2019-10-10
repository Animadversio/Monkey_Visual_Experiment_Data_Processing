%% Gaussian Fitting Statistics on the 2d tuning map 
pref_chan = 6;
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
avg_score = sum(score_mat,3)./cnt_mat;
results = autoGaussianSurf(PC2_ang, PC3_ang, avg_score, opts)
stat_str = sprintf("x0=%.1f, y0=%.1f, theta=%.1f, sgmx=%.1f, sgmy=%.1f, R^2=%.3f", results.x0, results.y0, results.angle, results.sigmax, results.sigmay, results.r2);
%
figure(8)
set(8, 'position', [541         358        1622         620])
subplot(1,2,1);imagesc(PC2_ang(:),PC3_ang(:),avg_score, [min(avg_score,[],'all'), max(avg_score,[],'all')]);axis image
title("Tuning map on PC2 3 subspace")

colorbar();
subplot(1,2,2);imagesc(PC2_ang(:),PC3_ang(:),results.G, [min(avg_score,[],'all'), max(avg_score,[],'all')]);axis image
title(["Fit Tuning map",stat_str])
colorbar();
%% %%%%%%%%%%%
%load(fullfile("\\storage1.ris.wustl.edu\crponce\Active\Data-Ephys-MAT","Beto64chan-02102019-003_formatted.mat"))
load(fullfile("\\storage1.ris.wustl.edu\crponce\Active\Data-Ephys-MAT","Beto64chan-03102019-002_formatted.mat"))
% meta,rasters,lfps,Trials
%% Loading code from "D:\Poncelab_Github\office-main\Project_Selectivity_Beto_loadRaw.m"
img_names = unique(Trials.imageName);
norm = 326; % 269 Day3 % 326 Day2 % 328 Day 1
ang_step = 18;
pref_chan = 5;
savepath = "C:\Users\ponce\OneDrive\Desktop\OneDrive_Binxu\OneDrive\PC_space_tuning\Exp2_chan06";%Exp3_chan05%Exp1_chan29
mkdir(savepath);
FitStats = {};
for pref_chan = 66
figure(6);clf
set(gcf, 'position', [131         145        1567         809]);
%%
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
avg_score = sum(score_mat,3)./cnt_mat;
opts.tilted = true;
[PC2_ang, PC3_ang] = meshgrid(18*(-5:5), 18*(-5:5));
figure,surf(PC2_ang,PC3_ang,avg_score)
results = autoGaussianSurf(PC2_ang, PC3_ang, avg_score, opts)
stat_str = sprintf("x0=%.1f, y0=%.1f, theta=%.1f,\n sgmx=%.1f, sgmy=%.1f, R^2=%.3f", results.x0, results.y0, results.angle, results.sigmax, results.sigmay, results.r2);
FitStats{pref_chan, 1} = results;

figure(3);clf
set(3, 'position', [541         358        1622         620])
subplot(1,2,1);
imagesc(PC2_ang(:),PC3_ang(:),avg_score, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
axis image
title("Tuning map on PC2 3 subspace")
ylabel("PC 2 degree")
xlabel("PC 3 degree")
shading flat
colorbar;
subplot(1,2,2);
imagesc(PC2_ang(:),PC3_ang(:),results.G, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
axis image
title(["Fit Tuning map",stat_str])
ylabel("PC 2 degree")
xlabel("PC 3 degree")
colorbar;

figure(6)
subplot(231)
imagesc(PC2_ang(:),PC3_ang(:),avg_score, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
imagesc(-90:18:90, -90:18:90, sum(score_mat,3)./cnt_mat)
ylabel("PC 2 degree")
xlabel("PC 3 degree")
title("Tuning map on PC2 3 subspace")
shading flat
axis image
colorbar
subplot(234)
imagesc(PC2_ang(:),PC3_ang(:),results.G, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
ylabel("PC 2 degree")
xlabel("PC 3 degree")
title(["Fit Tuning map",stat_str])
shading flat
axis image
colorbar
%%
score_mat = zeros(11,11,5);
cnt_mat = zeros(11,11);
for i = -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_PC49_%d_PC50_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters(pref_chan, :, img_idx);
        scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
        score_mat(i+6, j+6, 1:length(img_idx)) = scores;
    end
end
avg_score = sum(score_mat,3)./cnt_mat;
opts.tilted = true;
[PC2_ang, PC3_ang] = meshgrid(18*(-5:5), 18*(-5:5));
results = autoGaussianSurf(PC2_ang, PC3_ang, avg_score, opts)
stat_str = sprintf("x0=%.1f, y0=%.1f, theta=%.1f,\n sgmx=%.1f, sgmy=%.1f, R^2=%.3f", results.x0, results.y0, results.angle, results.sigmax, results.sigmay, results.r2);
FitStats{pref_chan, 2} = results;

figure(4);clf
set(4, 'position', [541         358        1622         620])
subplot(1,2,1);
imagesc(PC2_ang(:),PC3_ang(:),avg_score, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
axis image
ylabel("PC 49 degree")
xlabel("PC 50 degree")
title("Tuning map on PC49 50 subspace")
shading flat
colorbar;
subplot(1,2,2);
imagesc(PC2_ang(:),PC3_ang(:),results.G, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
axis image
title(["Fit Tuning map",stat_str])
ylabel("PC 2 degree")
xlabel("PC 3 degree")
colorbar;

figure(6)
subplot(232)
imagesc(PC2_ang(:),PC3_ang(:),avg_score, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
imagesc(-90:18:90, -90:18:90, sum(score_mat,3)./cnt_mat)
ylabel("PC 49 degree")
xlabel("PC 50 degree")
title("Tuning map on PC49 50 subspace")
shading flat
axis image
colorbar
subplot(235)
imagesc(PC2_ang(:),PC3_ang(:),results.G, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
ylabel("PC 49 degree")
xlabel("PC 50 degree")
title(["Fit Tuning map",stat_str])
shading flat
axis image
colorbar

%%
score_mat = zeros(11,11,5);
cnt_mat = zeros(11,11);
for i = -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_RND1_%d_RND2_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters(pref_chan, :, img_idx);
        scores = squeeze(mean(psths(1, 51:200, :)) - mean(psths(1, 1:50, :)) );
        score_mat(i+6, j+6, 1:length(img_idx)) = scores;
    end
end
avg_score = sum(score_mat,3)./cnt_mat;
opts.tilted = true;
[PC2_ang, PC3_ang] = meshgrid(18*(-5:5), 18*(-5:5));
results = autoGaussianSurf(PC2_ang, PC3_ang, avg_score, opts)
stat_str = sprintf("x0=%.1f, y0=%.1f, theta=%.1f\n sgmx=%.1f, sgmy=%.1f, R^2=%.3f", results.x0, results.y0, results.angle, results.sigmax, results.sigmay, results.r2);
FitStats{pref_chan, 3} = results;

figure(5);clf
set(5, 'position', [541         358        1622         620])
subplot(1,2,1);
imagesc(PC2_ang(:),PC3_ang(:),avg_score, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
axis image
ylabel("Rand vector 1 degree")
xlabel("Rand vector 2 degree")
title("Tuning map on Random Vector (outside first 50 PCs) subspace")
shading flat
colorbar;
subplot(1,2,2);
imagesc(PC2_ang(:),PC3_ang(:),results.G, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
axis image
title(["Fit Tuning map",stat_str])
ylabel("PC 2 degree")
xlabel("PC 3 degree")
colorbar;

figure(6)
subplot(233)
imagesc(PC2_ang(:),PC3_ang(:),avg_score, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
imagesc(-90:18:90, -90:18:90, sum(score_mat,3)./cnt_mat)
ylabel("Rand vector 1 degree")
xlabel("Rand vector 2 degree")
title("Tuning map on Random Vector (outside first 50 PCs) subspace")
shading flat
axis image
colorbar
subplot(236)
imagesc(PC2_ang(:),PC3_ang(:),results.G, [min(avg_score,[],'all'), max(avg_score,[],'all')]);
ylabel("Rand vector 1 degree")
xlabel("Rand vector 2 degree")
title(["Fit Tuning map",stat_str])
shading flat
axis image
colorbar
%%

saveas(6, fullfile(savepath, sprintf("chan%02d_PC_tune_cmp_fit.png", pref_chan)))
saveas(3, fullfile(savepath, sprintf("chan%02d_PC23_tune_fit.png", pref_chan)))
saveas(4, fullfile(savepath, sprintf("chan%02d_PC4950_tune_fit.png", pref_chan)))
saveas(5, fullfile(savepath, sprintf("chan%02d_RND_tune_fit.png", pref_chan)))
end
save(fullfile(savepath, "FitStats.mat"), 'FitStats')
