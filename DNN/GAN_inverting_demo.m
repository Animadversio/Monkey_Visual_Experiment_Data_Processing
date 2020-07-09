G = FC6Generator("matlabGANfc6.mat");
%% imresize()
gt_code = randn(1,4096,'single');
tar_img = G.visualize(gt_code);
%% Quick Demo
% out = G.dlforward(dlarray(randn(1,4096)));
code = dlarray(randn(1,4096));
code = gpuArray(code);
dltarget = dlarray(single(tar_img));
%%
tic
% iteration = 1;
averageGrad = [];
averageSqGrad = [];
for iteration = 1:100
[loss, grad] = dlfeval(@FitImgLoss, G, code, dltarget);
[code,averageGrad,averageSqGrad] = adamupdate(code,grad,averageGrad,averageSqGrad,iteration,0.02,0.8,0.98);
if mod(iteration,10)==0
    fprintf("%d %.3f\n",iteration,loss);
end
end
toc
%%

%%
imgfit = gather(G.visualize(code));
%%
figure(1);
subplot(221)
imshow(imgfit)
title("fit img")
subplot(222)
imshow(tar_img)
title("target img")
subplot(223)
imshow(tar_img - imgfit)
title("target - fit")
subplot(224)
imshow(imgfit - tar_img)
title("fit - target")
suptitle(compose("L1 distance %.3f / 255 per pixel corr %.3f",mean(abs(single(imgfit) - single(tar_img)),'all'),corr(single(imgfit(:)),single(tar_img(:)))))
%%
fit_code = single(gather(extractdata(code)));
%%
fprintf("ground truth norm: %.3f, fit norm %.3f, code correlation %.3f\n",norm(gt_code),norm(fit_code),corr(fit_code',gt_code'))
%% Spherical Linear interpolation of the fitted code and the ground truth code! 
slerp_imgs = G.visualize(3 * [cosd(0:10:90);sind(0:10:90)]' * [fit_code;gt_code]);
%%
figure(4);
montage(slerp_imgs,'Size',[1,size(slerp_imgs,4)])
title("SLERP interpolation between fit code(left) and ground truth code(right)")
%%
pasu_path = "S:\Stimuli\2019-Manifold\pasupathy-wg-f-4-ori";
savedir = "E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy";
imgnms = string(ls(pasu_path+"\*.jpg"));
%%
MAXSTEP = 300;
for imgi = 1:length(imgnms)
fprintf("Processing %d pasupathy image\n",imgi);
pasu_img = imread(fullfile(pasu_path, imgnms(imgi)));
rsz_img = imresize(pasu_img, [256, 256], 'bilinear', 'Antialiasing', true);
rsz_img = repmat(rsz_img, 1, 1, 3);
dltarget = dlarray(single(rsz_img));
code = dlarray(randn(1,4096));
code = gpuArray(code);
tic
averageGrad = [];
averageSqGrad = [];
for iteration = 1:MAXSTEP
[loss, grad] = dlfeval(@FitImgLoss, G, code, dltarget);
[code,averageGrad,averageSqGrad] = adamupdate(code,grad,averageGrad,averageSqGrad,iteration,0.02,0.8,0.98);
if mod(iteration,10)==0
    fprintf("%d %.3f\n",iteration,loss);
end
end
toc
img_fit = gather(G.visualize(code));
code_fit = gather(extractdata(code));
loss = extractdata(gather(loss));
save(fullfile(savedir, compose("%03d.mat", imgi)), "code_fit", 'loss')
figure(8);
montage({pasu_img,rsz_img,img_fit},'Size',[1,3])
title(compose("Original,  Resized,  Fit pasuspthy image %d",imgi))
saveas(8, fullfile(savedir, compose("%03d_cmp.jpg", imgi)))
imwrite(img_fit, fullfile(savedir, compose("%03d_fit.jpg", imgi)))
imwrite(img_fit, fullfile(savedir, compose("fit"+imgnms(imgi))))
end
%%
pasu_code = [];
for imgi = 1:length(imgnms)
% fprintf("Processing %d pasupathy image\n",imgi);
% pasu_img = imread(fullfile(pasu_path, imgnms(imgi)));
data = load(fullfile(savedir, compose("%03d.mat", imgi)), "code_fit", 'loss');
pasu_code = [pasu_code; data.code_fit];
end
save(fullfile(savedir,"pasu_fit_code.mat"), 'pasu_code', 'imgnms');

%%
rnd_cc = nan(1, 4096);
vec1 = randn(1,4096)';
vec2 = randn(4096,4096);
rnd_cc = corr(vec1, vec2);
figure;
hist(rnd_cc)
fprintf("mean %.4f std %.4f [0.1, 99.9] range [%.4f, %.4f]\n",mean(rnd_cc), std(rnd_cc), prctile(rnd_cc,0.1), prctile(rnd_cc,99.9))
%%
function [loss, grad] = FitImgLoss(G, code, tar_img)
out = G.dlforward(code);
loss = mean(abs(out(:,:,[3,2,1],:) - tar_img), 'all');
grad = dlgradient(loss, code);
end