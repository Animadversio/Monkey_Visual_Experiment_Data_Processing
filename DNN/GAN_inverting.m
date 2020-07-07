G = FC6Generator("matlabGANfc6.mat");
%% imresize()
gt_code = randn(1,4096,'single');
tar_img = G.visualize(gt_code);
%%
% out = G.dlforward(dlarray(randn(1,4096)));
code = dlarray(randn(1,4096));
code = gpuArray(code);
dltarget = dlarray(single(tar_img));
%
tic
% iteration = 1;
averageGrad = [];
averageSqGrad = [];
for iteration = 1:150
[loss, grad] = dlfeval(@FitImgLoss, G, code, dltarget);
[code,averageGrad,averageSqGrad] = adamupdate(code,grad,averageGrad,averageSqGrad,iteration,0.02,0.9,0.99);
if mod(iteration,10)==0
    fprintf("%d %.3f\n",iteration,loss);
end
end
toc
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
figure(3);
montage(slerp_imgs,'Size',[1,size(slerp_imgs,4)])
title("SLERP interpolation between fit code(left) and ground truth code(right)")
%%


function [loss, grad] = FitImgLoss(G, code, tar_img)
out = G.dlforward(code);
loss = mean(abs(out(:,:,[3,2,1],:) - tar_img), 'all');
grad = dlgradient(loss, code);
end