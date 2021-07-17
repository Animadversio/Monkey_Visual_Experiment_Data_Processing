function [code_fit, img_fit, loss] = GAN_invert_fun(G, target_img, MAXSTEPS, init_code)
if nargin <= 3
    init_code = [];
elseif nargin <= 2
    MAXSTEPS = 200;
end
% image resizing 
if size(target_img,1)==256 && size(target_img,2)==256 % resize to 256 256
rsz_img = target_img;
else
rsz_img = imresize(target_img, [256, 256], 'bilinear', 'Antialiasing', true);
end
if size(rsz_img, 3) == 1 % add RGB channel to gray image 
rsz_img = repmat(rsz_img, 1, 1, 3);
end
if max(rsz_img,[],'all') < 5 % rescale value to 255 range
rsz_img = clip(255.0 * rsz_img,0,255);
end
bsz = size(rsz_img,4); % batch processing is allowed.
% real fitting process
dltarget = dlarray(single(rsz_img));
if isempty(init_code)
code = dlarray(randn(bsz, 4096));
else
code = dlarray(init_code);
assert(size(init_code,1)==bsz, "initial code should match the batch size of target image")
end
code = gpuArray(code);
tic
averageGrad = [];
averageSqGrad = [];
for iteration = 1:MAXSTEPS
[loss, grad] = dlfeval(@FitImgLoss, G, code, dltarget);
[code,averageGrad,averageSqGrad] = adamupdate(code,grad,averageGrad,averageSqGrad,iteration,0.01,0.9,0.98);
if mod(iteration,10)==0
    fprintf("%d  %s\n",iteration,num2str(loss,'%.1f  '));
end
end
toc
img_fit = gather(G.visualize(code));
code_fit = extractdata(gather(code));
loss = CompImgLoss(G, code, dltarget); % compute image specific loss
loss = extractdata(gather(loss));
end 

function [loss, grad] = FitImgLoss(G, code, tar_img)
out = G.dlforward(code); % color channel doesn't get converted still 
loss = mean(abs(out(:,:,[3,2,1],:) - tar_img), [1,2,3]); % L1 loss
losssum = sum(loss);
grad = dlgradient(losssum, code);
end

function [loss] = CompImgLoss(G, code, tar_img)
out = G.dlforward(code); % color channel doesn't get converted still 
loss = squeeze(mean(abs(out(:,:,[3,2,1],:) - tar_img), [1,2,3])); % keep the 4th dimension as loss of each image. 
end