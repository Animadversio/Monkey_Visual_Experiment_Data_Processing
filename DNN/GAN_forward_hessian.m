%% GAN forward Hessian computation demo (translated from the Python code)
ref_code = randn(1, 4096);
tar_img = G.visualize(ref_code);
dltarget = dlarray(single(tar_img));
%%
ref_code = dlarray(ref_code);
%%
[vHv, hvp] = vHvoperator(G, ref_code, dltarget, randn(1,4096));
%%
HVPfun = @(x) HVPoperator_v(G, ref_code, dltarget, x);
%%
HVPfun(randn(4096,1))
%%
tic
[eigvecs,eigvals] = eigs(HVPfun, 4096, 1000, 'largestabs', 'IsFunctionSymmetric', 1);
toc % 77 sec for 100 eigvecs much slower than python backprop
% 530.1 sec for 1000 eigen vectors slower than the python version. 
%%
eigvals = diag(eigvals);
%%
figure;
plot(abs(eigvals))
set(gca, 'YScale', 'log')
%%
figure;
histogram(log10(abs(eigvals)),30)
set(gca, 'YScale', 'log')
%%
eig_ids = [1,2,5,10,15,20,30,40,50,60,70,80,90,100,150,200,250,300,400,500,600,700,800,900,1000];
img_all = [];
for eig_id = eig_ids
imglin = G.visualize(ref_code' + eigvecs(:, eig_id) * linspace(-50, 50, 11));
img_all = cat(4, img_all, imglin);
end
figure; 
montage(img_all,'size',[size(img_all,4) / 11,11])

%%
function hvp = HVPoperator_v(G, ref_code, dltarget, vector)
eps = norm_axis(ref_code,2) * 0.01;
code = ref_code' + vector * [eps , -eps] / norm_axis(vector,1);
[~, grad] = dlfeval(@ImgMetricLoss, G, code, dltarget);
hvp = double(extractdata((grad(:,1) - grad(:,2)) / 2 / eps));
end

function hvp = HVPoperator(G, ref_code, dltarget, vector)
eps = norm_axis(ref_code,2) * 0.01;
code = ref_code + [eps ; -eps] * vector / norm_axis(vector,2);
[loss, grad] = dlfeval(@ImgMetricLoss, G, code, dltarget);
hvp = (grad(1,:) - grad(2, :)) / 2 / eps;
end

function [vHv, hvp] = vHvoperator(G, ref_code, dltarget, vector)
eps = norm_axis(ref_code,2) * 0.01;
code = ref_code + [eps ; -eps] * vector / norm_axis(vector,2);
[loss, grad] = dlfeval(@ImgMetricLoss, G, code, dltarget);
hvp = (grad(1,:) - grad(2, :)) / 2 / eps;
vHv = mean(loss);
end

function [loss, grad] = ImgMetricLoss(G, code, tar_img)
out = G.dlforward(code); % color channel doesn't get converted still 
loss = mean(abs(out(:,:,[3,2,1],:) - tar_img), [1,2,3]);
loss_sum = sum(loss);
grad = dlgradient(loss_sum, code);
end