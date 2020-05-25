net = vgg16;
imgN = 121; B=10;
dlimg = randn(224,224,3,imgN);
feat_tsr = activations(net, dlimg, 'conv4_1');
%%
spmask = imgaussfilt(randn(size(feat_tsr,[1,2])), 5,'Padding' ,'symmetric');
figure;imagesc(spmask)
ft_tfm = randn(size(feat_tsr,3),1);
out = einsum(einsum(feat_tsr, spmask, 'ijkl,ij->kl'),ft_tfm,'kl,ki->li');
target = double(out + randn(size(out)) * 1);
%%
fprintf("Fa non-negative ALS fiting\n")
spN = prod(size(feat_tsr,[1,2]));
Fa_init = randn(spN,1);
Fb_init = randn(size(feat_tsr,3),1);
A = reshape(feat_tsr, [spN, size(feat_tsr,[3,4])]);
bcur = 0;
Facur = Fa_init;
Fbcur = Fb_init;
for k = 1:50
Xcur = double(einsum(A, Facur, 'ijk,il->kjl'));
% Fbcur_aug = regress(target, [Xcur, ones(size(Xcur,1),1)]);
% Fbcur_aug = lsqlin([Xcur, ones(size(Xcur,1),1)], target, [], [], [], [], ...
%     [zeros(1, length(Fbcur)), -inf], []);%, optimoptions('lsqlin','Algorithm','interior-point'));
% Fbcur = Fbcur_aug(1:end-1); bcur = Fbcur_aug(end);
[B_b,STATS_b] = lasso([Xcur, ones(size(Xcur,1),1)], target);
ci = round(size(B_b,2)/2); Fbcur = B_b(1:end-1, ci); bcur = B_b(end, ci);
Xcur = double(einsum(A, Fbcur, 'ijk,jl->kil'));
% Facur_aug = regress(target, [Xcur, ones(size(Xcur,1),1)]);
% Facur = Facur_aug(1:end-1); bcur = Facur_aug(end);
[B_a,STATS_a] = lasso([Xcur, ones(size(Xcur,1),1)], target);
ci = round(size(B_a,2)/2); Facur = B_a(1:end-1, ci); bcur = B_a(end, ci);

pred = einsum(einsum(A, Facur, 'ijk,il->ljk'), Fbcur, 'ljk,jl->k') + bcur;
res = target - pred;
fprintf("residue max %.1f norm %.1f\n", max(res), norm(res))
end
%%
figure;
subplot(131)
imagesc(reshape(Facur,size(feat_tsr,[1,2])))
subplot(132)
imagesc(spmask)
subplot(133)
scatter(Fbcur, ft_tfm)
% This is quite degenerate, need regularization