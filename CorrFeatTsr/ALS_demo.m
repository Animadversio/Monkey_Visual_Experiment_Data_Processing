%% A simple demo for alternative least square algorithm (no constraint) 
% Observation is more then any regressor
% If so the ALS fitting and convergence is fine! 
%% Set up synthetic data: Factorized linear model + noise
A = randn(10,15,100);
%A = randn(10,15,20);
Fa = randn(10,1);
Fb = randn(15,1);
out = einsum(einsum(A, Fa, 'ijk,il->jkl'), Fb, 'jk,ji->ki');
target = out + 0.5 * randn(size(out));
%% Main Algorithm for ALS
fprintf("Un constraint ALS fiting\n")
Fa_init = randn(10,1);
Fb_init = randn(15,1);
bcur = 0;
Facur = Fa_init;
Fbcur = Fb_init;
for k = 1:15
Xcur = einsum(A, Facur, 'ijk,il->kjl');
Fbcur_aug = regress(target, [Xcur, ones(size(Xcur,1),1)]);
Fbcur = Fbcur_aug(1:end-1); bcur = Fbcur_aug(end);
Xcur = einsum(A, Fbcur, 'ijk,jl->kil');
Facur_aug = regress(target, [Xcur, ones(size(Xcur,1),1)]);
Facur = Facur_aug(1:end-1); bcur = Facur_aug(end);
pred = einsum(einsum(A, Facur, 'ijk,il->ljk'), Fbcur, 'ljk,jl->k') + bcur;
res = target - pred;
fprintf("residue max %.1f norm %.1f\n", max(res), norm(res))
end
%% Main Algorithm for ALS (non negative for one side.)
fprintf("Fa non-negative ALS fiting\n")
Fa_init = randn(10,1);
Fb_init = randn(15,1);
bcur = 0;
Facur = Fa_init;
Fbcur = Fb_init;
for k = 1:30
Xcur = einsum(A, Facur, 'ijk,il->kjl');
% Fbcur_aug = regress(target, [Xcur, ones(size(Xcur,1),1)]);
Fbcur_aug = lsqlin([Xcur, ones(size(Xcur,1),1)], target, [], [], [], [], ...
    [zeros(1, length(Fbcur)), -inf], []);%, optimoptions('lsqlin','Algorithm','interior-point'));
Fbcur = Fbcur_aug(1:end-1); bcur = Fbcur_aug(end);
Xcur = einsum(A, Fbcur, 'ijk,jl->kil');
Facur_aug = regress(target, [Xcur, ones(size(Xcur,1),1)]);
Facur = Facur_aug(1:end-1); bcur = Facur_aug(end);
pred = einsum(einsum(A, Facur, 'ijk,il->ljk'), Fbcur, 'ljk,jl->k') + bcur;
res = target - pred;
fprintf("residue max %.1f norm %.1f\n", max(res), norm(res))
end
%% Visualize fitting
figure;set(gcf,'Position',[6         265        1547         426])
subplot(131)
scatter(pred, target);hold on
scatter(pred, out)
legend(["noisy Y","original Y"])
title("regressed output~noisy and original output")
subplot(132)
scatter(Facur, Fa)
title("Factor a ~ ground truth")
subplot(133)
scatter(Fbcur, Fb)
title("Factor b ~ ground truth")