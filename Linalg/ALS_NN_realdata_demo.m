%% Demo of the ALS factorized fitting method for real neural data and evolved images.
% the UI is improved then the ALS_demo etc. could be copied
% 
% # Result 
% regenerate neural scoring is easy. 
% But the mask is not easy to see. Adding Ridge regression will help make
% the mask smoother. 
net = vgg16;
%% Manifold Experiment
Expi = 3;
imgN = 121; B=10;
si=1;ui=1;Window=50:200;
% dlimg = randn(224,224,3,imgN);
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
score_grid = cellfun(@(psth) mean(psth(ui,Window,:),[2,3]), Stats(Expi).manif.psth{si});
%% Rearrange the score and corresponding image name
imgnm_vect = reshape(imgnm_grid, [], 1);
score_vect = reshape(score_grid, [], 1);
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(Stats(Expi).meta.stimuli, imgnm+".JPG")),[224,224]), ...
    imgnm_vect, 'UniformOutput', false);
%% Correlation Coefficient Clearly shows a spatial structure there
layername = 'conv4_3';
dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
tic
feat_tsr = activations(net, dlimg, layername);
toc
%% Setting Regularization typle for the 2 least squares
opt.regSp = 'nonneg'; %'L1' 'L2' 'nonneg'
opt.paramregSp = 1;
opt.regFt = 'L2'; %'L1' 'L2' 'nonneg'
opt.paramregFt = 1;
%% Try the factorized ALS fitting, 
fprintf("ALS fiting\n")
spN = prod(size(feat_tsr,[1,2]));
Fa_init = randn(spN,1);
Fb_init = randn(size(feat_tsr,3),1);
A = reshape(feat_tsr, [spN, size(feat_tsr,[3,4])]);
bcur = 0;
Facur = Fa_init;
Fbcur = Fb_init;
for k = 1:20
Xcur = double(einsum(A, Fbcur, 'ijk,jl->kil'));
switch opt.regSp
    case 'L2'
    Facur_aug = ridge(score_vect, Xcur,1,0); % ridge regression doesn't require add column of 1
    case 'L1'
    % [B_a,STATS_a] = lasso([Xcur, ones(size(Xcur,1),1)], target);
    % ci = round(size(B_a,2)/2); Facur = B_a(1:end-1, ci); bcur = B_a(end, ci);
    Facur_aug = lasso([ones(size(Xcur,1),1), Xcur], score_vect,'Lambda',opt.paramregSp);
    case 'nonneg'
    Facur_aug = lsqlin([ones(size(Xcur,1),1), Xcur], score_vect, [], [], [], [], ... % positive bound for spatial mask
        [zeros(1, length(Facur)), -inf], []);%, optimoptions('lsqlin','Algorithm','interior-point'));
    otherwise
    Facur_aug = regress(score_vect, [ones(size(Xcur,1),1), Xcur]);
end
% Facur = Facur_aug(2:end); bcur = Facur_aug(1); % ridge regression put the intercept at first column
Facur = Facur_aug(2:end)./max(Facur_aug); bcur = Facur_aug(1)/max(Facur_aug); % normalize the factors, make the spatial mask max to 1.

Xcur = double(einsum(A, Facur, 'ijk,il->kjl'));
switch opt.regFt
    case 'L2'
    Fbcur_aug = ridge(score_vect, Xcur, opt.paramregSp, 0); 
    case 'L1'
%     [B_b,STATS_b] = lasso([Xcur, ones(size(Xcur,1),1)], target);
%     ci = round(size(B_b,2)/2); Fbcur = B_b(1:end-1, ci); bcur = B_b(end, ci);
    Fbcur_aug = lasso([ones(size(Xcur,1),1), Xcur], score_vect,'Lambda',opt.paramregSp);
    case 'nonneg'
    Fbcur_aug = lsqlin([ones(size(Xcur,1),1), Xcur], score_vect, [], [], [], [], ...
        [zeros(1, length(Fbcur)), -inf], []);%, optimoptions('lsqlin','Algorithm','interior-point'));
    otherwise % none
    Fbcur_aug = regress(score_vect, [ones(size(Xcur,1),1), Xcur]);
end
Fbcur = Fbcur_aug(2:end); bcur = Fbcur_aug(1);

pred = einsum(einsum(A, Facur, 'ijk,il->ljk'), Fbcur, 'ljk,jl->k') + bcur;
res = score_vect - pred;
fprintf("residue max %.1f norm %.1f\n", max(res), norm(res))
end
%%
figure; 
subplot(121)
imagesc(reshape(pred,[11,11])) % recontruct the score is easy
subplot(122)
imagesc(reshape(Facur,size(feat_tsr,[1,2]))) % but the mask is not good looking?
