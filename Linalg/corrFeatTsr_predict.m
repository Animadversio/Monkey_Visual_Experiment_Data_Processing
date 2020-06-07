% corrFeatTsr predict
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%% Evol to Manif
Expi = 11; 
ui=1;si=1;
imgN=121; %B=10;
ExpType = "Manif";
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);

tmpfn = ls(fullfile(Stats(Expi).meta.stimuli, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");
suffix = "."+tmpparts{2};
%% Find, Load, Preprocess(resize) images and form an image tensor
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(Stats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
        imgnm_vect, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
layername = "conv4_3";
feat_tsr = activations(net, dlimg, layername);
%% Manif to Evol 
% R.cc_tsr(:,:,:,11)
ccWeight = R.cc_tsr(:,:,:,11);
pred_rsp = mean(feat_tsr .* ccWeight,[1,2,3]);
figure;
imagesc(reshape(pred_rsp,[11,11]))
%%
t_signif_tsr = (R.cc_tsr - cc_refM) ./ cc_refS;
%%
fi = 11;
ccWeight = R.cc_tsr;
ccWeight(t_signif_tsr < 5 & t_signif_tsr > -5) = 0;
pred_rsp = mean(feat_tsr .* ccWeight(:,:,:,fi),[1,2,3]);
figure;
imagesc(-90:18:90,-90:18:90,reshape(pred_rsp,[11,11]))
xlabel("PC3");ylabel("PC2");axis image;colorbar()
title("Use the significantly correlated channels for prediction")
%%
figure;
pred_rsp = mean(feat_tsr,[1,2,3]);
imagesc(reshape(pred_rsp,[11,11]))