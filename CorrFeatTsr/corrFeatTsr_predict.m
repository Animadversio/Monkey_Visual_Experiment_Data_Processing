% corrFeatTsr predict
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%% Evol to Manif
Expi = 11; 
ui=1; si=1;
imgN=121; %B=10;
ExpType = "Manif";
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
tmpfn = ls(fullfile(Stats(Expi).meta.stimuli, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");
suffix = "."+tmpparts{2};
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; % imgN by 200
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]

%% Find, Load, Preprocess(resize) images and form an image tensor
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(Stats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
        imgnm_vect, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
layername = "conv4_3";
feat_tsr = activations(net, dlimg, layername);
%%
ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
outfn = fullfile(ccmat_dir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername));
R = load(outfn,'cc_tsr', 'MFeat', 'StdFeat', 'wdw_vect');
%%
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
%% Manif to Evol 

%%
figure;
pred_rsp = mean(feat_tsr,[1,2,3]);
imagesc(reshape(pred_rsp,[11,11]))