%% Predict tuning across experiment

Animal="Alfa";
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(predsavedir,Animal+"_FeatTsrPredStats.mat"),'predStats')
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
prefchan_arr = arrayfun(@(S)S.units.pref_chan,EStats);
%%
srcExpi = 1;srcExpType = "Manif";
tgtExpi = 2;tgtExpType = "Manif";
tic
[cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, srcExpType, 1, "conv4_3");
[imgfullnm_vect, score_vect] = loadManifData(Stats, EStats, Animal, tgtExpType, 2, struct([]));
Wtsr = cc_tsr; Wtsr(abs(t_signif_tsr) < 5) = 0;
crosspred = LinPredictImgs(Wtsr, "conv4_3", imgfullnm_vect);
nlcrosspred = reshape(predStats(1).M2E.nlfunc(crosspred),size(crosspred));
toc
%%
nlpredcorr = corr(nlcrosspred(:,end), score_vect(:,end));
lpredcorr = corr(crosspred(:,end), score_vect(:,end));
nlpredrsquare = 1 - var(score_vect(:,end) - nlcrosspred(:,end)) / var(score_vect(:,end));
%%
figure(9);clf;
subtightplot(1,3,1,[0.02,0.02],0.08,0.03)
imagesc(reshape(score_vect(:,end),[11,11]));
axis image;colorbar();CLIM=caxis();
title("Real Response")
subtightplot(1,3,2,[0.02,0.02],0.08,0.03)
imagesc(reshape(nlcrosspred(:,end),[11,11]));
axis image;colorbar();caxis(CLIM);
title(compose("L-N model\ncc=%.3f r2=%.3f",nlpredcorr,nlpredrsquare))
subtightplot(1,3,3,[0.02,0.02],0.08,0.03)
imagesc(reshape(crosspred(:,end),[11,11]));
axis image;colorbar();
title(compose("L model on VGG %s\ncc=%.3f",layername,lpredcorr))
suptitle(compose("%s Manif Prediction from Manif Exp%d to Exp%d",Animal, srcExpi, tgtExpi))
% Alfa_cross_pred_Manifexp1to2.png
%%
function [pred_rsp] = LinPredictImgs(weight_tsr, layername, imgfullnms)
global net
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(imgnm)),[224,224]), imgfullnms, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[])); 
feat_tsr = activations(net, dlimg, layername,'MiniBatchSize',40);
pred_rsp = einsum(feat_tsr, weight_tsr, 'ijkl,ijkm->lm') / prod(size(feat_tsr, [1,2,3]));
% pred_rsp = mean(feat_tsr .* weight_tsr,[1,2,3]);
end

function [cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, ExpType, Expi, layername)
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
outfn = fullfile(savedir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername)); % LW means long window
load(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
end

function [imgfullnm_vect, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, flags)
ui=1;si=1;
assert(EStats(Expi).Animal == Animal && Stats(Expi).Animal == Animal)
fprintf("Processing %s Exp %d pref chan %d\n",ExpType,Expi,EStats(Expi).units.pref_chan)
if ExpType == "Manif"
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
stimpath = Stats(Expi).meta.stimuli;
% imgN=121; 
elseif ExpType == "Evol"
index_vect = cell2mat(EStats(Expi).evol.idx_seq');
imgnm_vect = EStats(Expi).imageName(index_vect); % reshape(imgnm_grid, [], 1);
stimpath = EStats(Expi).meta.stimuli;
end
imgN = length(imgnm_vect);
tmpfn = ls(fullfile(stimpath, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");suffix = "."+tmpparts{2};
imgfullnm_vect = cellfun(@(imgnm) fullfile(stimpath, imgnm+suffix),imgnm_vect);
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
if ExpType == "Manif"
psth_all = cellfun(@(psth) reshape(mean(psth(ui,:,:),[3]),1,1,[]), Stats(Expi).manif.psth{si}, ...
    'UniformOutput', false); % note there is trial averaging here. 
psth_all = reshape(cell2mat(psth_all),imgN,[]);
elseif ExpType == "Evol"
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; % imgN by 200
end
% This part could be abbrieviated. 
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
end