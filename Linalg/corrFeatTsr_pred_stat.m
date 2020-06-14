%% corrFeatTsr 
global net
net = vgg16;
%% 
%  use the Evolution to predict manifold experiments and vice versa.
%  Use the correlated voxels to predict 
%  Fit a model on Evolution experiment and test on manifold
%  
Animal="Alfa";
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
predsavedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict";
load(fullfile(predsavedir, Animal+"_FeatTsrPredStats.mat"),'predStats')
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
predStats = repmat(struct("E2M",struct(),"M2E",struct()),1,length(Stats));

%% This is the least square fit. For Poisson fit see the poisson_curv_fit demo
% ft = fittype( 'max(0, a*x - b)+c', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'max(0, a*(x-b))+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Lower = [0    -10000 0];
opts.Upper = [10000 10000 10000];
%%
EStats(19).meta.stimuli = "N:\Stimuli\2019-Manifold\alfa-191210a\backup_12_10_2019_13_07_57";
%% 
Animal="Alfa"; Expi = 11; layername = "fc6"; 
for Expi = 1:46
[imgfullnm_Evol, score_Evol] = loadManifData(Stats, EStats, Animal, "Evol", Expi, struct());
[imgfullnm_Manif, score_Manif] = loadManifData(Stats, EStats, Animal, "Manif", Expi, struct());
predStats(Expi).score_Evol = score_Evol;
predStats(Expi).score_Manif = score_Manif;
%% Manifold to Evol 
% layername = "conv4_3";%fi = 11;
tic
fprintf("Fitting Manif to Evol nonlinearity\n")
[cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, "Manif", Expi, layername); % correlation coefficient for manifold experiments
ccWeight = cc_tsr;
ccWeight(abs(t_signif_tsr) < 5) = 0;%threshold and get a weight tensor to do prediction
lfitscore_manif = LinPredictImgs(ccWeight, layername, imgfullnm_Manif);
% Fit a thresholding bias and a scaling at manifold exp
opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
        'Lower',[0, min(lfitscore_manif(:,end)), 0],...
        'Upper',[10000, max(lfitscore_manif(:,end)), max(score_Manif(:,end))]);
[relufit, gof, out] = fit(lfitscore_manif(:,end), score_Manif(:,end), ft, opts);
% fit the bias on 
NLpred = relufit(lfitscore_manif);
NLpred = reshape(NLpred, size(lfitscore_manif));
nlfitcorr = corr(NLpred(:,end), score_Manif(:,end));
lfitcorr = corr(lfitscore_manif(:,end), score_Manif(:,end));
toc
%%
tic
fprintf("Predicting Evolution experiments\n")
pred_evol_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_Evol);
NLpred_evol = relufit(pred_evol_rsp);
NLpred_evol = reshape(NLpred_evol, size(score_Evol));
nlpredcorr = corr(NLpred_evol(:,end),score_Evol(:,end));
lpredcorr = corr(pred_evol_rsp(:,end),score_Evol(:,end));
nlpredrsquare = 1 - var(NLpred_evol(:,end)-score_Evol(:,end)) / var(score_Evol(:,end));
toc
%%
predStats(Expi).M2E.nlfunc = relufit;
predStats(Expi).M2E.lfitscore = lfitscore_manif;
predStats(Expi).M2E.nlfitscore = NLpred;
predStats(Expi).M2E.nlfitcorr = nlfitcorr;
predStats(Expi).M2E.lfitcorr = lfitcorr;
predStats(Expi).M2E.gof = gof;
predStats(Expi).M2E.lpredscore = pred_evol_rsp;
predStats(Expi).M2E.nlpredscore = NLpred_evol;
predStats(Expi).M2E.nlpredcorr = nlpredcorr;
predStats(Expi).M2E.lpredcorr = lpredcorr;
predStats(Expi).M2E.nlpredrsquare = nlpredrsquare;
%%
%% Visualization
figure(1);clf;
subtightplot(1,3,1,[0.02,0.02],0.08,0.03)
imagesc(reshape(score_Manif(:,end),[11,11]));
axis image;colorbar();CLIM=caxis();
title("Real Response")
subtightplot(1,3,2,[0.02,0.02],0.08,0.03)
imagesc(reshape(NLpred(:,end),[11,11]));
axis image;colorbar();caxis(CLIM);
title(compose("L-N model\ncc=%.3f r2=%.3f",nlfitcorr,gof.rsquare))
subtightplot(1,3,3,[0.02,0.02],0.08,0.03)
imagesc(reshape(lfitscore_manif(:,end),[11,11]));
axis image;colorbar();
title(compose("L model on VGG %s\ncc=%.3f",layername,lfitcorr))
suptitle(compose("%s Manif Exp%d Goodness of NL Fit",Animal, Expi))

figure(5);clf;hold on
subtightplot(1,2,1,[0.03,0.07],0.08,0.07)
imagesc(10:10:190,1:size(score_Evol,1),NLpred_evol(:,1:end-5));% non-linearly predicted evolution score
title(compose("L-N model\n cc=%.3f r2=%.3f", nlpredcorr, nlpredrsquare))
subtightplot(1,2,2,[0.03,0.07],0.08,0.07)
imagesc(10:10:190,1:size(score_Evol,1),score_Evol(:,1:end-5));% neural score
title(compose("real psth (20ms mean rate)"))
suptitle(compose("%s Evol Exp%d (Manif2Evol) Prediction",Animal, Expi))

figure(4);clf
subtightplot(1,1,1,[0.02,0.02],0.08,0.03);hold on
plot(NLpred_evol(:,end)); % non-linearly predicted response
plot(score_Evol(:,end)); % Original score
plot(pred_evol_rsp(:,end)); % linearly predicted response
legend(["L-N Prediction", "Real FR", "L Prediction"])
xlim([0,size(score_Evol,1)])
suptitle(compose("%s Evol Exp%d (Manif2Evol) Prediction (50-200ms rate)\n L-N model cc=%.3f r2=%.3f, L model cc=%.3f",Animal, Expi, nlpredcorr, nlpredrsquare, lpredcorr))
saveas(1,fullfile(predsavedir,compose("%s_Exp%d_M2E_%s_ManifFit.jpg",Animal,Expi,layername)));
saveas(5,fullfile(predsavedir,compose("%s_Exp%d_M2E_%s_EvolPSTHPred.jpg",Animal,Expi,layername)));
saveas(4,fullfile(predsavedir,compose("%s_Exp%d_M2E_%s_EvolPred.jpg",Animal,Expi,layername)));

%% Evol to Manifold %%
tic
fprintf("Fitting Evol to Manif nonlinearity\n")
[cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, "Evol", Expi, layername); % correlation coefficient for manifold experiments
ccWeight = cc_tsr;
ccWeight(abs(t_signif_tsr) < 5) = 0;
lfitscore_evol = LinPredictImgs(ccWeight, layername, imgfullnm_Evol); % linearly predicted response
% Fit a thresholding bias and a scaling at manifold exp
opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
        'Lower',[0, min(lfitscore_evol(:,end)), 0],...
        'Upper',[10000, max(lfitscore_evol(:,end)), max(score_Evol(:,end))]);
[relufit, gof, out] = fit(lfitscore_evol(:,end), score_Evol(:,end), ft, opts); % fit the bias on evolution data
nlfitscore_evol = relufit(lfitscore_evol);
nlfitscore_evol = reshape(nlfitscore_evol, size(lfitscore_evol));
nlfitcorr = corr(nlfitscore_evol(:,end), score_Evol(:,end));
lfitcorr = corr(lfitscore_evol(:,end), score_Evol(:,end));
toc
%%
tic
fprintf("Predicting Manifold experiments\n")
pred_manif_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_Manif); % use the model to predict manifold exp
NLpred_manif = relufit(pred_manif_rsp);
NLpred_manif = reshape(NLpred_manif, size(score_Manif));
nlpredcorr = corr(NLpred_manif(:,end),score_Manif(:,end));
lpredcorr = corr(pred_manif_rsp(:,end),score_Manif(:,end));
nlpredrsquare = 1 - var(NLpred_manif(:,end)-score_Manif(:,end)) / var(score_Manif(:,end));
toc
%%
predStats(Expi).E2M.nlfunc = relufit;
predStats(Expi).E2M.lfitscore = lfitscore_evol;
predStats(Expi).E2M.nlfitscore = nlfitscore_evol;
predStats(Expi).E2M.nlfitcorr = nlfitcorr;
predStats(Expi).E2M.lfitcorr = lfitcorr;
predStats(Expi).E2M.gof = gof;
predStats(Expi).E2M.lpredscore = pred_manif_rsp;
predStats(Expi).E2M.nlpredscore = NLpred_manif;
predStats(Expi).E2M.nlpredcorr = nlpredcorr;
predStats(Expi).E2M.lpredcorr = lpredcorr;
predStats(Expi).E2M.nlpredrsquare = nlpredrsquare;
%%
%% Visualization
figure(1);clf;%set(1,'position',[40         187        1904         595])
subtightplot(1,3,1,[0.02,0.02],0.08,0.03)
imagesc(reshape(score_Manif(:,end),[11,11]));
axis image;colorbar();CLIM=caxis();
title("Real Response")
subtightplot(1,3,2,[0.02,0.02],0.08,0.03)
imagesc(reshape(NLpred_manif(:,end),[11,11]));
axis image;colorbar();caxis(CLIM);
title(compose("L-N model\ncc=%.3f r2=%.3f",nlpredcorr,nlpredrsquare))
subtightplot(1,3,3,[0.02,0.02],0.08,0.03)
imagesc(reshape(pred_manif_rsp(:,end),[11,11]));%lfitscore_manif there is a bug before
axis image;colorbar();
title(compose("L model on VGG %s\ncc=%.3f",layername,lpredcorr))
suptitle(compose("%s Manif Exp%d (pref chan %s) Prediction",Animal, Expi, prefchan_lab))

figure(5);clf;hold on;%set(5,'position',[718   199   513   766])
subtightplot(1,2,1,[0.03,0.07],0.08,0.07)
imagesc(10:10:190,1:size(score_Evol,1),nlfitscore_evol(:,1:end-5));
title(compose("L-N model\n cc=%.3f r2=%.3f", nlfitcorr, gof.rsquare))
subtightplot(1,2,2,[0.03,0.07],0.08,0.07)
imagesc(10:10:190,1:size(score_Evol,1),score_Evol(:,1:end-5))
title(compose("real psth (20ms mean rate)"))
suptitle(compose("%s Evol Exp%d (pref chan %s)\n(Manif2Evol) Goodness of NL Fit",Animal, Expi, prefchan_lab))

figure(4);clf;%set(4,'position',[332         297        1557         444])
subtightplot(1,1,1,[0.02,0.02],0.08,0.03);hold on
plot(nlfitscore_evol(:,end)); % fit 
plot(score_Evol(:,end)) % real neural score
plot(lfitscore_evol(:,end)) % linearly fit evolution response
legend(["L-N Prediction", "Real FR", "L Prediction"])
xlim([0,size(score_Evol,1)])
suptitle(compose("%s Evol Exp%d (pref chan %s) (Manif2Evol) Goodness of NL Fit (50-200ms rate)\n L-N model cc=%.3f r2=%.3f, L model cc=%.3f",Animal, Expi, prefchan_lab, nlfitcorr, gof.rsquare, lfitcorr))
saveas(1,fullfile(predsavedir,compose("%s_Exp%d_E2M_%s_ManifPred.jpg",Animal,Expi,layername)));
saveas(5,fullfile(predsavedir,compose("%s_Exp%d_E2M_%s_EvolPSTHFit.jpg",Animal,Expi,layername)));
saveas(4,fullfile(predsavedir,compose("%s_Exp%d_E2M_%s_EvolFit.jpg",Animal,Expi,layername)));
end
%%
save(fullfile(predsavedir,compose("%s_FeatTsrPredStats_%s.mat",Animal,layername)),'predStats')
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