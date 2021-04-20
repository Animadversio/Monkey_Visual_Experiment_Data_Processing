%% 
% Plan of this script is quite similar to pre_stat.m

global net
net = vgg16;

Animal="Alfa";
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
predsavedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict";
load(fullfile(predsavedir, Animal+"_FeatTsrPredStats.mat"),'predStats')
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
ft = fittype( 'max(0, a*(x-b))+c', 'independent', 'x', 'dependent', 'y' );
%%
diary("E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict\fitting_process.log")
predStats2 = repmat(struct("E2M",struct(),"M2E",struct()),1,length(Stats));
%%
Animal="Alfa"; Expi = 3; 
for Expi = 2:numel(Stats)
[imgfullnm_Evol, score_Evol] = loadPairedData(Stats, EStats, Animal, "Evol", Expi, struct()); % uniform interface to load spikes and data
[imgfullnm_Manif, score_Manif] = loadPairedData(Stats, EStats, Animal, "Manif", Expi, struct());
[imgfullnm_EvoRef, score_EvoRef] = loadPairedData(Stats, EStats, Animal, "EvoRef", Expi, struct());
%
[imgfullnm_Gabor, score_Gabor] = loadPairedData(Stats, EStats, Animal, "Gabor", Expi, struct());
[imgfullnm_Pasu, score_Pasu] = loadPairedData(Stats, EStats, Animal, "Pasu", Expi, struct());
% predStats(Expi).score_Evol = score_Evol;
% predStats(Expi).score_Manif = score_Manif;
%%
for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
fprintf("Using the correlation tensor from vgg16 %s layer\n", layername)
tic
fprintf("Fitting Manif to Evol nonlinearity\n")
[cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, "Manif", Expi, layername); % correlation coefficient for manifold experiments
ccWeight = cc_tsr;ccWeight(isnan(cc_tsr))=0; % there can be nans .....
thresh = prctile(t_signif_tsr(:,:,:,end),99,'all');
fprintf("Threshold %.3f (%d)\n",thresh,sum(t_signif_tsr(:,:,:,end) > thresh,'all'))
ccWeight(t_signif_tsr < thresh) = 0; %threshold and get a weight tensor to do prediction abs
lfitscore_manif = LinPredictImgs(ccWeight, layername, imgfullnm_Manif);
% Fit a thresholding bias and a scaling at manifold exp
% end the final col is the 50-200 range activation.
opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
        'Lower',[0, min(lfitscore_manif(:,end)), 0],...
        'Upper',[10000, max(lfitscore_manif(:,end)), max(score_Manif(:,end))]);
[relufit, gof, out] = fit(lfitscore_manif(:,end), score_Manif(:,end), ft, opts); 
% fit the bias on 
NLpred = relufit(lfitscore_manif);
NLpred = reshape(NLpred, size(lfitscore_manif));
[nlfitcorr,nlfitcorr_P] = corr(NLpred(:,end), score_Manif(:,end));
[lfitcorr,lfitcorr_P] = corr(lfitscore_manif(:,end), score_Manif(:,end));
fprintf("Fitting Manifold: Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e)\n",lfitcorr,lfitcorr_P,nlfitcorr,nlfitcorr_P)
toc
%%
tic
pred_EvoRef_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_EvoRef);
NLpred_EvoRef = relufit(pred_EvoRef_rsp);
NLpred_EvoRef = reshape(NLpred_EvoRef, size(score_EvoRef));
[nlpredcorr_evoref, nlpredcorr_P_evoref] = corr(NLpred_EvoRef(:,end),score_EvoRef(:,end));
[lpredcorr_evoref, lpredcorr_P_evoref] = corr(pred_EvoRef_rsp(:,end),score_EvoRef(:,end));
nlpredR2_evoref = 1 - var(NLpred_EvoRef(:,end)-score_EvoRef(:,end)) / var(score_EvoRef(:,end));
fprintf("Predicting EvoRef Images, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
    lpredcorr_evoref,lpredcorr_P_evoref,nlpredcorr_evoref,nlpredcorr_P_evoref,nlpredR2_evoref,numel(score_EvoRef(:,end)))
toc
%
if Stats(Expi).ref.didGabor
tic
pred_Gabor_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_Gabor);
NLpred_Gabor = relufit(pred_Gabor_rsp);
NLpred_Gabor = reshape(NLpred_Gabor, size(score_Gabor));
[nlpredcorr_gabor, nlpredcorr_P_gabor] = corr(NLpred_Gabor(:,end),score_Gabor(:,end));
[lpredcorr_gabor, lpredcorr_P_gabor] = corr(pred_Gabor_rsp(:,end),score_Gabor(:,end));
nlpredR2_gabor = 1 - var(NLpred_Gabor(:,end)-score_Gabor(:,end)) / var(score_Gabor(:,end));
fprintf("Predicting Gabor Images, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
    lpredcorr_gabor,lpredcorr_P_gabor,nlpredcorr_gabor,nlpredcorr_P_gabor,nlpredR2_gabor,numel(score_Gabor(:,end)))
toc
else
lpredcorr_gabor = nan;
lpredcorr_P_gabor = nan;
nlpredcorr_gabor = nan;
nlpredcorr_P_gabor = nan;
nlpredR2_gabor = nan;
pred_Gabor_rsp = [];
NLpred_Gabor = [];
end
%
if Stats(Expi).ref.didPasu
tic
pred_Pasu_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_Pasu);
NLpred_Pasu = relufit(pred_Pasu_rsp);
NLpred_Pasu = reshape(NLpred_Pasu, size(score_Pasu));
[nlpredcorr_pasu, nlpredcorr_P_pasu] = corr(NLpred_Pasu(:,end),score_Pasu(:,end));
[lpredcorr_pasu, lpredcorr_P_pasu] = corr(pred_Pasu_rsp(:,end),score_Pasu(:,end));
nlpredR2_pasu = 1 - var(NLpred_Pasu(:,end)-score_Pasu(:,end)) / var(score_Pasu(:,end));
fprintf("Predicting Pasu Images, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
    lpredcorr_pasu,lpredcorr_P_pasu,nlpredcorr_pasu,nlpredcorr_P_pasu,nlpredR2_pasu,numel(score_Pasu(:,end)))
toc
else
lpredcorr_pasu = nan;
lpredcorr_P_pasu = nan;
nlpredcorr_pasu = nan;
nlpredcorr_P_pasu = nan;
nlpredR2_pasu = nan;
pred_Pasu_rsp = [];
NLpred_Pasu = [];
end
tic
fprintf("Predicting Evolution experiments\n")
pred_evol_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_Evol);
NLpred_evol = relufit(pred_evol_rsp);
NLpred_evol = reshape(NLpred_evol, size(score_Evol));
[nlpredcorr, nlpredcorr_P] = corr(NLpred_evol(:,end),score_Evol(:,end));
[lpredcorr, lpredcorr_P] = corr(pred_evol_rsp(:,end),score_Evol(:,end));
nlpredR2 = 1 - var(NLpred_evol(:,end)-score_Evol(:,end)) / var(score_Evol(:,end));
fprintf("Predicting Evolution, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
	lpredcorr,lpredcorr_P,nlpredcorr,nlpredcorr_P,nlpredR2,numel(score_Evol(:,end)))
toc
predStats2(Expi).M2E.(layername).manif.origscore = score_Manif;
predStats2(Expi).M2E.(layername).manif.lfitscore = lfitscore_manif;
predStats2(Expi).M2E.(layername).manif.nlfitscore = NLpred;
predStats2(Expi).M2E.(layername).manif.nlfitcorr = nlfitcorr;
predStats2(Expi).M2E.(layername).manif.nlfitcorr_P = nlfitcorr_P;
predStats2(Expi).M2E.(layername).manif.lfitcorr = lfitcorr;
predStats2(Expi).M2E.(layername).manif.lfitcorr_P = lfitcorr_P;
predStats2(Expi).M2E.(layername).manif.nlfunc = relufit;
predStats2(Expi).M2E.(layername).manif.gof = gof;
predStats2(Expi).M2E.(layername).evol.lpredscore = pred_evol_rsp;
predStats2(Expi).M2E.(layername).evol.nlpredscore = NLpred_evol;
predStats2(Expi).M2E.(layername).evol.nlpredcorr = nlpredcorr;
predStats2(Expi).M2E.(layername).evol.nlpredcorr_P = nlpredcorr_P;
predStats2(Expi).M2E.(layername).evol.lpredcorr = lpredcorr;
predStats2(Expi).M2E.(layername).evol.lpredcorr_P = lpredcorr_P;
predStats2(Expi).M2E.(layername).evol.nlpredR2 = nlpredR2;
predStats2(Expi).M2E.(layername).evol.origscore = score_Evol;
predStats2(Expi).M2E.(layername).evoref.lpredscore = pred_EvoRef_rsp;
predStats2(Expi).M2E.(layername).evoref.nlpredscore = NLpred_EvoRef;
predStats2(Expi).M2E.(layername).evoref.nlpredcorr = nlpredcorr_evoref;
predStats2(Expi).M2E.(layername).evoref.nlpredcorr_P = nlpredcorr_P_evoref;
predStats2(Expi).M2E.(layername).evoref.lpredcorr = lpredcorr_evoref;
predStats2(Expi).M2E.(layername).evoref.lpredcorr_P = lpredcorr_P_evoref;
predStats2(Expi).M2E.(layername).evoref.nlpredR2 = nlpredR2_evoref;
predStats2(Expi).M2E.(layername).evoref.origscore = score_EvoRef;
predStats2(Expi).M2E.(layername).gabor.lpredscore = pred_Gabor_rsp;
predStats2(Expi).M2E.(layername).gabor.nlpredscore = NLpred_Gabor;
predStats2(Expi).M2E.(layername).gabor.nlpredcorr = nlpredcorr_gabor;
predStats2(Expi).M2E.(layername).gabor.nlpredcorr_P = nlpredcorr_P_gabor;
predStats2(Expi).M2E.(layername).gabor.lpredcorr = lpredcorr_gabor;
predStats2(Expi).M2E.(layername).gabor.lpredcorr_P = lpredcorr_P_gabor;
predStats2(Expi).M2E.(layername).gabor.nlpredR2 = nlpredR2_gabor;
predStats2(Expi).M2E.(layername).gabor.origscore = score_Gabor;
predStats2(Expi).M2E.(layername).pasu.lpredscore = pred_Pasu_rsp;
predStats2(Expi).M2E.(layername).pasu.nlpredscore = NLpred_Pasu;
predStats2(Expi).M2E.(layername).pasu.nlpredcorr = nlpredcorr_pasu;
predStats2(Expi).M2E.(layername).pasu.nlpredcorr_P = nlpredcorr_P_pasu;
predStats2(Expi).M2E.(layername).pasu.lpredcorr = lpredcorr_pasu;
predStats2(Expi).M2E.(layername).pasu.lpredcorr_P = lpredcorr_P_pasu;
predStats2(Expi).M2E.(layername).pasu.nlpredR2 = nlpredR2_pasu;
predStats2(Expi).M2E.(layername).pasu.origscore = score_Pasu;
end
end
diary off
%%
save(fullfile(MatStats_path,Animal+"_FeatTsrPredStats2.mat"),'predStats2')
%%
diary("E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict\fitting_process.log")
predStats2 = repmat(struct("E2M",struct(),"M2E",struct()),1,length(Stats));
%%
Animal="Alfa"; Expi = 3; 
for Expi = 2:numel(Stats)
[imgfullnm_Evol, score_Evol] = loadPairedData(Stats, EStats, Animal, "Evol", Expi, struct()); % uniform interface to load spikes and data
[imgfullnm_Manif, score_Manif] = loadPairedData(Stats, EStats, Animal, "Manif", Expi, struct());
[imgfullnm_EvoRef, score_EvoRef] = loadPairedData(Stats, EStats, Animal, "EvoRef", Expi, struct());
%
[imgfullnm_Gabor, score_Gabor] = loadPairedData(Stats, EStats, Animal, "Gabor", Expi, struct());
[imgfullnm_Pasu, score_Pasu] = loadPairedData(Stats, EStats, Animal, "Pasu", Expi, struct());
% predStats(Expi).score_Evol = score_Evol;
% predStats(Expi).score_Manif = score_Manif;
%%
for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
fprintf("Using the correlation tensor from vgg16 %s layer\n", layername)
tic
fprintf("Fitting Manif to Evol nonlinearity\n")
[cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, "Manif", Expi, layername); % correlation coefficient for manifold experiments
ccWeight = cc_tsr;ccWeight(isnan(cc_tsr))=0; % there can be nans .....
thresh = prctile(t_signif_tsr(:,:,:,end),99,'all');
fprintf("Threshold %.3f (%d)\n",thresh,sum(t_signif_tsr(:,:,:,end) > thresh,'all'))
ccWeight(t_signif_tsr < thresh) = 0; %threshold and get a weight tensor to do prediction abs
lfitscore_manif = LinPredictImgs(ccWeight, layername, imgfullnm_Manif);
% Fit a thresholding bias and a scaling at manifold exp
% end the final col is the 50-200 range activation.
opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
        'Lower',[0, min(lfitscore_manif(:,end)), 0],...
        'Upper',[10000, max(lfitscore_manif(:,end)), max(score_Manif(:,end))]);
[relufit, gof, out] = fit(lfitscore_manif(:,end), score_Manif(:,end), ft, opts); 
% fit the bias on 
NLpred = relufit(lfitscore_manif);
NLpred = reshape(NLpred, size(lfitscore_manif));
[nlfitcorr,nlfitcorr_P] = corr(NLpred(:,end), score_Manif(:,end));
[lfitcorr,lfitcorr_P] = corr(lfitscore_manif(:,end), score_Manif(:,end));
fprintf("Fitting Manifold: Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e)\n",lfitcorr,lfitcorr_P,nlfitcorr,nlfitcorr_P)
toc
%%
tic
pred_EvoRef_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_EvoRef);
NLpred_EvoRef = relufit(pred_EvoRef_rsp);
NLpred_EvoRef = reshape(NLpred_EvoRef, size(score_EvoRef));
[nlpredcorr_evoref, nlpredcorr_P_evoref] = corr(NLpred_EvoRef(:,end),score_EvoRef(:,end));
[lpredcorr_evoref, lpredcorr_P_evoref] = corr(pred_EvoRef_rsp(:,end),score_EvoRef(:,end));
nlpredR2_evoref = 1 - var(NLpred_EvoRef(:,end)-score_EvoRef(:,end)) / var(score_EvoRef(:,end));
fprintf("Predicting EvoRef Images, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
    lpredcorr_evoref,lpredcorr_P_evoref,nlpredcorr_evoref,nlpredcorr_P_evoref,nlpredR2_evoref,numel(score_EvoRef(:,end)))
toc
%
if Stats(Expi).ref.didGabor
tic
pred_Gabor_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_Gabor);
NLpred_Gabor = relufit(pred_Gabor_rsp);
NLpred_Gabor = reshape(NLpred_Gabor, size(score_Gabor));
[nlpredcorr_gabor, nlpredcorr_P_gabor] = corr(NLpred_Gabor(:,end),score_Gabor(:,end));
[lpredcorr_gabor, lpredcorr_P_gabor] = corr(pred_Gabor_rsp(:,end),score_Gabor(:,end));
nlpredR2_gabor = 1 - var(NLpred_Gabor(:,end)-score_Gabor(:,end)) / var(score_Gabor(:,end));
fprintf("Predicting Gabor Images, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
    lpredcorr_gabor,lpredcorr_P_gabor,nlpredcorr_gabor,nlpredcorr_P_gabor,nlpredR2_gabor,numel(score_Gabor(:,end)))
toc
else
lpredcorr_gabor = nan;
lpredcorr_P_gabor = nan;
nlpredcorr_gabor = nan;
nlpredcorr_P_gabor = nan;
nlpredR2_gabor = nan;
pred_Gabor_rsp = [];
NLpred_Gabor = [];
end
%
if Stats(Expi).ref.didPasu
tic
pred_Pasu_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_Pasu);
NLpred_Pasu = relufit(pred_Pasu_rsp);
NLpred_Pasu = reshape(NLpred_Pasu, size(score_Pasu));
[nlpredcorr_pasu, nlpredcorr_P_pasu] = corr(NLpred_Pasu(:,end),score_Pasu(:,end));
[lpredcorr_pasu, lpredcorr_P_pasu] = corr(pred_Pasu_rsp(:,end),score_Pasu(:,end));
nlpredR2_pasu = 1 - var(NLpred_Pasu(:,end)-score_Pasu(:,end)) / var(score_Pasu(:,end));
fprintf("Predicting Pasu Images, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
    lpredcorr_pasu,lpredcorr_P_pasu,nlpredcorr_pasu,nlpredcorr_P_pasu,nlpredR2_pasu,numel(score_Pasu(:,end)))
toc
else
lpredcorr_pasu = nan;
lpredcorr_P_pasu = nan;
nlpredcorr_pasu = nan;
nlpredcorr_P_pasu = nan;
nlpredR2_pasu = nan;
pred_Pasu_rsp = [];
NLpred_Pasu = [];
end
tic
fprintf("Predicting Evolution experiments\n")
pred_evol_rsp = LinPredictImgs(ccWeight, layername, imgfullnm_Evol);
NLpred_evol = relufit(pred_evol_rsp);
NLpred_evol = reshape(NLpred_evol, size(score_Evol));
[nlpredcorr, nlpredcorr_P] = corr(NLpred_evol(:,end),score_Evol(:,end));
[lpredcorr, lpredcorr_P] = corr(pred_evol_rsp(:,end),score_Evol(:,end));
nlpredR2 = 1 - var(NLpred_evol(:,end)-score_Evol(:,end)) / var(score_Evol(:,end));
fprintf("Predicting Evolution, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
	lpredcorr,lpredcorr_P,nlpredcorr,nlpredcorr_P,nlpredR2,numel(score_Evol(:,end)))
toc
predStats2(Expi).E2M.(layername).evol.lpredscore = pred_evol_rsp;
predStats2(Expi).E2M.(layername).evol.nlpredscore = NLpred_evol;
predStats2(Expi).E2M.(layername).evol.nlpredcorr = nlpredcorr;
predStats2(Expi).E2M.(layername).evol.nlpredcorr_P = nlpredcorr_P;
predStats2(Expi).E2M.(layername).evol.lpredcorr = lpredcorr;
predStats2(Expi).E2M.(layername).evol.lpredcorr_P = lpredcorr_P;
predStats2(Expi).E2M.(layername).evol.nlpredR2 = nlpredR2;

predStats2(Expi).E2M.(layername).manif.origscore = score_Manif;
predStats2(Expi).E2M.(layername).manif.lfitscore = lfitscore_manif;
predStats2(Expi).E2M.(layername).manif.nlfitscore = NLpred;
predStats2(Expi).E2M.(layername).manif.nlfitcorr = nlfitcorr;
predStats2(Expi).E2M.(layername).manif.nlfitcorr_P = nlfitcorr_P;
predStats2(Expi).E2M.(layername).manif.lfitcorr = lfitcorr;
predStats2(Expi).E2M.(layername).manif.lfitcorr_P = lfitcorr_P;
predStats2(Expi).E2M.(layername).manif.nlfunc = relufit;
predStats2(Expi).E2M.(layername).manif.gof = gof;

predStats2(Expi).E2M.(layername).evol.origscore = score_Evol;
predStats2(Expi).E2M.(layername).evoref.lpredscore = pred_EvoRef_rsp;
predStats2(Expi).E2M.(layername).evoref.nlpredscore = NLpred_EvoRef;
predStats2(Expi).E2M.(layername).evoref.nlpredcorr = nlpredcorr_evoref;
predStats2(Expi).E2M.(layername).evoref.nlpredcorr_P = nlpredcorr_P_evoref;
predStats2(Expi).E2M.(layername).evoref.lpredcorr = lpredcorr_evoref;
predStats2(Expi).E2M.(layername).evoref.lpredcorr_P = lpredcorr_P_evoref;
predStats2(Expi).E2M.(layername).evoref.nlpredR2 = nlpredR2_evoref;
predStats2(Expi).E2M.(layername).evoref.origscore = score_EvoRef;
predStats2(Expi).E2M.(layername).gabor.lpredscore = pred_Gabor_rsp;
predStats2(Expi).E2M.(layername).gabor.nlpredscore = NLpred_Gabor;
predStats2(Expi).E2M.(layername).gabor.nlpredcorr = nlpredcorr_gabor;
predStats2(Expi).E2M.(layername).gabor.nlpredcorr_P = nlpredcorr_P_gabor;
predStats2(Expi).E2M.(layername).gabor.lpredcorr = lpredcorr_gabor;
predStats2(Expi).E2M.(layername).gabor.lpredcorr_P = lpredcorr_P_gabor;
predStats2(Expi).E2M.(layername).gabor.nlpredR2 = nlpredR2_gabor;
predStats2(Expi).E2M.(layername).gabor.origscore = score_Gabor;
predStats2(Expi).E2M.(layername).pasu.lpredscore = pred_Pasu_rsp;
predStats2(Expi).E2M.(layername).pasu.nlpredscore = NLpred_Pasu;
predStats2(Expi).E2M.(layername).pasu.nlpredcorr = nlpredcorr_pasu;
predStats2(Expi).E2M.(layername).pasu.nlpredcorr_P = nlpredcorr_P_pasu;
predStats2(Expi).E2M.(layername).pasu.lpredcorr = lpredcorr_pasu;
predStats2(Expi).E2M.(layername).pasu.lpredcorr_P = lpredcorr_P_pasu;
predStats2(Expi).E2M.(layername).pasu.nlpredR2 = nlpredR2_pasu;
predStats2(Expi).E2M.(layername).pasu.origscore = score_Pasu;
end
end
diary off
%%



%% Util functions 
function [pred_rsp] = LinPredictImgs(weight_tsr, layername, imgfullnms)
global net
imgcol = cellfun(@(imgnm) imresize(to_rgb(imread(fullfile(imgnm))),[224,224]), imgfullnms, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[])); 
if size(dlimg,3)==1,dlimg = repmat(dlimg,1,1,3,1);end
feat_tsr = activations(net, dlimg, layername,'MiniBatchSize',40);
pred_rsp = einsum(feat_tsr, weight_tsr, 'ijkl,ijkm->lm') / prod(size(feat_tsr, [1,2,3]));
% pred_rsp = mean(feat_tsr .* weight_tsr,[1,2,3]);
end

function img = to_rgb(img)
if size(img,3)==1
   img = repmat(img,1,1,3);
end
end

function [cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, ExpType, Expi, layername)
% uniform interface to load the correlation tensor pre-computed with VGG activations. 
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
outfn = fullfile(savedir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername)); % LW means long window
load(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
end

function [imgfullnm_vect, score_vect] = loadPairedData(Stats, EStats, Animal, ExpType, Expi, flags)
% Uniform interface to load spikes and image path data, really well done!
% Maybe we should load some internal variance statistics to regularize the fit...
ui=1;si=1;
assert(EStats(Expi).Animal == Animal && Stats(Expi).Animal == Animal)
fprintf("Processing %s Exp %d pref chan %d\n",ExpType,Expi,EStats(Expi).units.pref_chan)
if ExpType == "Manif"
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
stimpath = Stats(Expi).meta.stimuli;
psth_all = cell2mat(cellfun(@(psth) mean(psth(ui,:,:),[3]), reshape(Stats(Expi).manif.psth{si},[],1),'Uni', 0)); ...);
% note there is trial averaging here!
% psth_all = reshape(cell2mat(psth_all),imgN,[]);
% imgN=121; 
elseif ExpType == "Evol"
index_vect = cell2mat(EStats(Expi).evol.idx_seq'); % concat the index into a vertical vect
imgnm_vect = EStats(Expi).imageName(index_vect); % reshape(imgnm_grid, [], 1);
stimpath = EStats(Expi).meta.stimuli;
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; 
elseif ExpType == "EvoRef"
imgnm_vect = EStats(Expi).ref.imgnm;
imgfullnm_vect = EStats(Expi).ref.impaths;
psth_all = cell2mat(cellfun(@(psth)mean(psth(ui,:,:),3),reshape(EStats(Expi).ref.psth_arr,[],1),'uni',0));
elseif ExpType == "Gabor"
if ~Stats(Expi).ref.didGabor,  
imgfullnm_vect = []; score_vect = []; return
else,
idx_vect = reshape(Stats(Expi).ref.gab_idx_grid,[],1);
validmsk = ~cellfun(@isempty,idx_vect);
idx_vect = idx_vect(validmsk);
imgnm_vect = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), idx_vect));
stimpath = "N:\Stimuli\2020-manifold-references";%Stats(Expi).meta.stimuli;
psth_all = cell2mat(cellfun(@(psth)mean(psth(ui,:,:),3),reshape(Stats(Expi).ref.gab_psths,[],1),'uni',0));
psth_all = psth_all(validmsk,:);
end
elseif ExpType == "Pasu"
if ~Stats(Expi).ref.didPasu, 
imgfullnm_vect = []; score_vect = []; return
else,
idx_vect = reshape(Stats(Expi).ref.pasu_idx_grid',[],1);
validmsk = ~cellfun(@isempty,idx_vect);
idx_vect = idx_vect(validmsk);
imgnm_vect = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), idx_vect));
stimpath = "N:\Stimuli\2020-manifold-references";%Stats(Expi).meta.stimuli;
psth_all = cell2mat(cellfun(@(psth)mean(psth(ui,:,:),3),reshape(Stats(Expi).ref.pasu_psths',[],1),'uni',0));
psth_all = psth_all(validmsk,:);
end
end
imgN = length(imgnm_vect);
if ~(ExpType == "EvoRef")
tmpfn = ls(fullfile(stimpath, imgnm_vect(1)+"*"));
if isempty(tmpfn), error("Stimpath not correct, examine!");end
tmpparts = split(tmpfn,".");suffix = "."+tmpparts{end};
imgfullnm_vect = cellfun(@(imgnm) fullfile(stimpath, imgnm+suffix),imgnm_vect);
end
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
% psth_all shape: imgN by 200
% This part could be abbrieviated. 
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
end