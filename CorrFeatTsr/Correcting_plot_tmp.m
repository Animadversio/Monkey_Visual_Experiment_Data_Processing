%%
Animal = "Beto";
predsavedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict";
load(fullfile(predsavedir,Animal+"_FeatTsrPredStats.mat"),'predStats')
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
for Expi = 1:45
prefchan_lab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
score_Evol = predStats(Expi).score_Evol;
score_Manif = predStats(Expi).score_Manif;
%% Redo the E2M plots for the cross prediction 
relufit = predStats(Expi).E2M.nlfunc;
lfitscore_evol = predStats(Expi).E2M.lfitscore;
nlfitscore_evol = predStats(Expi).E2M.nlfitscore;
nlfitcorr = predStats(Expi).E2M.nlfitcorr;
lfitcorr = predStats(Expi).E2M.lfitcorr;
gof = predStats(Expi).E2M.gof;
pred_manif_rsp = predStats(Expi).E2M.lpredscore;
NLpred_manif = predStats(Expi).E2M.nlpredscore;
nlpredcorr = predStats(Expi).E2M.nlpredcorr;
lpredcorr = predStats(Expi).E2M.lpredcorr;
nlpredrsquare = predStats(Expi).E2M.nlpredrsquare;
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