%% Plot the summary figure for predictability of manifold to Evolution and vice versa.
Animal = "Alfa";layername = "conv4-3";
predsavedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(predsavedir,Animal+"_FeatTsrPredStats.mat"),'predStats')
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
E2M_LinPredCorr = arrayfun(@(S)S.E2M.lpredcorr, predStats);
M2E_LinPredCorr = arrayfun(@(S)S.M2E.lpredcorr, predStats);
E2M_LNLPredCorr = arrayfun(@(S)S.E2M.nlpredcorr, predStats);
M2E_LNLPredCorr = arrayfun(@(S)S.M2E.nlpredcorr, predStats);
%% Pooled Area Statistics
figure(3);clf;hold on
ExpN = length(E2M_LinPredCorr);
scatter(ones(1,ExpN)+0.1*randn(1,ExpN),E2M_LinPredCorr);
scatter(1.5*ones(1,ExpN)+0.1*randn(1,ExpN),E2M_LNLPredCorr);
scatter(3*ones(1,ExpN)+0.1*randn(1,ExpN),M2E_LinPredCorr);
scatter(3.5*ones(1,ExpN)+0.1*randn(1,ExpN),M2E_LNLPredCorr);
legend(["E2M,LinearPred","E2M,LNPred","M2E,LinearPred","M2E,LNLPred"])
xticks([1.25,3.25]);xticklabels(["Evol->Manif","Manif->Evol"]);xlabel("Prediction Direction")
ylabel("Pearson Correlation")
title(compose("%s Cross Experiment Prediction Using Linear and L-N model \nwith Correlation Coefficient in VGG conv4-3",Animal))
%%
saveas(3,fullfile(predsavedir,compose("%s_conv4-3_pred_summary.png",Animal)))
savefig(3,fullfile(predsavedir,compose("%s_conv4-3_pred_summary.fig",Animal)))
%%
prefchan_arr = arrayfun(@(S)S.units.pref_chan,EStats);
V1msk = prefchan_arr <=48 & prefchan_arr>=33;
V4msk = prefchan_arr <=64 & prefchan_arr>=49;
ITmsk = prefchan_arr <=32 & prefchan_arr>=1;
%% Separate Area Statistics
figure(1);clf;hold on
scatter(1*ones(1,sum(V1msk))+0.1*randn(1,sum(V1msk)),E2M_LinPredCorr(V1msk));
scatter(1.5*ones(1,sum(V4msk))+0.1*randn(1,sum(V4msk)),E2M_LinPredCorr(V4msk));
scatter(2*ones(1,sum(ITmsk))+0.1*randn(1,sum(ITmsk)),E2M_LinPredCorr(ITmsk));
% scatter(3*ones(1,45)+0.1*randn(1,45),M2E_LinPredCorr);
% scatter(3.5*ones(1,45)+0.1*randn(1,45),M2E_LNLPredCorr);
scatter(3*ones(1,sum(V1msk))+0.1*randn(1,sum(V1msk)),M2E_LinPredCorr(V1msk));
scatter(3.5*ones(1,sum(V4msk))+0.1*randn(1,sum(V4msk)),M2E_LinPredCorr(V4msk));
scatter(4*ones(1,sum(ITmsk))+0.1*randn(1,sum(ITmsk)),M2E_LinPredCorr(ITmsk));

legend(["E2M,V1","E2M,V4","E2M,IT","M2E,V1","M2E,V4","M2E,IT"])
xticks([1,1.5,2,3,3.5,4]);xticklabels(["E2M,V1","E2M,V4","E2M,IT","M2E,V1","M2E,V4","M2E,IT"])
% xticks([1.5,3.5]);xticklabels(["Evol->Manif","Manif->Evol"])
ylabel("Pearson Correlation")
title(compose("%s Cross Experiment Prediction (each area) Using Linear model\nwith Correlation Coefficient in VGG conv4-3",Animal))
%%
saveas(1,fullfile(predsavedir,compose("%s_conv4-3_pred_area_summary.png",Animal)))
savefig(1,fullfile(predsavedir,compose("%s_conv4-3_pred_area_summary.fig",Animal)))
%%




