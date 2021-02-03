Animal = "Alfa";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path,Animal+"_FeatTsrPredStats2.mat"),'predStats2')
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
prefchan_arr = arrayfun(@(S)S.units.pref_chan,Stats);
ITmsk = prefchan_arr <=32;
V1msk = (prefchan_arr < 49) & (prefchan_arr>32);
V4msk = (prefchan_arr > 48);
Allmsk = ~isnan(prefchan_arr);
area_labs = ["V1","V4","IT","all"];area_msks = {V1msk, V4msk, ITmsk, Allmsk};
%%
diary(fullfile("O:\corrFeatTsr_Predict\summary\",compose("%s_E2M_PredStats.log",Animal)))
fprintf("Fitting on Evolution Deploy on Others:\n")
fprintf("Nonlinear fitting correlaiton:\n")
corr_summary_by_msk(predStats2,area_msks,area_labs,"evol","nlfitcorr","E2M")
fprintf("Prediction correlaiton by the nonlinear function:\n")
corr_summary_by_msk(predStats2,area_msks,area_labs,"manif","nlpredcorr","E2M")
corr_summary_by_msk(predStats2,area_msks,area_labs,"evoref","nlpredcorr","E2M")
corr_summary_by_msk(predStats2,area_msks,area_labs,"gabor","nlpredcorr","E2M")
corr_summary_by_msk(predStats2,area_msks,area_labs,"pasu","nlpredcorr","E2M")
diary off
%%
diary(fullfile("O:\corrFeatTsr_Predict\summary\",compose("%s_M2E_PredStats.log",Animal)))
fprintf("Fitting on Manifold Deploy on Others:\n")
fprintf("Nonlinear fitting correlaiton:\n")
corr_summary_by_msk(predStats2,area_msks,area_labs,"manif","nlfitcorr")
fprintf("Prediction correlaiton by the nonlinear function:\n")
corr_summary_by_msk(predStats2,area_msks,area_labs,"evol","nlpredcorr")
corr_summary_by_msk(predStats2,area_msks,area_labs,"evoref","nlpredcorr")
corr_summary_by_msk(predStats2,area_msks,area_labs,"gabor","nlpredcorr")
corr_summary_by_msk(predStats2,area_msks,area_labs,"pasu","nlpredcorr")
diary off
%%
nlpred_vec = ([predStats2(Expi).M2E.(layername).pasu.nlpredscore(:,end);predStats2(Expi).M2E.(layername).gabor.nlpredscore(:,end);predStats2(Expi).M2E.(layername).evoref.nlpredscore(:,end);]);
orig_vec = ([predStats2(Expi).M2E.(layername).pasu.origscore(:,end);predStats2(Expi).M2E.(layername).gabor.origscore(:,end);predStats2(Expi).M2E.(layername).evoref.origscore(:,end);]);
[cval,pval] = corr(orig_vec,nlpred_vec)
%% corr_summary_by_msk: function description
function [] = corr_summary_by_msk(predStats2,area_msks,area_labs,imagespace,corrvarnm,direction,areafirst)
if nargin <= 4, corrvarnm = "nlpredcorr"; end
if nargin <= 5, direction = "M2E"; end
if nargin <= 6, areafirst = true; end
fprintf("\nIn %s image space,\n",imagespace)
if areafirst
for L = 1:numel(area_labs)
    fprintf("Area %s\n",area_labs(L));
    msk = area_msks{L};
    for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
    fprintf("Layer %s ", layername)
    predcorr_P = arrayfun(@(Expi)single(predStats2(Expi).(direction).(layername).(imagespace).(corrvarnm+"_P")),1:numel(predStats2));
    predcorr = arrayfun(@(Expi)single(predStats2(Expi).(direction).(layername).(imagespace).(corrvarnm)),1:numel(predStats2));
    fprintf("Area %s Corr %.3f (%.1e) Prctile [%.3f %.3f] N=%d\n",area_labs(L),nanmedian(predcorr(msk)),nanmedian(predcorr_P(msk)),...
            prctile(predcorr(msk),20),prctile(predcorr(msk),80),sum(msk));
    end
end
else
for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
    predcorr_P = arrayfun(@(Expi)single(predStats2(Expi).(direction).(layername).(imagespace).(corrvarnm+"_P")),1:numel(predStats2));
    predcorr = arrayfun(@(Expi)single(predStats2(Expi).(direction).(layername).(imagespace).(corrvarnm)),1:numel(predStats2));
    fprintf("Layer %s Corr %.3f (%.1e) Prctile [%.3f %.3f]\n",layername,nanmedian(predcorr),nanmedian(predcorr_P),...
        prctile(predcorr,20),prctile(predcorr,80));
    for L = 1:numel(area_labs)
        msk = area_msks{L};
        fprintf("Area %s Corr %.3f (%.1e) Prctile [%.3f %.3f] N=%d\n",area_labs(L),nanmedian(predcorr(msk)),nanmedian(predcorr_P(msk)),...
            prctile(predcorr(msk),20),prctile(predcorr(msk),80),sum(msk));
    end
end
end 
end

% 
% fprintf("Natural image reference to Evolution\n")
% for L = 1:numel(area_labs)
%     fprintf("Area %s\n",area_labs(L));
%     msk = area_msks{L};
%     for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
%     fprintf("Layer %s ", layername)
%     evoref_predcorr_P = arrayfun(@(Expi)predStats2(Expi).M2E.(layername).evoref.nlpredcorr_P,1:numel(Stats));
%     evoref_predcorr = arrayfun(@(Expi)predStats2(Expi).M2E.(layername).evoref.nlpredcorr,1:numel(Stats));
%     fprintf("Area %s Corr %.3f (%.1e) Prctile [%.3f %.3f] N=%d\n",area_labs(L),nanmedian(evoref_predcorr(msk)),nanmedian(evoref_predcorr_P(msk)),...
%             prctile(evoref_predcorr(msk),20),prctile(evoref_predcorr(msk),80),sum(msk));
%     end
% end
% 
% 
% %%
% fprintf("Natural image reference to Evolution\n")
% for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
%     evoref_predcorr_P = arrayfun(@(Expi)predStats2(Expi).M2E.(layername).evoref.nlpredcorr_P,1:numel(Stats));
%     evoref_predcorr = arrayfun(@(Expi)predStats2(Expi).M2E.(layername).evoref.nlpredcorr,1:numel(Stats));
%     fprintf("Layer %s Corr %.3f (%.1e) Prctile [%.3f %.3f]\n",layername,nanmedian(evoref_predcorr),nanmedian(evoref_predcorr_P),...
%         prctile(evoref_predcorr,20),prctile(evoref_predcorr,80));
%     for L = 1:numel(area_labs)
%         msk = area_msks{L};
%         fprintf("Area %s Corr %.3f (%.1e) Prctile [%.3f %.3f] N=%d\n",area_labs(L),nanmedian(evoref_predcorr(msk)),nanmedian(evoref_predcorr_P(msk)),...
%             prctile(evoref_predcorr(msk),20),prctile(evoref_predcorr(msk),80),sum(msk));
%     end
% end
% fprintf("Generated images in evolution exp\n")
% for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
%     evol_predcorr_P = arrayfun(@(Expi)predStats2(Expi).M2E.(layername).evol.nlpredcorr_P,1:numel(Stats));
%     evol_predcorr = arrayfun(@(Expi)predStats2(Expi).M2E.(layername).evol.nlpredcorr,1:numel(Stats));
%     fprintf("Layer %s Corr %.3f (%.1e) Prctile [%.3f %.3f]\n",layername,nanmedian(evol_predcorr),nanmedian(evol_predcorr_P),...
%         prctile(evol_predcorr,20),prctile(evol_predcorr,80));
%     for L = 1:numel(area_labs)
%         msk = area_msks{L};
%         fprintf("Area %s Corr %.3f (%.1e) Prctile [%.3f %.3f] N=%d\n",area_labs(L),nanmedian(evol_predcorr(msk)),nanmedian(evol_predcorr_P(msk)),...
%             prctile(evol_predcorr(msk),20),prctile(evol_predcorr(msk),80),sum(msk));
%     end
% end
% fprintf("Manifold images \n")
% for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
%     manif_fitcorr_P = arrayfun(@(Expi)predStats2(Expi).M2E.(layername).manif.nlfitcorr_P,1:numel(Stats));
%     manif_fitcorr = arrayfun(@(Expi)predStats2(Expi).M2E.(layername).manif.nlfitcorr,1:numel(Stats));
%     fprintf("Layer %s Corr %.3f (%.1e) Prctile [%.3f %.3f]\n",layername,nanmedian(manif_fitcorr),nanmedian(manif_fitcorr_P),...
%         prctile(manif_fitcorr,20),prctile(manif_fitcorr,80));
%     for L = 1:numel(area_labs)
%         msk = area_msks{L};
%         fprintf("Area %s Corr %.3f (%.1e) Prctile [%.3f %.3f] N=%d\n",area_labs(L),nanmedian(manif_fitcorr(msk)),nanmedian(manif_fitcorr_P(msk)),...
%             prctile(manif_fitcorr(msk),20),prctile(manif_fitcorr(msk),80),sum(msk));
%     end
% end
% %%
% fprintf("Gabor patches \n")
% for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
%     gabor_predcorr_P = arrayfun(@(Expi)single(predStats2(Expi).M2E.(layername).gabor.nlpredcorr_P),1:numel(Stats));
%     gabor_predcorr = arrayfun(@(Expi)single(predStats2(Expi).M2E.(layername).gabor.nlpredcorr),1:numel(Stats));
%     fprintf("Layer %s Corr %.3f (%.1e) Prctile [%.3f %.3f]\n",layername,nanmedian(gabor_predcorr),nanmedian(gabor_predcorr_P),...
%         prctile(gabor_predcorr,20),prctile(gabor_predcorr,80));
%     for L = 1:numel(area_labs)
%         msk = area_msks{L};
%         fprintf("Area %s Corr %.3f (%.1e) Prctile [%.3f %.3f] N=%d\n",area_labs(L),nanmedian(gabor_predcorr(msk)),nanmedian(gabor_predcorr_P(msk)),...
%             prctile(gabor_predcorr(msk),20),prctile(gabor_predcorr(msk),80),sum(msk));
%     end
% end
% fprintf("Pasupathy patches \n")
% for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
%     pasu_predcorr_P = arrayfun(@(Expi)single(predStats2(Expi).M2E.(layername).pasu.nlpredcorr_P),1:numel(Stats));
%     pasu_predcorr = arrayfun(@(Expi)single(predStats2(Expi).M2E.(layername).pasu.nlpredcorr),1:numel(Stats));
%     fprintf("Layer %s Corr %.3f (%.1e) Prctile [%.3f %.3f]\n",layername,nanmedian(pasu_predcorr),nanmedian(pasu_predcorr_P),...
%         prctile(pasu_predcorr,20),prctile(pasu_predcorr,80));
%     for L = 1:numel(area_labs)
%         msk = area_msks{L};
%         fprintf("Area %s Corr %.3f (%.1e) Prctile [%.3f %.3f] N=%d\n",area_labs(L),nanmedian(pasu_predcorr(msk)),nanmedian(pasu_predcorr_P(msk)),...
%             prctile(pasu_predcorr(msk),20),prctile(pasu_predcorr(msk),80),sum(msk));
%     end
% end