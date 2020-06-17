%% 
%% Evolution Ref Signal to Noise Ratio
for Expi = 45:length(Stats)
refimgmean = cellfun(@(psth)mean(psth(1,window,:),'all'),EStats(Expi).ref.psth_arr);
refimgstd = cellfun(@(psth)std(mean(psth(1,window,:),2)),EStats(Expi).ref.psth_arr);
figure(3);
scatter(refimgmean,refimgstd)
ylabel("std");xlabel("mean")
title(compose("Exp %d pref chan %s", Expi, prefname_arr(Expi)));
pause
end
%%
ChanQualStats = repmat(struct(),length(Stats),1);
%%
for Expi = 1:length(Stats)
fprintf("Processing %s Evol Exp %d\n",Animal,Expi)
score_vect = cellfun(@(psth)squeeze(mean(psth(1,window,:),2)),EStats(Expi).ref.psth_arr,'UniformOutput',false);
basel_vect = cellfun(@(psth)squeeze(mean(psth(1,1:50,:),2)),EStats(Expi).ref.psth_arr,'UniformOutput',false);
groupsize = cellfun(@(scores) length(scores), score_vect);
idx_vect = arrayfun(@(L, idx) idx*ones(L,1), groupsize, [1:length(score_vect)]', 'UniformOutput', false);
idx_vect = cell2mat(idx_vect);
score_vect = cell2mat(score_vect);
basel_vect = cell2mat(basel_vect);
[P,ANOVATAB,STATS] = anova1(score_vect,idx_vect,'off');
ChanQualStats(Expi).evol.ref.F = ANOVATAB{2,5};
ChanQualStats(Expi).evol.ref.F_P = P;
[H,P,CI,STATS] = ttest2(score_vect,basel_vect);
ChanQualStats(Expi).evol.ref.T = STATS.tstat;
ChanQualStats(Expi).evol.ref.t_P = P;
end

%%
save(fullfile(MatStats_path,Animal+"_ChanQualStat.mat"),'ChanQualStats')
%%
figure;
subtightplot(1,2,1,0.05,0.08,0.05)
S1 = scatter(arrayfun(@(C)C.evol.ref.T,ChanQualStats), drivechan_tune_kappa, 25, ColorArr);
S1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp",1:length(Stats));
S1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("prefchan",prefname_arr);
ylabel("Manif Tuning Kappa");xlabel("t stat for Evol")
title(compose("Corr Coef %.3f",corr(tstat_arr,drivechan_tune_kappa)))
subtightplot(1,2,2,0.05,0.10,0.05)
S2 = scatter(arrayfun(@(C)C.evol.ref.F,ChanQualStats),  drivechan_tune_kappa, 25, ColorArr);
S2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp",1:length(Stats));
S2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("prefchan",prefname_arr);
msk = DAOA_arr < 2 & drivechan_tune_kappa <2 ;
xlabel("Delta Activation / Activation_0")
title(compose("Corr Coef %.3f (< 2 subregion Corr Coef %.3f)",...
    corr(DAOA_arr,drivechan_tune_kappa),corr(DAOA_arr(msk),drivechan_tune_kappa(msk))))
arrayfun(@(C)C.evol.ref.F,ChanQualStats)