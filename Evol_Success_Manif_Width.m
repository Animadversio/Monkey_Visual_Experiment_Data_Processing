%% Evolution Successfulness and Manifold Tuning Width
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
Kent_path = "E:\OneDrive - Washington University in St. Louis\Manif_SUHash\summary";
load(fullfile(Kent_path, Animal+"_Manif_Kent_Fit.mat"),"Kent_stats")

%% Evolution
EvolSuccStats = repmat(struct(),1,length(EStats));
for Expi = 1:length(Stats)
blockn = EStats(Expi).evol.block_n;
initpsths = cell2mat(reshape(EStats(Expi).evol.psth(2:3),1,1,[]));
endspsths = cell2mat(reshape(EStats(Expi).evol.psth(blockn-2:blockn-1),1,1,[]));
window = [51:200];
initacts = squeeze(mean(initpsths(1,window,:),2));
endsacts = squeeze(mean(endspsths(1,window,:),2));
[H,P,CI,STATS] = ttest2(endsacts, initacts);
EvolSuccStats(Expi).tstat = STATS.tstat;
EvolSuccStats(Expi).DAOA = (mean(endsacts) - mean(initacts)) / mean(initacts);
end
%% Get the Kappa fit for the Manifold data.
drivechan_tune_kappa = zeros(length(Stats),1);
for Expi = 1:length(Stats)
% ui = find(Stats(Expi).units.pref_chan_id ==
% EStats(Expi).units.pref_chan_id); % this can fail since the numbering
% could change from Evol to Manif
ui = EStats(Expi).units.unit_num_arr(EStats(Expi).units.pref_chan_id);
si = 1; % number of subspace, in case there are multiple subspace there.
kappa = Kent_stats(Expi).act_fit{ui,si}.coef(4);
drivechan_tune_kappa(Expi) = kappa;
end
%%
DAOA_arr = arrayfun(@(E)E.DAOA, EvolSuccStats)';
tstat_arr = arrayfun(@(E)E.tstat, EvolSuccStats)';
corr(DAOA_arr,  drivechan_tune_kappa) 
corr(tstat_arr, drivechan_tune_kappa) 
prefchan_arr = arrayfun(@(E) E.units.pref_chan,EStats);
prefname_arr = arrayfun(@(E) E.units.unit_name_arr(E.units.pref_chan_id),EStats);
%%
V1msk = prefchan_arr <=48 & prefchan_arr>=33;
V4msk = prefchan_arr <=64 & prefchan_arr>=49;
ITmsk = prefchan_arr <=32 & prefchan_arr>=1;
ColorArr = [V1msk',V4msk',ITmsk']; % R for V1, G for V4, B for IT
% ColorArr = 'b'; 
%%
figure(1);
subtightplot(1,2,1,0.05,0.08,0.05)
S1 = scatter(tstat_arr, drivechan_tune_kappa, 25, ColorArr);
S1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp",1:length(Stats));
S1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("prefchan",prefname_arr);
ylabel("Manif Tuning Kappa");xlabel("t stat for Evol")
title(compose("Corr Coef %.3f",corr(tstat_arr,drivechan_tune_kappa)))
subtightplot(1,2,2,0.05,0.10,0.05)
S2 = scatter(DAOA_arr,  drivechan_tune_kappa, 25, ColorArr);
S2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp",1:length(Stats));
S2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("prefchan",prefname_arr);
msk = DAOA_arr < 2 & drivechan_tune_kappa <2 ;
xlabel("Delta Activation / Activation_0")
title(compose("Corr Coef %.3f (< 2 subregion Corr Coef %.3f)",...
    corr(DAOA_arr,drivechan_tune_kappa),corr(DAOA_arr(msk),drivechan_tune_kappa(msk))))
suptitle(compose("%s Correlation between Evol successfulness and Manif tuning width\n(R,G,B for V1 V4 IT)",Animal))
%%
SuccCorr_dir = "E:\OneDrive - Washington University in St. Louis\Evol_Succ_Manif_Kappa";
saveas(2,fullfile(SuccCorr_dir, Animal+"_Evol_Succ_Manif_Kappa_corr.jpg"))
savefig(2,fullfile(SuccCorr_dir, Animal+"_Evol_Succ_Manif_Kappa_corr.fig"))
%%
saveas(2,fullfile(SuccCorr_dir, Animal+"_Evol_Succ_Manif_Kappa_corr_area.jpg"))

%% Evolution Ref Signal to Noise Ratio
for Expi = 1:length(Stats)
refimgmean = cellfun(@(psth)mean(psth(1,window,:),'all'),EStats(Expi).ref.psth_arr);
refimgstd = cellfun(@(psth)std(mean(psth(1,window,:),2)),EStats(Expi).ref.psth_arr);
figure(3);
scatter(refimgmean,refimgstd)
ylabel("std");xlabel("mean")
title(compose("Exp %d pref chan %s", Expi, prefname_arr(Expi)));
pause
end