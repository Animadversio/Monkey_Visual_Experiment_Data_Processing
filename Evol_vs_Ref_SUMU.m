% SU MU Ref images vs evolved images
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir,Animal+"_Evol_Stats.mat"), 'EStats');
%%
EStats(Expi).evol.psth
EStats(Expi).ref.psth_arr
%%
Scol = [];
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_Evol_Stats.mat"), 'EStats');
for Expi = 1:numel(EStats)
evol_mean = cellfun(@(P)mean(P(:,51:200,:),[2,3]),EStats(Expi).evol.psth(1:end-1));
ref_mean = cellfun(@(P)mean(P(:,51:200,:),[2,3]),EStats(Expi).ref.psth_arr);
% evol_bsl = cellfun(@(P)mean(P(:,1:50,:),[2,3]),EStats(Expi).evol.psth);
% ref_bsl = cellfun(@(P)mean(P(:,1:50,:),[2,3]),EStats(Expi).ref.psth_arr);
evol_psth_all = cat(3,EStats(Expi).evol.psth{:});
evol_bsl = mean(evol_psth_all(1,1:50,:),[2,3]);
ref_psth_all = cat(3,EStats(Expi).ref.psth_arr{:});
ref_bsl = mean(ref_psth_all(1,1:50,:),[2,3]);

unitN_in_prefchan = sum(EStats(Expi).units.spikeID(EStats(Expi).units.activ_msk) == EStats(Expi).units.pref_chan);
prefunit_i = EStats(Expi).units.unit_num_arr(EStats(Expi).units.pref_chan_id);
isSU = (unitN_in_prefchan > 1) && (prefunit_i==1);

S = struct();
S.Animal = Animal;
S.Expi = Expi;
S.isSU = isSU;
S.prefchan = EStats(Expi).units.pref_chan;
S.area = area_map(S.prefchan);
S.evol_max = max(evol_mean);
S.evol_min = min(evol_mean);
S.evol_bsl = evol_bsl;
S.ref_max = max(ref_mean);
S.ref_min = min(ref_mean);
S.ref_mean = mean(ref_mean);
S.ref_bsl = ref_bsl;
Scol = [Scol,S];
end
end
tab = struct2table(Scol);
disp(tab)
%%
outdir = "E:\OneDrive - Washington University in St. Louis\Evol_Ref_cmp\summary";
writetable(tab,fullfile(outdir,"evol_ref_SUMU_cmp.csv"))
%%

V4SUmsk = (tab.area == "V4") & (tab.isSU);
V4MUmsk = (tab.area == "V4") & (~tab.isSU);
ITSUmsk = (tab.area == "IT") & (tab.isSU);
ITMUmsk = (tab.area == "IT") & (~tab.isSU);
%%
paired_stripe_plot({tab.evol_max, tab.ref_max}, ["Evol", "NatRef"], ...
    {V4SUmsk,ITSUmsk}, ["V4 SU","IT SU"], 'MarkerEdgeAlpha', 0.65, 'LineWidth',2)
ylabel("evoked firing rate")
title("Max Activation for Evolved vs NatRef (SU)")
saveallform(outdir,"maxact_evol-ref_V4IT_SU")
%%
paired_stripe_plot({tab.evol_max, tab.ref_max}, ["Evol", "NatRef"], ...
    {V4MUmsk,ITMUmsk}, ["V4 MU","IT MU"], 'MarkerEdgeAlpha', 0.65, 'LineWidth',2)
ylabel("evoked firing rate")
title("Max Activation for Evolved vs NatRef (MU)")
saveallform(outdir,"maxact_evol-ref_V4IT_MU")
%%
paired_stripe_plot({tab.evol_max, tab.ref_mean}, ["Evol", "NatRef-avg"], ...
    {V4SUmsk,ITSUmsk}, ["V4 SU","IT SU"], 'MarkerEdgeAlpha', 0.65, 'LineWidth',2)
ylabel("evoked firing rate")
title("Max Activation for Evolved vs Average for NatRef (SU)")
saveallform(outdir,"maxact_evol-avg_ref_V4IT_SU")
%%
paired_stripe_plot({tab.evol_max, tab.ref_mean}, ["Evol", "NatRef-avg"], ...
    {V4MUmsk,ITMUmsk}, ["V4 MU","IT MU"], 'MarkerEdgeAlpha', 0.65, 'LineWidth',2)
ylabel("evoked firing rate")
title("Max Activation for Evolved vs Average for NatRef (MU)")
saveallform(outdir,"maxact_evol-avg_ref_V4IT_MU")
%%
mskcol = {V4SUmsk, V4MUmsk, ITSUmsk, ITMUmsk};
msklabs = ["V4-SU", "V4-MU", "IT-SU", "IT-MU"];
for mi = 1:4
msk = mskcol{mi};
fprintf(msklabs(mi)+" paired t test\n")
ttest2_print(tab.evol_max(msk), tab.ref_max(msk), "Evol", "NatRef-max", true);
ttest2_print(tab.evol_max(msk), tab.ref_mean(msk), "Evol", "NatRef-mean", true);
end
%%

% ylabel("normVUS bsl"); ylim([0,2*pi])
% title("SU-MU Non-Parametric Tuning Width Comparison")
% saveallform(nonpardir,"SUMU-normVUS_bsl_cmp_drv_sep")