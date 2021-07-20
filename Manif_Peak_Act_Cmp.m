%% Manif Peak Comp. Compare the peak activation in the few different spaces.
% This is for the preferred channel. 
Animal = "Both"; Set_Path;
mat_dir = "O:\Mat_Statistics";
%%
pasu_idx_vec = reshape(Stats(12).ref.pasu_idx_grid',[],1); % reshape into one row. But transpose to make similar object adjacent
pasu_val_msk = ~cellfun(@isempty,pasu_idx_vec);
assert(sum(pasu_val_msk)==186)
%%
figdir = "O:\Manif_Proto_Cmp";
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
ActStats = repmat(struct(), 1, numel(Stats));
for Expi = 1:numel(Stats)
fprintf("Processing Exp %d....\n",Expi)
ActStats(Expi).Animal = Animal;
ActStats(Expi).Expi = Expi;

% manifdir = Stats(Expi).meta.stimuli;
% evodir = EStats(Expi).meta.stimuli;
% evorefdir = fileparts(EStats(Expi).meta.stimuli);
% gabordir = "N:\Stimuli\2019-Manifold\gabor"; % N:\Stimuli\2019-Parameterized-Shapes\Gabors
% pasudir = "N:\Stimuli\2019-Manifold\pasupathy-wg-f-4-ori"; % "N:\Stimuli\2019-Parameterized-Shapes\Pasupathy"
% imgcol = {};
% imgnmcol = {};
% scorecol = [];
% semcol = [];
% stdcol = [];
% titstr_col = [];

bsl_VEC_ALL = [];
act_VEC_ALL = []; % Get activation of all the trials
% Manifold images. 
si=1; ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_sem_vec = cellfun(@(psth)sem(squeeze(mean(psth(ui,51:200,:),[1,2]))),reshape(Stats(Expi).manif.psth{si},[],1));
score_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).manif.psth{si},[],1),'uni',0));
score_act_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,50:200,:),[2])),reshape(Stats(Expi).manif.psth{si},[],1),'uni',0)); % single trial vector
bsl_VEC_ALL = [bsl_VEC_ALL; score_bsl_VEC];
act_VEC_ALL = [act_VEC_ALL; score_act_VEC];
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
ActStats(Expi).manif.maxScore = maxScore;
ActStats(Expi).manif.mean_vec = score_vec;
ActStats(Expi).manif.std_vec = score_std_vec;
ActStats(Expi).manif.sem_vec = score_sem_vec;

% if flag.doEvoRef
evoref_vec = cellfun(@(psth)mean(psth(1,51:200,:),'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_std_vec = cellfun(@(psth)std(mean(psth(1,51:200,:),[1,2]),1,'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_sem_vec = cellfun(@(psth)sem(squeeze(mean(psth(1,51:200,:),[1,2]))),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(1,1:45,:),[2])),reshape(EStats(Expi).ref.psth_arr,[],1),'uni',0));
evoref_act_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(1,50:200,:),[2])),reshape(EStats(Expi).ref.psth_arr,[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; evoref_bsl_VEC];
act_VEC_ALL = [act_VEC_ALL; evoref_act_VEC];
[evoref_sortScore,sortId] = sort(evoref_vec,'Descend');
[evoref_maxScore,evoref_maxId] = max(evoref_vec);
ActStats(Expi).evoref.maxScore = evoref_maxScore;
ActStats(Expi).evoref.mean_vec = evoref_vec;
ActStats(Expi).evoref.std_vec = evoref_std_vec;
ActStats(Expi).evoref.sem_vec = evoref_sem_vec;
evorefnm = EStats(Expi).ref.imgnm;
% end

% Pasupathy patches
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
nanmsk = ~pasu_val_msk|isnan(pasu_vec);
pasu_vec(nanmsk) = []; % isnan(pasu_vec) % get rid of non-existing pasupathy images. 
pasu_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
% pasu_sem_vec = cellfun(@(psth)sem(squeeze(mean(psth(ui,51:200,:),[1,2]))),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_sem_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all')/sqrt(size(psth,3)),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_std_vec(nanmsk) = [];
pasu_sem_vec(nanmsk) = [];
pasu_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).ref.pasu_psths',[],1),'uni',0));
pasu_act_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,50:200,:),[2])),reshape(Stats(Expi).ref.pasu_psths',[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; pasu_bsl_VEC];
act_VEC_ALL = [act_VEC_ALL; pasu_act_VEC];
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
ActStats(Expi).pasu.maxScore = pasu_maxScore;
ActStats(Expi).pasu.mean_vec = pasu_vec;
ActStats(Expi).pasu.std_vec = pasu_std_vec;
ActStats(Expi).pasu.sem_vec = pasu_sem_vec;
pasuimgnm = cellfun(@(idx)string(unique(Stats(Expi).imageName(idx))),reshape(Stats(Expi).ref.pasu_idx_grid',[],1),'uni',0);
pasuimgnm(nanmsk) = [];
pasuimgnm = string(pasuimgnm);
else
ActStats(Expi).pasu.maxScore = [];
ActStats(Expi).pasu.mean_vec = [];
ActStats(Expi).pasu.std_vec = [];
ActStats(Expi).pasu.sem_vec = [];
end
% Gabor patches
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
nanmsk = isnan(gab_vec);
gab_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
% gab_sem_vec = cellfun(@(psth)sem(squeeze(mean(psth(ui,51:200,:),[1,2]))),reshape(Stats(Expi).ref.gab_psths,[],1));
gab_sem_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all')/sqrt(size(psth,3)),reshape(Stats(Expi).ref.gab_psths,[],1));
gab_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).ref.gab_psths,[],1),'uni',0));
gab_act_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,50:200,:),[2])),reshape(Stats(Expi).ref.gab_psths,[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; gab_bsl_VEC];
act_VEC_ALL = [act_VEC_ALL; gab_act_VEC];
gaborimgnm = cellfun(@(idx)unique(Stats(Expi).imageName(idx)),reshape(Stats(Expi).ref.gab_idx_grid,[],1),'Uni',0);
gaborimgnm(nanmsk) = [];
gaborimgnm = string(gaborimgnm); 
gab_vec(nanmsk) = [];
gab_std_vec(nanmsk) = [];
gab_sem_vec(nanmsk) = [];
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
ActStats(Expi).gab.maxScore = gab_maxScore;
ActStats(Expi).gab.mean_vec = gab_vec;
ActStats(Expi).gab.std_vec = gab_std_vec;
ActStats(Expi).gab.sem_vec = gab_sem_vec;
else
ActStats(Expi).gab.maxScore = [];
ActStats(Expi).gab.mean_vec = [];
ActStats(Expi).gab.std_vec = [];
ActStats(Expi).gab.sem_vec = [];
end
ActStats(Expi).bsl_rate_mean = nanmean(bsl_VEC_ALL);
ActStats(Expi).bsl_rate_std = nanstd(bsl_VEC_ALL,1);
ActStats(Expi).bsl_rate_sem = sem(bsl_VEC_ALL,1);
ActStats(Expi).act_rate_mean = nanmean(act_VEC_ALL);
ActStats(Expi).act_rate_std = nanstd(act_VEC_ALL,1);
ActStats(Expi).act_rate_sem = sem(act_VEC_ALL,1);
for imspace = ["manif", "evoref", "pasu", "gab"]
ActStats(Expi).(imspace).mean_vec_Z = (ActStats(Expi).(imspace).mean_vec - ActStats(Expi).act_rate_mean) / ActStats(Expi).act_rate_std;
ActStats(Expi).(imspace).maxScore_Z = (ActStats(Expi).(imspace).maxScore - ActStats(Expi).act_rate_mean) / ActStats(Expi).act_rate_std;
end
% savenm = compose("%s_Exp%02d_pref%02d_prototypes",Animal,Expi,Stats(Expi).units.pref_chan);
end
save(fullfile(mat_dir, Animal+"_PeakActCmp.mat"), "ActStats")
end
%%
peaksumtab = table(); csr = 1;
for Animal = ["Alfa", "Beto"]
	load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
 	load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
 	load(fullfile(mat_dir, Animal+"_PeakActCmp.mat"), "ActStats")
 	for Expi = 1:numel(ActStats)
 		for imspace = ["manif", "evoref", "pasu", "gab"]
      if isempty(ActStats(Expi).(imspace).maxScore)
      peaksumtab.(imspace+"maxActZ")(csr) = nan;
 		peaksumtab.(imspace+"maxAct")(csr) = nan;
      else
 		peaksumtab.(imspace+"maxActZ")(csr) = ActStats(Expi).(imspace).maxScore_Z;
 		peaksumtab.(imspace+"maxAct")(csr) = ActStats(Expi).(imspace).maxScore;
      end
 		end
 		for varnm = ["bsl_rate_mean", "bsl_rate_std", "bsl_rate_sem", "act_rate_mean", "act_rate_std", "act_rate_sem"]
		peaksumtab.(varnm)(csr) = ActStats(Expi).(varnm);
 		end
 		csr = csr + 1;
 	end
end
%%
sumdir = "O:\Manif_NonParam\summary";
writetable(peaksumtab,fullfile(sumdir,"Both_PeakActStats.csv"))
%%
sumdir = "O:\Manif_NonParam\summary";
peaksumtab = readtable(fullfile(sumdir,"Both_PeakActStats.csv"));
pair_ttest_print(peaksumtab.manifmaxAct, peaksumtab.evorefmaxAct, "Manif", "EvoRef", "Max Activation (Sp/s)");
pair_ttest_print(peaksumtab.manifmaxActZ, peaksumtab.evorefmaxActZ, "Manif", "EvoRef", "Max Activation Zscored");
pair_ttest_print(peaksumtab.manifmaxAct, peaksumtab.pasumaxAct, "Manif", "Pasupathy", "Max Activation (Sp/s)");
pair_ttest_print(peaksumtab.manifmaxActZ, peaksumtab.pasumaxActZ, "Manif", "Pasupathy", "Max Activation Zscored");
pair_ttest_print(peaksumtab.manifmaxAct, peaksumtab.gabmaxAct, "Manif", "Gabor", "Max Activation (Sp/s)");
pair_ttest_print(peaksumtab.manifmaxActZ, peaksumtab.gabmaxActZ, "Manif", "Gabor", "Max Activation Zscored");
% [~,P,CI,TST] = ttest(peaksumtab.manifmaxAct, peaksumtab.evorefmaxAct);
% fprintf("Max Activation: Manif - EvoRef: %.3f(%.3f) - %.3f(%.3f). t=%.3f p=%.1e, df=%d")
function pair_ttest_print(vec1, vec2, var1lab, var2lab, cmpvar)
msk = ~(isnan(vec1) | isnan(vec2));
mean1 = mean(vec1(msk));
sem1 = sem(vec1(msk)); 
mean2 = mean(vec2(msk));
sem2 = sem(vec2(msk)); 
[~,P,CI,TST] = ttest(vec1, vec2);
fprintf("%s: %s - %s: %.3f(%.3f) - %.3f(%.3f). t=%.3f p=%.1e, df=%d\n",...
	cmpvar, var1lab, var2lab, mean1, sem1, mean2, sem2, TST.tstat, P, TST.df)
end

