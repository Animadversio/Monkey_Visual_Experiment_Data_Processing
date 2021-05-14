%% Create the summary plots for Compare evolution scores for Alfa/Beto.
%  See if the Evolution works
%% 
%% Plot the evolution trajectory comparison
Animal = "Both"; Set_Path;
%%
global figdir
figdir = "E:\OneDrive - Washington University in St. Louis\Evol_ReducDim\summary";
ExpType = "RDEvol";
Animal = "Both";
if strcmp(Animal,"Both") % load stats
A = load(fullfile(MatStats_path, "Alfa"+"_RDEvol_stats.mat"), 'RDStats');
B = load(fullfile(MatStats_path, "Beto"+"_RDEvol_stats.mat"), 'RDStats');
RDStats = [A.RDStats, B.RDStats];
else
load(fullfile(MatStats_path, Animal+"_RDEvol_stats.mat"), 'RDStats')
end

RDEvol_Stats = []; % struct list 
unitnum_arr = zeros(length(RDStats),1);
validmsk = ones(numel(RDStats), 1, 'logical'); % whether the exp should be excluded. (unconventional optimizer setup)
score_traces = cell(numel(RDStats), 2); 
zscore_traces = cell(numel(RDStats), 2); 
block_traces = cell(numel(RDStats), 2); % vector of block num for plotting
% ref_traces = cell(numel(RDStats), 2);
% ref_block_traces = cell(numel(RDStats), 2); % vector of block num for plotting
for Expi = 1:length(RDStats)
thread_n = RDStats(Expi).evol.thread_num;
prefchan_id = find((RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
chid = find((RDStats(Expi).units.unit_num_arr == RDStats(Expi).evol.unit_in_pref_chan(1)) ...
          & (RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
ui = find(prefchan_id==chid); % unit index in the evol.psth array
unitnum_arr(Expi) = RDStats(Expi).evol.unit_in_pref_chan(1);
prefchan_arr(Expi) = RDStats(Expi).evol.pref_chan(1);
% Collect scaler score of each trial and the block num into cell array
window = [51:200]; bslwdw = [1:50];
score_vec_col = cellfun(@(psth)squeeze(mean(psth(ui,window,:),[1,2]) - ...
                               mean(psth(ui,bslwdw,:),[1,2,3])),RDStats(Expi).evol.psth,'uni',0);
zscore_vec_col = zscore_cellarr(score_vec_col);
block_vec_col = cellfun(@(idx)RDStats(Expi).evol.block_arr(idx),...
                              RDStats(Expi).evol.idx_seq,'uni',0);

S = struct();
S.Animal = RDStats(Expi).Animal;
S.Expi = RDStats(Expi).Expi;
S.pref_chan = RDStats(Expi).evol.pref_chan(1);
S.pref_unit = RDStats(Expi).evol.unit_in_pref_chan(1);
% Generate stats on the cell array of scores. 
if all(RDStats(Expi).evol.optim_names == ["ZOHA Sphere lr euclid", "ZOHA Sphere lr euclid ReducDim"])
for thr_i = [1,2] % collect the vectorized scores into the summary stats
score_traces{Expi, thr_i} = cat(1,score_vec_col{thr_i,1:end-1});
zscore_traces{Expi, thr_i} = cat(1,zscore_vec_col{thr_i,1:end-1});
block_traces{Expi, thr_i} = cat(1,block_vec_col{thr_i,1:end-1});
% Collect Stats
end
midgen = round((size(score_vec_col,2) - 1)/2);
S = optimtraj_integral(score_vec_col, S);
S = Dprime_integral(score_vec_col,S);
S = score_cmp_stats(score_vec_col(:,2:3), "init23", S);
S = score_cmp_stats(score_vec_col(:,end-2:end-1), "last23", S);
S = score_cmp_stats(score_vec_col(:,midgen:midgen+1), "middle", S);
else
fprintf("Exp %02d, Skip Non standard optimnames",Expi)
disp(RDStats(Expi).evol.optim_names)
validmsk(Expi) = 0;
S = optimtraj_integral({}, S);
S = Dprime_integral({}, S);
S = score_cmp_stats({}, "init23", S);
S = score_cmp_stats({}, "last23", S);
S = score_cmp_stats({}, "middle", S);
end
RDEvol_Stats = [RDEvol_Stats, S];
end
%% Re-normalize the trajectories to show together.
score_C_m_trajs = {};
score_G_m_trajs = {};
score_C_s_trajs = {};
score_G_s_trajs = {};
block_trajs = {};
for Expi=1:size(score_traces,1)
    [score_C_m,score_C_s,blockvec] = sort_scoreblock(block_traces{Expi,1},score_traces{Expi,1});
    [score_G_m,score_G_s,blockvec] = sort_scoreblock(block_traces{Expi,2},score_traces{Expi,2});
    normmin = 0; % (score_C_m(1)+score_G_m(1))/2;%0; %
    normmax = max(max(score_C_m),max(score_G_m));
    scaling = normmax;%abs(normmax - normmin);
    score_C_m_norm = (score_C_m-normmin)/scaling;
    score_G_m_norm = (score_G_m-normmin)/scaling;
    score_C_m_trajs{Expi} = score_C_m_norm;
    score_G_m_trajs{Expi} = score_G_m_norm;
    score_C_s_trajs{Expi} = score_C_s/scaling;
    score_G_s_trajs{Expi} = score_G_s/scaling;
    block_trajs{Expi} = blockvec;
end
%%
%% Plot the trajectory comparison 
Corder = colororder;
h=figure;hold on;fignm=compose("%s_MaxNorm_scoreTraj_indiv_all", Animal);
set(h,'pos',[1000         315         765         663])
for Expi=1:numel(block_trajs)
    if isempty(block_trajs{Expi}), continue;end
    shadedErrorBar(block_trajs{Expi},score_C_m_trajs{Expi},score_C_s_trajs{Expi},'lineProps',{'Color',[Corder(2,:),0.6],'lineWidth',1},'patchSaturation',0.05)
    shadedErrorBar(block_trajs{Expi},score_G_m_trajs{Expi},score_G_s_trajs{Expi},'lineProps',{'Color',[Corder(1,:),0.6],'lineWidth',1},'patchSaturation',0.05)
end
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("%s All Sessions Evol Trajectory Comparison", Animal),'FontSize',14)
legend(["Full","50D"])
saveallform(figdir,fignm);
xlim([0,35]);fignm = fignm+"_Xlim";
saveallform(figdir,fignm);
%%
h=figure;hold on;fignm=compose("%s_MaxNorm_scoreTraj_avg_all", Animal);
set(h,'pos',[1000         315         765         663])
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{:}),cat(2,score_C_m_trajs{:}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{:}),cat(2,score_G_m_trajs{:}));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("%s All Sessions Averaged Optimization Trajectory",Animal))
legend(["Full","50D"])
saveallform(figdir,fignm);
xlim([0,35]);fignm = fignm+"_Xlim";
saveallform(figdir,fignm);
%%
V1msk = prefchan_arr <=48 & prefchan_arr>=33;
V4msk = prefchan_arr <=64 & prefchan_arr>=49;
ITmsk = prefchan_arr <=32 & prefchan_arr>=1;
msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
label_col = ["V1", "V4", "IT"];
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_avg_Area", Animal);
set(h,'pos',[285         388        1644         549])
T = tiledlayout(1,numel(msk_col),"pad",'compact');
for mski=1:3
nexttile(T,mski)
msk = msk_col{mski};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_C_m_trajs{msk}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_G_m_trajs{msk}));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory %s (n=%d)",label_col(mski),sum(msk)))
legend(["Full","50D"])
end
title(T, compose("%s Summary of mean evol trajectory for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for mski=1:3
   nexttile(T,mski)
   msk = msk_col{mski};
   xlim([1,prctile(cellfun(@numel,block_trajs(msk)),80)]);
end
saveallform(figdir,fignm+"_Xlim");
%%
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_indiv_Area", Animal);
set(h,'pos',[285         388        1644         549])
T = tiledlayout(1,numel(msk_col),"pad",'compact');
for mski=1:3
nexttile(T,mski)
msk = msk_col{mski};
for Expi = find(msk)
    shadedErrorBar(block_trajs{Expi},score_C_m_trajs{Expi},score_C_s_trajs{Expi},'lineProps',{'Color',[Corder(2,:),0.6],'lineWidth',1},'patchSaturation',0.08)
    shadedErrorBar(block_trajs{Expi},score_G_m_trajs{Expi},score_G_s_trajs{Expi},'lineProps',{'Color',[Corder(1,:),0.6],'lineWidth',1},'patchSaturation',0.08)
end
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory %s (n=%d)",label_col(mski),sum(msk)))
legend(["Full","50D"])
end
title(T, compose("%s Collect of mean evol trajectory for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for mski=1:3
   nexttile(T,mski)
   msk = msk_col{mski};
   xlim([1,prctile(cellfun(@numel,block_trajs(msk)),80)]);
end
saveallform(figdir,fignm+"_Xlim");
%% Collect stats 

%% Test on individual Session and Collect T stats on population
RDEvolTab = struct2table(RDEvol_Stats);
writetable(RDEvolTab, fullfile(figdir, Animal+"_RDEvol_trajcmp.csv"))
save(fullfile(figdir, Animal+"_RDEvol_summaryStats.mat"), "RDEvol_Stats")
%% Test
Alfamsk = (RDEvolTab.Animal=="Alfa");
Betomsk = (RDEvolTab.Animal=="Beto");
V1msk = (RDEvolTab.pref_chan<=48 & RDEvolTab.pref_chan>=33);
V4msk = (RDEvolTab.pref_chan>48);
ITmsk = (RDEvolTab.pref_chan<33);
%% Test on the aggregated mean value at population level with a T test. 
testProgression(RDEvolTab, "Dpr_int_norm", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
testProgression(RDEvolTab, "middle_cmp_dpr", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
testProgression(RDEvolTab, "last23_cmp_dpr", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
%%
testProgression(RDEvolTab, "last23_m_ratio", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
testProgression(RDEvolTab, "middle_m_ratio", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], "area", ...
    "Both Monk All Exp");
%%
h = stripe_plot(RDEvolTab, "last23_cmp_dpr", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], ...
                    "Both Monk All Exp", "area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_plot(RDEvolTab, "Dpr_int_norm", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], ...
                    "Both Monk All Exp", "area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_minor_plot(RDEvolTab, "Dpr_int_norm", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], ...
                   {Alfamsk, Betomsk}, ["Alfa", "Beto"], "Both Monk All Exp", "area_anim_sep", {[1,2],[2,3],[1,3]}, 'marker', 'MarkerEdgeAlpha',0.9);
%%
h = stripe_minor_plot(RDEvolTab, "traj_int_ratio", {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk}, ["V1","V4","IT"], ...
                   {Alfamsk, Betomsk}, ["Alfa", "Beto"], "Both Monk All Exp", "area_anim_sep", {[1,2],[2,3],[1,3]}, 'marker', 'MarkerEdgeAlpha',0.9);
%%
% score2cmp = score_vec_col(:,2:3);
% S = score_cmp_stats(score2cmp, "init23");
% S = score_cmp_stats(score_vec_col(:,end-2:end-1), "last23", S);
% Dprime_integral(score_vec_col)

function S = optimtraj_integral(score_traj2cmp, S)
% Integrate the area under the optimization trajectory.
% Return a structure. 
% 
% score_traj2cmp: cell array of scores, 2-by-blocknum. 
%      in each cell is a 1-by-N array of single trial image scores. 
% this function will clip the last generation for stability. 
% No need to clip the score trajectory before entering.
if nargin <= 1, S=struct(); end
if isempty(score_traj2cmp), % use this to generate nan default dict
    S.("traj_int") = nan(1, 2);
    S.("traj_int_norm") = nan(1, 2);
    S.("traj_int_ratio") = nan;
    S.("traj_incr_int") = nan(1, 2);
    S.("traj_incr_int_norm") = nan(1, 2);
    S.("traj_incr_int_ratio") = nan;
    return; 
end
score_m_arr = cellfun(@mean, score_traj2cmp(:,1:end-1));
traj_int = sum(score_m_arr,2)'; % AUC without subtracting 1st gen scores
traj_int_norm = mean(score_m_arr,2)';
traj_incr_int = sum(score_m_arr-score_m_arr(:,1),2)'; % AUC with subtracting 1st gen scores 
traj_incr_int_norm = mean(score_m_arr-score_m_arr(:,1),2)'; % Area for only increased activation.
S.("traj_int") = traj_int;
S.("traj_int_norm") = traj_int_norm;
S.("traj_int_ratio") = traj_int(2) / traj_int(1);
S.("traj_incr_int") = traj_incr_int;
S.("traj_incr_int_norm") = traj_incr_int_norm;
S.("traj_incr_int_ratio") = traj_incr_int(2) / traj_incr_int(1);
end

function S = Dprime_integral(score_traj2cmp, S)
% Integrate D prime between 2 threads. 
% 
% score_traj2cmp: cell array of scores, 2-by-blocknum. 
%      in each cell is a 1-by-N array of single trial image scores. 
% this function will clip the last generation for stability. 
% No need to clip the score trajectory before entering.
if nargin <= 1, S=struct(); end
if isempty(score_traj2cmp), % use this to generate nan default dict
    S.("Dpr_int") = nan;
    S.("Dpr_int_norm") = nan;
    return; 
end
Dpr_arr = [];
for blocki = 1:size(score_traj2cmp,2)-1 
Dpr_arr(blocki) = computeCohen_d(score_traj2cmp{1,blocki}, score_traj2cmp{2,blocki});
end
Dpr_int = sum(Dpr_arr);
Dpr_int_norm = sum(Dpr_arr)/numel(Dpr_arr);
S.("Dpr_int") = Dpr_int;
S.("Dpr_int_norm") = Dpr_int_norm;
end

function S = score_cmp_stats(score2cmp, prefix, S)
% score2cmp: a cell array of scores, 2 rows, each row is a thread / optimizer
% prefix: prefix to name the fields of the struct.
% S: Struct containing the stats 
%
% Example: 
% S = score_cmp_stats(score_vec_col(:,2:3), "init23");
% S = score_cmp_stats(score_vec_col(:,end-2:end-1), "last23");
if nargin == 1, prefix=""; end
if nargin <= 2, S=struct(); end
if isempty(score2cmp), % use this to generate nan default dict
    S.(prefix+"_cmp_t") = nan;
    S.(prefix+"_cmp_p") = nan;
    S.(prefix+"_cmp_dpr") = nan;
    S.(prefix+"_mean") = nan(1,2);
    S.(prefix+"_sem") = nan(1,2);
    S.(prefix+"_m_ratio") = nan;
return; 
end
score_thr1 = cat(1,score2cmp{1,:});
score_thr2 = cat(1,score2cmp{2,:});
[~,P,~,STS] = ttest2(score_thr1, score_thr2);
S.(prefix+"_cmp_t") = STS.tstat;
S.(prefix+"_cmp_p") = P;
S.(prefix+"_cmp_dpr") = computeCohen_d(score_thr1, score_thr2);
S.(prefix+"_mean") = [mean(score_thr1), mean(score_thr2)];
S.(prefix+"_sem") = [sem(score_thr1), sem(score_thr2)];
S.(prefix+"_m_ratio") = mean(score_thr2) / mean(score_thr1);
end

function [zscore_vec_col] = zscore_cellarr(score_vec_col)
score_all_vec = cat(1,score_vec_col{:});
scoreM = mean(score_all_vec);
scoreS = std(score_all_vec);
zscore_vec_col = cellfun(@(vec)(vec-scoreM) / scoreS, score_vec_col, 'uni', 0);
end

function [score_m,score_s,blockvec] = sort_scoreblock(blockarr,scorearr)
% sort an array of scores according to the block array labels. compute the
% mean and std for each block. 
% really useful function to summarize multiple evolution trajectories into
% a mean one. 
blockvec = min(blockarr):max(blockarr);
score_m = [];score_s = [];
for blocki = min(blockarr):max(blockarr)
    score_m(blocki) = mean(scorearr(blockarr==blocki));
    score_s(blocki) = sem(scorearr(blockarr==blocki));
end
end

function saveallform(figdir,fignm,h,sfxlist)
% Save a (current) figure with all suffices in a figdir. 
if nargin <=3, h=gcf; end
if nargin <=4, sfxlist = ["fig","pdf","png"]; end
for sfx = sfxlist
if strcmp(sfx, "fig")
   savefig(h,fullfile(figdir,fignm+"."+sfx))
else
   saveas(h,fullfile(figdir,fignm+"."+sfx))
end
end
end