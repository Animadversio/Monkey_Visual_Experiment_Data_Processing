%
% Analyze Optimization Trajectory in paired BigGAN and FC6 experiment. 
% Code pattern from Evol_Optimizer_CMAGA_cmp.m and Evol_RedDim_TrajsummaryPlot.m
%

mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_BigGAN_FC6_cmp"; 
figdir = fullfile(saveroot,"summary");
Animal = "Alfa"; Set_Path;
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
%%
score_traces = cell(numel(BFEStats),2);
block_traces = cell(numel(BFEStats),2);
wdw = [51:200]; bslwdw = [1:50];
for iTr = 1:numel(BFEStats)
ui = BFEStats(iTr).evol.unit_in_pref_chan(1);
activ_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:),[2])),BFEStats(iTr).evol.psth,"uni",0);
score_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:),[2]))-mean(psth(ui,bslwdw,:),[2,3]), BFEStats(iTr).evol.psth,"uni",0);
block_col = cellfun(@(idx) BFEStats(iTr).evol.block_arr(idx), BFEStats(iTr).evol.idx_seq,"uni",0);
assert(contains(BFEStats(iTr).evol.space_names{1},"fc6"))
assert(contains(BFEStats(iTr).evol.space_names{2},"BigGAN"))
for GANi = 1:2
score_traces{iTr,GANi} = cat(1,score_col{GANi,:});
% zscore_traces{iTr,GANi} = zscores_tsr(pref_chan_id, row_gen & thread_msks{threadi});
block_traces{iTr,GANi} = cat(1,block_col{GANi,:});
end
end
%%
%% Normalize code 
score_C_m_traj_col = [];
score_G_m_traj_col = [];
block_traj_col = [];
for iTr=1:size(score_traces,1)
    [score_C_m,score_C_s,blockvec] = sort_scoreblock(block_traces{iTr,1},score_traces{iTr,1});
    [score_G_m,score_G_s,blockvec] = sort_scoreblock(block_traces{iTr,2},score_traces{iTr,2});
    normmin = 0;%(score_C_m(1)+score_G_m(1))/2;
    normmax = max(score_C_m);%max(cat(2,score_C_m,score_G_m));
    scaling = abs(normmax - normmin);
    score_C_m_norm = (score_C_m-normmin)/scaling;
    score_G_m_norm = (score_G_m-normmin)/scaling;
    score_C_m_traj_col{iTr} = score_C_m_norm;
    score_G_m_traj_col{iTr} = score_G_m_norm;
    block_traj_col{iTr} = blockvec;
end
%%
Corder = colororder();
figure(1);clf;hold on;fignm=Animal+"_FC6MaxNorm_scoreTraj_ColAvg";
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_traj_col{:}),cat(2,score_C_m_traj_col{:}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_traj_col{:}),cat(2,score_G_m_traj_col{:}));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / FC6 Max each session")
title(Animal+compose(" Averaged Optimization Trajectory\n All Session"))
legend(["FC6","BigGAN"])
saveallform(figdir,fignm);
xlim([1,prctile(cellfun(@numel,block_traj_col),85)]);
fignm = fignm+"_Xlim";
saveallform(figdir,fignm);
%%
prefchan_arr = arrayfun(@(S)S.evol.pref_chan(1),BFEStats);
validmsk=true;
V1msk = prefchan_arr <=48 & prefchan_arr>=33;
V4msk = prefchan_arr <=64 & prefchan_arr>=49;
ITmsk = prefchan_arr <=32 & prefchan_arr>=1;
msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
label_col = ["V1", "V4", "IT"];
h=figure(2); fignm=compose("%s_FC6MaxNorm_scoreTraj_avg_Area", Animal);
set(h,'pos',[285         388        1644         549])
T = tiledlayout(1,numel(msk_col),"pad",'compact');
for mski=1:3
msk = msk_col{mski};expidx = find(msk);
if sum(msk)==0,continue;end
nexttile(T,mski)
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_traj_col{expidx}),cat(2,score_C_m_traj_col{expidx}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_traj_col{expidx}),cat(2,score_G_m_traj_col{expidx}));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / FC6 Max each session")
title(compose("Averaged Optimization Trajectory %s (n=%d)",label_col(mski),sum(msk)))
legend(["FC6","BigGAN"])
end
title(T, compose("%s Summary of mean evol trajectory for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for mski=1:3
   msk = msk_col{mski};
   if sum(msk)==0,continue;end
   nexttile(T,mski)
   xlim([1,prctile(cellfun(@numel,block_traj_col(msk)),80)]);
end
saveallform(figdir,fignm+"_Xlim");