

figdir = "O:\Evol_ReducDim\summary";
outdir = "O:\Manuscript_Manifold\Figure3\RedDimEffectProg";
prefchan_arr = arrayfun(@(S)S.evol.pref_chan(1),RDStats');
area_arr = arrayfun(@area_map,prefchan_arr);
anim_arr = [RDStats.Animal]';
%% creat masks to filter array
V1msk = area_arr=="V1";
V4msk = area_arr=="V4";
ITmsk = area_arr=="IT";
Alfamsk = anim_arr=="Alfa";
Betomsk = anim_arr=="Beto";
anysucsmsk = any(RDEvolTab.t_p_succ<0.01,2);
allsucsmsk = all(RDEvolTab.t_p_succ<0.01,2);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Older more verbose API Trajectory Summary plot: 
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
%% Plot trajectory comparison for each individual session.
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

%% Plot trajectory comparison averaging all sessions.
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

%% Plot trajectory comparison averaging sessions with pref chan in V1, V4, IT
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

%% Plot trajectory comparison averaging sessions with pref chan in V1, V4, IT
msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
label_col = ["V1", "V4", "IT"];
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_avg_Area_movmean", Animal);
set(h,'pos',[285         388        1644         549])
T = tiledlayout(1,numel(msk_col),"pad",'compact');
for mski=1:3
nexttile(T,mski)
msk = msk_col{mski};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_C_m_trajs{msk}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_G_m_trajs{msk}));
score_C_col_mov_m = movmean(score_C_col_m,3);
score_G_col_mov_m = movmean(score_G_col_m,3);
score_C_col_mov_s = movmean(score_C_col_s,3);
score_G_col_mov_s = movmean(score_G_col_s,3);
shadedErrorBar(block_colvec,score_C_col_mov_m,score_C_col_mov_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_mov_m,score_G_col_mov_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory %s (n=%d)",label_col(mski),sum(msk)))
legend(["Full","50D"])
end
title(T, compose("%s Summary of mean evol trajectory for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for mski=1:3 % relimit x axis by 80 percentile of block number
   nexttile(T,mski)
   msk = msk_col{mski};
   xlim([1,prctile(cellfun(@numel,block_trajs(msk)),80)]);
end
saveallform(figdir,fignm+"_Xlim");

%% Plot trajectory comparison averaging sessions with pref chan in V1, V4, IT in A B (Final version)
msk = validmsk&anysucsmsk;
msk_col = {V1msk&msk, V4msk&msk, ITmsk&msk};
anim_msks = {Alfamsk&msk, Betomsk&msk};
label_col = ["V1", "V4", "IT"];
anim_col = ["Alfa","Beto"];
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_avg_Area_Anim_movmean_anysucs", Animal);
set(h,'pos',[125   258   935   715])
T = tiledlayout(2,numel(msk_col),"pad",'compact',"tilespac",'compact');
for animi=1:2
for mski=1:3
nexttile(T,mski+3*(animi-1))
msk = msk_col{mski} & anim_msks{animi};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_C_m_trajs{msk}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_G_m_trajs{msk}));
score_C_col_mov_m = movmean(score_C_col_m,3);
score_G_col_mov_m = movmean(score_G_col_m,3);
score_C_col_mov_s = movmean(score_C_col_s,3);
score_G_col_mov_s = movmean(score_G_col_s,3);
shadedErrorBar(block_colvec,score_C_col_mov_m,score_C_col_mov_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_mov_m,score_G_col_mov_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory\n%s %s (n=%d)",anim_col(animi),label_col(mski),sum(msk)))
legend(["Full","50D"])
end
end
title(T, compose("%s Summary of mean evol trajectory (3 block movmean) for each area",Animal), 'FontSize',16)
saveallform(figdir,fignm);
for animi=1:2
for mski=1:3
    nexttile(T,mski+3*(animi-1))
    msk = msk_col{mski} & anim_msks{animi};
    xlim([1,prctile(cellfun(@numel,block_trajs(msk)),80)]);
end
end
saveallform(figdir,fignm+"_Xlim");

%% Area average with individual curves plotted on it.
% msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
msk = validmsk&anysucsmsk;
msk_col = {V1msk&msk, V4msk&msk, ITmsk&msk};
label_col = ["V1", "V4", "IT"];
h=figure; fignm=compose("%s_MaxNorm_scoreTraj_avg_ws_indiv_Area_anysucs", Animal);
set(h,'pos',[285         388        1644         549])
T = tiledlayout(1,numel(msk_col),"pad",'compact');
for mski=1:3
nexttile(T,mski);hold on
msk = msk_col{mski};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_C_m_trajs{msk}));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_trajs{msk}),cat(2,score_G_m_trajs{msk}));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:),'LineWidth',1.5})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:),'LineWidth',1.5})
for Expi = find(msk)' % plot individual experiments' optim traj in the mask. 
    plot(block_trajs{Expi},movmean(score_C_m_trajs{Expi},3),'Color',[Corder(2,:),0.4],'lineWidth',0.5)
    plot(block_trajs{Expi},movmean(score_G_m_trajs{Expi},3),'Color',[Corder(1,:),0.4],'lineWidth',0.5)
end
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


%%
% score2cmp = score_vec_col(:,2:3);
% S = score_cmp_stats(score2cmp, "init23");
% S = score_cmp_stats(score_vec_col(:,end-2:end-1), "last23", S);
% Dprime_integral(score_vec_col)
