%
% Analyze Optimization Trajectory in paired BigGAN and FC6 experiment. 
% Code pattern from Evol_Optimizer_CMAGA_cmp.m and Evol_RedDim_TrajsummaryPlot.m
% Updated to modern pattern @ July 12th. 
mat_dir = "O:\Mat_Statistics"; 
saveroot = "O:\Evol_BigGAN_FC6_cmp"; 
outdir = "O:\ThesisProposal\BigGAN";
figdir = fullfile(saveroot,"summary");
%
Animal = "Both"; Set_Path;
if Animal == "Both"
A = load(fullfile(mat_dir, "Alfa" + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats');
B = load(fullfile(mat_dir, "Beto" + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats');
BFEStats = [A.BFEStats;B.BFEStats];
clear A B
else
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
end
%% Collect the score each block 
score_traces = cell(numel(BFEStats),2);
block_traces = cell(numel(BFEStats),2);
wdw = [51:200]; bslwdw = [1:50];
for iTr = 1:numel(BFEStats)
% ui = BFEStats(iTr).evol.unit_in_pref_chan(1);
ui=1;
activ_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:),[2])),BFEStats(iTr).evol.psth,"uni",0); % cell arr, vector of score per gen. 
score_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:),[2]))-mean(psth(ui,bslwdw,:),[2,3]), BFEStats(iTr).evol.psth,"uni",0);
block_col = cellfun(@(idx) BFEStats(iTr).evol.block_arr(idx), BFEStats(iTr).evol.idx_seq,"uni",0); % cell arr, vector of block id per gen. 
assert(contains(BFEStats(iTr).evol.space_names{1},"fc6"))
assert(contains(BFEStats(iTr).evol.space_names{2},"BigGAN"))
for GANi = 1:2
score_traces{iTr,GANi} = cat(1,score_col{GANi,1:end-1});
% zscore_traces{iTr,GANi} = zscores_tsr(pref_chan_id, row_gen & thread_msks{threadi});
block_traces{iTr,GANi} = cat(1,block_col{GANi,1:end-1});
end
end
%% Newer API verions
[score_m_traj_col, block_traj_col, score_m_traj_extrap_col, block_traj_extrap_col, ...
   extrap_mask_col, sucsmsk_end, sucsmsk_max, tval_end_arr, pval_end_arr, tval_max_arr, pval_max_arr] =...
   optim_traj_process(block_traces, score_traces, ["FC6", "BigGAN"], "max1", 35);

%%
outdir = "O:\ThesisProposal\BigGAN";
prefchan_arr = arrayfun(@(S)S.evol.pref_chan(1),BFEStats);
area_arr = arrayfun(@area_map,prefchan_arr);
anim_arr = arrayfun(@(S)S.Animal,BFEStats);
V4msk = area_arr=="V4";
ITmsk = area_arr=="IT";
Amsk = anim_arr=="Alfa";
Bmsk = anim_arr=="Beto";

%% statisitcs of succeed for individual GAN, BigGAN and FC6
GANstr = ["FC6","BigGAN"];
GANsucsmsks = {pval_end_arr(:,1)<0.01, pval_end_arr(:,2)<0.01}; % end vs init 
GANsucsmsks = {pval_max_arr(:,1)<0.01, pval_max_arr(:,2)<0.01}; % max vs init 
msk_col = {Amsk,Bmsk,V4msk,ITmsk};
label_col = ["Alfa","Beto","V4","IT"];
for GANi = 1:numel(GANstr)
GANsucsmsk = GANsucsmsks{GANi};
fprintf("Using GAN %s\n",GANstr(GANi))
fprintf("Success experiments/Total %d/%d  (%.1f%%)\n",...
    sum(GANsucsmsk),numel(GANsucsmsk),100*sum(GANsucsmsk)/numel(GANsucsmsk))
for mi = 1:numel(msk_col)
    msk = msk_col{mi};
    msklabel = label_col(mi);
    fprintf("%s Success experiments/Total %d/%d (%.1f%%)\n",msklabel,...
        sum(msk&GANsucsmsk),sum(msk),sum(msk&GANsucsmsk)/sum(msk)*100)
end
end
%% basic statisitcs of succeed
sucsmsk = sucsmsk_max;
msk_col = {Amsk,Bmsk,V4msk,ITmsk};
label_col = ["Alfa","Beto","V4","IT"];
fprintf("Success experiments/Total %d/%d  (%.1f%%)\n",...
    sum(sucsmsk),numel(sucsmsk),100*sum(sucsmsk)/numel(sucsmsk))
for mi = 1:numel(msk_col)
    msk = msk_col{mi};
    msklabel = label_col(mi);
    fprintf("%s Success experiments/Total %d/%d (%.1f%%)\n",msklabel,...
        sum(msk&sucsmsk),sum(msk),sum(msk&sucsmsk)/sum(msk)*100)
end
%%
figh = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
	{V4msk&sucsmsk_end,ITmsk&sucsmsk_end}, ["V4","IT"], {}, [], ["FC6", "BigGAN"], false);
%% Final version
[score_m_traj_col, block_traj_col, score_m_traj_extrap_col, block_traj_extrap_col, ...
   extrap_mask_col, sucsmsk_end, sucsmsk_max, tval_end_arr, pval_end_arr, tval_max_arr, pval_max_arr] =...
   optim_traj_process(block_traces, score_traces, ["FC6", "BigGAN"], "max1", 35);
[figh,T] = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
   {V4msk&sucsmsk_max,ITmsk&sucsmsk_max}, ["V4","IT"], {Amsk,Bmsk}, ["Alfa","Beto"], ["FC6", "BigGAN"], true);
saveallform([outdir,figdir],"Both_area_sep_FC6maxnorm_pop_w_indiv",figh)
%% Final Published version
[score_m_traj_col, block_traj_col, score_m_traj_extrap_col, block_traj_extrap_col, ...
   extrap_mask_col, sucsmsk_end, sucsmsk_max, tval_end_arr, pval_end_arr, tval_max_arr, pval_max_arr] =...
   optim_traj_process(block_traces, score_traces, ["FC6", "BigGAN"], "max12", 35);
[figh,T] = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
   {V4msk&sucsmsk_max,ITmsk&sucsmsk_max}, ["V4","IT"], {Amsk,Bmsk}, ["Alfa","Beto"], ["FC6", "BigGAN"], true);
for i=1:4,nexttile(T,i);ylim([0,1]);end
saveallform([outdir,figdir],"Both_area_sep_maxbothnorm_pop_w_indiv",figh)
% for i=1:4,nexttile(T,i);ylim([0,1]);end

%%
figh = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
   {V4msk&sucsmsk_max,ITmsk&sucsmsk_max}, ["V4","IT"], {Amsk,Bmsk}, ["Alfa","Beto"], ["FC6", "BigGAN"], true);
saveallform([outdir,figdir],"Both_area_sep_maxnorm_pop_w_indiv",figh)
%
figh = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
   {V4msk&sucsmsk_max,ITmsk&sucsmsk_max}, ["V4","IT"], {Amsk,Bmsk}, ["Alfa","Beto"], ["FC6", "BigGAN"], false);
saveallform([outdir,figdir],"Both_area_sep_maxnorm_pop",figh)
%% Collect stats 
% normalize by max of the two trajs
% normalize by min or init gen the two trajs
normscore_mat = nan(size(score_traces));
for iTr=1:size(score_traces,1)
    [score_FC_m,score_FC_s,blockvec] = sort_scoreblock(block_traces{iTr,1},score_traces{iTr,1});
    [score_BG_m,score_BG_s,blockvec] = sort_scoreblock(block_traces{iTr,2},score_traces{iTr,2});
%     normalizer = min([mean(score_FC_m(2:3)),mean(score_BG_m(2:3))]);
    normalizer = max(score_FC_m);
%     normalizer = max([score_FC_m,score_BG_m]);
%     normalizer = mean(score_FC_m(end-1:end));
    normscore_mat(iTr,1) = max(score_FC_m)/normalizer;
    normscore_mat(iTr,2) = max(score_BG_m)/normalizer;
%     normscore_mat(iTr,1) = mean(score_FC_m(end-1:end))/normalizer;
%     normscore_mat(iTr,2) = mean(score_BG_m(end-1:end))/normalizer;
end
%%
cc_col = {normscore_mat(V4msk&sucsmsk_max,1),
          normscore_mat(V4msk&sucsmsk_max,2),
          normscore_mat(ITmsk&sucsmsk_max,1),
          normscore_mat(ITmsk&sucsmsk_max,2),};
label_arr = ["V4 FC6", "V4 BigGAN",...
             "IT FC6", "IT BigGAN"];
figh=figure(2);clf
violinplot_cell(cc_col, label_arr,'showData',true,'GroupOrder',cellstr(label_arr))
%%
cc_col = {normscore_mat(V4msk&sucsmsk_max&Amsk,2),
          normscore_mat(V4msk&sucsmsk_max&Bmsk,2),
          normscore_mat(ITmsk&sucsmsk_max&Amsk,2),
          normscore_mat(ITmsk&sucsmsk_max&Bmsk,2),};
label_arr = ["V4 BigGAN A", "V4 BigGAN B", ...
             "IT BigGAN A", "IT BigGAN B", ];
figh=figure(3); clf;
viol_arr = violinplot_cell(cc_col, label_arr, ...
           'showData',true,'GroupOrder',cellstr(label_arr));
hline(1.0); ylabel("Max BigGAN Activ Normalized By Max FC6")
title("Activation Ratio Comparison")
saveallform([outdir,figdir], "BigGAN_area_anim_cmp_FC6norm", figh);
%%
ttest2_print(normscore_mat(V4msk&sucsmsk_max&Amsk,2),normscore_mat(ITmsk&sucsmsk_max&Amsk,2),'V4 A','IT A')
ttest2_print(normscore_mat(V4msk&sucsmsk_max&Bmsk,2),normscore_mat(ITmsk&sucsmsk_max&Bmsk,2),'V4 B','IT B')
ttest2_print(normscore_mat(V4msk&sucsmsk_max,2),normscore_mat(ITmsk&sucsmsk_max,2),'V4 AB','IT AB')
%% Older API version
%% Normalize score by max 
score_C_m_traj_col = [];
score_G_m_traj_col = [];
score_C_traj_extrap_col = [];
score_G_traj_extrap_col = [];
block_traj_col = [];
blockN_extrap = 50;
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
    extrap_val_C = mean(score_C_m_norm(end));
    extrap_val_G = mean(score_G_m_norm(end));
    if numel(blockvec) < blockN_extrap
        score_C_traj_extrap_col{iTr} = [score_C_m_norm,extrap_val_C * ones(1,blockN_extrap-numel(score_C_m_norm))];
        score_G_traj_extrap_col{iTr} = [score_G_m_norm,extrap_val_G * ones(1,blockN_extrap-numel(score_G_m_norm))];
        block_traj_extrap_col{iTr} = 1:blockN_extrap;
    else
        score_C_traj_extrap_col{iTr} = score_C_m_norm;
        score_G_traj_extrap_col{iTr} = score_G_m_norm;
        block_traj_extrap_col{iTr} = blockvec;
    end
end
%% Create an Success Mask
GANstr = ["FC6","BigGAN"];
tval_end_arr = []; pval_end_arr = [];
tval_max_arr = []; pval_max_arr = [];
for iTr = 1:size(score_traces, 1)
for GANi=1:2
fprintf(GANstr(GANi)+" ")
blockN = max(block_traces{iTr,GANi});
[tval,pval] = ttest2_print(score_traces{iTr,GANi}(any(block_traces{iTr,GANi}==[1,2],2)),...
    score_traces{iTr,GANi}(any(block_traces{iTr,GANi}==[blockN-1,blockN],2)),"init12","last12");
tval_end_arr(iTr,GANi) = tval;
pval_end_arr(iTr,GANi) = pval;
[score_m,score_s,blockvec] = sort_scoreblock(block_traces{iTr,GANi},...
                score_traces{iTr,GANi});
[~,maxN] = max(score_m);
[tval,pval] = ttest2_print(score_traces{iTr,GANi}(any(block_traces{iTr,GANi}==[1,2],2)),...
    score_traces{iTr,GANi}(any(block_traces{iTr,GANi}==[maxN-1,maxN],2)),"init12","max12");
tval_max_arr(iTr,GANi) = tval;
pval_max_arr(iTr,GANi) = pval;
end
end
% masks of experiments that any one of it succeed
sucsmsk_end = any((tval_end_arr<0)&(pval_end_arr<0.001), 2);
sucsmsk_max = any((tval_max_arr<0)&(pval_max_arr<0.001), 2);
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