%% Evol_Converg_Speed

Animal = "Both"; Set_Path;
mat_dir = "O:\Mat_Statistics";

for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
tic
EvolTrajStat = repmat(struct(),1,numel(EStats));
for Expi = 1:numel(EStats)
    ui = EStats(Expi).evol.unit_in_pref_chan;
    act_col = cellfun(@(psth)squeeze(mean(psth(ui,51:200,:))), EStats(Expi).evol.psth, 'uni',0);
    bsl_col = cellfun(@(psth)squeeze(mean(psth(ui,1:50,:))), EStats(Expi).evol.psth, 'uni',0);
    gen_col = arrayfun(@(i)i*ones(size(EStats(Expi).evol.psth{i},3),1),1:numel(EStats(Expi).evol.psth),'uni',0);
    gen_vec = cell2mat(gen_col');
    act_vec = cell2mat(reshape(act_col,[],1));
    bsl_vec = cell2mat(reshape(bsl_col,[],1));
    bsl_mean = mean(bsl_vec);
    bsl_std = std(bsl_vec);
    bsl_sem = sem(bsl_vec);
    act_traj_mean = cellfun(@mean, act_col);
    act_traj_std = cellfun(@std, act_col);
    act_traj_sem = cellfun(@sem, act_col);
    gprMdl = fitrgp(gen_vec, act_vec);
    gpr_fit = gprMdl.predict([1:max(gen_vec)-1]');
    init_act = gpr_fit(1);
    fina_act = gpr_fit(end);
    [max_act, maxid] = max(gpr_fit);
    bnd = [1, maxid]; %max(gen_vec)-1];
    step50 = fzero(@(x)gprMdl.predict(x)-((max_act-init_act) * 0.5 + init_act), bnd);
    step63 = fzero(@(x)gprMdl.predict(x)-((max_act-init_act) * 0.6321 + init_act), bnd);
    step85 = fzero(@(x)gprMdl.predict(x)-((max_act-init_act) * 0.85 + init_act), bnd);
    EvolTrajStat(Expi).step50 = step50;
    EvolTrajStat(Expi).step63 = step63;
    EvolTrajStat(Expi).step85 = step85;
    EvolTrajStat(Expi).init_act = init_act;
    EvolTrajStat(Expi).fina_act = fina_act;
    EvolTrajStat(Expi).max_act = max_act;
    EvolTrajStat(Expi).act_vec = act_vec;
    EvolTrajStat(Expi).bsl_vec = bsl_vec;
    EvolTrajStat(Expi).gprMdl = gprMdl;
    EvolTrajStat(Expi).gpr_fit = gpr_fit;
    EvolTrajStat(Expi).act_traj_mean = act_traj_mean;
    EvolTrajStat(Expi).act_traj_sem = act_traj_sem;
    EvolTrajStat(Expi).act_traj_std = act_traj_std;
    [~,T_P,~,TST] = ttest2(act_col{1}, act_col{maxid});
    EvolTrajStat(Expi).t_initmax = TST.tstat;
    EvolTrajStat(Expi).t_p_initmax = T_P;
    [~,T_P,~,TST] = ttest2(act_col{1}, act_col{end-1});
    EvolTrajStat(Expi).t_initend = TST.tstat;
    EvolTrajStat(Expi).t_p_initend = T_P;
    EvolTrajStat(Expi).DAOA_initmax = (max_act - init_act) / init_act;
    EvolTrajStat(Expi).DAOA_initend = (fina_act - init_act) / init_act;
%     figure(4);clf;
%     plot([1:max(gen_vec) - 1]',gpr_fit);hold on 
%     scatter(gen_vec, act_vec);
%     title(Expi)
%     box off
    toc
%     break
end
%%
save(fullfile(mat_dir,Animal+"_EvolTrajStats.mat"),'EvolTrajStat')
%%
end 
%%
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir,Animal+"_EvolTrajStats.mat"),'EvolTrajStat')
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
StatTab = repmat(struct(),1,numel(EStats));
for Expi = 1:numel(EStats)
	StatTab(Expi).Expi = Expi;
	StatTab(Expi).Animal = Animal;
	StatTab(Expi).pref_chan = EStats(Expi).evol.pref_chan;
	StatTab(Expi).unit_in_pref_chan = EStats(Expi).evol.unit_in_pref_chan;
	StatTab(Expi).imgsize = EStats(Expi).evol.imgsize;
	StatTab(Expi).imgpos = EStats(Expi).evol.imgpos;
    StatTab(Expi).step50 = EvolTrajStat(Expi).step50;
    StatTab(Exp&mski).step63 = EvolTrajStat&msk(Expi).step63;
    StatTab(Expi).step85 = EvolTrajStat(Expi)&msk.step85;
    StatTab(Expi).init_act = EvolTrajStat(Expi)&msk.init_act;
    StatTab(Expi).fina_act = EvolTrajStat(Expi)&msk.fina_act;
    StatTab(Expi).t_initmax = EvolTrajStat(Expi).t_initmax;
	StatTab(Expi).t_p_initmax = E&mskvolTrajStat(Expi).t_p_initmax;
	StatTab(Expi).t_initend = EvolTrajStat(Expi).t_initend;
	StatTab(Expi).t_p_initend = EvolTrajStat(Expi).t_p_initend;
	StatTab(Expi).DAOA_initmax = EvolTrajStat(Expi).DAOA_initmax;
	StatTab(Expi).DAOA_initend = EvolTrajStat(Expi).DAOA_initend;
end
StatTab = struct2table(StatTab);
%%
writetable(StatTab, fullfile(mat_dir, Animal+"_EvolTrajStats.csv"))
end
%% Reload and merge the stats for evolution trajectories. 
tabA = readtable(fullfile(mat_dir, "Alfa"+"_EvolTrajStats.csv"));
tabB = readtable(fullfile(mat_dir, "Beto"+"_EvolTrajStats.csv"));
StatTab = cat(1, tabA, tabB);
writetable(StatTab, fullfile(mat_dir, "Both"+"_EvolTrajStats.csv"));
%%
msk = (StatTab.t_p_initmax<1E-2);
V1msk = StatTab.pref_chan <=48 & StatTab.pref_chan >= 33;
V4msk = StatTab.pref_chan <=64 & StatTab.pref_chan >= 49;
ITmsk = StatTab.pref_chan <=32 & StatTab.pref_chan >= 1;
% ttest2(StatTab.step63(V4msk), StatTab.step63(V1msk));
[~,P,CI,TST] = ttest2(StatTab.step63(V4msk&msk), StatTab.step63(V1msk&msk));
fprintf("V4 - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
[~,P,CI,TST] = ttest2(StatTab.step63(ITmsk&msk), StatTab.step63(V1msk&msk));
fprintf("IT - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
[~,P,CI,TST] = ttest2(StatTab.step63(ITmsk&msk), StatTab.step63(V4msk&msk));
fprintf("IT - V4: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
%%
[~,P,CI,TST] = ttest2(StatTab.step50(V4msk&msk), StatTab.step50(V1msk&msk));
fprintf("V4 - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
[~,P,CI,TST] = ttest2(StatTab.step50(ITmsk&msk), StatTab.step50(V1msk&msk));
fprintf("IT - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
[~,P,CI,TST] = ttest2(StatTab.step50(ITmsk&msk), StatTab.step50(V4msk&msk));
fprintf("IT - V4: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
%%
global sumdir
sumdir = "O:\EvolTraj_Cmp\summary";
Alfamsk = StatTab.Animal == "Alfa";
Betomsk = StatTab.Animal == "Beto";
h = stripe_plot(StatTab, "step63", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1", "V4", "IT"], ...
    "all exps", "area_cmp", {[3,1],[2,1],[3,2]});
%%
h = stripe_minor_plot(StatTab, "step63", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1", "V4", "IT"], ...
    {Alfamsk&msk, Betomsk&msk}, ["Alfa", "Beto"], "all exps", "area_anim_cmp", {[3,1],[2,1],[3,2]});
%%
msk = (StatTab.t_p_initmax<1E-2);
area_prog_cmp(StatTab, "step63", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1", "V4", "IT"], {[3,1],[2,1],[3,2]})
area_prog_cmp(StatTab, "step85", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1", "V4", "IT"], {[3,1],[2,1],[3,2]})
area_prog_cmp(StatTab, "step50", {V1msk&msk, V4msk&msk, ITmsk&msk}, ["V1", "V4", "IT"], {[3,1],[2,1],[3,2]})
%% Plot example trajectories 
Animal = "Alfa";
A = load(fullfile(mat_dir,Animal+"_EvolTrajStats.mat"),'EvolTrajStat')
Animal = "Beto";
B = load(fullfile(mat_dir,Animal+"_EvolTrajStats.mat"),'EvolTrajStat')
%% Creat masks
sucsmsk = (StatTab.t_p_initmax<1E-2);
V1msk = StatTab.pref_chan <=48 & StatTab.pref_chan >= 33;
V4msk = StatTab.pref_chan <=64 & StatTab.pref_chan >= 49;
ITmsk = StatTab.pref_chan <=32 & StatTab.pref_chan >= 1;
Alfamsk = StatTab.Animal == "Alfa";
Betomsk = StatTab.Animal == "Beto";
%%
figdir = "O:\EvolTraj_Cmp\summary";
EvolTrajStat = [A.EvolTrajStat,B.EvolTrajStat];
traj_col = arrayfun(@(E)E.act_traj_mean, EvolTrajStat,'uni',0);
summarize_trajs_tile(traj_col,"max",1,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim");
summarize_trajs_tile(traj_col,"max",3,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_movmean");
%%
summarize_trajs_tile(traj_col,"max",3,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    {},["Both"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_movmean");
%%
summarize_trajs_tile(traj_col,"max",1,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_sucs");
summarize_trajs_tile(traj_col,"max",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_sucs_movmean");
%% Merge trajs in the same tile for an monk 
summarize_trajs_merge(traj_col,"max",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_Anim_sucs_movmean_merge");
summarize_trajs_merge(traj_col,"max",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {},["Both"],figdir,"Both_MaxNorm_scoreTraj_avg_Area_sucs_movmean_merge");
%%
figure(14);AlignAxisLimits(get(14,'Child'))
%%
summarize_trajs_merge(traj_col,"movmean_rng",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_SmthRngNorm_scoreTraj_avg_Area_Anim_sucs_movmean_merge");
summarize_trajs_merge(traj_col,"movmean_rng",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {},["Both"],figdir,"Both_SmthRngNorm_scoreTraj_avg_Area_sucs_movmean_merge");
%%
summarize_trajs_merge(traj_col,"rng",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {Alfamsk,Betomsk},["Alfa","Beto"],figdir,"Both_RngNorm_scoreTraj_avg_Area_Anim_sucs_movmean_merge");
summarize_trajs_merge(traj_col,"rng",3,{V1msk&sucsmsk,V4msk&sucsmsk,ITmsk&sucsmsk},["V1","V4","IT"],...
    {},["Both"],figdir,"Both_RngNorm_scoreTraj_avg_Area_sucs_movmean_merge");
function [h,score_m,score_s,blockvec] = summarize_trajs_tile(traj_col,norm_mode,movmean_N,major_msk,major_lab,minor_msk,minor_lab,figdir,fignm)
% msk_col = {V1msk&validmsk, V4msk&validmsk, ITmsk&validmsk};
% anim_msks = {Alfamsk&validmsk, Betomsk&validmsk};
% label_col = ["V1", "V4", "IT"];
% anim_col = ["Alfa","Beto"];
% msk on the experimental session level. 
if nargin<2, norm_mode="max";end
if nargin<3, movmean_N=3;end
if nargin<4, major_msk={ones(size(traj_col,1),1,'logical')}; major_lab=["all"]; end
if nargin<6 || isempty(minor_msk), 
    minor_msk={ones(size(traj_col,1),1,'logical')}; 
end
if nargin<7 || isempty(minor_lab), minor_lab = ["all"]; end 
if nargin<8, figdir=""; end
if nargin<9, fignm=compose("%sNorm_scoreTraj_avg_Area_Anim_movmean",norm_mode); end
Corder = colororder();
% Normalize the trajs 
if strcmp(norm_mode,"max") % normalize each traj to the max in that exp. 
normed_traj_col = cellfun(@(traj) traj/max(traj),traj_col,'uni',0);
elseif strcmp(norm_mode,"rng") % normalize each traj, min to 0, max to 1
normed_traj_col = cellfun(@(traj) (traj-min(traj))/(max(traj)-min(traj)),traj_col,'uni',0);
elseif strcmp(norm_mode,"movmean_rng") % normalize each traj to min max of smooth traj
normed_traj_col = [];
for i = 1:numel(traj_col)
smooth_traj = movmean(traj_col{i},3);
MIN = min(smooth_traj); 
MAX = max(smooth_traj);
normed_traj_col{i} = (traj_col{i}-MIN)/(MAX-MIN);
end
normed_traj_col = reshape(normed_traj_col, size(traj_col));
else
normed_traj_col = cellfun(@(traj) traj/max(traj),traj_col,'uni',0);
end
% block cell array
block_col = cellfun(@(traj)1:numel(traj),traj_col,'uni',0);
% [score_m,score_s,blockvec] = sort_scoreblock(block_col, normed_traj_col)
nrow = numel(minor_msk); ncol = numel(major_msk);
h=figure; %fignm=compose("%s_MaxNorm_scoreTraj_avg_Area_Anim_movmean", Animal);
set(h,'pos',[125,   258,   310*ncol, 100+310*nrow])
T = tiledlayout(nrow,ncol,"pad",'compact',"tilespac",'compact');
for animi=1:nrow
for mski=1:ncol
nexttile(T,mski+ncol*(animi-1))
msk = major_msk{mski} & minor_msk{animi};
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_col(msk),normed_traj_col(msk));
score_C_col_mov_m = movmean(score_C_col_m,movmean_N);
score_C_col_mov_s = movmean(score_C_col_s,movmean_N);
shadedErrorBar(block_colvec,score_C_col_mov_m,score_C_col_mov_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory\n%s %s (n=%d)",minor_lab(animi),major_lab(mski),sum(msk)))
% legend(["Full","50D"])
end
end
title(T, compose("Summary of mean evol trajectory (%d block movmean) for each area",movmean_N), 'FontSize',16)
saveallform(figdir,fignm);
for animi=1:nrow
for mski=1:ncol
    nexttile(T,mski+ncol*(animi-1))
    msk = major_msk{mski} & minor_msk{animi};
    xlim([1,prctile(cellfun(@numel,block_col(msk)),80)]);
end
end
saveallform(figdir,fignm+"_Xlim");
end

function area_prog_cmp(StatTab, varnm, msks, labels, pairs)
fprintf("Comparison of %s among areas\n",varnm)
fullvec = [];
idxvec = [];
for i = 1:numel(msks)
varvec = StatTab.(varnm)(msks{i});
fprintf("%s: %.1f(%.1f,n=%d)  ",labels(i),mean(varvec),sem(varvec),numel(varvec));
fullvec = [fullvec; varvec];
idxvec = [idxvec; i * ones(numel(varvec), 1)];
end
fprintf("\n")
for pi = 1:numel(pairs)
i = pairs{pi}(1); j = pairs{pi}(2);
[~,P,CI,TST] = ttest2(StatTab.(varnm)(msks{i}), StatTab.(varnm)(msks{j}));
fprintf("%s - %s: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",labels(i),labels(j),P,TST.tstat,TST.df,CI(1),CI(2));
end
[cval, pval] = corr(idxvec, fullvec, 'Type','Spearman');
fprintf("Spearman corrlation of index and %s: %.3f (p=%.1e,n=%d)\n",varnm,cval,pval,numel(fullvec))
[F_p,F_tbl] = anova1(fullvec, idxvec, 'off');
Fval = F_tbl{2,5}; F_df = F_tbl{4,3};
anovastr = compose("ANOVA F=%.3f p=%.1e(df=%d,%d)\n",Fval,F_p,F_tbl{2,3},F_tbl{3,3});
lm = fitlm(idxvec, fullvec);
lmstr = compose("Linear Regres %s = %.3f + %s * %.3f \n Intercept %.3f+-%.3f, Slope %.3f+-%.3f\n Slope!=0: T=%.1f P=%.1e\n Fit Rsquare=%.3f\n",...
                varnm, lm.Coefficients.Estimate(1), "area", lm.Coefficients.Estimate(2),...
                lm.Coefficients.Estimate(1), lm.Coefficients.SE(1), ...
                lm.Coefficients.Estimate(2), lm.Coefficients.SE(2), ...
                lm.Coefficients.tStat(2), lm.Coefficients.pValue(2), ...
                lm.Rsquared.Ordinary);
fprintf(anovastr +lmstr+"\n")
% [~,P,CI,TST] = ttest2(StatTab.step50(ITmsk&msk), StatTab.step50(V1msk&msk));
% fprintf("IT - V1: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
% [~,P,CI,TST] = ttest2(StatTab.step50(ITmsk&msk), StatTab.step50(V4msk&msk));
% fprintf("IT - V4: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
end

function [gpr_fit, xlinsp, AOC, gprMdl] = GPRfitting(X, Y)
gprMdl = fitrgp(X, Y, 'SigmaLowerBound', 5E-3); % Edit this lower bound if encoutering fitting problem! 
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
AOC = trapz(xlinsp, gpr_fit);  % integrate under the gaussian process fitting curve. 
% Plot curve on a fig
% plot(xlinsp, gpr_fit, vargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend
end
function h = stripe_plot(StatsTab_sum, statname, masks, labels, titstr, savestr, Tpairs)
if nargin<7, Tpairs = {}; end
global sumdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
xjitter = 0.15*randn(numel(statscol{i}),1);
scatter(i+xjitter, statscol{i})
% medval = median(StatsTab_sum.(statname)(masks{i}));
% vline(medval,"-.",labels(i)+" "+num2str(medval))
end
xticks(1:numel(masks));xticklabels(labels)
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
Ns = cellfun(@sum,masks);
ylabel(statname)
title(title_str)
legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmp.png", statname, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmp.pdf", statname, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_stripcmp.fig", statname, savestr)))
end

function h = stripe_minor_plot(StatsTab_sum, statname, masks, labels, minormsks, minorlabels, titstr, savestr, Tpairs)
% Stripe plot comparing a scaler data w.r.t. 2 variables, signified by masks, labels  and  minormsks, minorlabels
%  masks, labels are plotted in different columns; minormsks, minorlabels are plotted in same columne by different color.
if nargin<7, Tpairs = {}; end
marker = 'o*x^v';
global sumdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
labelcol = []; Ns = []; Stds = []; Means = [];
for i = 1:numel(masks)
    for j = 1:numel(minormsks)
    M = masks{i}&minormsks{j};
    xjitter = 0.15*randn(sum(M),1);
    scatter(i+xjitter, StatsTab_sum.(statname)(M),'Marker',marker(j))
    labelcol = [labelcol, labels(i)+minorlabels(j)];
    Ns = [Ns, sum(M)];
    Stds = [Stds, std(StatsTab_sum.(statname)(M))];
    Means = [Means, mean(StatsTab_sum.(statname)(M))];
    end
end
xticks(1:numel(masks));xticklabels(labels)
ylabel(statname)
statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
FStat = anova_cells(statscol); % F stats is across the major variables/ msks
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e df=%d)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P,FStat.STATS.df);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
title(title_str)
% Ns = cellfun(@sum,masks);
% legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
legend(compose("%s(%d) %.2f (%.2f)",labelcol',Ns',Means',Stds'),'Location','best') % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmpmarker.png", statname, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmpmarker.pdf", statname, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_stripcmpmarker.fig", statname, savestr)))
end
