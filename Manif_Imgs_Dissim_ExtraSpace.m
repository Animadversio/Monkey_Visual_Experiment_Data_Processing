global summarydir
summarydir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning\summary";
%%
Animal = "Beto";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
% load(fullfile(mat_dir, "gab_imdist.mat"),'gab_imdist')
% load(fullfile(mat_dir, "pasu_imdist.mat"),'pasu_imdist')
% load(fullfile(mat_dir, Animal+"_EvoRef_ImDist.mat"),"EvoRefImDistStat")
load(fullfile(mat_dir, Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
D = torchImDist();
%%
ManifImDistStat_Ext = struct();
metname = "squ";
for Expi = 1:10
for si=1:numel(Stats(Expi).manif.psth)
    tic
imnm_grid = string(cellfun(@(idx)Stats(Expi).imageName{idx(1)},Stats(Expi).manif.idx_grid{si},"Uni",0));
manif_imgrid = arrayfun(@(nm)imread(fullfile(Stats(Expi).meta.stimuli,nm+".jpg")),imnm_grid,"Uni",0);
manif_img_tsr = cell2mat(reshape(manif_imgrid,1,1,1,[]));
D.select_metric("squeeze");
ManifImDistStat_Ext(Expi).squ{si} = D.distmat_B(manif_img_tsr);
toc
end
end
%%
save(fullfile(mat_dir, Animal+"_Manif_ImDist_Ext.mat"),"ManifImDistStat_Ext")
%%
%
tab_all = [];
%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
summarydir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning\summary";
for Animal = ["Alfa","Beto"]
% load(fullfile(mat_dir, "gab_imdist.mat"),'gab_imdist')
% load(fullfile(mat_dir, "pasu_imdist.mat"),'pasu_imdist')
% load(fullfile(mat_dir, Animal+"_EvoRef_ImDist.mat"),"EvoRefImDistStat")
load(fullfile(mat_dir, Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
manifAUCWid = [];
metname = "squ";
csri = 1;
for Expi = 1:numel(Stats)
for si=1:numel(Stats(Expi).manif.psth)
manifAUCWid(csri).Expi = Expi;
manifAUCWid(csri).Animal = Animal;
manifAUCWid(csri).pref_chan = EStats(Expi).units.pref_chan;
manifAUCWid(csri).pref_unit = EStats(Expi).evol.unit_in_pref_chan;
if EStats(Expi).units.pref_chan < 33
manifAUCWid(csri).area = "IT";
elseif EStats(Expi).units.pref_chan < 49
manifAUCWid(csri).area = "V1";  
else
manifAUCWid(csri).area = "V4"; 
end
ui=EStats(Expi).evol.unit_in_pref_chan;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_sem_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all')/sqrt(size(psth,3)),reshape(Stats(Expi).manif.psth{si},[],1));
score_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).manif.psth{si},[],1),'uni',0)); % single trial
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
manifAUCWid(csri).maxScore = maxScore;
% for mi = 1:numel(metric_list) % Loop through different distance metrics
% 	metname = metric_list(mi);
if si==1
dist_vec = ManifImDistStat(Expi).(metname)(:,maxId);
else
dist_vec = ManifImDistStat_Ext(Expi).(metname){si}(:,maxId);
end
[gpr_fit, ~, AUC, gprMdl] = GPRfitting(dist_vec, score_vec);
[~,~,AUC_linsmth,smth_Y] = BinSmoothing(dist_vec, score_vec);
AUC_raw = trapz_raw(dist_vec, score_vec);
manifAUCWid(csri).AUC_linsmth = AUC_linsmth; 
manifAUCWid(csri).AUC_raw = AUC_raw; 
manifAUCWid(csri).AUC = AUC; 
manifAUCWid(csri).normAUC_linsmth = AUC_linsmth / maxScore; 
manifAUCWid(csri).normAUC_raw = AUC_raw / maxScore; 
manifAUCWid(csri).normAUC = AUC / maxScore; 


actmap_col = cellfun(@(psth)squeeze(mean(psth(ui,51:200,:),2)),...
    reshape(Stats(Expi).manif.psth{si},[],1),'uni',0);
anovaStats = anova_cells(actmap_col);
manifAUCWid(csri).F = anovaStats.F;
manifAUCWid(csri).F_P = anovaStats.F_P;
manifAUCWid(csri).F_df = anovaStats.STATS.df;
manifAUCWid(csri).FSTATS = anovaStats.STATS;
csri = csri + 1;
end
end
tab = struct2table(manifAUCWid);
writetable(tab,fullfile(summarydir,Animal+"_ManifAUC_Tab_squ_ExtSp.csv"));
tab_all = [tab_all;tab];
end
%%
writetable(tab_all, fullfile(summarydir,"Both"+"_ManifAUC_Tab_squ_ExtSp.csv"));
%%
AUCtab = tab_all;
validmsk = ~ (AUCtab.Animal=="Alfa" & AUCtab.Expi==10);
Fmsk = AUCtab.F_P < 1E-3;
Alfamsk = (AUCtab.Animal=="Alfa");
Betomsk = (AUCtab.Animal=="Beto");
V1msk = AUCtab.area=="V1";%(AUCtab.chan<=48 & AUCtab.chan>=33);
V4msk = AUCtab.area=="V4";%(AUCtab.chan>48);
ITmsk = AUCtab.area=="IT";%(AUCtab.chan<33);
%%
msk = Fmsk & validmsk;
diary(fullfile(summarydir,'normAUC_WidthProg_Stats.log'))
fprintf("Test areal progression of statistics normAUC_bsl :\n")
test_progression(AUCtab, "normAUC", {V1msk & msk, V4msk & msk, ITmsk & msk});
diary off
%%
stripe_minor_plot(AUCtab, "normAUC", {V1msk & msk, V4msk & msk, ITmsk & msk},...
    ["V1","V4","IT"], {Alfamsk, Betomsk}, ["Alfa", "Beto"], "all signif ANOVA P<1E-3",...
    "area_anim_sep", {[3,2],[2,1],[3,1]})

function lm = test_progression(tab,varnm,msks,labels)
var_vec = []; idx_vec = [];
for i = 1:numel(msks)
msk = msks{i};
var_arr = tab.(varnm)(msk);
idx_arr = i*ones(numel(var_arr),1);
var_vec = [var_vec; var_arr];
idx_vec = [idx_vec; idx_arr];
end
[cval,pval] = corr(idx_vec, var_vec, 'Type', 'Spearman');
fprintf("Spearman Correlation %.3f(%.1e) df=%d\n",cval,pval,numel(var_vec)-1)
lm = fitlm(idx_vec, var_vec);
fprintf("Linear Regression t=%.3f(p=%.1e) R2=%.3f\n",...
    lm.Coefficients.tStat(2),lm.Coefficients.pValue(2),lm.Rsquared.Adjusted)
disp(lm)
end

function [gpr_fit, xlinsp, AUC, gprMdl] = GPRfitting(X, Y)
gprMdl = fitrgp(X, Y, 'SigmaLowerBound', 5E-3); % Edit this lower bound if encoutering fitting problem! 
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
AUC = trapz(xlinsp, gpr_fit);
% plot(xlinsp, gpr_fit, vargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend
end

function [r,m,b] = Linearfitting(X, Y, varargin) % Util function to add a linear reg line to a scatter
if nargin==2, varargin={};end
[r,m,b] = regression(reshape(X,1,[]), reshape(Y,1,[]));
% xmin = min(X); xmax = max(X);
% p = plot([xmin,xmax],[xmin,xmax].*m+b,varargin{:});
% set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

function [smth_fit, xlinsp, AUC, smthcrv] = BinSmoothing(X, Y, interpmethod)
% defaultly use linear interpolation between smoothed data points 
if nargin==2, varargin={}; interpmethod = 'linear'; 
elseif nargin==3, varargin={}; end
smthcrv = smooth(X,Y);
xlinsp = linspace(0,max(X),100);
[uniqX, iX, iUniq] = unique(X); % get rid of redundancy in X if there is any.
smth_fit = interp1(X(iX),smthcrv(iX),xlinsp,interpmethod,'extrap'); %'spline'
AUC = trapz(xlinsp, smth_fit);  % integrate under the gaussian process fitting curve. 
% % Plot curve on a fig
% scatter(X, smthcrv, 'DisplayName', 'Data Smooth')
% plot(xlinsp, smth_fit,'Disp',strcat('smooth+interp ',interpmethod), varargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend
end

function [AUC_raw] = trapz_raw(X,Y)
[sorted_X, sortid] = sort(X);
AUC_raw = trapz(sorted_X, Y(sortid));
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