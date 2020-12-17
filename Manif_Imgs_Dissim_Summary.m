%% Manif_Imgs_Dissim_Summary
%  Structure is borrowed from Manif_Imgs_Dissim_Merge, and it's more about collecting stats and plot.
%  This script is extremely well structured! 
Animal = "Beto";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir,"gab_imdist.mat"),'gab_imdist')
load(fullfile(mat_dir,"pasu_imdist.mat"),'pasu_imdist')
load(fullfile(mat_dir,Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";

%% Pre step, get the valid mask for pasupathy patches
if Animal == "Alfa"
pasu_idx_vec = reshape(Stats(1).ref.pasu_idx_grid',[],1); % 12% reshape into one row. But transpose to make similar object adjacent
elseif Animal == "Beto"
pasu_idx_vec = reshape(Stats(12).ref.pasu_idx_grid',[],1);
end
pasu_val_msk = ~cellfun(@isempty,pasu_idx_vec);

%%
% RadTuneStats = repmat(struct(), 1, numel(Stats));
%% Taking in all ImDist structure and COMPUTE Tuning Statistics w.r.t. it.
bestcorr = 0;
metric_list = ["squ","SSIM","L2","FC6","FC6_corr"];
label_list = ["LPIPS (SqueezeNet)", "SSIM", "L2", "FC6 (L2)", "FC6 (1 - corr)"];
for Expi=42:numel(Stats)
% Manifold images. 
fprintf("Processing %s Exp %d\n", Animal,Expi)
RadTuneStats(Expi).Animal = Animal; 
RadTuneStats(Expi).Expi = Expi;
si=1; ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
RadTuneStats(Expi).manif.maxScore = maxScore;
for mi = 1:numel(metric_list) % Loop through different distance metrics
	metname = metric_list(mi);
	[gpr_fit, ~, AOC, gprMdl] = GPRfitting(ManifImDistStat(Expi).(metname)(:,maxId), score_vec); % Gaussian Process Smoothing or Fitting
	[R2, slope, intercept] = Linearfitting(ManifImDistStat(Expi).(metname)(:,maxId), score_vec);
	% titstr{mi} = { compose("Manif: pear %.3f spear %.3f",
	RadTuneStats(Expi).manif.AOC.(metname) = AOC; 
	RadTuneStats(Expi).manif.gprMdl.(metname) = gprMdl;
	RadTuneStats(Expi).manif.linR2.(metname) = R2;
	RadTuneStats(Expi).manif.lincc.(metname) = [slope, intercept];
	RadTuneStats(Expi).manif.corr.(metname) = corr(ManifImDistStat(Expi).(metname)(:,maxId),score_vec);
	RadTuneStats(Expi).manif.corr_sp.(metname) = corr(ManifImDistStat(Expi).(metname)(:,maxId),score_vec,'Type','Spearman');
end
% Pasupathy patches
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(~pasu_val_msk) = []; % isnan(pasu_vec) % get rid of non-existing pasupathy images. 
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
RadTuneStats(Expi).pasu.maxScore = pasu_maxScore; 
for mi = 1:numel(metric_list)
	metname = metric_list(mi);
	[gpr_fit, ~, AOC, gprMdl] = GPRfitting(pasu_imdist.(metname)(:,pasu_maxId), pasu_vec);
	[R2, slope, intercept] = Linearfitting(pasu_imdist.(metname)(:,pasu_maxId), pasu_vec);
	% titstr{mi}{end+1} = compose("Pasu: pear %.3f spear %.3f",
	RadTuneStats(Expi).pasu.AOC.(metname) = AOC; 
	RadTuneStats(Expi).pasu.gprMdl.(metname) = gprMdl;
	RadTuneStats(Expi).pasu.linR2.(metname) = R2;
	RadTuneStats(Expi).pasu.lincc.(metname) = [slope, intercept];
	RadTuneStats(Expi).pasu.corr.(metname) = corr(pasu_imdist.(metname)(:,pasu_maxId),pasu_vec,'Rows','complete');
	RadTuneStats(Expi).pasu.corr_sp.(metname) = corr(pasu_imdist.(metname)(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete');
end	
legend(["Manifold","Manifold","Pasupathy","Pasupathy"])
end
% Gabor patches
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
RadTuneStats(Expi).gabor.maxScore = gab_maxScore; 
for mi = 1:numel(metric_list)
    metname = metric_list(mi);
    [gpr_fit, ~, AOC, gprMdl] = GPRfitting(gab_imdist.(metname)(:,gab_maxId), gab_vec);% % Plot the Interpolation of it
	[R2, slope, intercept] = Linearfitting(gab_imdist.(metname)(:,gab_maxId), gab_vec);
    % titstr{mi}{end+1} = compose("Gabor: pear %.3f spear %.3f",
	RadTuneStats(Expi).gabor.AOC.(metname) = AOC; 
	RadTuneStats(Expi).gabor.gprMdl.(metname) = gprMdl;
	RadTuneStats(Expi).gabor.linR2.(metname) = R2;
	RadTuneStats(Expi).gabor.lincc.(metname) = [slope, intercept];
	RadTuneStats(Expi).gabor.corr.(metname) = corr(gab_imdist.(metname)(:,gab_maxId),gab_vec,'Rows','complete');
	RadTuneStats(Expi).gabor.corr_sp.(metname) = corr(gab_imdist.(metname)(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete');
end
end
% saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.pdf",Animal,Expi,Stats(Expi).units.pref_chan)))
end
%%
save(fullfile(mat_dir,Animal+"_ManifRadiTuneStats.mat"),'RadTuneStats')
%% Collect the Hierarchical Stats in Mat file into a csv table for summarizing
summarydir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning\summary";
tab_all = table();
metname = "squ"; % pick the stats you like
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifRadiTuneStats.mat"),'RadTuneStats')
tab = struct();
for Expi = 1:numel(RadTuneStats)
tab(Expi).Expi = Expi;
tab(Expi).Animal = Animal;
tab(Expi).normAOC_mani = RadTuneStats(Expi).manif.AOC.(metname) / RadTuneStats(Expi).manif.maxScore;
tab(Expi).AOC_mani = RadTuneStats(Expi).manif.AOC.(metname);
tab(Expi).corr_mani = RadTuneStats(Expi).manif.corr.(metname);
if ~isempty(RadTuneStats(Expi).pasu)
tab(Expi).normAOC_pasu = RadTuneStats(Expi).pasu.AOC.(metname) / RadTuneStats(Expi).pasu.maxScore;
tab(Expi).AOC_pasu = RadTuneStats(Expi).pasu.AOC.(metname); 
tab(Expi).corr_pasu = RadTuneStats(Expi).pasu.corr.(metname); 
else
tab(Expi).AOC_pasu = nan; tab(Expi).normAOC_pasu = nan; tab(Expi).corr_pasu = nan;
end
if ~isempty(RadTuneStats(Expi).gabor)
tab(Expi).normAOC_gab = RadTuneStats(Expi).gabor.AOC.(metname) / RadTuneStats(Expi).gabor.maxScore;
tab(Expi).AOC_gab = RadTuneStats(Expi).gabor.AOC.(metname); 
tab(Expi).corr_gab = RadTuneStats(Expi).gabor.corr.(metname); 
else
tab(Expi).AOC_gab = nan; tab(Expi).normAOC_gab = nan; tab(Expi).corr_gab = nan;
end
% tab = [tab; [normAOC, AOC, normAOC_pasu, AOC_pasu, normAOC_gab, AOC_gab]];
end
tab = struct2table(tab);
writetable(tab, fullfile(summarydir,Animal+"_RadialTuningStatsTab.csv"))
tab_all = [tab_all; tab];
end
writetable(tab_all, fullfile(summarydir,"Both"+"_RadialTuningStatsTab.csv"))
%%
tab_all = readtable(fullfile(summarydir,"Both"+"_RadialTuningStatsTab.csv"));
tab = tab_all;
%%
Animal = "Both";
h=figure; set(h,'pos',[1000, 413, 640, 570]); hold on
scatter(tab.normAOC_mani, tab.corr_mani)
scatter(tab.normAOC_pasu, tab.corr_pasu)
scatter(tab.normAOC_gab, tab.corr_gab)
plot([tab.normAOC_mani,tab.normAOC_pasu]', [tab.corr_mani, tab.corr_pasu]', 'Color', [0,0,0,0.2],'HandleVisibility','off')
plot([tab.normAOC_mani,tab.normAOC_gab]', [tab.corr_mani, tab.corr_gab]', 'Color', [0,0,1,0.2],'HandleVisibility','off')
legend(["GAN Manifold Space", "Pasupathy", "Gabor"])
xlabel("Normalized AUC"); ylabel("Correlation Coefficient") 
title(["Manifold Exps Tuning AUC ~ Correlation "+Animal,"Squeeze Net LPIPS metric"])
savefig(h,fullfile(summarydir, Animal+"_AUC_Corr_scatt_squ.fig")) 
saveas(h,fullfile(summarydir, Animal+"_AUC_Corr_scatt_squ.png"))
saveas(h,fullfile(summarydir, Animal+"_AUC_Corr_scatt_squ.pdf"))
%%
h2=figure; set(h2,'pos',[1000, 413, 640, 570]); hold on
histogram(tab.normAOC_mani, 10, 'FaceAlpha', 0.4)%, tab.corr_mani)
histogram(tab.normAOC_pasu, 10, 'FaceAlpha', 0.4)%, tab.corr_pasu)
histogram(tab.normAOC_gab, 10, 'FaceAlpha', 0.4)%, tab.corr_gab)
% plot([tab.normAOC_mani,tab.normAOC_pasu]', [tab.corr_mani, tab.corr_pasu]', 'Color', [0,0,0,0.3],'HandleVisibility','off')
legend(["GAN Manifold Space", "Pasupathy", "Gabor"])
xlabel("Normalized AUC"); %ylabel("Correlation Coefficient") 
title(["Manifold Exps Density of Tuning AUC "+Animal,"Squeeze Net LPIPS metric"])
savefig(h2,fullfile(summarydir, Animal+"_AUC_hist_squ.fig")) 
saveas(h2,fullfile(summarydir, Animal+"_AUC_hist_squ.png"))
saveas(h2,fullfile(summarydir, Animal+"_AUC_hist_squ.pdf"))
%%
h3 = figure; hold on; set(h3,'pos',[1000,366,438,650]);
jit = randn(size(tab,1),1) * 0.1;
scatter(jit + 1, tab.normAOC_mani);
scatter(jit + 2, tab.normAOC_pasu);
scatter(jit + 3, tab.normAOC_gab);
plot([jit + [1:3]]', [tab.normAOC_mani,tab.normAOC_pasu,tab.normAOC_gab]', 'Color', [0,0,0,0.2]);
ylabel("Normalized AUC");xlabel("Tuning Space")
title(["Manifold Exps Tuning AUC "+Animal,"Squeeze Net LPIPS metric"])
xticks([1,2,3]); xticklabels(["GAN Manifold", "Pasupathy", "Gabor"])
savefig(h3, fullfile(summarydir, Animal+"_AUC_scatt_squ.fig")) 
saveas(h3, fullfile(summarydir, Animal+"_AUC_scatt_squ.png"))
saveas(h3, fullfile(summarydir, Animal+"_AUC_scatt_squ.pdf"))
%%
% Linearfitting(1:10,-[1:10]+0.1*randn(1,10))
% Util functions to compute Gaussian Process fitting and linear fitting.
function [gpr_fit, xlinsp, AOC, gprMdl] = GPRfitting(X, Y)
gprMdl = fitrgp(X, Y, 'SigmaLowerBound', 5E-3); % Edit this lower bound if encoutering fitting problem! 
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
AOC = trapz(xlinsp, gpr_fit);
% plot(xlinsp, gpr_fit, vargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend
end
function [r,m,b] = Linearfitting(X, Y, varargin) % Util function to add a linear reg line to a scatter
if nargin==2, varargin={};end
[r,m,b] = regression(X, Y);
% xmin = min(X); xmax = max(X);
% p = plot([xmin,xmax],[xmin,xmax].*m+b,varargin{:});
% set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end