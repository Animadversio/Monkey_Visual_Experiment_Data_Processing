%% Manif_Imgs_Dissim_Summary
%  Structure is borrowed from Manif_Imgs_Dissim_Merge, and it's more about collecting stats and plot.
%  This script is extremely well structured and concise!
%  Used to generate panels for figure 2 E  
%  Emphasize comparison across stimuli spaces.
%  Add the evoref images to
Animal = "Alfa";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir,"gab_imdist.mat"),'gab_imdist')
load(fullfile(mat_dir,"pasu_imdist.mat"),'pasu_imdist')
load(fullfile(mat_dir,Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
load(fullfile(mat_dir,Animal+"_EvoRef_ImDist.mat"),"EvoRefImDistStat") % newly added, image distances for evolved images
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
%% Taking in all ImDist structure and COMPUTE Tuning Statistics in all spaces w.r.t. it.
metric_list = ["squ","SSIM","L2","FC6","FC6_corr"];
label_list = ["LPIPS (SqueezeNet)", "SSIM", "L2", "FC6 (L2)", "FC6 (1 - corr)"];
for Expi=1:numel(Stats)
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
	dist_vec = ManifImDistStat(Expi).(metname)(:,maxId);
	[gpr_fit, ~, AUC, gprMdl] = GPRfitting(dist_vec, score_vec); % Gaussian Process Smoothing or Fitting
	[R2, slope, intercept] = Linearfitting(dist_vec, score_vec);% Linear fitting
	[~,~,AUC_linsmth,smth_Y] = BinSmoothing(dist_vec, score_vec);
	AUC_raw = trapz_raw(dist_vec, score_vec);
	% titstr{mi} = { compose("Manif: pear %.3f spear %.3f",
	RadTuneStats(Expi).manif.AUC_linsmth.(metname) = AUC_linsmth; 
	RadTuneStats(Expi).manif.AUC_raw.(metname) = AUC_raw; 
	RadTuneStats(Expi).manif.AUC.(metname) = AUC; 
	RadTuneStats(Expi).manif.gprMdl.(metname) = gprMdl;
	RadTuneStats(Expi).manif.linR2.(metname) = R2;
	RadTuneStats(Expi).manif.lincc.(metname) = [slope, intercept];
	[cval,pval] = corr(dist_vec,score_vec);
	RadTuneStats(Expi).manif.corr.(metname) = cval;
	RadTuneStats(Expi).manif.corr_P.(metname) = pval;
	[cval,pval] = corr(dist_vec,score_vec,'Type','Spearman');
	RadTuneStats(Expi).manif.corr_sp.(metname) = cval;
	RadTuneStats(Expi).manif.corr_sp_P.(metname) = pval;
end
% Evolution Reference images
evoref_vec = cellfun(@(psth)mean(psth(1,51:200,:),'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
[sortScore,sortId] = sort(evoref_vec,'Descend');
[maxScore,maxId] = max(evoref_vec);
RadTuneStats(Expi).evoref.maxScore = maxScore;
for mi = 1:numel(metric_list) % Loop through different distance metrics
	metname = metric_list(mi);
	evoref_dist_vec = EvoRefImDistStat(Expi).(metname)(:,maxId);
	[gpr_fit, ~, AUC, gprMdl] = GPRfitting(evoref_dist_vec, evoref_vec); % Gaussian Process Smoothing or Fitting
	[R2, slope, intercept] = Linearfitting(evoref_dist_vec, evoref_vec);% Linear fitting
	[~,~,AUC_linsmth,smth_Y] = BinSmoothing(evoref_dist_vec, evoref_vec);
	AUC_raw = trapz_raw(evoref_dist_vec, evoref_vec);
	% titstr{mi} = { compose("Manif: pear %.3f spear %.3f",
	RadTuneStats(Expi).evoref.AUC_linsmth.(metname) = AUC_linsmth; 
	RadTuneStats(Expi).evoref.AUC_raw.(metname) = AUC_raw; 
	RadTuneStats(Expi).evoref.AUC.(metname) = AUC; 
	RadTuneStats(Expi).evoref.gprMdl.(metname) = gprMdl;
	RadTuneStats(Expi).evoref.gpr_fit.(metname) = gpr_fit;
	RadTuneStats(Expi).evoref.linR2.(metname) = R2;
	RadTuneStats(Expi).evoref.lincc.(metname) = [slope, intercept];
	[cval,pval] = corr(evoref_dist_vec,evoref_vec);
	RadTuneStats(Expi).evoref.corr.(metname) = cval;
	RadTuneStats(Expi).evoref.corr_P.(metname) = pval;
	[cval,pval] = corr(evoref_dist_vec,evoref_vec,'Type','Spearman');
	RadTuneStats(Expi).evoref.corr_sp.(metname) = cval;
	RadTuneStats(Expi).evoref.corr_sp_P.(metname) = pval;
    fprintf("%d %s corr %.3f(%.1e) ",Expi,metname,RadTuneStats(Expi).evoref.corr.(metname),RadTuneStats(Expi).evoref.corr_P.(metname))
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
	pasu_dist_vec = pasu_imdist.(metname)(:,pasu_maxId);
	[gpr_fit, ~, AUC, gprMdl] = GPRfitting(pasu_dist_vec, pasu_vec);
	[R2, slope, intercept] = Linearfitting(pasu_dist_vec, pasu_vec);% Linear fitting
	[~,~,AUC_linsmth,smth_Y] = BinSmoothing(pasu_dist_vec, pasu_vec);
	AUC_raw = trapz_raw(pasu_dist_vec, pasu_vec);
	% titstr{mi}{end+1} = compose("Pasu: pear %.3f spear %.3f",
	RadTuneStats(Expi).pasu.AUC_linsmth.(metname) = AUC_linsmth; 
	RadTuneStats(Expi).pasu.AUC_raw.(metname) = AUC_raw; 
	RadTuneStats(Expi).pasu.AUC.(metname) = AUC; 
	RadTuneStats(Expi).pasu.gprMdl.(metname) = gprMdl;
	RadTuneStats(Expi).pasu.linR2.(metname) = R2;
	RadTuneStats(Expi).pasu.lincc.(metname) = [slope, intercept];
	[cval,pval] = corr(pasu_dist_vec,pasu_vec,'Rows','complete');
	RadTuneStats(Expi).pasu.corr.(metname) = cval;
	RadTuneStats(Expi).pasu.corr_P.(metname) = pval;
	[cval,pval] = corr(pasu_dist_vec,pasu_vec,'Type','Spearman','Rows','complete');
	RadTuneStats(Expi).pasu.corr_sp.(metname) = cval;
	RadTuneStats(Expi).pasu.corr_sp_P.(metname) = pval;
end	
end
% Gabor patches
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
RadTuneStats(Expi).gabor.maxScore = gab_maxScore; 
for mi = 1:numel(metric_list)
    metname = metric_list(mi);
    gab_dist_vec = gab_imdist.(metname)(:,gab_maxId);
    [gpr_fit, ~, AUC, gprMdl] = GPRfitting(gab_dist_vec, gab_vec);% % Plot the Interpolation of it
	[R2, slope, intercept] = Linearfitting(gab_dist_vec, gab_vec);% Linear fitting
	[~,~,AUC_linsmth,smth_Y] = BinSmoothing(gab_dist_vec, gab_vec);
	AUC_raw = trapz_raw(gab_dist_vec, gab_vec);
    % titstr{mi}{end+1} = compose("Gabor: pear %.3f spear %.3f",
	RadTuneStats(Expi).gabor.AUC_linsmth.(metname) = AUC_linsmth; 
	RadTuneStats(Expi).gabor.AUC_raw.(metname) = AUC_raw; 
	RadTuneStats(Expi).gabor.AUC.(metname) = AUC; 
	RadTuneStats(Expi).gabor.gprMdl.(metname) = gprMdl;
	RadTuneStats(Expi).gabor.linR2.(metname) = R2;
	RadTuneStats(Expi).gabor.lincc.(metname) = [slope, intercept];
	[cval,pval] = corr(gab_dist_vec,gab_vec,'Rows','complete');
	RadTuneStats(Expi).gabor.corr.(metname) = cval;
	RadTuneStats(Expi).gabor.corr_P.(metname) = pval;
	[cval,pval] = corr(gab_dist_vec,gab_vec,'Type','Spearman','Rows','complete');
	RadTuneStats(Expi).gabor.corr_sp.(metname) = cval;
	RadTuneStats(Expi).gabor.corr_sp_P.(metname) = pval;
end
end
% saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.pdf",Animal,Expi,Stats(Expi).units.pref_chan)))
end
%%
save(fullfile(mat_dir,Animal+"_ManifRadiTuneStats.mat"),'RadTuneStats')
%%

%% Collect the Hierarchical Stats Struct in Mat into a csv table for summarizing and ploting
global summarydir
summarydir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning\summary";
tab_all = table();
metname = "squ"; % pick the stats you like
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_ManifRadiTuneStats.mat"),'RadTuneStats')
load(fullfile(mat_dir,Animal+"_Evol_stats.mat"),'EStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
tab = struct();
for Expi = 1:numel(RadTuneStats)
tab(Expi).Expi = Expi;
tab(Expi).Animal = Animal;
tab(Expi).pref_chan = EStats(Expi).units.pref_chan;
tab(Expi).pref_unit = EStats(Expi).evol.unit_in_pref_chan;
if EStats(Expi).units.pref_chan < 33
tab(Expi).area = "IT";
elseif EStats(Expi).units.pref_chan < 49
tab(Expi).area = "V1";  
else
tab(Expi).area = "V4"; 
end
tab(Expi).is_evoref_gab = all(contains(EStats(Expi).ref.imgnm,"gab")); % Boolean, if all ref images are gabor.
tab(Expi).peak_evoref = RadTuneStats(Expi).evoref.maxScore;
tab(Expi).normAUC_evoref = RadTuneStats(Expi).evoref.AUC.(metname) / RadTuneStats(Expi).evoref.maxScore;
tab(Expi).AUC_evoref = RadTuneStats(Expi).evoref.AUC.(metname);
tab(Expi).normAUC_raw_evoref = RadTuneStats(Expi).evoref.AUC_raw.(metname) / RadTuneStats(Expi).evoref.maxScore;
tab(Expi).AUC_raw_evoref = RadTuneStats(Expi).evoref.AUC_raw.(metname);
tab(Expi).normAUC_lin_evoref = RadTuneStats(Expi).evoref.AUC_linsmth.(metname) / RadTuneStats(Expi).evoref.maxScore;
tab(Expi).AUC_lin_evoref = RadTuneStats(Expi).evoref.AUC_linsmth.(metname);
tab(Expi).corr_evoref = RadTuneStats(Expi).evoref.corr.(metname);
tab(Expi).corr_P_evoref = RadTuneStats(Expi).evoref.corr_P.(metname);

tab(Expi).peak_mani = RadTuneStats(Expi).manif.maxScore;
tab(Expi).normAUC_mani = RadTuneStats(Expi).manif.AUC.(metname) / RadTuneStats(Expi).manif.maxScore;
tab(Expi).AUC_mani = RadTuneStats(Expi).manif.AUC.(metname);
tab(Expi).normAUC_raw_mani = RadTuneStats(Expi).manif.AUC_raw.(metname) / RadTuneStats(Expi).manif.maxScore;
tab(Expi).AUC_raw_mani = RadTuneStats(Expi).manif.AUC_raw.(metname);
tab(Expi).normAUC_lin_mani = RadTuneStats(Expi).manif.AUC_linsmth.(metname) / RadTuneStats(Expi).manif.maxScore;
tab(Expi).AUC_lin_mani = RadTuneStats(Expi).manif.AUC_linsmth.(metname);
tab(Expi).corr_mani = RadTuneStats(Expi).manif.corr.(metname);
tab(Expi).corr_P_mani = RadTuneStats(Expi).manif.corr_P.(metname);
if ~isempty(RadTuneStats(Expi).pasu)
tab(Expi).peak_pasu = RadTuneStats(Expi).pasu.maxScore;
tab(Expi).normAUC_pasu = RadTuneStats(Expi).pasu.AUC.(metname) / RadTuneStats(Expi).pasu.maxScore;
tab(Expi).AUC_pasu = RadTuneStats(Expi).pasu.AUC.(metname); 
tab(Expi).normAUC_raw_pasu = RadTuneStats(Expi).pasu.AUC_raw.(metname) / RadTuneStats(Expi).pasu.maxScore;
tab(Expi).AUC_raw_pasu = RadTuneStats(Expi).pasu.AUC_raw.(metname); 
tab(Expi).normAUC_lin_pasu = RadTuneStats(Expi).pasu.AUC_linsmth.(metname) / RadTuneStats(Expi).pasu.maxScore;
tab(Expi).AUC_lin_pasu = RadTuneStats(Expi).pasu.AUC_linsmth.(metname); 
tab(Expi).corr_pasu = RadTuneStats(Expi).pasu.corr.(metname); 
tab(Expi).corr_P_pasu = RadTuneStats(Expi).pasu.corr_P.(metname); 
else
for varnm = ["peak_pasu", "normAUC_pasu", "AUC_pasu", "normAUC_raw_pasu", "AUC_raw_pasu", "normAUC_lin_pasu", "AUC_lin_pasu", "corr_pasu", "corr_P_pasu"]
	tab(Expi).(varnm) = nan;
end
end
if ~isempty(RadTuneStats(Expi).gabor)
tab(Expi).peak_gab = RadTuneStats(Expi).gabor.maxScore;
tab(Expi).normAUC_gab = RadTuneStats(Expi).gabor.AUC.(metname) / RadTuneStats(Expi).gabor.maxScore;
tab(Expi).AUC_gab = RadTuneStats(Expi).gabor.AUC.(metname); 
tab(Expi).normAUC_raw_gab = RadTuneStats(Expi).gabor.AUC_raw.(metname) / RadTuneStats(Expi).gabor.maxScore;
tab(Expi).AUC_raw_gab = RadTuneStats(Expi).gabor.AUC_raw.(metname); 
tab(Expi).normAUC_lin_gab = RadTuneStats(Expi).gabor.AUC_linsmth.(metname) / RadTuneStats(Expi).gabor.maxScore;
tab(Expi).AUC_lin_gab = RadTuneStats(Expi).gabor.AUC_linsmth.(metname); 
tab(Expi).corr_gab = RadTuneStats(Expi).gabor.corr.(metname); 
tab(Expi).corr_P_gab = RadTuneStats(Expi).gabor.corr_P.(metname); 
else
for varnm = ["peak_gab", "normAUC_gab", "AUC_gab", "normAUC_raw_gab", "AUC_raw_gab", "normAUC_lin_gab", "AUC_lin_gab", "corr_gab", "corr_P_gab"]
	tab(Expi).(varnm) = nan;
end
end
% iCh = EStats(Expi).units.pref_chan_id; % Note this is not correct! 
% iCh in manifold and in evolution may not match.
% Add F statistics to filter Manifold
iCh = find(Stats(Expi).units.spikeID==EStats(Expi).units.pref_chan &...
           Stats(Expi).units.unit_num_arr==EStats(Expi).evol.unit_in_pref_chan);
si = 1;
actmap_col = cellfun(@(A)...
    squeeze(A(iCh,:)),MapVarStats(Expi).manif.act_col{si},'uni',0);
anovaStats = anova_cells(actmap_col);
tab(Expi).F = anovaStats.F;
tab(Expi).F_P = anovaStats.F_P;
tab(Expi).F_df = anovaStats.STATS.df;
tab(Expi).FSTATS = anovaStats.STATS;
end
tab = struct2table(tab);
writetable(tab, fullfile(summarydir,Animal+"_RadialTuningStatsTab_"+metname+".csv"))
tab_all = [tab_all; tab]; % Cat the table of 2 monkeys together 
end
writetable(tab_all, fullfile(summarydir,"Both"+"_RadialTuningStatsTab_"+metname+".csv"))
%% Add F statistics to filter 
% for ri = 1:size(tab_all,1)
% anovaStats_exp = cellfun(@anova_cells,actmap_col);
% end
%% Load the pre-computed tables. 
summarydir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning\summary";
tab_all = readtable(fullfile(summarydir,"Both"+"_RadialTuningStatsTab_squ.csv"));
tab = tab_all;
%% Plot Summary figure from this 
Animal = "Both";
h=figure; set(h,'pos',[1000, 413, 640, 570]); hold on
scatter(tab.normAUC_mani, tab.corr_mani)
scatter(tab.normAUC_pasu, tab.corr_pasu)
scatter(tab.normAUC_gab, tab.corr_gab)
plot([tab.normAUC_mani,tab.normAUC_pasu]', [tab.corr_mani, tab.corr_pasu]', 'Color', [0,0,0,0.2],'HandleVisibility','off')
plot([tab.normAUC_mani,tab.normAUC_gab]', [tab.corr_mani, tab.corr_gab]', 'Color', [0,0,1,0.2],'HandleVisibility','off')
legend(["GAN Manifold Space", "Pasupathy", "Gabor"])
xlabel("Normalized AUC"); ylabel("Correlation Coefficient") 
title(["Manifold Exps Tuning AUC ~ Correlation "+Animal,"Squeeze Net LPIPS metric"])
savefig(h,fullfile(summarydir, Animal+"_AUC_Corr_scatt_squ.fig")) 
saveas(h,fullfile(summarydir, Animal+"_AUC_Corr_scatt_squ.png"))
saveas(h,fullfile(summarydir, Animal+"_AUC_Corr_scatt_squ.pdf"))
%%
h2=figure; set(h2,'pos',[1000, 413, 640, 570]); hold on
histogram(tab.normAUC_mani, 10, 'FaceAlpha', 0.4)%, tab.corr_mani)
histogram(tab.normAUC_pasu, 10, 'FaceAlpha', 0.4)%, tab.corr_pasu)
histogram(tab.normAUC_gab, 10, 'FaceAlpha', 0.4)%, tab.corr_gab)
% plot([tab.normAUC_mani,tab.normAUC_pasu]', [tab.corr_mani, tab.corr_pasu]', 'Color', [0,0,0,0.3],'HandleVisibility','off')
legend(["GAN Manifold Space", "Pasupathy", "Gabor"])
xlabel("Normalized AUC"); % ylabel("Correlation Coefficient") 
title(["Manifold Exps Density of Tuning AUC "+Animal,"Squeeze Net LPIPS metric"])
savefig(h2,fullfile(summarydir, Animal+"_AUC_hist_squ.fig")) 
saveas(h2,fullfile(summarydir, Animal+"_AUC_hist_squ.png"))
saveas(h2,fullfile(summarydir, Animal+"_AUC_hist_squ.pdf"))
%%
h3 = figure; hold on; set(h3,'pos',[1000,366,438,650]);
jit = randn(size(tab,1),1) * 0.1;
scatter(jit + 1, tab.normAUC_mani);
scatter(jit + 2, tab.normAUC_pasu);
scatter(jit + 3, tab.normAUC_gab);
plot([jit + [1:3]]', [tab.normAUC_mani,tab.normAUC_pasu,tab.normAUC_gab]', 'Color', [0,0,0,0.2]);
ylabel("Normalized AUC");xlabel("Tuning Space")
title(["Manifold Exps Tuning AUC "+Animal,"Squeeze Net LPIPS metric"])
xticks([1,2,3]); xticklabels(["GAN Manifold", "Pasupathy", "Gabor"])
savefig(h3, fullfile(summarydir, Animal+"_AUC_scatt_squ.fig")) 
saveas(h3, fullfile(summarydir, Animal+"_AUC_scatt_squ.png"))
saveas(h3, fullfile(summarydir, Animal+"_AUC_scatt_squ.pdf"))
%%
Animal = "Both";
h4 = figure; hold on; set(h4,'pos',[816   312   622   684]);
Alfamsk = (tab.Animal == "Alfa"); Betomsk = (tab.Animal == "Beto");
spacelist = ["evoref", "mani", "pasu", "gab"];
Cord = colororder;
jit = randn(size(tab,1),1) * 0.1;
for spi = 1:numel(spacelist)
space = spacelist(spi);
scatter(jit(Alfamsk) + spi, tab.("normAUC_"+space)(Alfamsk), 36, Cord(1,:),'*');
scatter(jit(Betomsk) + spi, tab.("normAUC_"+space)(Betomsk), 36, Cord(1,:),'o');
end
plot([jit + [1:4]]', [tab.normAUC_evoref,tab.normAUC_mani,tab.normAUC_pasu,tab.normAUC_gab]', 'Color', [0,0,0,0.2]);
ylabel("Normalized AUC");xlabel("Tuning Space")
title(["Manifold Exps Tuning AUC "+Animal,"Squeeze Net LPIPS metric"])
xticks([1,2,3,4]); xticklabels(["Natural", "GAN Manifold", "Pasupathy", "Gabor"])
savefig(h4, fullfile(summarydir, Animal+"_AUC_scatt_squ_nat.fig")) 
saveas(h4, fullfile(summarydir, Animal+"_AUC_scatt_squ_nat.png"))
saveas(h4, fullfile(summarydir, Animal+"_AUC_scatt_squ_nat.pdf"))
%%
h4 = figure; hold on; set(h4,'pos',[816   312   622   684]);
Alfamsk = (tab.Animal == "Alfa"); Betomsk = (tab.Animal == "Beto");
spacelist = ["evoref", "mani", "pasu", "gab"];
Cord = colororder;
jit = randn(size(tab,1),1) * 0.1;
for spi = 1:numel(spacelist)
space = spacelist(spi);
scatter(jit(Alfamsk) + spi, tab.("peak_"+space)(Alfamsk), 36, Cord(1,:),'*');
scatter(jit(Betomsk) + spi, tab.("peak_"+space)(Betomsk), 36, Cord(1,:),'o');
end
plot([jit + [1:4]]', [tab.peak_evoref,tab.peak_mani,tab.peak_pasu,tab.peak_gab]', 'Color', [0,0,0,0.2]);
ylabel("Peak Activation");xlabel("Tuning Space")
title(["Manifold Exps Peak Firing Rate "+Animal,"Squeeze Net LPIPS metric"])
xticks([1,2,3,4]); xticklabels(["Natural", "GAN Manifold", "Pasupathy", "Gabor"])
savefig(h4, fullfile(summarydir, Animal+"_peak_scatt_squ_nat.fig")) 
saveas(h4, fullfile(summarydir, Animal+"_peak_scatt_squ_nat.png"))
saveas(h4, fullfile(summarydir, Animal+"_peak_scatt_squ_nat.pdf"))
%%
Animal = "Both";
gabrefmsk = logical(tab.is_evoref_gab);
h4 = figure; hold on; set(h4,'pos',[816   312   622   684]);
Alfamsk = (tab.Animal == "Alfa"); Betomsk = (tab.Animal == "Beto");
spacelist = ["evoref", "mani", "pasu", "gab"];
Cord = colororder;
jit = randn(size(tab,1),1) * 0.1;
for spi = 1:numel(spacelist)
space = spacelist(spi);
if space == "evoref", 
scatter(jit(Alfamsk & gabrefmsk) + spi, tab.("normAUC_"+space)(Alfamsk & gabrefmsk), 10, Cord(1,:),'*');
scatter(jit(Betomsk & gabrefmsk) + spi, tab.("normAUC_"+space)(Betomsk & gabrefmsk), 10, Cord(1,:),'o');
scatter(jit(Alfamsk & ~gabrefmsk) + spi, tab.("normAUC_"+space)(Alfamsk & ~gabrefmsk), 36, Cord(1,:),'*');
scatter(jit(Betomsk & ~gabrefmsk) + spi, tab.("normAUC_"+space)(Betomsk & ~gabrefmsk), 36, Cord(1,:),'o');
else
scatter(jit(Alfamsk) + spi, tab.("normAUC_"+space)(Alfamsk), 36, Cord(1,:),'*');
scatter(jit(Betomsk) + spi, tab.("normAUC_"+space)(Betomsk), 36, Cord(1,:),'o');
end
end
plot([jit + [1:4]]', [tab.normAUC_evoref,tab.normAUC_mani,tab.normAUC_pasu,tab.normAUC_gab]', 'Color', [0,0,0,0.2]);
ylabel("Normalized AUC");xlabel("Tuning Space")
title(["Manifold Exps Tuning AUC "+Animal,"Squeeze Net LPIPS metric","* Alfa o Beto, small symbol gabor reference"])
xticks([1,2,3,4]); xticklabels(["Natural", "GAN Manifold", "Pasupathy", "Gabor"])
savefig(h4, fullfile(summarydir, Animal+"_AUC_scatt_squ_nat_gab.fig")) 
saveas(h4, fullfile(summarydir, Animal+"_AUC_scatt_squ_nat_gab.png"))
saveas(h4, fullfile(summarydir, Animal+"_AUC_scatt_squ_nat_gab.pdf"))
%%
h4 = figure; hold on; set(h4,'pos',[816   312   622   684]);
gabrefmsk = logical(tab.is_evoref_gab);
Alfamsk = (tab.Animal == "Alfa"); Betomsk = (tab.Animal == "Beto");
spacelist = ["evoref", "mani", "pasu", "gab"];
Cord = colororder;
jit = randn(size(tab,1),1) * 0.1;
for spi = 1:numel(spacelist)
space = spacelist(spi);
if space == "evoref", 
scatter(jit(Alfamsk & gabrefmsk) + spi, tab.("peak_"+space)(Alfamsk & gabrefmsk), 10, Cord(1,:),'*');
scatter(jit(Betomsk & gabrefmsk) + spi, tab.("peak_"+space)(Betomsk & gabrefmsk), 10, Cord(1,:),'o');
scatter(jit(Alfamsk & ~gabrefmsk) + spi, tab.("peak_"+space)(Alfamsk & ~gabrefmsk), 36, Cord(1,:),'*');
scatter(jit(Betomsk & ~gabrefmsk) + spi, tab.("peak_"+space)(Betomsk & ~gabrefmsk), 36, Cord(1,:),'o');
else
scatter(jit(Alfamsk) + spi, tab.("peak_"+space)(Alfamsk), 36, Cord(1,:),'*');
scatter(jit(Betomsk) + spi, tab.("peak_"+space)(Betomsk), 36, Cord(1,:),'o');
end
end
plot([jit + [1:4]]', [tab.peak_evoref,tab.peak_mani,tab.peak_pasu,tab.peak_gab]', 'Color', [0,0,0,0.2]);
ylabel("Peak Activation");xlabel("Tuning Space")
title(["Manifold Exps Peak Firing Rate "+Animal,"Squeeze Net LPIPS metric","* Alfa o Beto, small symbol gabor reference"])
xticks([1,2,3,4]); xticklabels(["Natural", "GAN Manifold", "Pasupathy", "Gabor"])
savefig(h4, fullfile(summarydir, Animal+"_peak_scatt_squ_nat_gab.fig")) 
saveas(h4, fullfile(summarydir, Animal+"_peak_scatt_squ_nat_gab.png"))
saveas(h4, fullfile(summarydir, Animal+"_peak_scatt_squ_nat_gab.pdf"))
%%
tab_nog = tab;
gabrefmsk = logical(tab.is_evoref_gab);
for varname = ["peak_evoref", "normAUC_evoref", "AUC_evoref", "corr_evoref"]
	tab_nog.(varname)(gabrefmsk) = nan;
end
%%
Cord = colororder;
plotJitScatterVar("AUC_", "UnNormalized Area Under Curve", {'m', Cord(1,:), Cord(3,:), 'g'}, tab, summarydir)
plotJitScatterVar("normAUC_", "Normalized Area Under Curve", {'m', Cord(1,:), Cord(3,:), 'g'}, tab, summarydir)
%%
plotJitScatterVar("AUC_lin_", "UnNormalized Area Under Curve(lin smooth)", {'m', Cord(1,:), Cord(3,:), 'g'}, tab, summarydir)
plotJitScatterVar("normAUC_lin_", "Normalized Area Under Curve(lin smooth)", {'m', Cord(1,:), Cord(3,:), 'g'}, tab, summarydir)
plotJitScatterVar("AUC_raw_", "UnNormalized Area Under Curve(raw trapz)", {'m', Cord(1,:), Cord(3,:), 'g'}, tab, summarydir)
plotJitScatterVar("normAUC_raw_", "Normalized Area Under Curve(raw trapz)", {'m', Cord(1,:), Cord(3,:), 'g'}, tab, summarydir)

%%

msk = (tab.F_P<0.001) & ~((tab.Animal=="Alfa") & (tab.Expi==10)) ; % maybe add F criterion to get rid of flat curves....
V1msk = tab.area == "V1";
V4msk = tab.area == "V4";
ITmsk = tab.area == "IT";
Alfamsk = (tab.Animal == "Alfa"); 
Betomsk = (tab.Animal == "Beto");
h = stripe_plot(tab_all, "normAUC_mani", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver PC23)", "anim_sep");
h = stripe_plot(tab_all, "normAUC_mani", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver PC23)", "area_sep",{[3,1],[2,1],[3,2]});
h = stripe_minor_plot(tab_all, "normAUC_mani", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], "All Exp (driver PC23)", "area_anim", {[3,1],[2,1],[3,2]});
%%
Fmsk = (tab.F_P<0.001) & msk; 
h = stripe_plot(tab_all, "normAUC_mani", {Fmsk&Alfamsk, Fmsk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver PC23 P<0.001)", "anim_sep_Fsig");
h = stripe_plot(tab_all, "normAUC_mani", {Fmsk&V1msk, Fmsk&V4msk, Fmsk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver PC23 P<0.001)", "area_sep_Fsig",{[3,1],[2,1],[3,2]});
h = stripe_minor_plot(tab_all, "normAUC_mani", {Fmsk&V1msk, Fmsk&V4msk, Fmsk&ITmsk}, ["V1","V4","IT"],...
    {Fmsk&Alfamsk, Fmsk&Betomsk}, ["Alfa","Beto"], "All Exp (driver PC23 P<0.001)", "area_anim_Fsig", {[3,1],[2,1],[3,2]});
%%

%%
h = stripe_plot(tab_all, "normAUC_evoref", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver)", "_anim_sep");
h = stripe_plot(tab_all, "normAUC_evoref", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver)", "_area_sep",{[3,1],[2,1],[3,2]});
h = stripe_plot(tab_all, "normAUC_pasu", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver)", "_anim_sep");
h = stripe_plot(tab_all, "normAUC_pasu", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver)", "_area_sep",{[3,1],[2,1],[3,2]});
h = stripe_plot(tab_all, "normAUC_gab", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver)", "_anim_sep");
h = stripe_plot(tab_all, "normAUC_gab", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver)", "_area_sep",{[3,1],[2,1],[3,2]});

%% Plotting routines (some borrowed from Manif_Fit_summary)
function plotJitScatterVar(varprefix, label, colors, tab, summarydir)
Animal = "Both";
h4 = figure; hold on; set(h4,'pos',[816   312   622   684]);
gabrefmsk = logical(tab.is_evoref_gab);
Alfamsk = (tab.Animal == "Alfa"); Betomsk = (tab.Animal == "Beto");
spacelist = ["evoref", "mani", "pasu", "gab"];
Cord = colororder;
jit = randn(size(tab,1),1) * 0.1;
for spi = 1:numel(spacelist)
space = spacelist(spi);
if space == "evoref"
scatter(jit(Alfamsk & gabrefmsk) + spi, tab.(varprefix+space)(Alfamsk & gabrefmsk), 10, colors{spi},'*');
scatter(jit(Betomsk & gabrefmsk) + spi, tab.(varprefix+space)(Betomsk & gabrefmsk), 10, colors{spi},'o');
scatter(jit(Alfamsk & ~gabrefmsk) + spi, tab.(varprefix+space)(Alfamsk & ~gabrefmsk), 36, colors{spi},'*');
scatter(jit(Betomsk & ~gabrefmsk) + spi, tab.(varprefix+space)(Betomsk & ~gabrefmsk), 36, colors{spi},'o');
else
scatter(jit(Alfamsk) + spi, tab.(varprefix+space)(Alfamsk), 36, colors{spi},'*');
scatter(jit(Betomsk) + spi, tab.(varprefix+space)(Betomsk), 36, colors{spi},'o');
end
end
plot([jit + [1:4]]', [tab.(varprefix+"evoref"),tab.(varprefix+"mani"),tab.(varprefix+"pasu"),tab.(varprefix+"gab")]', 'Color', [0,0,0,0.2]);
ylabel(label);xlabel("Tuning Space")
title(["Manifold Exps "+label+" "+Animal,"Squeeze Net LPIPS metric","* Alfa o Beto, small symbol gabor reference"])
xticks([1,2,3,4]); xticklabels(["Natural", "GAN Manifold", "Pasupathy", "Gabor"])
savefnm = compose("%s_%sscatt_squ_nat_gab",Animal,varprefix);
savefig(h4, fullfile(summarydir, savefnm+".fig")) 
saveas(h4, fullfile(summarydir, savefnm+".png"))
saveas(h4, fullfile(summarydir, savefnm+".pdf"))
end

function h = stripe_plot(tab, statname, masks, labels, titstr, savestr, Tpairs)
if nargin<7, Tpairs = {}; end
global summarydir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
statscol = cellfun(@(M)tab.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
xjitter = 0.15*randn(numel(statscol{i}),1);
scatter(i+xjitter, statscol{i})
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
% xlabel(statname) % "d'( Dirichlet Energy vs Shuffling Ctrl )"
title(title_str)
legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
saveas(h,fullfile(summarydir,compose("%s_%s_stripcmp.png", statname, savestr)))
saveas(h,fullfile(summarydir,compose("%s_%s_stripcmp.pdf", statname, savestr)))
savefig(h,fullfile(summarydir,compose("%s_%s_stripcmp.fig", statname, savestr)))
end

function h = stripe_minor_plot(tab, statname, masks, labels, minormsks, minorlabels, titstr, savestr, Tpairs)
if nargin<7, Tpairs = {}; end
marker = 'o*x^v';
global summarydir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
labelcol = []; Ns = [];
for i = 1:numel(masks)
    for j = 1:numel(minormsks)
    M = masks{i}&minormsks{j};
    xjitter = 0.15*randn(sum(M),1);
    scatter(i+xjitter, tab.(statname)(M),'Marker',marker(j))
    labelcol = [labelcol, labels(i)+minorlabels(j)];
    Ns = [Ns, sum(M)];
    end
end
xticks(1:numel(masks));xticklabels(labels)
ylabel(statname)
statscol = cellfun(@(M)tab.(statname)(M), masks,'uni',0);
FStat = anova_cells(statscol);
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
legend(compose("%s(%d)",labelcol',Ns'),'Location','best') % ["Driver", "Non-Drivers"]
saveas(h,fullfile(summarydir,compose("%s_%s_stripcmpmarker.png", statname, savestr)))
saveas(h,fullfile(summarydir,compose("%s_%s_stripcmpmarker.pdf", statname, savestr)))
savefig(h,fullfile(summarydir,compose("%s_%s_stripcmpmarker.fig", statname, savestr)))
end

%% SummaryStats 
%% compare_summary: function description
function [statstr] = xspace_cmp_summary(arg)
	statstr = "";
end

%% Util functions to compute Gaussian Process fitting and linear fitting and compute AUC
% Linearfitting(1:10,-[1:10]+0.1*randn(1,10))
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