%% Manif_Imgs_Dissim_Merge
% Load the relevant statistics
Animal = "Alfa";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir,"gab_imdist.mat"),'gab_imdist')
load(fullfile(mat_dir,"pasu_imdist.mat"),'pasu_imdist')
load(fullfile(mat_dir,Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
%% Get the invalid mask for Pasupathy images.
% pasu_nm_grid = cellfun(@(idx)unique(Stats(12).imageName(idx)),Stats(12).ref.pasu_idx_grid,'Uni',0);
if Animal == "Alfa"
pasu_idx_vec = reshape(Stats(1).ref.pasu_idx_grid',[],1); % 12% reshape into one row. But transpose to make similar object adjacent
elseif Animal == "Beto"
pasu_idx_vec = reshape(Stats(12).ref.pasu_idx_grid',[],1);
end
pasu_val_msk = ~cellfun(@isempty,pasu_idx_vec);
%% Lengthy full version, more readable. 
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
Expi = 35;
bestcorr = 0;
Cord = colororder;
for Expi=11:12%34:numel(Stats)
figure(21);
T=tiledlayout(1,4,'TileSpacing','compact');
title(T,compose("%s Exp %d prefchan %d",Animal,Expi,Stats(Expi).units.pref_chan),'FontSize',16)
% Manifold images. 
si=1; ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
nexttile(1);hold on
scatter(ManifImDistStat(Expi).squ(:,maxId), score_vec)
plotGPRfitting(ManifImDistStat(Expi).squ(:,maxId), score_vec)% Gaussian Process Smoothing or Fitting
xlabel("Perceptual Similarity(SqueezeNet)");ylabel("Neural Activation")
titstr1 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).squ(:,maxId),score_vec),corr(ManifImDistStat(Expi).squ(:,maxId),score_vec,'Type','Spearman'))};
nexttile(2);hold on
scatter(ManifImDistStat(Expi).SSIM(:,maxId),score_vec)
plotGPRfitting(ManifImDistStat(Expi).SSIM(:,maxId),score_vec)
xlabel("SSIM")
nexttile(3);hold on
scatter(ManifImDistStat(Expi).L2(:,maxId),score_vec)
plotGPRfitting(ManifImDistStat(Expi).L2(:,maxId),score_vec)
xlabel("L2")
nexttile(4);hold on
scatter(ManifImDistStat(Expi).FC6(:,maxId),score_vec)
plotGPRfitting(ManifImDistStat(Expi).FC6(:,maxId),score_vec)
titstr4 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec),corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec,'Type','Spearman'))};
xlabel("FC6 (1-corr)")
% Pasupathy patches
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(~pasu_val_msk) = []; % isnan(pasu_vec)
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
nexttile(1);
scatter(pasu_imdist.squ(:,pasu_maxId),pasu_vec)
plotGPRfitting(pasu_imdist.squ(:,pasu_maxId), pasu_vec)
titstr1{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Rows','complete'),corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
nexttile(2);
scatter(pasu_imdist.SSIM(:,pasu_maxId),pasu_vec)
plotGPRfitting(pasu_imdist.SSIM(:,pasu_maxId),pasu_vec)
nexttile(3);
scatter(pasu_imdist.L2(:,pasu_maxId),pasu_vec)
plotGPRfitting(pasu_imdist.L2(:,pasu_maxId),pasu_vec)
nexttile(4);
scatter(pasu_imdist.FC6(:,pasu_maxId),pasu_vec)
plotGPRfitting(pasu_imdist.FC6(:,pasu_maxId),pasu_vec)
legend(["Manifold","Manifold","Pasupathy","Pasupathy"])
titstr4{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Rows','complete'),corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
end
% Gabor patches
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
nexttile(1);hold on
scatter(gab_imdist.squ(:,gab_maxId), gab_vec,'g')
plotGPRfitting(gab_imdist.squ(:,gab_maxId), gab_vec,{'g'})% % Plot the Interpolation of it
titstr1{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Rows','complete'),corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
nexttile(2);hold on
scatter(gab_imdist.SSIM(:,gab_maxId),gab_vec,'g')
plotGPRfitting(gab_imdist.SSIM(:,gab_maxId),gab_vec,{'g'})
nexttile(3);hold on
scatter(gab_imdist.L2(:,gab_maxId),gab_vec,'g')
plotGPRfitting(gab_imdist.L2(:,gab_maxId),gab_vec,{'g'})
nexttile(4);hold on
scatter(gab_imdist.FC6(:,gab_maxId),gab_vec,'g')
plotGPRfitting(gab_imdist.FC6(:,gab_maxId),gab_vec,{'g'})
legend(["Manifold","Manifold","Pasupathy","Pasupathy","Gabor","Gabor"])
titstr4{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Rows','complete'),corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
end
nexttile(1);title(titstr1)
nexttile(4);title(titstr4)
% saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.png",Animal,Expi,Stats(Expi).units.pref_chan)))
% savefig(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.fig",Animal,Expi,Stats(Expi).units.pref_chan)))
end


%% newer compact version of code.
bestcorr = 0;
Cord = colororder;
metric_list = ["squ","SSIM","L2","FC6","FC6_corr"];
label_list = ["LPIPS (SqueezeNet)", "SSIM", "L2", "FC6 (L2)", "FC6 (1 - corr)"];
for Expi=42:42%numel(Stats)
titstr = cell(1,numel(metric_list));
figure(21);
T=tiledlayout(1,numel(metric_list),'TileSpacing','compact','Padding','compact');
title(T,compose("%s Exp %d prefchan %d",Animal,Expi,Stats(Expi).units.pref_chan),'FontSize',16)
% Manifold images. 
si=1; ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
for mi = 1:numel(metric_list) % Loop through different distance metrics
	nexttile(mi);  metname = metric_list(mi); hold on
    scatter(ManifImDistStat(Expi).(metname)(:,maxId), score_vec)
	plotGPRfitting(ManifImDistStat(Expi).(metname)(:,maxId), score_vec)% Gaussian Process Smoothing or Fitting
	titstr{mi} = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).(metname)(:,maxId),score_vec),...
				corr(ManifImDistStat(Expi).(metname)(:,maxId),score_vec,'Type','Spearman'))};
	xlabel(label_list(mi));
end	
% Pasupathy patches
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(~pasu_val_msk) = []; % isnan(pasu_vec) % get rid of non-existing pasupathy images. 
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
for mi = 1:numel(metric_list)
	nexttile(mi);  metname = metric_list(mi);
    scatter(pasu_imdist.(metname)(:,pasu_maxId),pasu_vec)
	plotGPRfitting(pasu_imdist.(metname)(:,pasu_maxId), pasu_vec)
	titstr{mi}{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.(metname)(:,pasu_maxId),pasu_vec,'Rows','complete'),...
				corr(pasu_imdist.(metname)(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
end	
legend(["Manifold","Manifold","Pasupathy","Pasupathy"])
end
% Gabor patches
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
for mi = 1:numel(metric_list)
    nexttile(mi);  metname = metric_list(mi);
    scatter(gab_imdist.(metname)(:,gab_maxId), gab_vec,'g')
    plotGPRfitting(gab_imdist.(metname)(:,gab_maxId), gab_vec,{'g'})% % Plot the Interpolation of it
    titstr{mi}{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.(metname)(:,gab_maxId),gab_vec,'Rows','complete'),...
                                    corr(gab_imdist.(metname)(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
end
legend(["Manifold","Manifold","Pasupathy","Pasupathy","Gabor","Gabor"])
end
for mi = 1:numel(metric_list)
nexttile(mi);title(titstr{mi})
end
nexttile(1); ylabel("Neural Activation")
saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.png",Animal,Expi,Stats(Expi).units.pref_chan)))
savefig(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.fig",Animal,Expi,Stats(Expi).units.pref_chan)))
saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.pdf",Animal,Expi,Stats(Expi).units.pref_chan)))
end
%%



%%
function plotGPRfitting(X, Y, vargin)
if nargin == 2, vargin={};end
gprMdl = fitrgp(X, Y, 'SigmaLowerBound', 1E-3); % This lower bound is to fix a problem! 
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit, vargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend
end
