%% Manif_Imgs_Dissim_Merge
Animal = "Beto";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir,"gab_imdist.mat"),'gab_imdist')
load(fullfile(mat_dir,"pasu_imdist.mat"),'pasu_imdist')
load(fullfile(mat_dir,Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
%%
% pasu_nm_grid = cellfun(@(idx)unique(Stats(12).imageName(idx)),Stats(12).ref.pasu_idx_grid,'Uni',0);
pasu_idx_vec = reshape(Stats(12).ref.pasu_idx_grid',[],1); % reshape into one row. But transpose to make similar object adjacent
pasu_val_msk = ~cellfun(@isempty,pasu_idx_vec);
%%
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
%%
metric_list = ["squ","SSIM","L2","FC6"]

%%
saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.png",Animal,Expi,Stats(Expi).units.pref_chan)))
saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.pdf",Animal,Expi,Stats(Expi).units.pref_chan)))
savefig(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.fig",Animal,Expi,Stats(Expi).units.pref_chan)))

%%
function plotGPRfitting(X, Y, vargin)
if nargin == 2, vargin={};end
gprMdl = fitrgp(X, Y);
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit, vargin{:})
end
