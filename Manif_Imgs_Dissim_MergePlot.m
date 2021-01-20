%% Manif_Imgs_Dissim_Merge
%  This code is to plot the radial tuning curves for several different spaces together 
% Load the relevant statistics
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
% for Animal = ["Beto", "Alfa"]%
Animal = "Beto";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir, "gab_imdist.mat"),'gab_imdist')
load(fullfile(mat_dir, "pasu_imdist.mat"),'pasu_imdist')
load(fullfile(mat_dir, Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
load(fullfile(mat_dir, Animal+"_EvoRef_ImDist.mat"),"EvoRefImDistStat")
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%% Get the invalid mask for Pasupathy images.
% pasu_nm_grid = cellfun(@(idx)unique(Stats(12).imageName(idx)),Stats(12).ref.pasu_idx_grid,'Uni',0);
if Animal == "Alfa"
pasu_idx_vec = reshape(Stats(1).ref.pasu_idx_grid',[],1); % 12% reshape into one row. But transpose to make similar object adjacent
elseif Animal == "Beto"
pasu_idx_vec = reshape(Stats(12).ref.pasu_idx_grid',[],1);
end
pasu_val_msk = ~cellfun(@isempty,pasu_idx_vec);

%% newer compact version of code
flag.doEvoRef = 1;
flag.doError = 1;
bestcorr = 0;
Cord = colororder;
metric_list = ["squ","SSIM","L2","FC6","FC6_corr"];
label_list = ["LPIPS (SqueezeNet)", "SSIM", "L2", "FC6 (L2)", "FC6 (1 - corr)"];
for Expi=11%1:numel(Stats)
titstr = cell(1,numel(metric_list));
figure(21);
T=tiledlayout(1,numel(metric_list),'TileSpacing','compact','Padding','compact');
title(T,compose("%s Exp %d prefchan %d",Animal,Expi,Stats(Expi).units.pref_chan),'FontSize',16)
bsl_VEC_ALL = [];
% Manifold images. 
si=1; ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_sem_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all')/sqrt(size(psth,3)),reshape(Stats(Expi).manif.psth{si},[],1));
score_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).manif.psth{si},[],1),'uni',0)); % single trial
bsl_VEC_ALL = [bsl_VEC_ALL; score_bsl_VEC];
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
for mi = 1:numel(metric_list) % Loop through different distance metrics
	nexttile(mi);  metname = metric_list(mi); hold on
	distmat = ManifImDistStat(Expi).(metname);
    scatter(distmat(:,maxId), score_vec, 36, Cord(1,:))
    if flag.doError, errorbar(distmat(:,maxId), score_vec, score_sem_vec,'Color',Cord(1,:),'LineStyle','none','LineWidth',0.25); end 
	plotGPRfitting(distmat(:,maxId), score_vec, mean(score_std_vec), {'Color',Cord(1,:)})% Gaussian Process Smoothing or Fitting
	[cc_s,pval_s] = corr(distmat(:,maxId),score_vec,'Type','Spearman');
	[cc_p,pval_p] = corr(distmat(:,maxId),score_vec);
	titstr{mi} = { compose("Manif: pear %.3f (%.1e) spear %.3f (%.1e)",cc_p,pval_p,cc_s,pval_s)};
	xlabel(label_list(mi));
end	

if flag.doEvoRef
evoref_vec = cellfun(@(psth)mean(psth(1,51:200,:),'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_std_vec = cellfun(@(psth)std(mean(psth(1,51:200,:),[1,2]),1,'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_sem_vec = cellfun(@(psth)std(mean(psth(1,51:200,:),[1,2]),1,'all')/sqrt(size(psth,3)),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(1,1:45,:),[2])),reshape(EStats(Expi).ref.psth_arr,[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; evoref_bsl_VEC];
[evoref_sortScore,sortId] = sort(evoref_vec,'Descend');
[evoref_maxScore,evoref_maxId] = max(evoref_vec);
for mi = 1:numel(metric_list)
    nexttile(mi);  metname = metric_list(mi);
	distmat = EvoRefImDistStat(Expi).(metname);
    scatter(distmat(:,evoref_maxId), evoref_vec,'m')
    if flag.doError, errorbar(distmat(:,evoref_maxId), evoref_vec,evoref_sem_vec,'m','LineStyle','none','LineWidth',0.25); end 
    plotGPRfitting(distmat(:,evoref_maxId), evoref_vec, mean(evoref_std_vec), {'m'})% % Plot the Interpolation of it
    [cc_s,pval_s] = corr(distmat(:,evoref_maxId),evoref_vec,'Type','Spearman','Rows','complete');
	[cc_p,pval_p] = corr(distmat(:,evoref_maxId),evoref_vec,'Rows','complete');    	
    titstr{mi}{end+1} = compose("Natural Ref: pear %.3f (%.1e) spear %.3f (%.1e)",cc_p,pval_p,cc_s,pval_s);
end
legend(["Manifold","Manifold","NatRef","NatRef","Pasupathy","Pasupathy","Gabor","Gabor"])
end
% Pasupathy patches
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(~pasu_val_msk) = []; % isnan(pasu_vec) % get rid of non-existing pasupathy images. 
pasu_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_sem_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all')/sqrt(size(psth,3)),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_std_vec(~pasu_val_msk) = [];
pasu_sem_vec(~pasu_val_msk) = [];
pasu_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).ref.pasu_psths',[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; pasu_bsl_VEC];
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
for mi = 1:numel(metric_list)
	nexttile(mi);  metname = metric_list(mi);
    scatter(pasu_imdist.(metname)(:,pasu_maxId),pasu_vec,36,Cord(3,:))
    if flag.doError, errorbar(pasu_imdist.(metname)(:,pasu_maxId),pasu_vec,pasu_sem_vec,'Color',Cord(3,:),'LineStyle','none','LineWidth',0.25); end 
	plotGPRfitting(pasu_imdist.(metname)(:,pasu_maxId), pasu_vec, nanmean(pasu_std_vec), {'Color',Cord(3,:)})
	[cc_s,pval_s] = corr(pasu_imdist.(metname)(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete');
	[cc_p,pval_p] = corr(pasu_imdist.(metname)(:,pasu_maxId),pasu_vec,'Rows','complete');		
	titstr{mi}{end+1} = compose("Pasu: pear %.3f (%.1e) spear %.3f (%.1e)",cc_p,pval_p,cc_s,pval_s);
end	
legend(["Manifold","Manifold","NatRef","NatRef","Pasupathy","Pasupathy"])
end
% Gabor patches
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
gab_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
gab_sem_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all')/sqrt(size(psth,3)),reshape(Stats(Expi).ref.gab_psths,[],1));
gab_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).ref.gab_psths,[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; gab_bsl_VEC];
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
for mi = 1:numel(metric_list)
    nexttile(mi);  metname = metric_list(mi);
    scatter(gab_imdist.(metname)(:,gab_maxId), gab_vec,'g')
    if flag.doError, errorbar(gab_imdist.(metname)(:,gab_maxId), gab_vec,gab_sem_vec,'g','LineStyle','none','LineWidth',0.25); end 
    plotGPRfitting(gab_imdist.(metname)(:,gab_maxId), gab_vec, nanmean(gab_std_vec), {'g'})% % Plot the Interpolation of it
    [cc_s,pval_s] = corr(gab_imdist.(metname)(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete');
	[cc_p,pval_p] = corr(gab_imdist.(metname)(:,gab_maxId),gab_vec,'Rows','complete');    	
    titstr{mi}{end+1} = compose("Gabor: pear %.3f (%.1e) spear %.3f (%.1e)",cc_p,pval_p,cc_s,pval_s);
end
legend(["Manifold","Manifold","NatRef","NatRef","Pasupathy","Pasupathy","Gabor","Gabor"])
end

bsl_rate = nanmean(bsl_VEC_ALL);
bsl_rate_std = nanstd(bsl_VEC_ALL,1);
bsl_rate_sem = sem(bsl_VEC_ALL,1);
for mi = 1:numel(metric_list) % add title and baseline firing rate 
nexttile(mi);title(titstr{mi})
hline(bsl_rate,'k-');hline(bsl_rate+[-bsl_rate_sem,bsl_rate_sem],'-.')
YLIM = ylim();ylim(max(0,YLIM)); % no need to show negative firing rate. 
end
nexttile(1); ylabel("Neural Activation")
% Saving part. 
savefnm = compose("%s_Exp%02d_pref%02d_peak",Animal,Expi,Stats(Expi).units.pref_chan); 
if flag.doEvoRef, savefnm = savefnm + "_nat";end
if flag.doError, savefnm = savefnm + "_err";end
saveas(21,fullfile(figdir, savefnm + ".png"))
savefig(21,fullfile(figdir, savefnm + ".fig"))
saveas(21,fullfile(figdir, savefnm + ".pdf"))
end
% end
%%

%%
function plotGPRfitting(X, Y, Sigma, vargin)
if nargin == 3, vargin={};end
if nargin == 2, Sigma=std(Y)/sqrt(2);end
gprMdl = fitrgp(X, Y, 'SigmaLowerBound', 5E-1,'Standardize',true,'Sigma',Sigma); % This lower bound is to fix a problem! 
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit, vargin{:}) %,'HandleVisibility','off' % adding these arguments will remove it from legend
% Spline interpolation is not good....
% ypred = makima(X,Y,xlinsp);
% plot(xlinsp, ypred, vargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend

% Smoothing spline can be really bad. 
% options = fitoptions('Method','Smooth','SmoothingParam',0.8);
% [f,gof,out] = fit(X,Y,'smooth',options);
% ypred = feval(f,xlinsp);
% plot(xlinsp, ypred, vargin{:})
end

%%% Obsolete full version, more readable but verbose. 
% figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
% Expi = 35;
% bestcorr = 0;
% Cord = colororder;
% for Expi=11:12%34:numel(Stats)
% figure(21);
% T=tiledlayout(1,4,'TileSpacing','compact');
% title(T,compose("%s Exp %d prefchan %d",Animal,Expi,Stats(Expi).units.pref_chan),'FontSize',16)
% % Manifold images. 
% si=1; ui=1;
% score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
% [sortScore,sortId]=sort(score_vec,'Descend');
% [maxScore,maxId]=max(score_vec);
% nexttile(1);hold on
% scatter(ManifImDistStat(Expi).squ(:,maxId), score_vec)
% plotGPRfitting(ManifImDistStat(Expi).squ(:,maxId), score_vec)% Gaussian Process Smoothing or Fitting
% xlabel("Perceptual Similarity(SqueezeNet)");ylabel("Neural Activation")
% titstr1 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).squ(:,maxId),score_vec),corr(ManifImDistStat(Expi).squ(:,maxId),score_vec,'Type','Spearman'))};
% nexttile(2);hold on
% scatter(ManifImDistStat(Expi).SSIM(:,maxId),score_vec)
% plotGPRfitting(ManifImDistStat(Expi).SSIM(:,maxId),score_vec)
% xlabel("SSIM")
% nexttile(3);hold on
% scatter(ManifImDistStat(Expi).L2(:,maxId),score_vec)
% plotGPRfitting(ManifImDistStat(Expi).L2(:,maxId),score_vec)
% xlabel("L2")
% nexttile(4);hold on
% scatter(ManifImDistStat(Expi).FC6(:,maxId),score_vec)
% plotGPRfitting(ManifImDistStat(Expi).FC6(:,maxId),score_vec)
% titstr4 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec),corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec,'Type','Spearman'))};
% xlabel("FC6 (1-corr)")
% % Pasupathy patches
% if Stats(Expi).ref.didPasu
% pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
% pasu_vec(~pasu_val_msk) = []; % isnan(pasu_vec)
% [pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
% [pasu_maxScore,pasu_maxId]=max(pasu_vec);
% nexttile(1);
% scatter(pasu_imdist.squ(:,pasu_maxId),pasu_vec)
% plotGPRfitting(pasu_imdist.squ(:,pasu_maxId), pasu_vec)
% titstr1{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Rows','complete'),corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
% nexttile(2);
% scatter(pasu_imdist.SSIM(:,pasu_maxId),pasu_vec)
% plotGPRfitting(pasu_imdist.SSIM(:,pasu_maxId),pasu_vec)
% nexttile(3);
% scatter(pasu_imdist.L2(:,pasu_maxId),pasu_vec)
% plotGPRfitting(pasu_imdist.L2(:,pasu_maxId),pasu_vec)
% nexttile(4);
% scatter(pasu_imdist.FC6(:,pasu_maxId),pasu_vec)
% plotGPRfitting(pasu_imdist.FC6(:,pasu_maxId),pasu_vec)
% legend(["Manifold","Manifold","Pasupathy","Pasupathy"])
% titstr4{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Rows','complete'),corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
% end
% % Gabor patches
% if Stats(Expi).ref.didGabor
% gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
% [gab_sortScore,sortId]=sort(gab_vec,'Descend');
% [gab_maxScore,gab_maxId]=max(gab_vec);
% nexttile(1);hold on
% scatter(gab_imdist.squ(:,gab_maxId), gab_vec,'g')
% plotGPRfitting(gab_imdist.squ(:,gab_maxId), gab_vec,{'g'})% % Plot the Interpolation of it
% titstr1{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Rows','complete'),corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
% nexttile(2);hold on
% scatter(gab_imdist.SSIM(:,gab_maxId),gab_vec,'g')
% plotGPRfitting(gab_imdist.SSIM(:,gab_maxId),gab_vec,{'g'})
% nexttile(3);hold on
% scatter(gab_imdist.L2(:,gab_maxId),gab_vec,'g')
% plotGPRfitting(gab_imdist.L2(:,gab_maxId),gab_vec,{'g'})
% nexttile(4);hold on
% scatter(gab_imdist.FC6(:,gab_maxId),gab_vec,'g')
% plotGPRfitting(gab_imdist.FC6(:,gab_maxId),gab_vec,{'g'})
% legend(["Manifold","Manifold","Pasupathy","Pasupathy","Gabor","Gabor"])
% titstr4{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Rows','complete'),corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
% end
% nexttile(1);title(titstr1)
% nexttile(4);title(titstr4)
% % saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.png",Animal,Expi,Stats(Expi).units.pref_chan)))
% % savefig(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.fig",Animal,Expi,Stats(Expi).units.pref_chan)))
% end