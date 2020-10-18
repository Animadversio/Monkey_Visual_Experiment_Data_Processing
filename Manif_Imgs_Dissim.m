%% Manif Image Dissimilarity
%  First half of the code is on computing the image dissimilarity of
%  different spaces. Gabor Pasupathy, Manifold 
%  How neural firing change w.r.t. image dissimilarity
py.torch.cuda.empty_cache()
D = torchImDist();
%%
net = alexnet;
%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%% Pasu Distance metric
pasu_path = "N:\Stimuli\2019-Manifold\pasupathy-wg-f-4-ori";
imlist = dir(fullfile(pasu_path,"*.jpg"));
imlist = struct2cell(imlist);
imgnms = string(cellstr(imlist(1,:)))';
imgi = 1;
img = imread(fullfile(pasu_path, imgnms{imgi}));
imgs_all = arrayfun(@(nm)imread(fullfile(pasu_path, nm)),imgnms,'Uni',0);
rsz_imgs_all = cellfun(@(img)repmat(imresize(img,[256,256]),1,1,3),imgs_all,'Uni',0);
rsz_imgs_tsr = cell2mat(shiftdim(rsz_imgs_all,-3));
%%
pasu_imdist = struct();
%%
tic
D = D.select_metric("squeeze");
pasu_imdist.squ = D.distmat_B(rsz_imgs_tsr);
toc
L2dist = pdist(reshape(rsz_imgs_tsr,[],size(rsz_imgs_tsr,4))'/255.0);
pasu_imdist.L2 = squareform(L2dist);% D = D.select_metric("L2");%D.distmat(rsz_imgs_tsr);
toc% 
% Compute FC6 corr and L2 
acts = net.activations(rsz_imgs_tsr,"fc6");
FC6dist = squareform(pdist(squeeze(acts)','euclidean'));
pasu_imdist.FC6 = FC6dist;
FC6dist_corr = squareform(pdist(squeeze(acts)','correlation'));
pasu_imdist.FC6_corr = FC6dist_corr;
toc
tic
D = D.select_metric("SSIM");
pasu_imdist.SSIM = D.distmat(rsz_imgs_tsr);
toc% 119sec
%%
save(fullfile(mat_dir,"pasu_imdist.mat"),'pasu_imdist')
%%
%% Gabor distance
gab_nm_grid = cellfun(@(idx)unique(Stats(40).imageName(idx)),Stats(40).ref.gab_idx_grid,'Uni',1);
gab_nm_grid = reshape(gab_nm_grid,[],1); % reshape into one row. But transpose to make similar object adjacent
gab_nms = string(gab_nm_grid);
%%
gab_path = "N:\Stimuli\2019-Manifold\gabor";
gab_imgs_all = arrayfun(@(nm)imread(fullfile(gab_path, nm+".bmp")),gab_nms,'Uni',0);
gab_rsz_imgs_all = cellfun(@(img)repmat(imresize(img,[256,256]),1,1,3),gab_imgs_all,'Uni',0);
gab_rsz_imgs_tsr = cell2mat(shiftdim(gab_rsz_imgs_all,-3));
%%
gab_imdist = struct();
tic
D=D.select_metric("squeeze");
gab_imdist.squ = D.distmat(gab_rsz_imgs_tsr);
toc
L2dist = pdist(reshape(gab_rsz_imgs_tsr,[],size(gab_rsz_imgs_tsr,4))'/255.0);
gab_imdist.L2 = squareform(L2dist);% D = D.select_metric("L2");%D.distmat(rsz_imgs_tsr);
toc% 
tic
D = D.select_metric("SSIM");
gab_imdist.SSIM = D.distmat(gab_rsz_imgs_tsr);
toc% 119sec
acts = net.activations(gab_rsz_imgs_tsr,"fc6");
FC6dist_corr = squareform(pdist(squeeze(acts)','correlation'));
gab_imdist.FC6_corr = FC6dist_corr;
FC6dist = squareform(pdist(squeeze(acts)','euclidean'));
gab_imdist.FC6 = FC6dist;
toc
%%
save(fullfile(mat_dir,"gab_imdist.mat"),'gab_imdist')
%%
Expi=1;si=1;
imnm_grid = string(cellfun(@(idx)Stats(Expi).imageName{idx(1)},Stats(Expi).manif.idx_grid{si},'UniformOutput',false));
manif_imgrid = arrayfun(@(nm)imread(fullfile(Stats(Expi).meta.stimuli,nm+".jpg")),imnm_grid,"Uni",0);
manif_img_tsr = cell2mat(reshape(manif_imgrid,1,1,1,[]));
%% Demo
tic
D.select_metric("squeeze");
manif_imdist.squ = D.distmat(manif_img_tsr);
toc
% D.select_metric("SSIM");
% manif_imdist.SSIM = D.distmat(manif_img_tsr);
% toc% 119sec
% D.select_metric("L2");
% manif_imdist.L2 = D.distmat(manif_img_tsr);
% toc% 
%%
save("manif_imdist.mat",'manif_imdist')
%% Compute the image distance function of all Images on manifold
ManifImDistStat = repmat(struct(),1,numel(Stats));
for Expi=1:numel(Stats)
tic;
fprintf("Processing exp %d\n",Expi)
si=1;
imnm_grid = string(cellfun(@(idx)Stats(Expi).imageName{idx(1)},Stats(Expi).manif.idx_grid{si},'Uni',0));
manif_imgrid = arrayfun(@(nm)imread(fullfile(Stats(Expi).meta.stimuli,nm+".jpg")),imnm_grid,"Uni",0);
manif_img_tsr = cell2mat(reshape(manif_imgrid,1,1,1,[]));
toc
D = D.select_metric("squeeze");
ManifImDistStat(Expi).squ = D.distmat_B(manif_img_tsr); % squeeze net. 
toc
L2dist = pdist(reshape(manif_img_tsr,[],size(manif_img_tsr,4))'/255.0);
ManifImDistStat(Expi).L2 = squareform(L2dist); % L2 distance. 
toc
acts = net.activations(manif_img_tsr,"fc6"); 
FC6dist_corr = squareform(pdist(squeeze(acts)','correlation'));
ManifImDistStat(Expi).FC6_corr = FC6dist_corr; % FC6 corr
FC6dist = squareform(pdist(squeeze(acts)','euclidean'));
ManifImDistStat(Expi).FC6 = FC6dist; % FC6 L2 distance
toc
tic
D = D.select_metric("SSIM");
ManifImDistStat(Expi).SSIM = D.distmat(manif_img_tsr); % SSIM distance
toc% 119sec
% D = D.select_metric("SSIM");
% ManifImDistStat(Expi).SSIM = D.distmat(manif_img_tsr);
% toc% 119sec
% D = D.select_metric("L2");
% ManifImDistStat(Expi).L2 = D.distmat(manif_img_tsr);
% toc% 
end
%%
save(fullfile(mat_dir,Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
%%
for Expi=1:numel(Stats)
tic;
fprintf("Processing exp %d\n",Expi)
si=1;
imnm_grid = string(cellfun(@(idx)Stats(Expi).imageName{idx(1)},Stats(Expi).manif.idx_grid{si},'Uni',0));
manif_imgrid = arrayfun(@(nm)imread(fullfile(Stats(Expi).meta.stimuli,nm+".jpg")),imnm_grid,"Uni",0);
manif_img_tsr = cell2mat(reshape(manif_imgrid,1,1,1,[]));
acts = net.activations(manif_img_tsr,"fc6");
FC6dist = squareform(pdist(squeeze(acts)','correlation'));
ManifImDistStat(Expi).FC6 = FC6dist;
toc;
end
%%
save(fullfile(mat_dir,Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
%%
pasu_nm_grid = cellfun(@(idx)unique(Stats(1).imageName(idx)),Stats(1).ref.pasu_idx_grid,'Uni',0);
pasu_nm_grid = reshape(pasu_nm_grid',[],1); % reshape into one row. But transpose to make similar object adjacent
pasu_val_msk = ~cellfun(@isempty,pasu_nm_grid);
% pasu_nm_grid(cellfun(@isempty,pasu_nm_grid)) = []; % get rid of empty entries
% pasu_nms = string(pasu_nm_grid);
%% Load image distance matrices
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
%% Compute Summary Statistics and Collect stats into a table.
%  the relationship between the image metric and the firing. 
%  Summarize it in a large csv
Expi = 35;
nanstrct = struct('Manif_squ_pear',nan,'Manif_squ_spear',nan,'Manif_FC6_pear',nan,'Manif_FC6_spear',nan,'Manif_squ_pear_max',nan,'Manif_squ_spear_max',nan,'Manif_FC6_pear_max',nan,'Manif_FC6_spear_max',nan,'Pasu_squ_pear',nan,'Pasu_squ_spear',nan,'Pasu_FC6_pear',nan,'Pasu_FC6_spear',nan,'Pasu_squ_pear_max',nan,'Pasu_squ_spear_max',nan,'Pasu_FC6_pear_max',nan,'Pasu_FC6_spear_max',nan,'Gabor_squ_pear',nan,'Gabor_squ_spear',nan,'Gabor_FC6_pear',nan,'Gabor_FC6_spear',nan,'Gabor_squ_pear_max',nan,'Gabor_squ_spear_max',nan,'Gabor_FC6_pear_max',nan,'Gabor_FC6_spear_max',nan,"Manif_squ_reg_b",nan(2,1),"Manif_squ_reg_Rsq",nan,"Manif_squ_reg_F",nan,"Manif_squ_reg_P",nan,"Manif_FC6_reg_b",nan(2,1),"Manif_FC6_reg_Rsq",nan,"Manif_FC6_reg_F",nan,"Manif_FC6_reg_P",nan,"Pasu_squ_reg_b",nan(2,1),"Pasu_squ_reg_Rsq",nan,"Pasu_squ_reg_F",nan,"Pasu_squ_reg_P",nan,"Pasu_FC6_reg_b",nan(2,1),"Pasu_FC6_reg_Rsq",nan,"Pasu_FC6_reg_F",nan,"Pasu_FC6_reg_P",nan,"Gabor_squ_reg_b",nan(2,1),"Gabor_squ_reg_Rsq",nan,"Gabor_squ_reg_F",nan,"Gabor_squ_reg_P",nan,"Gabor_FC6_reg_b",nan(2,1),"Gabor_FC6_reg_Rsq",nan,"Gabor_FC6_reg_F",nan,"Gabor_FC6_reg_P",nan);
metCorrStats = repmat(nanstrct,1,numel(Stats));
for Expi=1:numel(Stats) % basic stats for exp
metCorrStats(Expi).Animal = Animal;
metCorrStats(Expi).Expi = Expi;
metCorrStats(Expi).area = chan2area(Stats(Expi).units.pref_chan);
metCorrStats(Expi).prefchan = Stats(Expi).units.pref_chan;
metCorrStats(Expi).prefchan_ui = EStats(Expi).evol.unit_in_pref_chan; % unit in pref chan
metCorrStats(Expi).imgpos = EStats(Expi).evol.imgpos; % position of image
metCorrStats(Expi).imgsize = EStats(Expi).evol.imgsize; % position of image
end
for Expi=1:numel(Stats) % Stats for manifold tuning specfically
si=1;ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
metCorrStats(Expi).Manif_rng = [max(score_vec),min(score_vec)];
metCorrStats(Expi).Manif_squ_pear = corr(ManifImDistStat(Expi).squ(:,maxId),score_vec);
metCorrStats(Expi).Manif_squ_spear = corr(ManifImDistStat(Expi).squ(:,maxId),score_vec,'Type','Spearman');
metCorrStats(Expi).Manif_FC6_pear = corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec);
metCorrStats(Expi).Manif_FC6_spear = corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec,'Type','Spearman');

metCorrStats(Expi).Manif_squ_pear_max = min(corr(ManifImDistStat(Expi).squ,score_vec));
metCorrStats(Expi).Manif_squ_spear_max = min(corr(ManifImDistStat(Expi).squ,score_vec,'Type','Spearman'));
metCorrStats(Expi).Manif_FC6_pear_max = min(corr(ManifImDistStat(Expi).FC6,score_vec));
metCorrStats(Expi).Manif_FC6_spear_max = min(corr(ManifImDistStat(Expi).FC6,score_vec,'Type','Spearman'));
% Regression coefficient
[b,~,stats] = regress_(ManifImDistStat(Expi).squ(:,maxId), score_vec);
metCorrStats(Expi).Manif_squ_reg_b = b; 
metCorrStats(Expi).Manif_squ_reg_Rsq = stats(1);
metCorrStats(Expi).Manif_squ_reg_F = stats(2);
metCorrStats(Expi).Manif_squ_reg_P = stats(3);
[b,~,stats] = regress_(ManifImDistStat(Expi).FC6(:,maxId), score_vec);
metCorrStats(Expi).Manif_FC6_reg_b = b;
metCorrStats(Expi).Manif_FC6_reg_Rsq = stats(1);
metCorrStats(Expi).Manif_FC6_reg_F = stats(2);
metCorrStats(Expi).Manif_FC6_reg_P = stats(3);
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(~pasu_val_msk) = [];
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
metCorrStats(Expi).Pasu_rng = [max(pasu_vec),min(pasu_vec)];
metCorrStats(Expi).Pasu_incomp = sum(isnan(pasu_vec))>2;
metCorrStats(Expi).Pasu_squ_pear = corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Rows','complete');
metCorrStats(Expi).Pasu_squ_spear = corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete');
metCorrStats(Expi).Pasu_FC6_pear = corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Rows','complete');
metCorrStats(Expi).Pasu_FC6_spear = corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete');

metCorrStats(Expi).Pasu_squ_pear_max = min(corr(pasu_imdist.squ,pasu_vec,'Rows','complete'));
metCorrStats(Expi).Pasu_squ_spear_max = min(corr(pasu_imdist.squ,pasu_vec,'Type','Spearman','Rows','complete'));
metCorrStats(Expi).Pasu_FC6_pear_max = min(corr(pasu_imdist.FC6,pasu_vec,'Rows','complete'));
metCorrStats(Expi).Pasu_FC6_spear_max = min(corr(pasu_imdist.FC6,pasu_vec,'Type','Spearman','Rows','complete'));
% Regression coefficient
[b,~,stats] = regress_(pasu_imdist.squ(:,pasu_maxId),pasu_vec);
metCorrStats(Expi).Pasu_squ_reg_b = b; 
metCorrStats(Expi).Pasu_squ_reg_Rsq = stats(1);
metCorrStats(Expi).Pasu_squ_reg_F = stats(2);
metCorrStats(Expi).Pasu_squ_reg_P = stats(3);
[b,~,stats] = regress_(pasu_imdist.FC6(:,pasu_maxId),pasu_vec);
metCorrStats(Expi).Pasu_FC6_reg_b = b;
metCorrStats(Expi).Pasu_FC6_reg_Rsq = stats(1);
metCorrStats(Expi).Pasu_FC6_reg_F = stats(2);
metCorrStats(Expi).Pasu_FC6_reg_P = stats(3);
end
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
metCorrStats(Expi).Gabor_rng = [max(gab_vec),min(gab_vec)];
metCorrStats(Expi).Gabor_incomp = sum(isnan(gab_vec))>2;
metCorrStats(Expi).Gabor_squ_pear = corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Rows','complete');
metCorrStats(Expi).Gabor_squ_spear = corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete');
metCorrStats(Expi).Gabor_FC6_pear = corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Rows','complete');
metCorrStats(Expi).Gabor_FC6_spear = corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete');

metCorrStats(Expi).Gabor_squ_pear_max = min(corr(gab_imdist.squ,gab_vec,'Rows','complete'));
metCorrStats(Expi).Gabor_squ_spear_max = min(corr(gab_imdist.squ,gab_vec,'Type','Spearman','Rows','complete'));
metCorrStats(Expi).Gabor_FC6_pear_max = min(corr(gab_imdist.FC6,gab_vec,'Rows','complete'));
metCorrStats(Expi).Gabor_FC6_spear_max = min(corr(gab_imdist.FC6,gab_vec,'Type','Spearman','Rows','complete'));
% Regression coefficient
[b,~,stats] = regress_(gab_imdist.squ(:,gab_maxId),gab_vec);
metCorrStats(Expi).Gabor_squ_reg_b = b; 
metCorrStats(Expi).Gabor_squ_reg_Rsq = stats(1);
metCorrStats(Expi).Gabor_squ_reg_F = stats(2);
metCorrStats(Expi).Gabor_squ_reg_P = stats(3);
[b,~,stats] = regress_(gab_imdist.FC6(:,gab_maxId),gab_vec);
metCorrStats(Expi).Gabor_FC6_reg_b = b;
metCorrStats(Expi).Gabor_FC6_reg_Rsq = stats(1);
metCorrStats(Expi).Gabor_FC6_reg_F = stats(2);
metCorrStats(Expi).Gabor_FC6_reg_P = stats(3);
end
end
%%
rspCorrTab = struct2table(metCorrStats);
writetable(rspCorrTab,fullfile(figdir,Animal+"_ImMetricCorrTab.csv"));

%% Visualize the trend of the image similarity and firing
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
Expi = 35;
bestcorr = 0;
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
% Gaussian Process Smoothing
xlinsp = linspace(0, max(ManifImDistStat(Expi).squ(:,maxId)), 50);
gprMdl = fitrgp(ManifImDistStat(Expi).squ(:,maxId), score_vec);
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit)
xlabel("Perceptual Similarity(SqueezeNet)");ylabel("Neural Activation")
titstr1 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).squ(:,maxId),score_vec),corr(ManifImDistStat(Expi).squ(:,maxId),score_vec,'Type','Spearman'))};
% nexttile(2)
% scatter(ManifImDistStat(Expi).SSIM(:,maxId),score_vec)
% xlabel("SSIM")
% nexttile(3)
% scatter(ManifImDistStat(Expi).L2(:,maxId),score_vec)
% xlabel("L2")
nexttile(4)
scatter(ManifImDistStat(Expi).FC6(:,maxId),score_vec)
titstr4 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec),corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec,'Type','Spearman'))};
xlabel("FC6 (1-corr)")
% Pasupathy patches
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(~pasu_val_msk) = []; % isnan(pasu_vec)
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
nexttile(1);hold on
scatter(pasu_imdist.squ(:,pasu_maxId),pasu_vec)
% Gaussian Process Smoothing
xlinsp = linspace(0,max(pasu_imdist.squ(:,pasu_maxId)),100);
gprMdl = fitrgp(pasu_imdist.squ(:,pasu_maxId), pasu_vec);
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit)
titstr1{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Rows','complete'),corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
% nexttile(2);hold on
% scatter(pasu_imdist.SSIM(:,pasu_maxId),pasu_vec)
% nexttile(3);hold on
% scatter(pasu_imdist.L2(:,pasu_maxId),pasu_vec)
nexttile(4);hold on
scatter(pasu_imdist.FC6(:,pasu_maxId),pasu_vec)
legend(["Manifold","Pasupathy"])
titstr4{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Rows','complete'),corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
end
% Gabor patches
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
nexttile(1);hold on
scatter(gab_imdist.squ(:,gab_maxId), gab_vec,'g')
% Plot the Interpolation of it
xlinsp = linspace(0, max(gab_imdist.squ(:,gab_maxId)),100);
gprMdl = fitrgp(gab_imdist.squ(:,gab_maxId), gab_vec);
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit,'g')
titstr1{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Rows','complete'),corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
% nexttile(2);hold on
% scatter(gab_imdist.SSIM(:,gab_maxId),gab_vec,'g')
% nexttile(3);hold on
% scatter(gab_imdist.L2(:,gab_maxId),gab_vec,'g')
nexttile(4);hold on
scatter(gab_imdist.FC6(:,gab_maxId),gab_vec,'g')
legend(["Manifold","Pasupathy","Gabor"])
titstr4{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Rows','complete'),corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
end
nexttile(1);title(titstr1)
nexttile(4);title(titstr4)
% saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.png",Animal,Expi,Stats(Expi).units.pref_chan)))
% savefig(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.fig",Animal,Expi,Stats(Expi).units.pref_chan)))
end
% title(T,compose("%s Exp %d prefchan %d Pasupathy",Animal,Expi,Stats(Expi).units.pref_chan))
%% Plot the same thing with maximum correlation stimuli. 
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
Expi = 35;
bestcorr = 1;
for Expi=34:numel(Stats)
si=1;ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
if bestcorr
CC_pos = corr(score_vec,ManifImDistStat(Expi).squ,'type','Spearman','Rows','complete');
[maxCC,maxId] = min(CC_pos);
end
figure(22);set(22,'pos',[0   286  2000  544])
T=tiledlayout(1,4,'TileSpacing','compact');
title(T,compose("%s Exp %d prefchan %d",Animal,Expi,Stats(Expi).units.pref_chan),'FontSize',16)

nexttile(1) %subplot(141);
scatter(ManifImDistStat(Expi).squ(:,maxId),score_vec)
xlabel("Perceptual Similarity(SqueezeNet)");ylabel("Neural Activation")
titstr1 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).squ(:,maxId),score_vec),corr(ManifImDistStat(Expi).squ(:,maxId),score_vec,'Type','Spearman'))};
% nexttile %subplot(142);
% scatter(ManifImDistStat(Expi).SSIM(:,maxId),score_vec)
% xlabel("SSIM")
% nexttile %subplot(143);
% scatter(ManifImDistStat(Expi).L2(:,maxId),score_vec)
% xlabel("L2")
nexttile(4) %subplot(144);
scatter(ManifImDistStat(Expi).FC6(:,maxId),score_vec)
titstr4 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec),corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec,'Type','Spearman'))};
xlabel("FC6 (1-corr)")
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(~pasu_val_msk) = [];
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
if bestcorr
CC_pos = corr(pasu_vec,pasu_imdist.squ,'type','Spearman');
[~,pasu_maxId] = min(CC_pos);
end
nexttile(1);hold on%subplot(141);
scatter(pasu_imdist.squ(:,pasu_maxId),pasu_vec)
titstr1{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Rows','complete'),corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
% nexttile(2);hold on%subplot(142);
% scatter(pasu_imdist.SSIM(:,pasu_maxId),pasu_vec)
% nexttile(3);hold on%subplot(143);
% scatter(pasu_imdist.L2(:,pasu_maxId),pasu_vec)
nexttile(4);hold on%subplot(144);
scatter(pasu_imdist.FC6(:,pasu_maxId),pasu_vec)
legend(["Manifold","Pasupathy"])
titstr4{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Rows','complete'),corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
end
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
if bestcorr
CC_pos = corr(gab_vec,gab_imdist.squ,'type','Spearman');
[~,gab_maxId] = min(CC_pos);
end
nexttile(1);hold on
scatter(gab_imdist.squ(:,gab_maxId),gab_vec,50,[0,1,0])
titstr1{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Rows','complete'),corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
% nexttile(2);hold on
% scatter(gab_imdist.SSIM(:,gab_maxId),gab_vec,50,[0,1,0])
% nexttile(3);hold on
% scatter(gab_imdist.L2(:,gab_maxId),gab_vec,50,[0,1,0])
nexttile(4);hold on
scatter(gab_imdist.FC6(:,gab_maxId),gab_vec,50,[0,1,0])
legend(["Manifold","Pasupathy","Gabor"])
titstr4{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Rows','complete'),corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
end
nexttile(1);title(titstr1)
nexttile(4);title(titstr4)
saveas(22,fullfile(figdir,compose("%s_Exp%02d_pref%02d_bestcorr.png",Animal,Expi,Stats(Expi).units.pref_chan)))
savefig(22,fullfile(figdir,compose("%s_Exp%02d_pref%02d_bestcorr.fig",Animal,Expi,Stats(Expi).units.pref_chan)))
end

%%
% [b,bint,stats] = regress_(ManifImDistStat(Expi).squ(:,maxId),score_vec)

function area_arr = chan2area(chans)
area_arr = strings(size(chans));
area_arr(chans<=32)="IT";
area_arr(33<=chans & chans<=48)="V1";
area_arr(49<=chans)="V4";
end
function [b,bint,stats] = regress_(dist, act)
[b,bint,~,~,stats] = regress(act, [dist, ones(numel(dist), 1)]);
end

