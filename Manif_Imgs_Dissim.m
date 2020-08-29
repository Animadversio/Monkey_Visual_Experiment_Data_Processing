%% Manif Image Dissimilarity
py.torch.cuda.empty_cache()
D = torchImDist();
%%
Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
net = alexnet;
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
pasu_imdist.squ = DistMat;
%%
tic
D.select_metric();
pasu_imdist.squ = D.distmat(rsz_imgs_tsr);
toc
D.select_metric("SSIM");
pasu_imdist.SSIM = D.distmat(rsz_imgs_tsr);
toc% 119sec
D.select_metric("L2");
pasu_imdist.L2 = D.distmat(rsz_imgs_tsr);
toc% 
%%
acts = net.activations(rsz_imgs_tsr,"fc6");
FC6dist = squareform(pdist(squeeze(acts)','correlation'));
pasu_imdist.FC6 = FC6dist;
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
% D=D.select_metric("SSIM");
% gab_imdist.SSIM = D.distmat(gab_rsz_imgs_tsr);
% toc% 119sec
% D=D.select_metric("L2");
% gab_imdist.L2 = D.distmat(gab_rsz_imgs_tsr);
gab_imdist.L2 = gab_imdist.squ;
gab_imdist.SSIM = gab_imdist.squ;
toc% 
acts = net.activations(gab_rsz_imgs_tsr,"fc6");
FC6dist = squareform(pdist(squeeze(acts)','correlation'));
gab_imdist.FC6 = FC6dist;
toc
%%
acts = net.activations(rsz_imgs_tsr,"fc6");
FC6dist = squareform(pdist(squeeze(acts)','correlation'));
pasu_imdist.FC6 = FC6dist;
%%
save(fullfile(mat_dir,"gab_imdist.mat"),'gab_imdist')
%%
Expi=1;si=1;
imnm_grid = string(cellfun(@(idx)Stats(Expi).imageName{idx(1)},Stats(Expi).manif.idx_grid{si},'UniformOutput',false));
manif_imgrid = arrayfun(@(nm)imread(fullfile(Stats(Expi).meta.stimuli,nm+".jpg")),imnm_grid,"Uni",0);
manif_img_tsr = cell2mat(reshape(manif_imgrid,1,1,1,[]));
%%
tic
D.select_metric("squeeze");
manif_imdist.squ = D.distmat(manif_img_tsr);
toc
D.select_metric("SSIM");
manif_imdist.SSIM = D.distmat(manif_img_tsr);
toc% 119sec
D.select_metric("L2");
manif_imdist.L2 = D.distmat(manif_img_tsr);
toc% 
%%
save("manif_imdist.mat",'manif_imdist')
%% Get the image distance function of all Images on manifold
ManifImDistStat = repmat(struct(),1,numel(Stats));
for Expi=1:numel(Stats)
tic;
fprintf("Processing exp %d\n",Expi)
si=1;
imnm_grid = string(cellfun(@(idx)Stats(Expi).imageName{idx(1)},Stats(Expi).manif.idx_grid{si},'Uni',0));
manif_imgrid = arrayfun(@(nm)imread(fullfile(Stats(Expi).meta.stimuli,nm+".jpg")),imnm_grid,"Uni",0);
manif_img_tsr = cell2mat(reshape(manif_imgrid,1,1,1,[]));
%%
toc
D = D.select_metric("squeeze");
ManifImDistStat(Expi).squ = D.distmat(manif_img_tsr);
toc
D = D.select_metric("SSIM");
ManifImDistStat(Expi).SSIM = D.distmat(manif_img_tsr);
toc% 119sec
D = D.select_metric("L2");
ManifImDistStat(Expi).L2 = D.distmat(manif_img_tsr);
toc% 
end
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
end
%%
save(fullfile(mat_dir,Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
%%
pasu_nm_grid = cellfun(@(idx)unique(Stats(Expi).imageName(idx)),Stats(Expi).ref.pasu_idx_grid,'Uni',0);
pasu_nm_grid = reshape(pasu_nm_grid',[],1); % reshape into one row. But transpose to make similar object adjacent
pasu_nm_grid(cellfun(@isempty,pasu_nm_grid)) = []; % get rid of empty entries
pasu_nms = string(pasu_nm_grid);
%%
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
Expi = 35;
bestcorr = 0;
for Expi=1:45
si=1;ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
figure(21);
T=tiledlayout(1,4,'TileSpacing','compact');
title(T,compose("%s Exp %d prefchan %d",Animal,Expi,Stats(Expi).units.pref_chan),'FontSize',16)

nexttile%subplot(141);
scatter(ManifImDistStat(Expi).squ(:,maxId),score_vec)
xlabel("Perceptual Similarity(SqueezeNet)");ylabel("Neural Activation")
titstr1 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).squ(:,maxId),score_vec),corr(ManifImDistStat(Expi).squ(:,maxId),score_vec,'Type','Spearman'))};
nexttile%subplot(142);
scatter(ManifImDistStat(Expi).SSIM(:,maxId),score_vec)
xlabel("SSIM")
nexttile%subplot(143);
scatter(ManifImDistStat(Expi).L2(:,maxId),score_vec)
xlabel("L2")
nexttile%subplot(144);
scatter(ManifImDistStat(Expi).FC6(:,maxId),score_vec)
titstr4 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec),corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec,'Type','Spearman'))};
xlabel("FC6 (1-corr)")
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(isnan(pasu_vec)) = [];
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
nexttile(1);hold on%subplot(141);
scatter(pasu_imdist.squ(:,pasu_maxId),pasu_vec)
titstr1{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec),corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Type','Spearman'));
nexttile(2);hold on%subplot(142);
scatter(pasu_imdist.SSIM(:,pasu_maxId),pasu_vec)
nexttile(3);hold on%subplot(143);
scatter(pasu_imdist.L2(:,pasu_maxId),pasu_vec)
nexttile(4);hold on%subplot(144);
scatter(pasu_imdist.FC6(:,pasu_maxId),pasu_vec)
legend(["Manifold","Pasupathy"])
titstr4{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec),corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Type','Spearman'));
end

if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
nexttile(1);hold on
scatter(gab_imdist.squ(:,gab_maxId),gab_vec,'g')
titstr1{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.squ(:,gab_maxId),gab_vec),corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Type','Spearman'));
nexttile(2);hold on
scatter(gab_imdist.SSIM(:,gab_maxId),gab_vec,'g')
nexttile(3);hold on
scatter(gab_imdist.L2(:,gab_maxId),gab_vec,'g')
nexttile(4);hold on
scatter(gab_imdist.FC6(:,gab_maxId),gab_vec,'g')
legend(["Manifold","Pasupathy","Gabor"])
titstr4{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.FC6(:,gab_maxId),gab_vec),corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Type','Spearman'));
end
nexttile(1);title(titstr1)
nexttile(4);title(titstr4)
saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.png",Animal,Expi,Stats(Expi).units.pref_chan)))
savefig(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.fig",Animal,Expi,Stats(Expi).units.pref_chan)))
end
% title(T,compose("%s Exp %d prefchan %d Pasupathy",Animal,Expi,Stats(Expi).units.pref_chan))
%%
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
Expi = 35;
bestcorr = 1;
for Expi=11:45
si=1;ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
if bestcorr
CC_pos = corr(score_vec,ManifImDistStat(Expi).squ,'type','Spearman');
[maxCC,maxId] = min(CC_pos);
end
figure(22);set(22,'pos',[0   286  2000  544])
T=tiledlayout(1,4,'TileSpacing','compact');
title(T,compose("%s Exp %d prefchan %d",Animal,Expi,Stats(Expi).units.pref_chan),'FontSize',16)

nexttile%subplot(141);
scatter(ManifImDistStat(Expi).squ(:,maxId),score_vec)
xlabel("Perceptual Similarity(SqueezeNet)");ylabel("Neural Activation")
titstr1 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).squ(:,maxId),score_vec),corr(ManifImDistStat(Expi).squ(:,maxId),score_vec,'Type','Spearman'))};
nexttile%subplot(142);
scatter(ManifImDistStat(Expi).SSIM(:,maxId),score_vec)
xlabel("SSIM")
nexttile%subplot(143);
scatter(ManifImDistStat(Expi).L2(:,maxId),score_vec)
xlabel("L2")
nexttile%subplot(144);
scatter(ManifImDistStat(Expi).FC6(:,maxId),score_vec)
titstr4 = { compose("Manif: pear %.3f spear %.3f",corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec),corr(ManifImDistStat(Expi).FC6(:,maxId),score_vec,'Type','Spearman'))};
xlabel("FC6 (1-corr)")
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(isnan(pasu_vec)) = [];
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
if bestcorr
CC_pos = corr(pasu_vec,pasu_imdist.squ,'type','Spearman');
[~,pasu_maxId] = min(CC_pos);
end
nexttile(1);hold on%subplot(141);
scatter(pasu_imdist.squ(:,pasu_maxId),pasu_vec)
titstr1{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec),corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Type','Spearman'));
nexttile(2);hold on%subplot(142);
scatter(pasu_imdist.SSIM(:,pasu_maxId),pasu_vec)
nexttile(3);hold on%subplot(143);
scatter(pasu_imdist.L2(:,pasu_maxId),pasu_vec)
nexttile(4);hold on%subplot(144);
scatter(pasu_imdist.FC6(:,pasu_maxId),pasu_vec)
legend(["Manifold","Pasupathy"])
titstr4{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec),corr(pasu_imdist.FC6(:,pasu_maxId),pasu_vec,'Type','Spearman'));
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
titstr1{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.squ(:,gab_maxId),gab_vec),corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Type','Spearman'));
nexttile(2);hold on
scatter(gab_imdist.SSIM(:,gab_maxId),gab_vec,50,[0,1,0])
nexttile(3);hold on
scatter(gab_imdist.L2(:,gab_maxId),gab_vec,50,[0,1,0])
nexttile(4);hold on
scatter(gab_imdist.FC6(:,gab_maxId),gab_vec,50,[0,1,0])
legend(["Manifold","Pasupathy","Gabor"])
titstr4{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.FC6(:,gab_maxId),gab_vec),corr(gab_imdist.FC6(:,gab_maxId),gab_vec,'Type','Spearman'));
end
nexttile(1);title(titstr1)
nexttile(4);title(titstr4)
saveas(22,fullfile(figdir,compose("%s_Exp%02d_pref%02d_bestcorr.png",Animal,Expi,Stats(Expi).units.pref_chan)))
savefig(22,fullfile(figdir,compose("%s_Exp%02d_pref%02d_bestcorr.fig",Animal,Expi,Stats(Expi).units.pref_chan)))
end
%% Summary Statistics


