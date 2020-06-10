net = vgg16;
%%
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
%%
Animal="Alfa"; 
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
Animal="Alfa"; ExpType="Manif"; 
flags = struct("batch",50,"online_compute",0,"load_all_img",1,"shuffleN",100);
wdw_vect = [[1, 20] + 10 * [0:18]'; [1,50]+[0:50:150]'; [51,200]];
% Expi = 11; 
for Expi = 2:46
T00 = tic;
[imgfn, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, flags);
toc(T00)
totl_vox_num = [];
corr_vox_num = zeros(24,[]);
corr_vox_prct = zeros(24,[]);
med_pos_cc = zeros(24,[]);
med_neg_cc = zeros(24,[]);
mean_pos_t = zeros(24,[]);
mean_neg_t = zeros(24,[]);
%%
layernames = ["fc8", "fc7", "fc6", "conv5_3", "conv4_3", "conv3_3"];
for iLayer = 1:length(layernames)
layername = layernames(iLayer);
[cc_tsr, MFeat, StdFeat, cc_refM, cc_refS] = corrFeatTsr_func(imgfn, score_vect, net, layername, flags);
outfn = fullfile(savedir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername)); % LW means long window
save(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
%%
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
nfeat = prod(size(t_signif_tsr,[1,2,3]));
totl_vox_num(iLayer) = nfeat;
% View the number of signicantly correlated voxels as a function of time window.
fprintf("Correlation layer %s (vox # %d)\n",layername,nfeat)
fid = fopen(fullfile(hier_savedir,Animal+"_Manif_all_VGG.log"),'w+');
fprintf(fid,"Correlation layer %s (vox # %d)\n",layername,nfeat);
fi = 24; wdw = wdw_vect(fi,:);%ci=6; 
for fi = 1:24
wdw = wdw_vect(fi,:);
tcol = t_signif_tsr(:,:,:,fi);%yi,xi,
cc_tmp = cc_tsr(:,:,:,fi);
signif_n = sum(tcol > 5 | tcol<-5,'all');
fprintf("Firing rate in [%d, %d] ms Signif corr voxel num %d (%.1f), pos corr median %.3f (t=%.2f), neg corr median %.3f (t=%.2f)\n",...
    wdw(1),wdw(2),signif_n,signif_n/nfeat*100,median(cc_tmp(tcol>5)),mean(tcol(tcol>5)),...
    median(cc_tmp(tcol<-5)),mean(tcol(tcol<-5)))
fprintf(fid, "Firing rate in [%d, %d] ms Signif corr voxel num %d (%.1f), pos corr median %.3f (t=%.2f), neg corr median %.3f (t=%.2f)\n",...
    wdw(1),wdw(2),signif_n,signif_n/nfeat*100,median(cc_tmp(tcol>5)),mean(tcol(tcol>5)),...
    median(cc_tmp(tcol<-5)),mean(tcol(tcol<-5)));
corr_vox_num(fi,iLayer) = signif_n;
med_pos_cc(fi,iLayer) = median(cc_tmp(tcol>5));
med_neg_cc(fi,iLayer) = median(cc_tmp(tcol<-5));
mean_pos_t(fi,iLayer) = mean(tcol(tcol>5));
mean_neg_t(fi,iLayer) = mean(tcol(tcol<-5));
end
fclose(fid);
% save(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
%     "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t");
toc(T00)
end
save(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t");
end
%%
%
Animal="Alfa"; ExpType="Evol"; 
flags = struct("batch",60,"online_compute",1,"load_all_img",1,"shuffleN",100);
wdw_vect = [[1, 20] + 10 * [0:18]'; [1,50]+[0:50:150]'; [51,200]];
% Expi = 11; 
for Expi = 1:46
T00 = tic;
[imgfn, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, flags);
toc(T00)
totl_vox_num = [];
corr_vox_num = zeros(24,[]);
corr_vox_prct = zeros(24,[]);
med_pos_cc = zeros(24,[]);
med_neg_cc = zeros(24,[]);
mean_pos_t = zeros(24,[]);
mean_neg_t = zeros(24,[]);
%%
layernames = ["fc8", "fc7", "fc6", "conv5_3", "conv4_3", "conv3_3"];
for iLayer = 1:length(layernames)
layername = layernames(iLayer);
[cc_tsr, MFeat, StdFeat, cc_refM, cc_refS] = corrFeatTsr_func(imgfn, score_vect, net, layername, flags);
outfn = fullfile(savedir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername)); % LW means long window
save(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
%%
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
nfeat = prod(size(t_signif_tsr,[1,2,3]));
totl_vox_num(iLayer) = nfeat;
% View the number of signicantly correlated voxels as a function of time window.
fid = fopen(fullfile(hier_savedir,Animal + "_Evol_all_VGG.log"),'w+');
fprintf("Correlation layer %s (vox # %d)\n",layername,nfeat)
fprintf(fid,"Correlation layer %s (vox # %d)\n",layername,nfeat);
fi = 24; wdw = wdw_vect(fi,:);%ci=6; 
for fi = 1:24
wdw = wdw_vect(fi,:);
tcol = t_signif_tsr(:,:,:,fi);%yi,xi,
cc_tmp = cc_tsr(:,:,:,fi);
signif_n = sum(tcol > 5 | tcol<-5,'all');
fprintf("Firing rate in [%d, %d] ms Signif corr voxel num %d (%.1f), pos corr median %.3f (t=%.2f), neg corr median %.3f (t=%.2f)\n",...
    wdw(1),wdw(2),signif_n,signif_n/nfeat*100,median(cc_tmp(tcol>5)),mean(tcol(tcol>5)),...
    median(cc_tmp(tcol<-5)),mean(tcol(tcol<-5)))
fprintf(fid, "Firing rate in [%d, %d] ms Signif corr voxel num %d (%.1f), pos corr median %.3f (t=%.2f), neg corr median %.3f (t=%.2f)\n",...
    wdw(1),wdw(2),signif_n,signif_n/nfeat*100,median(cc_tmp(tcol>5)),mean(tcol(tcol>5)),...
    median(cc_tmp(tcol<-5)),mean(tcol(tcol<-5)));
corr_vox_num(fi,iLayer) = signif_n;
med_pos_cc(fi,iLayer) = median(cc_tmp(tcol>5));
med_neg_cc(fi,iLayer) = median(cc_tmp(tcol<-5));
mean_pos_t(fi,iLayer) = mean(tcol(tcol>5));
mean_neg_t(fi,iLayer) = mean(tcol(tcol<-5));
end
fclose(fid);
% save(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
%     "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t");
toc(T00)
end
save(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t");
end
%%

%%
function dlimg = loadimges(images)
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(imgnm)),[224,224]), images, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[])); 
end

function [imgfullnm_vect, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, flags)
ui=1;si=1;
assert(EStats(Expi).Animal == Animal && Stats(Expi).Animal == Animal)
fprintf("Processing %s Exp %d pref chan %d\n",ExpType,Expi,EStats(Expi).units.pref_chan)
if ExpType == "Manif"
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
stimpath = Stats(Expi).meta.stimuli;
% imgN=121; 
elseif ExpType == "Evol"
index_vect = cell2mat(EStats(Expi).evol.idx_seq');
imgnm_vect = EStats(Expi).imageName(index_vect); % reshape(imgnm_grid, [], 1);
stimpath = EStats(Expi).meta.stimuli;
end
imgN = length(imgnm_vect);
tmpfn = ls(fullfile(stimpath, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");suffix = "."+tmpparts{2};
imgfullnm_vect = cellfun(@(imgnm) fullfile(stimpath, imgnm+suffix),imgnm_vect);
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
if ExpType == "Manif"
psth_all = cellfun(@(psth) reshape(mean(psth(ui,:,:),[3]),1,1,[]), Stats(Expi).manif.psth{si}, ...
    'UniformOutput', false); % note there is trial averaging here. 
psth_all = reshape(cell2mat(psth_all),imgN,[]);
elseif ExpType == "Evol"
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; % imgN by 200
end
% This part could be abbrieviated. 
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
end