%% Create the summary plots for Compare evolution scores for Beto.
%  See if the Evolution works
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
%% 
Animal = "Beto"; ExpType = "RDEvol"; 
flags = struct("batch",40,"online_compute",1,"load_all_img",1,"shuffleN",100);
wdw_vect = [[1, 20] + 10 * [0:18]'; [1,50]+[0:50:150]'; [51,200]];
for Expi = 1:36
T0 = tic;
[imgfn, score_vect] = loadData(RDStats, "Beto", ExpType, Expi);
totl_vox_num = [];
corr_vox_num = zeros(24,[]);
corr_vox_prct = zeros(24,[]);
med_pos_cc = zeros(24,[]);
med_neg_cc = zeros(24,[]);
mean_pos_t = zeros(24,[]);
mean_neg_t = zeros(24,[]);
layernames = ["fc8", "fc7", "fc6", "conv5_3", "conv4_3", "conv3_3"];
for iLayer = 1:length(layernames)
layername = layernames(iLayer);fprintf("Correlate %s layer\n",layername)
[cc_tsr, MFeat, StdFeat, cc_refM, cc_refS] = corrFeatTsr_func(imgfn, score_vect, net, layername, flags);
outfn = fullfile(savedir, compose("%s_%s_Exp%d_%s_both.mat",Animal,ExpType,Expi,layername)); % LW means long window
save(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
nfeat = prod(size(t_signif_tsr,[1,2,3]));
totl_vox_num(iLayer) = nfeat;
for fi = 1:24
wdw = wdw_vect(fi,:);
tcol = t_signif_tsr(:,:,:,fi);%yi,xi,
cc_tmp = cc_tsr(:,:,:,fi);
signif_n = sum(tcol > 5 | tcol<-5,'all');
corr_vox_num(fi,iLayer) = signif_n;
med_pos_cc(fi,iLayer) = median(cc_tmp(tcol>5));
med_neg_cc(fi,iLayer) = median(cc_tmp(tcol<-5));
mean_pos_t(fi,iLayer) = mean(tcol(tcol>5));
mean_neg_t(fi,iLayer) = mean(tcol(tcol<-5));
end
fprintf("Finished %s (%.1f sec)\n",layername,toc(T0))
end
save(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16_both.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t");
end
%
for Expi = 1:36
[imgfn, score_vect] = loadData(RDStats, "Beto", ExpType, Expi, 1);
totl_vox_num = [];
corr_vox_num = zeros(24,[]);
corr_vox_prct = zeros(24,[]);
med_pos_cc = zeros(24,[]);
med_neg_cc = zeros(24,[]);
mean_pos_t = zeros(24,[]);
mean_neg_t = zeros(24,[]);
layernames = ["fc8", "fc7", "fc6", "conv5_3", "conv4_3", "conv3_3"];
for iLayer = 1:length(layernames)
layername = layernames(iLayer);fprintf("Correlate %s layer\n",layername)
[cc_tsr, MFeat, StdFeat, cc_refM, cc_refS] = corrFeatTsr_func(imgfn, score_vect, net, layername, flags);
outfn = fullfile(savedir, compose("%s_%s_Exp%d_%s_thread1.mat",Animal,ExpType,Expi,layername)); % LW means long window
save(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
nfeat = prod(size(t_signif_tsr,[1,2,3]));
totl_vox_num(iLayer) = nfeat;
for fi = 1:24
wdw = wdw_vect(fi,:);
tcol = t_signif_tsr(:,:,:,fi);%yi,xi,
cc_tmp = cc_tsr(:,:,:,fi);
signif_n = sum(tcol > 5 | tcol<-5,'all');
corr_vox_num(fi,iLayer) = signif_n;
med_pos_cc(fi,iLayer) = median(cc_tmp(tcol>5));
med_neg_cc(fi,iLayer) = median(cc_tmp(tcol<-5));
mean_pos_t(fi,iLayer) = mean(tcol(tcol>5));
mean_neg_t(fi,iLayer) = mean(tcol(tcol<-5));
end
end
save(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16_thread1.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t");
end
%%
% [imgfullnm_vect, score_vect] = loadData(RDStats, Animal, ExpType, Expi, [1,2]);
%%
function [imgfullnm_vect, score_vect] = loadData(RDStats, Animal, ExpType, Expi, threads)
if nargin == 4
threads = 1:RDStats(Expi).evol.thread_num; % specify which thread to extract.
end
assert(RDStats(Expi).Animal == Animal && RDStats(Expi).Animal == Animal)
fprintf("Processing %s Exp %d pref chan %d (thread %s)\n",ExpType,Expi,RDStats(Expi).units.pref_chan,num2str(threads))
% ui=1;si=1;
prefchan_id = find((RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
chid = find((RDStats(Expi).units.unit_num_arr == 1) & (RDStats(Expi).units.spikeID == RDStats(Expi).evol.pref_chan(1)));
ui = find(prefchan_id==chid);

if ExpType == "Manif"
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
stimpath = Stats(Expi).meta.stimuli;
% imgN=121; 
elseif ExpType == "Evol"
index_vect = cell2mat(EStats(Expi).evol.idx_seq');
imgnm_vect = EStats(Expi).imageName(index_vect); % reshape(imgnm_grid, [], 1);
stimpath = EStats(Expi).meta.stimuli;
elseif ExpType == "RDEvol"
% index_vect1 = cell2mat(RDStats(Expi).evol.idx_seq(1,:)');
% index_vect2 = cell2mat(RDStats(Expi).evol.idx_seq(2,:)');
% Reshape concatenate different threads vertically. and concatenate the idx accordingly. 
index_vect = cell2mat(reshape(RDStats(Expi).evol.idx_seq(threads,:)',[],1)) ;
imgnm_vect = RDStats(Expi).imageName(index_vect); % [index_vect1;index_vect2]
stimpath = RDStats(Expi).meta.stimuli;
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
elseif ExpType == "RDEvol"
% psth_all1 = squeeze(cell2mat(reshape(RDStats(Expi).evol.psth(1,:),1,1,[])))';
% psth_all2 = squeeze(cell2mat(reshape(RDStats(Expi).evol.psth(2,:),1,1,[])))';
% psth_all  = cat(1,psth_all1,psth_all2);
psth_ui = cellfun(@(psth)psth(ui,:,:), RDStats(Expi).evol.psth, 'UniformOutput', false); % size, thread_n by block_n
psth_all = squeeze(cell2mat(reshape(psth_ui(threads,:)',1,1,[])))'; % reshape gets all the 
end
% This part could be abbrieviated. 
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
end