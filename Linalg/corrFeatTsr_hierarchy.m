net = vgg16;
%% corrFeatTsr shuffle test 
% The code shuffle the trials 
% For Manifold experiment, the image number could be handled in a single
% batch, so don't need online batch computation of correlation value. 
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
online_compute = true;
%% Load image name and firing data
ExpType = "Evol";
Expi = 11;
fprintf("Processing Evol Exp %d pref chan %d\n",Expi,EStats(Expi).units.pref_chan)
index_vect = cell2mat(EStats(Expi).evol.idx_seq');
imgN=length(index_vect); 
imgnm_vect = EStats(Expi).imageName(index_vect); % reshape(imgnm_grid, [], 1);
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; % imgN by 200
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
% Get Image suffix (Assume all evolved images use the same suffix)
tmpfn = ls(fullfile(EStats(Expi).meta.stimuli, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");
suffix = "."+tmpparts{2}; % suffix = ".bmp";
%% Significantly correlated "pixels" across hierachy
T0 = tic;
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(EStats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
        imgnm_vect, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
toc(T0); 
%% Load image name and firing data
ExpType = "Manif";
Expi = 11;
si=1;ui=1;
fprintf("Processing Manif Exp %d pref chan %d\n",Expi,EStats(Expi).units.pref_chan)
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
imgN=length(imgnm_vect); 
psth_all = cellfun(@(psth) reshape(mean(psth(ui,:,:),[3]),1,1,[]), Stats(Expi).manif.psth{si}, ...
    'UniformOutput', false);
psth_all = reshape(cell2mat(psth_all),imgN,[]);
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
% Get Image suffix (Assume all evolved images use the same suffix)
tmpfn = ls(fullfile(Stats(Expi).meta.stimuli, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");
suffix = "."+tmpparts{2}; % suffix = ".bmp";
%% Significantly correlated "pixels" across hierachy
T0 = tic;
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(strrep(Stats(Expi).meta.stimuli,"N:\","S:\"), imgnm+suffix)),[224,224]), ...
        imgnm_vect, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
toc(T0); 
%% Prepare to collect statistics across the hierarchy.
corr_vox_num = zeros(24,[]);
med_pos_cc = zeros(24,[]);
med_neg_cc = zeros(24,[]);
mean_pos_t = zeros(24,[]);
mean_neg_t = zeros(24,[]);
layernames = ["fc8", "fc7", "fc6", "conv5_3", "conv4_3", "conv3_3"];%,"conv2_2"]; 
%% 
for iLayer = 1:length(layernames)
layername = layernames(iLayer);% 
% Shuffle the thing first and then compute
Bsz = 60;
shuffleN = 100;
score_shuffle = [score_vect]; % put the real scores at first place in the trials, to compute correlation together. 
for i = 1:shuffleN
    score_shuffle = [score_shuffle, score_vect(randperm(imgN),:)];
end
%
if online_compute % batch computation of correlation coefficient tensor. 
%% initialize the "sufficient statistics"
dummy = activations(net, zeros(224,224,3), layername);
ft_shape = size(dummy); % size(feat_tsr,[1,2,3]);
nfeat = prod(ft_shape); % number of feature predictors in total
ntpnt = size(score_vect,2);
SSqFeat = zeros(nfeat,1);
SFeat = zeros(nfeat,1);
SSqrsp = zeros(ntpnt*(shuffleN+1),1);
Srsp = zeros(ntpnt*(shuffleN+1),1);
InnProd = zeros(nfeat, ntpnt*(shuffleN+1));
% Batch by batch
curN = 0;
csr = 1; 
T0_all = tic;
while curN < imgN
    T0 = tic;
    csr = curN + 1; 
    csr_end = min(imgN, csr + Bsz -1);
    rsp = score_shuffle(csr:csr_end,:); % change this into 
%     imgcol = cellfun(@(imgnm) imresize(imread(fullfile(EStats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
%         imgnm_vect(csr:csr_end), 'UniformOutput', false);
    T1 = toc(T0);
%     dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
    feat_tsr_tmp = activations(net, dlimg(:,:,:,csr:csr_end), layername);
    T2 = toc(T0);
    % 
    feat_tsr_tmp = reshape(feat_tsr_tmp, nfeat, []); % [nfeat, batchsize]
    SSqFeat = SSqFeat + sum(feat_tsr_tmp.^2, 2);
    SFeat = SFeat + sum(feat_tsr_tmp, 2);
    SSqrsp = SSqrsp + sum(rsp.^2, 1)'; % summation over samples
    Srsp = Srsp + sum(rsp, 1)';
    InnProd = InnProd + feat_tsr_tmp*rsp;
    curN = csr_end;
    T3 = toc(T0);
    fprintf("Latencies: load img %.1f CNN proc %.1f compute innprod %.1f\n", T1, T2, T3);
end
toc(T0_all)
% Compute mean, and std from the running stats. 
Mrsp = Srsp / curN;
MFeat = SFeat / curN;
Stdrsp = sqrt(SSqrsp / curN - Mrsp.^2);
StdFeat = sqrt(SSqFeat / curN - MFeat.^2);
cc_tsr_ref = (InnProd./curN - MFeat * Mrsp') ./ (StdFeat * Stdrsp');
cc_tsr_ref = reshape(cc_tsr_ref, [ft_shape, ntpnt, shuffleN+1]);
cc_tsr = cc_tsr_ref(:,:,:,:,1);
cc_tsr_ref = cc_tsr_ref(:,:,:,:,2:end);
MFeat = reshape(MFeat, ft_shape);
StdFeat = reshape(StdFeat, ft_shape);
else
T0 = tic;
T1 = toc(T0);
feat_tsr = activations(net, dlimg, layername, 'MiniBatchSize', Bsz); % higher in hierachy need smaller batch 
T2 = toc(T0);
%
nfeat = prod(size(feat_tsr, [1,2,3]));
ntpnt = size(score_vect, 2);
tic
cc_tsr_ref = corr(reshape(feat_tsr,nfeat,imgN)', score_shuffle);
cc_tsr_ref = single(reshape(cc_tsr_ref, [size(feat_tsr, [1,2,3]), ntpnt, shuffleN + 1]));
cc_tsr = cc_tsr_ref(:,:,:,:,1);
cc_tsr_ref = cc_tsr_ref(:,:,:,:,2:end);
MFeat = mean(feat_tsr,4);
StdFeat = std(feat_tsr,0,4);
T3 = toc(T0);
fprintf("Latencies: load img %.1f CNN proc %.1f compute innprod %.1f\n", T1, T2, T3);
end
% Compute the t-value tensor
cc_refM = mean(cc_tsr_ref,5);
cc_refS = std(cc_tsr_ref,0,5);
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
% View the number of signicantly correlated voxels as a function of time window.
fprintf("Correlation layer %s (vox # %d)\n",layername,nfeat)
yi=1; xi=1; fi = 24; wdw = wdw_vect(fi,:);%ci=6; 
for fi = 1:24
wdw = wdw_vect(fi,:);
tcol = t_signif_tsr(:,:,:,fi);%yi,xi,
cc_tmp = cc_tsr(:,:,:,fi);
signif_n = sum(tcol > 5 | tcol<-5,'all');
fprintf("Firing rate in [%d, %d] ms Signif corr voxel num %d (%.1f), pos corr median %.3f (t=%.2f), neg corr median %.3f (t=%.2f)\n",...
    wdw(1),wdw(2),signif_n,signif_n/nfeat*100,median(cc_tmp(tcol>5)),mean(tcol(tcol>5)),...
    median(cc_tmp(tcol<-5)),mean(tcol(tcol<-5)))
corr_vox_num(fi,iLayer) = signif_n;
med_pos_cc(fi,iLayer) = median(cc_tmp(tcol>5));
med_neg_cc(fi,iLayer) = median(cc_tmp(tcol<-5));
mean_pos_t(fi,iLayer) = mean(tcol(tcol>5));
mean_neg_t(fi,iLayer) = mean(tcol(tcol<-5));
end
end
% signif_n = 0;
% tcol = [];
% for ci = 1:size(cc_tsr_ref,3)
% nulldist = squeeze(cc_tsr_ref(1,1,ci,fi,:));
% % pthresh = prctile(nulldist,[1,99]);
% pthresh = mean(nulldist) + 5 * std(nulldist) * [-1, 1]; % 5 sigma 1.5E-12
% if cc_tsr(yi,xi,ci,fi) < pthresh(2) && cc_tsr(yi,xi,ci,fi)  > pthresh(1)
%     continue
% else
%     signif_n = signif_n +1;
%     tval = (cc_tsr(yi,xi,ci,fi) - mean(nulldist)) / std(nulldist);
%     tcol = [tcol, tval];
% end
% end
%% 
figure(11);
signif_n = 0;
yi=1;xi=1;fi = 11; wdw = wdw_vect(fi,:);%ci=6; 
for ci = 1:size(cc_tsr_ref,3)
nulldist = squeeze(cc_tsr_ref(1,1,ci,fi,:));
% pthresh = prctile(nulldist,[1,99]);
pthresh = mean(nulldist) + 5 * std(nulldist) * [-1, 1]; % 5 sigma 1.5E-12
if cc_tsr(yi,xi,ci,fi) < pthresh(2) && cc_tsr(yi,xi,ci,fi)  > pthresh(1)
    continue
else
    signif_n = signif_n +1;
    pval = erfc(abs(cc_tsr(yi,xi,ci,fi) - mean(nulldist)) / std(nulldist));
end
hist(squeeze(cc_tsr_ref(yi,xi,ci,fi,:)))
line([1, 1]*cc_tsr(yi,xi,ci,fi),ylim()) 
xlabel("correlation coefficient")
title(compose("%s corrcoef VS shuffled corrcoef distribution\n x=%d,y=%d,chan=%d,fi=%d ([%d,%d] ms)\np Val=%.1E",layername,xi,yi,ci,fi, wdw(1), wdw(2),pval))
pause
end
fprintf("Significantly correlated voxel num %d\n",signif_n)

%%
corr_vox_num = [];
corr_vox_prct = [];
%% significantly correlated voxel number
layernames = ["fc8", "fc7", "fc6", "conv5-3", "conv4-3", "conv3-3", "conv2-2"];
totl_vox_num = [1000, 4096, 4096, 100352, 401408, 802816, 1605632];
corr_vox_num = reshape(corr_vox_num,24,[]);
corr_vox_prct = corr_vox_num ./ totl_vox_num;
med_pos_cc = reshape(med_pos_cc,24,[]);
med_neg_cc = reshape(med_neg_cc,24,[]);
%%
figure(2);clf;hold on 
subplot(131)
plot(1:24,corr_vox_prct(1:24,:))
legend(layernames, 'Location',"Best")
title("Correlated Voxel Percents")
subplot(132)
plot(1:24,med_pos_cc(1:24,:))
title("Median Positive Correlated Coefficient")
subplot(133)
plot(1:24,med_neg_cc(1:24,:))
title("Median Negative Correlated Coefficient")
%%
figure(2);clf;hold on 

suptitle(compose("Beto Evol Exp11 VGG16 Correlation Percent and cc Distribution"))
subplot(131)
plot([1:19,NaN,20:23,NaN,24],[corr_vox_prct(1:19,:);nan(1,7);corr_vox_prct(20:23,:);nan(1,7);corr_vox_prct(24,:)], "-o")
legend(layernames, 'Location',"Best")
title("Correlated Voxel Percents")
subplot(132)
plot([1:19,NaN,20:23,NaN,24],[med_pos_cc(1:19,:);nan(1,7);med_pos_cc(20:23,:);nan(1,7);med_pos_cc(24,:)], "-o")
legend(layernames, 'Location',"Best")
title("Median Positive Correlated Coefficient")
subplot(133)
plot([1:19,NaN,20:23,NaN,24],[med_neg_cc(1:19,:);nan(1,7);med_neg_cc(20:23,:);nan(1,7);med_neg_cc(24,:)], "-o")
legend(layernames, 'Location',"Best")
title("Median Negative Correlated Coefficient")
saveas(fullfile(savedir,compose("Beto_Evol_Exp11_VGG16_cc.jpg")))
%%
figure(2);clf;hold on 

suptitle(compose("Beto Evol Exp11 VGG16 Correlation Percent and cc Distribution"))
subplot(131)
plot([1:19,NaN,20:23,NaN,24],[corr_vox_prct(1:19,:);nan(1,6);corr_vox_prct(20:23,:);nan(1,6);corr_vox_prct(24,:)], "-o")
legend(layernames, 'Location',"Best")
title("Correlated Voxel Percents")
subplot(132)
plot([1:19,NaN,20:23,NaN,24],[med_pos_cc(1:19,:);nan(1,6);med_pos_cc(20:23,:);nan(1,6);med_pos_cc(24,:)], "-o")
legend(layernames, 'Location',"Best")
title("Median Positive Correlated Coefficient")
subplot(133)
plot([1:19,NaN,20:23,NaN,24],[med_neg_cc(1:19,:);nan(1,6);med_neg_cc(20:23,:);nan(1,6);med_neg_cc(24,:)], "-o")
legend(layernames, 'Location',"Best")
title("Median Negative Correlated Coefficient")
saveas(2, fullfile(savedir,compose("Beto_Manif_Exp11_VGG16_cc.jpg")))
%%
savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
save(fullfile(savedir, compose("Beto_Evol_Exp11_VGG16.mat")),"layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc");
%%
save(fullfile(savedir, compose("Beto_Manif_Exp11_VGG16.mat")),"layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc");