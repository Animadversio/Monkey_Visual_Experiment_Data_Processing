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

online_compute = false;
%% Load image name and firing data
ExpType = "Evol";
Expi = 11;
fprintf("Processing Evol Exp %d pref chan %d\n",Expi,EStats(Expi).units.pref_chan)
index_vect = cell2mat(EStats(Expi).evol.idx_seq');
imgN=length(index_vect); 
imgnm_vect = EStats(Expi).imageName(index_vect); % reshape(imgnm_grid, [], 1);
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; % imgN by 200
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
% Get Image suffix (Assume all evolved images use the same suffix)
tmpfn = ls(fullfile(EStats(Expi).meta.stimuli, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");
suffix = "."+tmpparts{2};% suffix = ".bmp";

%% Batch computation of correlation coefficient. 
%% Need to have the scores first
layername = "conv4_3";
imgN=length(index_vect); Bsz=128;% make this smaller if the layer is deep! 
score_shuffle = [];
shuffleN = 100;
for i = 1:shuffleN
    score_shuffle = [score_shuffle, score_vect(randperm(imgN),:)];
end

if online_compute % batch computation of correlation coefficient tensor. 
%% initialize the "sufficient statistics"
dummy = activations(net, zeros(224,224,3), layername);
ft_shape = size(dummy); % size(feat_tsr,[1,2,3]);
nfeat = prod(ft_shape); % number of feature predictors in total
ntpnt = size(score_vect,2);
SSqFeat = zeros(nfeat,1,'single');
SFeat = zeros(nfeat,1,'single');
SSqrsp = zeros(ntpnt*shuffleN,1,'single');
Srsp = zeros(ntpnt*shuffleN,1,'single');
InnProd = zeros(nfeat, ntpnt*shuffleN,'single');
%% batch computation of correlation coefficient tensor.
%% Batch by batch
curN = 0;
csr = 1; 
T0_all = tic;
while curN < imgN
    T0 = tic;
    csr = curN + 1; 
    csr_end = min(imgN, csr + Bsz -1);
    rsp = score_shuffle(csr:csr_end,:); % change this into 
    imgcol = cellfun(@(imgnm) imresize(imread(fullfile(EStats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
        imgnm_vect(csr:csr_end), 'UniformOutput', false);
    T1 = toc(T0);
    dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
    feat_tsr_tmp = activations(net, dlimg, layername);
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
cc_tsr = (InnProd./curN - MFeat * Mrsp') ./ (StdFeat * Stdrsp');
cc_tsr = single(reshape(cc_tsr, [ft_shape, ntpnt, shuffleN]));
MFeat = reshape(MFeat, ft_shape);
StdFeat = reshape(StdFeat, ft_shape);

else % Commpute feature tensor over 3000+ images directly using activations. (faster then doing loop)
T0 = tic;
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(EStats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
        imgnm_vect, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
T1 = toc(T0);
feat_tsr = activations(net, dlimg, layername, 'MiniBatchSize', Bsz);
T2 = toc(T0);
% shuffleN = 100;
% score_shuffle = [];
% for i = 1:shuffleN
%     score_shuffle = [score_shuffle, score_vect(randperm(imgN),:)];
% end
cc_tsr_ref = corr(reshape(feat_tsr,nfeat,imgN)', score_shuffle);
cc_tsr_ref = single(reshape(cc_tsr_ref, [size(feat_tsr, [1,2,3]), ntpnt, shuffleN]));
MFeat = mean(feat_tsr,4);
StdFeat = std(feat_tsr,0,4);
T3 = toc(T0);
fprintf("Latencies: load img %.1f CNN proc %.1f compute innprod %.1f\n", T1, T2, T3);

end
%% Load the real correlation coefficients
ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
outfn = fullfile(ccmat_dir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername));
R = load(outfn,'cc_tsr', 'MFeat', 'StdFeat', 'wdw_vect');

%% Correlation Coefficient Clearly shows a spatial structure there
figure(19);
corr_tsr_L1 = squeeze(mean(abs(cc_tsr(:,:,:,:)),3)); % H, W, timefr
CLIM = prctile(cc_tsr, [2.5, 98], 'all'); 
for fi=1:size(cc_tsr,4)
wdw = wdw_vect(fi,:);
imagesc(corr_tsr_L1(:,:,fi));axis image
title(sprintf("%s Manif Exp %d Pref chan %d\nMean CorreCoef in of VGG16 %s feature\n with [%d,%d] ms firing rate", ...
    Animal, Expi, EStats(Expi).units.pref_chan, strrep(layername,"_",'-'), wdw(1), wdw(2)))
caxis(CLIM);colorbar()
pause(0.05)
% saveas(19,fullfile(savedir, compose("%s_Manif_Exp%d_%s_wdw%d.png",Animal,Expi,layername,fi)))
end

%% Directly computing feature tensor 
tic
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(EStats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
        imgnm_vect, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
toc
feat_tsr_tmp = activations(net, dlimg, layername);
toc % takes 103sec in total, but 37 sec in CNN processing, 60+ sec in loading images.
%%
T0 = tic;
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(EStats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
        imgnm_vect, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
%% Commpute feature tensor over 3000+ images directly using activations. (faster then doing loop)
T0 = tic;
% dummy = activations(net, zeros(224,224,3), layername);
% ft_shape = size(dummy); %size(feat_tsr,[1,2,3]);
% nfeat = prod(ft_shape); 
T1 = toc(T0);
feat_tsr = activations(net, dlimg, layername);
T2 = toc(T0);
tic
shuffleN = 100;
score_shuffle = [];
for i = 1:shuffleN
    score_shuffle = [score_shuffle, score_vect(randperm(imgN),:)];
end
cc_tsr = corr(reshape(feat_tsr,nfeat,imgN)', score_shuffle);
cc_tsr = single(reshape(cc_tsr, [size(feat_tsr, [1,2,3]), ntpnt, shuffleN]));
MFeat = mean(feat_tsr,4);
StdFeat = std(feat_tsr,0,4);
toc
T3 = toc(T0);
fprintf("Latencies: load img %.1f CNN proc %.1f compute innprod %.1f\n", T1, T2, T3);
%%
cc_tsr_ref = cc_tsr;
%%
figure(20);
corr_tsr_L1 = squeeze(mean(abs(cc_tsr_ref),3)); % H, W, timefr, rep
CLIM = prctile(corr_tsr_L1, [2.5, 98], 'all'); 
for fi = 11
for ti = 1:size(cc_tsr, 5)
wdw = wdw_vect(fi,:);
imagesc(corr_tsr_L1(:,:,fi,ti));axis image
title(sprintf("%s Manif Exp %d Pref chan %d\nMean CorreCoef in of VGG16 %s feature\n with [%d,%d] ms firing rate", ...
    Animal, Expi, EStats(Expi).units.pref_chan, strrep(layername,"_",'-'), wdw(1), wdw(2)))
caxis(CLIM);colorbar()
pause(0.1)
% saveas(19,fullfile(savedir, compose("%s_Manif_Exp%d_%s_wdw%d.png",Animal,Expi,layername,fi)))
end
end
%%
figure(1);
signif_n = 0;
xi=27; yi=5; fi = 11;wdw = wdw_vect(fi,:);%ci=6; 
for ci = 1:512
nulldist = squeeze(cc_tsr_ref(yi,xi,ci,fi,:));
% pthresh = prctile(nulldist,[1,99]);
pthresh = mean(nulldist) + 5 * std(nulldist) * [-1, 1]; % 5 sigma 1.5E-12
if R.cc_tsr(yi,xi,ci,fi) < pthresh(2) && R.cc_tsr(yi,xi,ci,fi)  > pthresh(1)
    continue
else
    signif_n = signif_n +1;
    pval = erfc(abs(R.cc_tsr(yi,xi,ci,fi) - mean(nulldist)) / std(nulldist));
end
hist(squeeze(cc_tsr_ref(yi,xi,ci,fi,:)))
line([1, 1]*R.cc_tsr(yi,xi,ci,fi),ylim()) 
xlabel("correlation coefficient")
title(compose("Conv4-3 corrcoef VS shuffled corrcoef distribution\n x=%d,y=%d,chan=%d,fi=%d ([%d,%d] ms)\np Val=%.1E",xi,yi,ci,fi, wdw(1), wdw(2),pval))
% pause
end
fprintf("Significantly correlated voxel num %d\n",signif_n)
%% Significance tensor
cc_refM = mean(cc_tsr_ref,5);
cc_refS = std(cc_tsr_ref,0,5);
t_signif_tsr = (R.cc_tsr - cc_refM) ./ cc_refS;
%% thresholding t value to create a mask
t_signif_tsr(t_signif_tsr < 5 & t_signif_tsr > -5) =0;
%% Visualize the t mask. 
t_sum_tsr = mean(abs(t_signif_tsr(:,:,:,:)),3);
CLIM = prctile(t_sum_tsr,[2.5,98],'all');
figure;
for fi=1:size(t_signif_tsr,4)
wdw = wdw_vect(fi,:);
imagesc(t_sum_tsr(:,:,fi));axis image
title(sprintf("%s Manif Exp %d Pref chan %d\nMean CorreCoef in of VGG16 %s feature\n with [%d,%d] ms firing rate", ...
    Animal, Expi, EStats(Expi).units.pref_chan, strrep(layername,"_",'-'), wdw(1), wdw(2)))
caxis(CLIM);colorbar()
pause(0.15)
% saveas(19,fullfile(savedir, compose("%s_Manif_Exp%d_%s_wdw%d.png",Animal,Expi,layername,fi)))
end 
%%