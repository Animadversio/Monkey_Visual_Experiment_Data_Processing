addpath CorrFeatTsr
%% untrained Corr Feature Tensor
Animal = "Alfa";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
predsavedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict";
% load(fullfile(predsavedir, Animal+"_FeatTsrPredStats.mat"),'predStats')
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
net = vgg16;
net_un = getUntrainedNet('vgg16');
%%
imgpath = "S:\Stimuli\2019-Selectivity\gabors_41\gab_ori_30.0_1.0.bmp";
img = imresize(imread(imgpath),[224,224]);
img = repmat(img,1,1,3) / 255.0;
rawacts = squeeze(activations(net_un,img,'fc6'));
acts = squeeze(activations(net,img,'fc6'));
figure;scatter(acts(:),rawacts(:))

%%
Animal = "Alfa";ExpType = "Evol";Expi = 29;
[imgfn, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, struct());
%%
layername = "conv3_3";
flags = struct("batch",50,"online_compute",1,"load_all_img",1,"shuffleN",50);
wdw_vect = [[1, 20] + 10 * [0:18]'; [1,50]+[0:50:150]'; [51,200]];
[cc_tsr, MFeat, StdFeat, cc_refM, cc_refS] = corrFeatTsr_func(imgfn, score_vect, net_un, layername, flags);
t_signif_tsr = (cc_tsr - cc_refM)./cc_refS;
%%
thresh="both";sum_method = "L1";
if thresh == "both"
maskTsr = abs(t_signif_tsr)>=5; % Bi sided thresholding 
elseif thresh == "pos"
maskTsr = t_signif_tsr>=5; % positively correlated thresholding
elseif thresh == "neg"
maskTsr = t_signif_tsr<=-5; % negatively correlated thresholding
end
plotTsr = cc_tsr;
plotTsr(~maskTsr) = 0; 
if sum_method=="L1"
L1plotTsr = squeeze(mean(abs(plotTsr),3));
elseif sum_method=="L1signif"
L1plotTsr = squeeze(sum(abs(plotTsr),3)./(sum(maskTsr,3)));
elseif sum_method=="max"
L1plotTsr = squeeze(max(abs(plotTsr),[],3));
end

CLIM_arr = zeros(24,2);
CLIM_arr(1:19,:) = CLIM_arr(1:19,:) + prctile(L1plotTsr(:,:,1:19),[2,98],'all')' + [0, 1E-4];
CLIM_arr(20:23,:) = CLIM_arr(20:23,:) + prctile(L1plotTsr(:,:,20:23),[2,98],'all')'+ [0, 1E-4];
CLIM_arr(24,:) = CLIM_arr(24,:) + prctile(L1plotTsr(:,:,24),[2,98],'all')'+ [0, 1E-4];
%
CLIM_arr(isnan(CLIM_arr)) = 0; % put 0 in the nan place (if there is no activated voxel. then prctile will be all nan)

% if doSave
% v = VideoWriter(fullfile(savedir,compose("%s_%s_Exp%d_%s_cc_%s-%s.mp4",Animal,ExpType,Expi,layername,thresh,sum_method)));
% v.FrameRate = 2;open(v);
% end
figure(13);%(12);
% imagesc(mean(abs(cc_tsr(:,:,:,end)),3));colorbar()
for fi = 1:24
imagesc(L1plotTsr(:,:,fi));axis image;caxis(CLIM_arr(fi,:));colorbar();
title(compose("%s %s Exp%d Corr Coef with untrained VGG16 %s\n rate in [%d,%d]ms\n Channel compressed with %s thresholded %s",...
    Animal,ExpType,Expi,layername,wdw_vect(fi,1),wdw_vect(fi,2),thresh,sum_method))
pause(0.2)
% if doSave,writeVideo(v,getframe(13));end
drawnow;
end
% if doSave,close(v);end
%%
load(fullfile(MatStats_path, "Alfa_Evol_ccFtMask.mat"))
%%
trainedccMsk = ccMskStats(Expi).(ExpType).(layername).(sum_method);
untrainedccMsk = L1plotTsr(:,:,end);
ccbtwnet = corr(trainedccMsk(:), untrainedccMsk(:));% correlation of masks between trained and untrained network
%%
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