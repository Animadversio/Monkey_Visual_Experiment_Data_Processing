%% 
% Plan of this script is quite similar to pre_stat.m

global net
net = vgg16;
global ft
ft = fittype( 'max(0, a*(x-b))+c', 'independent', 'x', 'dependent', 'y' );
%%
Animal="Beto";
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
predsavedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict";
% load(fullfile(predsavedir, Animal+"_FeatTsrPredStats.mat"),'predStats')
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')

%%
load(fullfile(MatStats_path,Animal+"_FeatTsrPredStats2.mat"),'predStats2')
diary(fullfile(predsavedir,"Beto_fitting_process.log"))
% predStats2 = repmat(struct("E2M",struct(),"M2E",struct()),1,length(Stats));
%%
Animal="Beto"; Expi = 3; 
%%
for Expi = 21:numel(Stats)
[imgfullnm_Evol, score_Evol] = loadPairedData(Stats, EStats, Animal, "Evol", Expi, struct()); % uniform interface to load spikes and data
[imgfullnm_Manif, score_Manif] = loadPairedData(Stats, EStats, Animal, "Manif", Expi, struct());
[imgfullnm_EvoRef, score_EvoRef] = loadPairedData(Stats, EStats, Animal, "EvoRef", Expi, struct());
[imgfullnm_Gabor, score_Gabor] = loadPairedData(Stats, EStats, Animal, "Gabor", Expi, struct());
[imgfullnm_Pasu, score_Pasu] = loadPairedData(Stats, EStats, Animal, "Pasu", Expi, struct());
%%
for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
fprintf("Using the correlation tensor from vgg16 %s layer\n", layername)
fprintf("Fitting Manif to Evol nonlinearity\n")
[cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, "Manif", Expi, layername); % correlation coefficient for manifold experiments
ccWeight = cc_tsr;ccWeight(isnan(cc_tsr))=0; % there can be nans .....
thresh = prctile(t_signif_tsr(:,:,:,end),99,'all');
fprintf("Threshold %.3f (%d)\n",thresh,sum(t_signif_tsr(:,:,:,end) > thresh,'all'))
ccWeight(t_signif_tsr < thresh) = 0; %threshold and get a weight tensor to do prediction abs
[S_manif,relufit] = fitrelu_get_stats(score_Manif, imgfullnm_Manif, "Manif", ccWeight, layername);
predStats2(Expi).M2E.(layername).manif = S_manif;
%%
[S_evol] = pred_get_stats(score_Evol, imgfullnm_Evol, "Evol", ccWeight, layername, relufit);
predStats2(Expi).M2E.(layername).evol = S_evol;
%%
[S_evoref] = pred_get_stats(score_EvoRef, imgfullnm_EvoRef, "EvoRef", ccWeight, layername, relufit);
predStats2(Expi).M2E.(layername).evoref = S_evoref;
%%
[S_gabor] = pred_get_stats(score_Gabor, imgfullnm_Gabor, "Gabor", ccWeight, layername, relufit);
predStats2(Expi).M2E.(layername).gabor = S_gabor;
[S_pasu] = pred_get_stats(score_Pasu, imgfullnm_Pasu, "Pasu", ccWeight, layername, relufit);
predStats2(Expi).M2E.(layername).pasu = S_pasu;
end
end
diary off


%%
save(fullfile(MatStats_path,Animal+"_FeatTsrPredStats2.mat"),'predStats2')
diary off
%%
save(fullfile(MatStats_path,Animal+"_FeatTsrPredStats2.mat"),'predStats2')
%%
diary("E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict\fitting_process.log")
% predStats2 = repmat(struct("E2M",struct(),"M2E",struct()),1,length(Stats));
%%
Animal="Beto"; %Expi = 3; 
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
predsavedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Predict";
% load(fullfile(predsavedir, Animal+"_FeatTsrPredStats.mat"),'predStats')
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
% load(fullfile(MatStats_path,Animal+"_FeatTsrPredStats2.mat"),'predStats2')
predStats2 = repmat(struct("E2M",struct(),"M2E",struct()),1,length(Stats));
%%
for Expi = 1:numel(Stats)
[imgfullnm_Evol, score_Evol] = loadPairedData(Stats, EStats, Animal, "Evol", Expi, struct()); % uniform interface to load spikes and data
[imgfullnm_Manif, score_Manif] = loadPairedData(Stats, EStats, Animal, "Manif", Expi, struct());
[imgfullnm_EvoRef, score_EvoRef] = loadPairedData(Stats, EStats, Animal, "EvoRef", Expi, struct());
[imgfullnm_Gabor, score_Gabor] = loadPairedData(Stats, EStats, Animal, "Gabor", Expi, struct());
[imgfullnm_Pasu, score_Pasu] = loadPairedData(Stats, EStats, Animal, "Pasu", Expi, struct());
%%
for layername = ["conv3_3","conv4_3","conv5_3","fc6","fc7","fc8"] %
fprintf("Using the correlation tensor from vgg16 %s layer\n", layername)
fprintf("Fitting Manif to Evol nonlinearity\n")
[cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, "Manif", Expi, layername); % correlation coefficient for manifold experiments
ccWeight = cc_tsr;ccWeight(isnan(cc_tsr))=0; % there can be nans .....
thresh = prctile(t_signif_tsr(:,:,:,end),99,'all');
fprintf("Threshold %.3f (%d)\n",thresh,sum(t_signif_tsr(:,:,:,end) > thresh,'all'))
ccWeight(t_signif_tsr < thresh) = 0; %threshold and get a weight tensor to do prediction abs

[S_evol,relufit] = fitrelu_get_stats(score_Evol, imgfullnm_Evol, "Evol", ccWeight, layername);
predStats2(Expi).E2M.(layername).evol = S_evol;
%%
[S_manif] = pred_get_stats(score_Manif, imgfullnm_Manif, "Manif", ccWeight, layername, relufit);
predStats2(Expi).E2M.(layername).manif = S_manif;
[S_evoref] = pred_get_stats(score_EvoRef, imgfullnm_EvoRef, "EvoRef", ccWeight, layername, relufit);
predStats2(Expi).E2M.(layername).evoref = S_evoref;
[S_gabor] = pred_get_stats(score_Gabor, imgfullnm_Gabor, "Gabor", ccWeight, layername, relufit);
predStats2(Expi).E2M.(layername).gabor = S_gabor;
[S_pasu] = pred_get_stats(score_Pasu, imgfullnm_Pasu, "Pasu", ccWeight, layername, relufit);
predStats2(Expi).E2M.(layername).pasu = S_pasu;
end
end
diary off
%%
save(fullfile(MatStats_path,Animal+"_FeatTsrPredStats2.mat"),'predStats2')


%% Util functions 
function [pred_rsp] = LinPredictImgs(weight_tsr, layername, imgfullnms)
global net
imgcol = cellfun(@(imgnm) imresize(to_rgb(imread(fullfile(imgnm))),[224,224]), imgfullnms, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[])); 
if size(dlimg,3)==1,dlimg = repmat(dlimg,1,1,3,1);end
feat_tsr = activations(net, dlimg, layername,'MiniBatchSize',40);
pred_rsp = einsum(feat_tsr, weight_tsr, 'ijkl,ijkm->lm') / prod(size(feat_tsr, [1,2,3]));
% pred_rsp = mean(feat_tsr .* weight_tsr,[1,2,3]);
end

%% pred_get_stats: function description
function [S] = pred_get_stats(scores_mat, imageNames, imglabel, ccWeight, layername, relufit)
if ~(isempty(scores_mat) || isempty(imageNames))
tic
pred_rsp = LinPredictImgs(ccWeight, layername, imageNames);
NLpred_rsp = relufit(pred_rsp);
NLpred_rsp = reshape(NLpred_rsp, size(scores_mat));
[nlpredcorr, nlpredcorr_P] = corr(NLpred_rsp(:,end),scores_mat(:,end));
[lpredcorr, lpredcorr_P] = corr(pred_rsp(:,end),scores_mat(:,end));
nlpredR2 = 1 - var(NLpred_rsp(:,end)-scores_mat(:,end)) / var(scores_mat(:,end));
fprintf("Predicting %s Images, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
    imglabel,lpredcorr,lpredcorr_P,nlpredcorr,nlpredcorr_P,nlpredR2,numel(scores_mat(:,end)))
toc
S.lpredscore = pred_rsp;
S.nlpredscore = NLpred_rsp;
S.nlpredcorr = nlpredcorr;
S.nlpredcorr_P = nlpredcorr_P;
S.lpredcorr = lpredcorr;
S.lpredcorr_P = lpredcorr_P;
S.nlpredR2 = nlpredR2;
S.origscore = scores_mat;
else
S.lpredscore = [];
S.nlpredscore = [];
S.nlpredcorr = nan;
S.nlpredcorr_P = nan;
S.lpredcorr = nan;
S.lpredcorr_P = nan;
S.nlpredR2 = nan;
S.origscore = [];
end
end

function [S,relufit] = fitrelu_get_stats(scores_mat, imageNames, imglabel, ccWeight, layername)
global ft % relu function
tic
lfitscore = LinPredictImgs(ccWeight, layername, imageNames);
opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
        'Lower',[0, min(lfitscore(:,end)), 0],...
        'Upper',[10000, max(lfitscore(:,end)), max(scores_mat(:,end))]);
[relufit, gof, out] = fit(lfitscore(:,end), scores_mat(:,end), ft, opts); 
% fit the bias on 
NLpred = relufit(lfitscore);
NLpred = reshape(NLpred, size(lfitscore));

[nlfitcorr, nlfitcorr_P] = corr(NLpred(:,end),scores_mat(:,end));
[lfitcorr, lfitcorr_P] = corr(lfitscore(:,end),scores_mat(:,end));
nlpredR2 = 1 - var(NLpred(:,end)-scores_mat(:,end)) / var(scores_mat(:,end));
fprintf("Fitting %s Images, Linear corr %.3f(%.1e), Nonlinear corr %.3f(%.1e) R2 %.3f n=%d\n",...
    imglabel,lfitcorr,lfitcorr_P,nlfitcorr,nlfitcorr_P,nlpredR2,numel(scores_mat(:,end)))
toc
S.origscore = scores_mat;
S.lfitscore = lfitscore;
S.nlfitscore = NLpred;
S.nlfitcorr = nlfitcorr;
S.nlfitcorr_P = nlfitcorr_P;
S.lfitcorr = lfitcorr;
S.lfitcorr_P = lfitcorr_P;
S.nlfitR2 = nlpredR2;
S.nlfunc = relufit;
S.gof = gof;
end

function img = to_rgb(img)
if size(img,3)==1
   img = repmat(img,1,1,3);
end
end

function [cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, ExpType, Expi, layername)
% uniform interface to load the correlation tensor pre-computed with VGG activations. 
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
outfn = fullfile(savedir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername)); % LW means long window
load(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
end

function [imgfullnm_vect, score_vect] = loadPairedData(Stats, EStats, Animal, ExpType, Expi, flags)
% Uniform interface to load spikes and image path data, really well done!
% Maybe we should load some internal variance statistics to regularize the fit...
ui=1;si=1;
assert(EStats(Expi).Animal == Animal && Stats(Expi).Animal == Animal)
fprintf("Processing %s Exp %d pref chan %d\n",ExpType,Expi,EStats(Expi).units.pref_chan)
if ExpType == "Manif"
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
stimpath = Stats(Expi).meta.stimuli;
psth_all = cell2mat(cellfun(@(psth) mean(psth(ui,:,:),[3]), reshape(Stats(Expi).manif.psth{si},[],1),'Uni', 0)); ...);
% note there is trial averaging here!
% psth_all = reshape(cell2mat(psth_all),imgN,[]);
% imgN=121; 
elseif ExpType == "Evol"
index_vect = cell2mat(EStats(Expi).evol.idx_seq'); % concat the index into a vertical vect
imgnm_vect = EStats(Expi).imageName(index_vect); % reshape(imgnm_grid, [], 1);
stimpath = EStats(Expi).meta.stimuli;
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; 
elseif ExpType == "EvoRef"
imgnm_vect = EStats(Expi).ref.imgnm;
imgfullnm_vect = EStats(Expi).ref.impaths;
psth_all = cell2mat(cellfun(@(psth)mean(psth(ui,:,:),3),reshape(EStats(Expi).ref.psth_arr,[],1),'uni',0));
elseif ExpType == "Gabor"
if ~Stats(Expi).ref.didGabor,  
imgfullnm_vect = []; score_vect = []; return
else,
idx_vect = reshape(Stats(Expi).ref.gab_idx_grid,[],1);
validmsk = ~cellfun(@isempty,idx_vect);
idx_vect = idx_vect(validmsk);
imgnm_vect = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), idx_vect));
stimpath = "N:\Stimuli\2020-manifold-references";%Stats(Expi).meta.stimuli;
psth_all = cell2mat(cellfun(@(psth)mean(psth(ui,:,:),3),reshape(Stats(Expi).ref.gab_psths,[],1),'uni',0));
psth_all = psth_all(validmsk,:);
end
elseif ExpType == "Pasu"
if ~Stats(Expi).ref.didPasu, 
imgfullnm_vect = []; score_vect = []; return
else,
idx_vect = reshape(Stats(Expi).ref.pasu_idx_grid',[],1);
validmsk = ~cellfun(@isempty,idx_vect);
idx_vect = idx_vect(validmsk);
imgnm_vect = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), idx_vect));
stimpath = "N:\Stimuli\2020-manifold-references";%Stats(Expi).meta.stimuli;
psth_all = cell2mat(cellfun(@(psth)mean(psth(ui,:,:),3),reshape(Stats(Expi).ref.pasu_psths',[],1),'uni',0));
psth_all = psth_all(validmsk,:);
end
end
imgN = length(imgnm_vect);
if ~(ExpType == "EvoRef")
% if contains(stimpath,"N:\"),stimpath = strrep(stimpath,"N:\","S:\");end
% if contains(stimpath,"\\storage1.ris.wustl.edu\crponce\Active\"),
%     stimpath = strrep(stimpath,"\\storage1.ris.wustl.edu\crponce\Active\","S:\");end
tmpfn = ls(fullfile(stimpath, imgnm_vect(1)+"*"));
if isempty(tmpfn), error("Stimpath not correct, examine!");end
tmpparts = split(tmpfn,".");suffix = "."+tmpparts{end};
imgfullnm_vect = cellfun(@(imgnm) fullfile(stimpath, imgnm+suffix), imgnm_vect);
end
% imgfullnm_vect = arrayfun(@(impath) strrep(impath,"N:\","S:\"), imgfullnm_vect);
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
% psth_all shape: imgN by 200
% This part could be abbrieviated. 
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
end