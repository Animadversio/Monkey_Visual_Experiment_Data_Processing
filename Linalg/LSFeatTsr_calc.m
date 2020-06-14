%% LeastSquare FeatTsr
Animal="Alfa"; 
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%

[imgfullnm_vect, score_Manif] = loadManifData(Stats, EStats, Animal, "Manif", 3, struct());
%%
feat_tsr = FeatTsrImgs("conv5_3", imgfullnm_vect);
%% Compute least square manually
% Gram = einsum(feat_tsr, feat_tsr, 'ijkl,ijml->ijkm')
tic
Gram_pix = zeros([size(feat_tsr,[1,2,3]),size(feat_tsr,3)]);
for i = 1:size(feat_tsr,1)
    for j = 1:size(feat_tsr,2)
        Gram_pix(i,j,:,:) = squeeze(feat_tsr(i,j,:,:)) * squeeze(feat_tsr(i,j,:,:))';
    end
end
toc
%%
InnProd = einsum(feat_tsr, score_Manif, 'ijkl,lm->ijkm');
%%
pixWtsr = zeros([size(feat_tsr,[1,2,3])]);
rsqmat = zeros([size(feat_tsr,[1,2])]);
xi = 10; yi = 10; fi = 24;
Lambda = 100000000;
tic
for xi = 1:size(feat_tsr,1)
    for yi = 1:size(feat_tsr,2)
    W = (squeeze(Gram_pix(xi,yi,:,:)) + Lambda * eye(512))\ squeeze(InnProd(xi,yi,:,fi));
    predfr = squeeze(feat_tsr(xi,yi,:,:))'*W;
    pixWtsr(xi,yi,:) = W;
    rsqmat(xi,yi) = 1-var(score_Manif(:,end) - predfr) / var(score_Manif(:,end));
    end
end
toc
%%
MSEMap = zeros([size(feat_tsr,[1,2])]);
LassoWTsr = zeros(size(feat_tsr,[1,2,3]));
xi=18;yi=17;
tic
for xi = 1:size(feat_tsr,1)
    for yi = 1:size(feat_tsr,2)
    [B,LassoStats] = lasso(squeeze(feat_tsr(xi,yi,:,:))',score_Manif(:,end),'Lambda',10);
    MSEMap(xi,yi)=LassoStats.MSE;
    LassoWTsr(xi,yi,:)=B;
    end
    toc
end
toc
LassoRsqMap = 1- MSEMap./var(score_Manif(:,end));
%
figure;
imagesc(LassoRsqMap)%rsqmat
%%
xi = 10;yi = 10;fi = 24;
figure;
% imagesc(mean(abs(InnProd(3:end-2,3:end-2,:,fi)),3))
imagesc(mean(abs(InnProd(:,:,:,fi)),3));axis image
%%
function [feat_tsr] = FeatTsrImgs(layername, imgfullnms)
global net
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(imgnm)),[224,224]), imgfullnms, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[])); 
feat_tsr = activations(net, dlimg, layername,'MiniBatchSize',40);
% pred_rsp = einsum(feat_tsr, weight_tsr, 'ijkl,ijkm->lm') / prod(size(feat_tsr, [1,2,3]));
% pred_rsp = mean(feat_tsr .* weight_tsr,[1,2,3]);
end

function [cc_tsr, t_signif_tsr] = load_cc_tsr(Animal, ExpType, Expi, layername)
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
outfn = fullfile(savedir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername)); % LW means long window
load(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
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