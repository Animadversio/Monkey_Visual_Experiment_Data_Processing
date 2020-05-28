net = vgg16;
%%
Animal = "Alfa";
MatStats_path = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\CNNFeatCorr";
%% Manifold Experiment
imgN=121; %B=10;
for Expi = 2:45
fprintf("Processing Manif Exp %d pref chan %d\n",Expi,EStats(Expi).units.pref_chan)
si=1;ui=1;%Window=50:200;
% dlimg = randn(224,224,3,imgN);
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
% Only one single time slice
% score_grid = cellfun(@(psth) mean(psth(ui,Window,:),[2,3]), Stats(Expi).manif.psth{si});
% score_vect = reshape(score_grid, imgN, []);
% Below add more time slices into it
psth_all = cellfun(@(psth) reshape(mean(psth(ui,:,:),[3]),1,1,[]), Stats(Expi).manif.psth{si}, ...
    'UniformOutput', false);
psth_all = reshape(cell2mat(psth_all),imgN,[]);
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window average 
score_vect = score_vect(:,1:10:end); % subsample to decrease redunancy
wdw_vect=[1, 20] + 10 * [0:18]';
%%
tmpfn = ls(fullfile(Stats(Expi).meta.stimuli, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");
suffix = "."+tmpparts{2};
%% Rearrange the score and corresponding image name
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(Stats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
        imgnm_vect, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[])); % put all the images along the 3rd dim
%% Correlation Coefficient Clearly shows a spatial structure there
savedir = fullfile(result_dir,compose("%s_Manif_Exp%d",Animal,Expi));
mkdir(savedir);
for layername = ["conv3_1", "conv4_3", "conv5_3"]
% layername = 'conv4_3';
T0 = tic;
% dummy = activations(net, zeros(224,224,3), layername);
% ft_shape = size(dummy); %size(feat_tsr,[1,2,3]);
% nfeat = prod(ft_shape); 
T1 = toc(T0);
feat_tsr = activations(net, dlimg, layername);
T2 = toc(T0);
corr_tsr = corr(reshape(feat_tsr,[],imgN)', score_vect);
corr_tsr = reshape(corr_tsr, [size(feat_tsr, [1,2,3]), size(score_vect, 2)]);
MFeat = mean(feat_tsr,4);
StdFeat = std(feat_tsr,0,4);
T3 = toc(T0);
fprintf("Latencies: load img %.1f CNN proc %.1f compute innprod %.1f\n", T1, T2, T3);
%
% figure;
% imagesc(mean(abs(corr_tsr),3));axis image
% title(sprintf("Mean CorreCoef in with VGG16 %s feature", layername))
% colorbar()
outfn = fullfile(result_dir, compose("%s_Manif_Exp%d_%s.mat",Animal,Expi,layername)); % LW means long window
save(outfn,'corr_tsr','MFeat','StdFeat','wdw_vect');
%% Correlation Coefficient Clearly shows a spatial structure there
figure(19);
corr_tsr_L1 = squeeze(mean(abs(corr_tsr(:,:,:,:)),3)); % H, W, timefr
CLIM = prctile(corr_tsr_L1, [2.5, 98], 'all'); 
for fi=1:size(corr_tsr,4)
wdw = wdw_vect(fi,:);
imagesc(corr_tsr_L1(:,:,fi));axis image
title(sprintf("%s Manif Exp %d Pref chan %d\nMean CorreCoef in of VGG16 %s feature\n with [%d,%d] ms firing rate", ...
    Animal, Expi, EStats(Expi).units.pref_chan, strrep(layername,"_",'-'), wdw(1), wdw(2)))
caxis(CLIM);colorbar()
saveas(19,fullfile(savedir, compose("%s_Manif_Exp%d_%s_wdw%d.png",Animal,Expi,layername,fi)))
end
end

end
%%

%%
Animal = "Beto";
MatStats_path = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
result_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\CNNFeatCorr";
result_dir = "S:\CNNFeatCorr";
mkdir(result_dir)
%%
if Animal == "Alfa"
EStats(1).meta.stimuli = "N:\Stimuli\2019-Manifold\alfa-191119a\backup_11_19_2019_11_58_11";
EStats(2).meta.stimuli = "N:\Stimuli\2019-Manifold\alfa-191120a\backup_11_20_2019_15_21_40";
EStats(7).meta.stimuli = "N:\Stimuli\2019-Manifold\alfa-191125a\backup_11_25_2019_12_43_48";
EStats(10).meta.stimuli = "N:\Stimuli\2019-Manifold\alfa-191127a\backup_11_27_2019_12_34_04";
EStats(19).meta.stimuli = "N:\Stimuli\2019-Manifold\alfa-191210a\full_evolution\backup_12_10_2019_13_07_57";
EStats(27).meta.stimuli = "N:\Stimuli\2019-Manifold\alfa-200409a\2020-04-09-10-16-14";
EStats(28).meta.stimuli = "N:\Stimuli\2019-Manifold\alfa-200409b\2020-04-09-10-54-58";
% assert(exist(strrep(EStats(Expi).meta.stimuli,"2019-06-Evolutions","2019-Manifold"),'dir')>0)
% EStats(Expi).meta.stimuli = strrep(EStats(Expi).meta.stimuli,"2019-06-Evolutions","2019-Manifold");
end
%%
for Expi = 1:length(EStats)
    EStats(Expi).meta.stimuli = strrep(EStats(Expi).meta.stimuli,"\\storage1.ris.wustl.edu\crponce\Active\","N:\");
    if exist(strrep(EStats(Expi).meta.stimuli,"2019-06-Evolutions","2019-Manifold"),'dir') ==0
       disp(Expi)
       disp(EStats(Expi).meta.stimuli)
       keyboard
    end
end
%%
for Expi = 1:length(Stats)
    Stats(Expi).meta.stimuli = strrep(Stats(Expi).meta.stimuli,"\\storage1.ris.wustl.edu\crponce\Active\","N:\");
    if exist(strrep(Stats(Expi).meta.stimuli,"2019-06-Evolutions","2019-Manifold"),'dir') == 0
       disp(Expi)
       disp(Stats(Expi).meta.stimuli)
       keyboard
    end
end
%% Evolution Experiment, far more images (>*30) more noisy! 
%  Harder to find the pattern
ExpType = "Evol";
for Expi = 33:length(EStats) % Starts 14:16
%%
fprintf("Processing Evol Exp %d pref chan %d\n",Expi,EStats(Expi).units.pref_chan)
% Formulate data 
index_vect = cell2mat(EStats(Expi).evol.idx_seq');
imgnm_vect = EStats(Expi).imageName(index_vect);%reshape(imgnm_grid, [], 1);
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; % imgN by 200
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window average 
score_vect = score_vect(:,1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % subsample to decrease redunancy
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]];
% score_vect = mean(psth_all(:,80:200)); 
imgN=length(index_vect); Bsz=128;
%%  online calculation of mean and Sq to compute Correlation
savedir = fullfile(result_dir,compose("%s_Evol_Exp%d",Animal,Expi));
mkdir(savedir);
for layername = ["conv1_2","conv2_2","conv3_1", "conv4_3", "conv5_3"]
% layername = 'conv4_3';
dummy = activations(net, zeros(224,224,3), layername);
ft_shape = size(dummy); %size(feat_tsr,[1,2,3]);
nfeat = prod(ft_shape); 
ntpnt = size(score_vect,2);
SSqFeat = zeros(nfeat,1);
SFeat = zeros(nfeat,1);
SSqrsp = zeros(ntpnt,1);
Srsp = zeros(ntpnt,1);
InnProd = zeros(nfeat, ntpnt);
%% Assume all images use the same suffix
% suffix = ".bmp";
tmpfn = ls(fullfile(EStats(Expi).meta.stimuli, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");
suffix = "."+tmpparts{2};
%%
curN = 0;
csr = 1; 
T0_all = tic;
while curN < imgN
    T0 = tic;
    csr = curN + 1; 
    csr_end = min(imgN, csr + Bsz -1);
    rsp = score_vect(csr:csr_end,:); 
    imgcol = cellfun(@(imgnm) imresize(imread(fullfile(EStats(Expi).meta.stimuli, imgnm+suffix)),[224,224]), ...
        imgnm_vect(csr:csr_end), 'UniformOutput', false);
    T1 = toc(T0);
    dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
    feat_tsr_tmp = activations(net, dlimg, layername);
    T2 = toc(T0);
    feat_tsr_tmp = reshape(feat_tsr_tmp, nfeat, []);
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
cc_tsr = reshape(cc_tsr, [ft_shape, ntpnt]);
MFeat = reshape(MFeat, ft_shape);
StdFeat = reshape(StdFeat, ft_shape);
% Save to dist
% outfn = fullfile(result_dir, compose("%s_Evol_Exp%d_%s.mat",Animal,Expi,layername));
% save(outfn,'cc_tsr','MFeat','StdFeat');
outfn = fullfile(result_dir, compose("%s_Evol_Exp%d_%s.mat",Animal,Expi,layername)); % LW means long window
save(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect');

%% Correlation Coefficient Clearly shows a spatial structure there
figure(17);
for fi=1:size(cc_tsr,4)
wdw = wdw_vect(fi,:);
imagesc(mean(abs(cc_tsr(:,:,:,fi)),3));axis image
title(sprintf("Exp %d Pref chan %d\nMean CorreCoef in of VGG16 %s feature\n with [%d,%d] ms firing rate", ...
    Expi, EStats(Expi).units.pref_chan, layername, wdw(1), wdw(2)))
colorbar()
saveas(17,fullfile(savedir, compose("%s_Evol_Exp%d_%s_wdw%d.png",Animal,Expi,layername,fi)))
end

end

end
%%
