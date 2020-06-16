%% correlation tensor for GAN layer instead of cnn layers.
global G 
G = FC6Generator("matlabGANfc6.mat");
%% Animal specific data.
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
% MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
load(fullfile(MatStats_path, compose("%s_Manif_RFstats.mat", Animal)), 'RFStats')
%%
ExpType = "Manif";
Expi = 11;
si=1;ui=1;
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
imgN = size(imgnm_vect,1);
% Below add more time slices into it
psth_all = cellfun(@(psth) reshape(mean(psth(ui,:,:),[3]),1,1,[]), Stats(Expi).manif.psth{si}, ...
    'UniformOutput', false);
psth_all = reshape(cell2mat(psth_all),imgN,[]);
% The score(firing rate at different time slices) the time window info in
% recorded in `wdw_vect`
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window average 
score_vect = score_vect(:,1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % subsample to decrease redunancy
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]];
%%
basis_path = fullfile(EStats(Expi).meta.stimuli,"PC_imgs","PC_vector_data.npz");
f = py.numpy.load(basis_path);
PC_Vec = f.get('PC_vecs').double;
sphere_norm = f.get('sphere_norm').double;
f.close();
% Get the mat containing all the codes of the last generation. To see
% whether we should inverse PC1 
matfns = string(ls(fullfile(EStats(Expi).meta.stimuli,"*.mat")));
code_tmp = load(fullfile(EStats(Expi).meta.stimuli,matfns(end)));
proj_coord = mean(code_tmp.codes,1) * PC_Vec';
if proj_coord(1)>0
    basis = PC_Vec(1:3,:);
else
    fprintf("The evolution direction is inverse to the PC1 direction of PCA. Inverse PC1 as basis\n")
    basis = [-1,1,1]' .* PC_Vec(1:3,:);% Note the final PC may need to reverse! not always the same dir!
end
clear code_tmp matfns
%% Prepare the vectors for
[theta, phi] = meshgrid(-90:18:90,-90:18:90);
anglemap = cat(3, cosd(theta).* cosd(phi), sind(theta).* cosd(phi), sind(phi));
codes_all = sphere_norm * einsum(anglemap, basis, 'ijk,kl->ijl');
%%
figure(3);
imshow(G.visualize(squeeze(codes_all(6,6,:))))
%%
layername = 'relu_conv3_1';%'relu_deconv2';%'conv3_1';
acts_tsr = G.activations(reshape(codes_all,[],4096),layername);
ft_shape = size(acts_tsr,[1,2,3]);
nfeat = prod(ft_shape);
ntpnt = size(score_vect,2);
Gcc_tsr = corr(reshape(acts_tsr, [], imgN)', score_vect);
Gcc_tsr = reshape(Gcc_tsr, [ft_shape, ntpnt]);
%%

%%
shuffleN = 100;
score_shuffle = [];
for i = 1:shuffleN
    score_shuffle = [score_shuffle, score_vect(randperm(imgN),:)];
end
Gcc_tsr_ref = corr(reshape(acts_tsr,nfeat,imgN)', score_shuffle);
Gcc_tsr_ref = single(reshape(Gcc_tsr_ref, [ft_shape, ntpnt, shuffleN]));
Gcc_refM = mean(Gcc_tsr_ref,5);
Gcc_refS = std(Gcc_tsr_ref,0,5);
t_signif_tsr = (Gcc_tsr - Gcc_refM) ./ Gcc_refS;
%%
figure(1);
for fi=1:size(Gcc_tsr,4)
wdw = wdw_vect(fi,:);
imagesc(mean(abs(Gcc_tsr(:,:,:,fi)),3)); axis image
title(sprintf("Exp %d Pref chan %d\nMean CorreCoef in of fc6 GAN %s feature\n with [%d,%d] ms firing rate", ...
    Expi, EStats(Expi).units.pref_chan, layername, wdw(1), wdw(2)))
colorbar()
pause
% saveas(17,fullfile(savedir, compose("%s_Evol_Exp%d_%s_wdw%d.png",Animal,Expi,layername,fi)))
end
%%
figure(4);fi=11;
for ci=1:size(Gcc_tsr,3)
wdw = wdw_vect(fi,:);
imagesc(Gcc_tsr(:,:,ci,fi)); axis image
title(sprintf("Exp %d Pref chan %d\nMean CorreCoef in of fc6 GAN %s feature %d\n with [%d,%d] ms firing rate", ...
    Expi, EStats(Expi).units.pref_chan, layername, ci, wdw(1), wdw(2)))
colorbar()
pause(0.3)
% saveas(17,fullfile(savedir, compose("%s_Evol_Exp%d_%s_wdw%d.png",Animal,Expi,layername,fi)))
end
%% Plot the mean t value for thresholded tensor
threshold = 5;
t_thresh = t_signif_tsr;
t_thresh(abs(t_signif_tsr) < threshold) = 0;
%%
savedir = "E:\OneDrive - Washington University in St. Louis\corrGANTsr";
v = VideoWriter(fullfile(savedir,compose("%s_%s_Exp%d_fc6GAN%s_T_thresh%.1f.avi",Animal,ExpType,Expi,layername,threshold)));
v.FrameRate = 2;open(v);
plot_tsr = mean(abs(t_thresh),3);
CLIM = prctile(plot_tsr, [2,98], 'all'); 
figure(1);
for fi=1:size(Gcc_tsr,4)
wdw = wdw_vect(fi,:);
imagesc(plot_tsr(:,:,:,fi)); caxis(CLIM);axis image
title(sprintf("Exp %d Pref chan %d\nMean CorreCoef in of fc6 GAN %s feature\n with [%d,%d] ms firing rate", ...
    Expi, EStats(Expi).units.pref_chan, layername, wdw(1), wdw(2)))
colorbar()
Fs = getframe(1);
writeVideo(v,Fs);
pause(0.1)
% saveas(17,fullfile(savedir, compose("%s_Evol_Exp%d_%s_wdw%d.png",Animal,Expi,layername,fi)))
end
close(v);
%%
threshold = 1;
plot_tsr = Gcc_tsr;
plot_tsr(abs(t_signif_tsr) < threshold) = 0;
plot_tsr = mean(abs(plot_tsr),3);
%% Plot the mean correlation value for the thresholded tensor.
v = VideoWriter(fullfile(savedir,compose("%s_%s_Exp%d_fc6GAN%s_cc_thresh%.1f.avi",Animal,ExpType,Expi,layername,threshold)));
v.FrameRate = 2;open(v);
CLIM = prctile(plot_tsr, [2,98], 'all'); 
figure(3);
for fi=1:size(Gcc_tsr,4)
wdw = wdw_vect(fi,:);
imagesc(plot_tsr(:,:,1,fi)); caxis(CLIM);axis image
title(sprintf("Exp %d Pref chan %d\nMean CorreCoef in of fc6 GAN %s feature\n with [%d,%d] ms firing rate", ...
    Expi, EStats(Expi).units.pref_chan, layername, wdw(1), wdw(2)))
colorbar()
Fs = getframe(3);
writeVideo(v,Fs);
pause(0.1)
% saveas(17,fullfile(savedir, compose("%s_Evol_Exp%d_%s_wdw%d.png",Animal,Expi,layername,fi)))
end
close(v);