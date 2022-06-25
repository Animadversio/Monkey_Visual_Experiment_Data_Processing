RFStat;% = RFS_col(iRF);
maskS;%

%% Compute the Alpha masks 
iSz = 1; % choose size 1! this is kind of arbitrary. 
FLATMAP_VAL = 0.1;
THR_RATIO = 0.3;%0.500
alphamask = {};
alphamask_thr = {};
alphamask_bin = {};
for ich = 1:numel(maskS)
predmap = maskS(ich).pred_rfmat{iSz};
peak = max(predmap,[],'all');
predRNG = (max(predmap(:)) - min(predmap(:)));
if predRNG > 1E-2 % dynamic range of predicted value is not zero
cval = 0.606 * peak;
alphamask{end+1} = min(predmap / cval,1);
thrval = THR_RATIO * peak;
alphamask_thr{end+1} = max((predmap - thrval) / (peak - thrval), 0);
alphamask_bin{end+1} = double((predmap - thrval)>0);
else % the predicted map constant thus it is not informative! just assume flat prior
alphamask{end+1} = FLATMAP_VAL * ones(size(predmap));
alphamask_thr{end+1} = FLATMAP_VAL * ones(size(predmap));
alphamask_bin{end+1} = FLATMAP_VAL * ones(size(predmap));
end
end
alphamasks = cat(3, alphamask{:}); % cat all the units' masks along the 
alphamasks_thr = cat(3, alphamask_thr{:});
alphamasks_bin = cat(3, alphamask_bin{:});

%%
iCh = find((RFStat.meta.spikeID == S.units.pref_chan) &...
    (RFStat.meta.unitID == S.units.unit_in_pref_chan));
% iCh = findAlignUnits
%% Crop (interpolate) out the image version of the mask.
pixelnum = imgpix;%256

imgpix = S.stim.imsize_pix;
imgsize = S.stim.imsize_deg;
imgpos = S.stim.impos;
imgXlim = [-0.5,0.5]*imgsize + imgpos(1);
imgYlim = [-0.5,0.5]*imgsize + imgpos(2);
imgXq = linspace(imgXlim(1),imgXlim(2),pixelnum);
imgYq = linspace(imgYlim(1),imgYlim(2),pixelnum);
[XXq, YYq] = meshgrid(imgXq, imgYq);
alpha_msk.ceil = alphamasks(:,:,iCh);
alpha_msk.thr = alphamasks_thr(:,:,iCh);
alpha_msk.bin = alphamasks_bin(:,:,iCh);
msktypelist = string(fieldnames(alpha_msk))';
for msktype = msktypelist
    img_alpha_msk.(msktype) = griddata(maskS(1).XX,maskS(1).YY,...
                    alpha_msk.(msktype),XXq,YYq,'cubic');
end
%%

%%
figure;imagesc(maskS(iCh).pred_rfmat{3});axis image
figure;imagesc(img_alpha_msk.thr);axis image
%% load and format images. 
imgtsr = cellfun(@imread, S.stim.imgfps, 'uni',0);
imgtsr = cat(4, imgtsr{:});
%% Resize to the pix size
imgrsztsr = imresize(imgtsr,[imgpix,imgpix]);
%% masking the image tensor
imgrsztsr_msk = uint8(double(imgrsztsr) .* img_alpha_msk.thr);
%%
figure;
montage(imgrsztsr_msk, 'Size', [10 6], 'BorderSize', 2)
%% Loading in the evolution images.
FC6lastgen = BFEStats.imageName(BFEStats.evol.idx_seq{1,end-1});
BGlastgen = BFEStats.imageName(BFEStats.evol.idx_seq{2,end-1});
evol_imgfps = string(cat(1,FC6lastgen,BGlastgen));
evol_imgfps = fullfile(BFEStats.meta.stimuli,evol_imgfps+".bmp");
evolimgtsr = cellfun(@imread, evol_imgfps, 'uni',0);
evolimgtsr = cat(4, evolimgtsr{:});
% masking the image tensor
evolimgrsztsr = imresize(evolimgtsr,[imgpix,imgpix]);
evolimgrsztsr_msk = uint8(double(evolimgrsztsr) .* img_alpha_msk.thr);
%%
figure;
montage(evolimgrsztsr_msk, 'Size', [10 6], 'BorderSize', 2)
%%
D = torchImDist("squeeze");
% D=D.load_net();
%%
distMat = D.distmat_B(imgrsztsr_msk);
%%
figh3 = figure('pos',[680   408   590   570]);
imagesc(distMat);axis image;colorbar()
title("LPIPS-squ distance matrix under RFmask")
saveas(figh3, fullfile(S.meta.figdir,...
    compose("distmat_LPIPS_squ_RFmask_%s.png",S.units.unit_name_arr(iCh))))

%%
distMat_evol = D.distmat_B(cat(4,imgrsztsr_msk,evolimgrsztsr_msk));
%%
figh4 = figure('pos',[680   408   590   570]);
imagesc(distMat_evol);axis image;colorbar()
title("LPIPS-squ distance matrix under RFmask")
saveas(figh4, fullfile(S.meta.figdir,...
    compose("distmat_w_evol_LPIPS_squ_RFmask_%s.png",S.units.unit_name_arr(iCh))))
%%
distvec_mean = mean(distMat_evol(1:60,61:90),2);
distvec_sem = sem(distMat_evol(1:60,61:90),2);
distvec_min = min(distMat_evol(1:60,61:90),[],2);
distvec_mean_BG = mean(distMat_evol(1:60,91:end),2);
distvec_sem_BG = sem(distMat_evol(1:60,91:end),2);
distvec_min_BG = min(distMat_evol(1:60,91:end),[],2);
respvec = S.resp.meanvec_pref;
respvec_sem = S.resp.semvec_pref;
%%
figh5 = figure('pos',[680   531   500   440]);
errorbar(distvec_mean, respvec, ...
    respvec_sem, respvec_sem, distvec_sem, distvec_sem, 'o')
[rho,pval] = corr(distvec_mean, respvec');
title(compose("Distance to FC6 prototypes vs firing rate\ncorr %.3f (P=%.1e)",...
    rho,pval))
ylabel("response firing rate")
xlabel("Distance to prototypes")
saveas(figh5, fullfile(S.meta.figdir,...
    compose("dist-scatter_LPIPS_squ_RFmask_FC6evol_%s.png",S.units.unit_name_arr(iCh))))
%%
figh5 = figure('pos',[680   531   500   440]);
% scatter(meandistvec, respvec)
errorbar(distvec_mean_BG, respvec, ...
    respvec_sem, respvec_sem, distvec_sem_BG, distvec_sem_BG, 'o')
[rho,pval] = corr(distvec_mean_BG, respvec');
title(compose("Distance to BigGAN prototypes vs firing rate\ncorr %.3f (P=%.1e)",...
    rho,pval))
ylabel("response firing rate")
xlabel("Distance to prototypes")
saveas(figh5, fullfile(S.meta.figdir,...
    compose("dist-scatter_LPIPS_squ_RFmask_BGevol_%s.png",S.units.unit_name_arr(iCh))))
%%
figh6 = figure('pos',[680   531   500   440]);
% scatter(meandistvec, respvec)
errorbar(distvec_min_BG, respvec, respvec_sem, respvec_sem, 'o')
[rho,pval] = corr(distvec_min_BG, respvec');
title(compose("Min Distance to BigGAN prototypes vs firing rate\ncorr %.3f (P=%.1e)",...
    rho,pval))
ylabel("response firing rate")
xlabel("Distance to prototypes")
saveas(figh6, fullfile(S.meta.figdir,...
    compose("dist-scatter_LPIPS_squ_min_RFmask_BGevol_%s.png",S.units.unit_name_arr(iCh))))
%%
figh6 = figure('pos',[680   531   500   440]);
% scatter(meandistvec, respvec)
errorbar(distvec_min, respvec, respvec_sem, respvec_sem, 'o')
[rho,pval] = corr(distvec_min, respvec');
title(compose("Min Distance to FC6 prototypes vs firing rate\ncorr %.3f (P=%.1e)",...
    rho,pval))
ylabel("response firing rate")
xlabel("Distance to prototypes")
saveas(figh6, fullfile(S.meta.figdir,...
    compose("dist-scatter_LPIPS_squ_min_RFmask_FC6evol_%s.png",S.units.unit_name_arr(iCh))))
