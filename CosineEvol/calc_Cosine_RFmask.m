

dataroot = "O:\Evol_Cosine";
area_colormap = containers.Map({'V1','V4','IT'},{[1,0,0],[0,1,0],[0,0,1]});

figh = 11; figh2 = 12; figh3=13; figmtg = 14;

for iCS = 1:numel(CStats)
CStat = CStats(iCS);
RFdir = fullfile(dataroot,compose("%s-%s-RF", ...
    datestr(CStat.meta.expday,"yyyy-mm-dd"), CStat.Animal));
load(fullfile(RFdir,'RFStat.mat'),'RFStat');
load(fullfile(RFdir,'maskStat.mat'),'maskS');
assert(RFStat.meta.datetime == CStat.meta.expday)

%% get configs of Cosine Evol 
Fmsk = CStat.targ.maskVec{1};
imgsize = CStat.evol.imgsize{1};
imgpos = CStat.evol.imgpos{1};
imgXlim = [-0.5,0.5]*imgsize + imgpos(1);
imgYlim = [-0.5,0.5]*imgsize + imgpos(2);

[targareamsk, targarea] = parse_mode2maskVec(CStat.targ.score_mode, ...
    CStat.meta.array_layout, CStat.meta.spikeID);
targlist = find(targareamsk & Fmsk);
targAct = CStat.targ.targetActVec{1}(targlist);
% Find the units that are involved in the objective. Idx according to the
% RF experiments
RFchlist = findAlignUnits(targlist, CStat, RFStat);
assert(numel(targAct)==numel(RFchlist))
% Find units that has significant RFs experiments
% Note that, some units don't have a good RF estimate but it's still
% selective. For those we need to assign uniform prior over the image. 
TFvalidlist = find((RFStat.stats.bestT_P<1E-3) & (RFStat.stats.F_P<1E-5));
%% Plot overall contours 
% 
% figure(figh);clf;set(figh,'pos',[684   327   480   490]);hold on
% for ich = RFchlist'
%     chan_num = RFStat.unit.chan_num_arr(ich);
%     area = area_map(chan_num,CStat.meta.array_layout); 
%     clr = area_colormap(area);
%     cval = 0.606 * max(maskS(ich).pred_rfmat{1},[],'all');
%     contour(maskS(1).XX,maskS(1).YY,maskS(ich).pred_rfmat{1}, [cval,cval],'linecolor',[clr])
% end
% rectangle('pos',[imgXlim(1),imgYlim(1),imgsize,imgsize],'linesty','-.')
% title(compose("RF of Optim Objective Units\n%s",CStat.meta.explabel),'interpr','none')
% axis equal
% saveallform(CStat.meta.figdir,"RFcontour_target_fullview",figh)
% %% Plot RF contours for all channel with valid RF stats
% figure(figh2);clf;set(figh2,'pos',[684   327   480   490]);hold on
% for ich = TFvalidlist'
%     chan_num = RFStat.unit.chan_num_arr(ich);
%     area = area_map(chan_num,CStat.meta.array_layout); 
%     clr = area_colormap(area);
%     cval = 0.606 * max(maskS(ich).pred_rfmat{1},[],'all');
%     contour(maskS(1).XX,maskS(1).YY,maskS(ich).pred_rfmat{1}, [cval,cval],'linecolor',[clr])
% end
% rectangle('pos',[imgXlim(1),imgYlim(1),imgsize,imgsize],'linesty','-.')
% title(compose("RF of All Valid Units\n%s",CStat.meta.explabel),'interpr','none')
% axis equal
% saveallform(CStat.meta.figdir,"RFcontour_all_fullview",figh2)

%% Collect all Alpha masks 
iSz = 1;
FLATMAP_VAL = 0.1;
alphamask = {};
alphamask_thr = {};
alphamask_bin = {};
for ich = RFchlist'
predmap = maskS(ich).pred_rfmat{iSz};
peak = max(predmap,[],'all');
predRNG = (max(predmap(:)) - min(predmap(:)));
if predRNG > 1E-2
cval = 0.606 * peak;
alphamask{end+1} = min(predmap / cval,1);
thrval = 0.500 * peak;
alphamask_thr{end+1} = max((predmap - thrval) / (peak - thrval), 0);
alphamask_bin{end+1} = double((predmap - thrval)>0);
else % the predicted map is not informative! just assume flat prior
alphamask{end+1} = FLATMAP_VAL * ones(size(predmap));
alphamask_thr{end+1} = FLATMAP_VAL * ones(size(predmap));
alphamask_bin{end+1} = FLATMAP_VAL * ones(size(predmap));
end
end
alphamasks = cat(3, alphamask{:});
alphamasks_thr = cat(3, alphamask_thr{:});
alphamasks_bin = cat(3, alphamask_bin{:});
%% Average the RF masks collected.
% Weighted average of the alpha masks save to the struct
alpha_msk = struct();
alpha_msk.unif_wt = mean(alphamasks,3);
alpha_msk.unif_wt_thr = mean(alphamasks_thr,3);
alpha_msk.unif_wt_bin = mean(alphamasks_bin,3);
% Weighted average using the amplitude of target activation pattern.
% (negative target matters)
weights = abs(targAct);
alpha_msk.abs_wt = sum(alphamasks .* reshape(weights,1,1,[]),3) / sum(weights);
alpha_msk.abs_wt_thr = sum(alphamasks_thr .* reshape(weights,1,1,[]),3) / sum(weights);
alpha_msk.abs_wt_bin = sum(alphamasks_bin .* reshape(weights,1,1,[]),3) / sum(weights);
% Weighted average using the rectified target activation pattern. 
% (negative target doesn't count)
weights = max(targAct,0);
alpha_msk.rect_wt = sum(alphamasks .* reshape(weights,1,1,[]),3) / sum(weights);
alpha_msk.rect_wt_thr = sum(alphamasks_thr .* reshape(weights,1,1,[]),3) / sum(weights);
alpha_msk.rect_wt_bin = sum(alphamasks_bin .* reshape(weights,1,1,[]),3) / sum(weights);

msktypelist = string(fieldnames(alpha_msk))';

imgXq = linspace(imgXlim(1),imgXlim(2),256);
imgYq = linspace(imgYlim(1),imgYlim(2),256);
[XXq, YYq] = meshgrid(imgXq, imgYq);
for msktype = msktypelist
    img_alpha_msk.(msktype) = griddata(maskS(1).XX,maskS(1).YY,...
                    alpha_msk.(msktype),XXq,YYq,'cubic');
end

%%
for msktype = msktypelist
    % plot the RF mask with the full visual field [-8, 8] deg
    figure(figh3);set(figh3,'pos',[684   327   480   490])
    imagesc(maskS(1).Xq, maskS(1).Yq, alpha_msk.(msktype))
    set(gca,'Ydir','normal');axis image; 
    colormap('gray'); colorbar();
    title(compose("RF of Optim Objective Units MaskType %s\n%s",msktype,CStat.meta.explabel),'interpr','none')
    saveallform(CStat.meta.figdir,"RF_mean_mask_"+msktype+"_fullview",figh3,["png"])
    % plot the RF mask with the image frame. 
    figure(figh3);
    imagesc(imgXq,imgYq,img_alpha_msk.(msktype))
    set(gca,'Ydir','normal');axis image; 
    colormap('gray'); colorbar();
    title(compose("RF of Optim Objective Units MaskType %s\n%s",msktype,CStat.meta.explabel),'interpr','none')
    saveallform(CStat.meta.figdir,"RF_mean_mask_"+msktype+"_imgview",figh3,["png"])
end

%% Produce Alpha mask of images. Use 0.606 as vmax, and set everything above 0.606 as 1.

% Note up down of the array need to be flipped. 
% the mask is saved by smaller Y on top. This is not a problem for imagesc,
% but it is for saving images.
alpha_msk_flip = struct();
for msktype = msktypelist
popmean_peak = max(img_alpha_msk.(msktype),[],'all');
meanmask_clip = min(1, img_alpha_msk.(msktype) / (0.80 * popmean_peak));
meanmask_clip_flip = flipud(meanmask_clip);
imwrite(zeros(256,256),fullfile(CStat.meta.figdir,"alpha_mask_"+msktype+".png"),...
        'Alpha', (1 - meanmask_clip_flip))
alpha_msk_flip.(msktype) = meanmask_clip_flip;
end

RFchname = RFStat.unit.unit_name_arr(RFchlist);
save(fullfile(CStat.meta.figdir,"Evol_RF_alpha_mask.mat"),...
            'img_alpha_msk','alpha_msk','alpha_msk_flip',...
            'alphamasks','alphamasks_thr','alphamasks_bin',...
            'RFchlist','RFchname')

%% Use the mask to produce masked images 
msktype = "rect_wt_thr";
visualize_PCCosine_imageEvol_wMask(CStat, alpha_msk_flip, msktype, figmtg)
visualize_PCCosine_imageEvol(CStat,figmtg+1,figmtg+1)
end

%%
[XXq, YYq] = meshgrid(imgXq, imgYq);
imgalphamsk_wt = griddata(maskS(1).XX,maskS(1).YY,alpha_msk_wt_mean,XXq,YYq);

figure()
imagesc(XXq,YYq,imgalphamsk_wt)
set(gca,'Ydir','normal');axis image; colormap('gray')
%%
%% Debugging part
% figure(9);
% for i = 10%1:numel(alphamask)
%     ich = RFchlist(i);
%     predmap = maskS(ich).pred_rfmat{iSz};
%     subplot(1,4,1)
%     imagesc(predmap)
%     axis equal tight
%     subplot(1,4,2)
%     imagesc(alphamask{i})
%     axis equal tight
%     subplot(1,4,3)
%     imagesc(alphamask_thr{i})
%     axis equal tight
%     subplot(1,4,4)
%     imagesc(alphamask_bin{i})
%     axis equal tight
%     sgtitle([i, RFStat.unit.unit_name_arr(ich)])
%     pause
% end

%%

function visualize_PCCosine_imageEvol_wMask(CStat,alpha_msk_flip, msktype, figmtg, figmtg2)
imgalphamsk = alpha_msk_flip.(msktype);

best_img_col = {};
mask_img_col = {};
best_imnam_col = strings();
gen_idx_seq = CStat.stim.gen_idx_seq;
for geni = 1:numel(gen_idx_seq)
    [maxscore,maxidx] = max(CStat.score.offline_vec(gen_idx_seq{geni}));
    imgnm = CStat.imageName(gen_idx_seq{geni}(maxidx));
    best_imnam_col(geni) = fullfile(CStat.meta.stimuli, imgnm+".bmp");
    best_img_col{geni} = imread(best_imnam_col(geni));
    mask_img_col{geni} = single(best_img_col{geni}) / 255.0 .* imgalphamsk;
end
scores_rec = cellfun(@(idx)CStat.score.offline_vec(idx),gen_idx_seq,'uni',0);
meanscores = cellfun(@mean, scores_rec);
% Use end-1 to get rid of final generation. the score of last generation is
% more variable and could distort the colormap with its outlier values. 
[scoreframe_imgs, Clim] = score_frame_image_arr(mask_img_col(1:end-1),...
                                        meanscores(1:end-1));

figure(figmtg);
montage(mask_img_col(1:end-1))
set(figmtg,'pos',[511    72   970   900])
title(compose("%s  Mask type %s", CStat.meta.explabel, msktype),'interpreter','none','FontSize',12)
saveallform(CStat.meta.figdir,"Image_Evol_per_gen_wRFmsk_"+msktype,figmtg,["jpg","pdf"])

figmtg = figure(figmtg);
montage(scoreframe_imgs)
caxis([Clim]); cb = colorbar(); cb.Label.Interpreter = 'none';
cb.Label.String = CStat.targ.score_mode{1}; cb.Label.FontSize = 12;
set(figmtg,'pos',[511    72   1030   900])
title(compose("%s  Mask type %s", CStat.meta.explabel, msktype),'interpreter','none','FontSize',12)
saveallform(CStat.meta.figdir,"Image_Evol_per_gen_score_framed_wRFmsk_"+msktype,figmtg,["jpg","pdf"])
end