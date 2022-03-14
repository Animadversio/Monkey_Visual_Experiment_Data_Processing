%% Find most aligned natural images with a PC Cosine Exp 
figh5 = 15;figmtg3 = 16;
for iCS = 1:numel(CosStats)
CStat = CosStats(iCS);
RFdir = fullfile(dataroot,compose("%s-%s-RF", ...
    datestr(CStat.meta.expday,"yyyy-mm-dd"), CStat.Animal));
seldir = fullfile(dataroot,compose("%s-%s-refRepr", ...
    datestr(CStat.meta.expday,"yyyy-mm-dd"), CStat.Animal));
load(fullfile(RFdir,'RFStat.mat'),'RFStat');
load(fullfile(RFdir,'maskStat.mat'),'maskS');
load(fullfile(seldir,'ReprStat.mat'),'ReprStat');
load(fullfile(CStat.meta.figdir,"Evol_RF_alpha_mask.mat"),...
            'img_alpha_msk','alpha_msk','alpha_msk_flip',...
            'alphamasks','alphamasks_thr','alphamasks_bin',...
            'RFchlist','RFchname')
assert(RFStat.meta.datetime == CStat.meta.expday)
%%
respmat = []; % online matrix of 
for ich = 1:numel(ReprStat.units.spikeID)
if ReprStat.units.unit_num_arr(ich) == 0, continue;end
chan = ReprStat.units.spikeID(ich);
unit = ReprStat.units.unit_num_arr(ich);
respmat = [respmat; CStat.targ.repr_D.responseTensor(chan,unit+1,:)];
end
respmat = squeeze(respmat);
size(ReprStat.resp.meanMat)
size(respmat)
offline_respmat = ReprStat.resp.meanMat(ReprStat.units.unit_num_arr~=0,:);
offline_evkmat = offline_respmat - ReprStat.resp.bslmean(ReprStat.units.unit_num_arr~=0,:);
on_offline_corr = corr(respmat(:), offline_evkmat(:), 'rows','complete');
xcorrmat = (corr(offline_evkmat', respmat'));
selcorrmat = (corr(offline_evkmat', offline_evkmat'));
xcorr_chan = diag(xcorrmat);
fprintf("Correlation of online-offline resp %.3f\n",on_offline_corr)
fprintf("single chan correlation %.3f+-%.3f\n",nanmean(xcorr_chan),nanstd(xcorr_chan))
%%
Fmsk = CStat.targ.maskVec{1};
[targareamsk, targarea] = parse_mode2maskVec(CStat.targ.score_mode, ...
    CStat.meta.array_layout, CStat.meta.spikeID);
targlist = find(targareamsk & Fmsk);
targAct = CStat.targ.targetActVec{1}(targlist);
%%
Reprchlist = findAlignUnits(targlist, CStat, ReprStat);
assert(numel(targAct)==numel(Reprchlist))
Reprmsk = zeros(size(ReprStat.units.spikeID),'logical');
Reprmsk(Reprchlist) = 1;
targAct_sel = nan(size(ReprStat.units.spikeID));
targAct_sel(Reprchlist) = targAct; 
% targAct_sel has same length as units in selectivity exp. 
% use to compute similarity with the representation in selectivity experiments.
%% Find the natural image with most similar representation
% mean and std computed from offline 
meanVec= mean(ReprStat.resp.meanMat - ReprStat.resp.bslmean,2); 
stdVec = std(ReprStat.resp.meanMat,1,2);
normReprMat = (ReprStat.resp.meanMat - ReprStat.resp.bslmean - meanVec)./stdVec;

msktype = "rect_wt_thr";
scores_cosine = scorePopulationVec_direction(ReprStat.resp.meanMat - ReprStat.resp.bslmean,...
    targAct_sel, meanVec, stdVec, Reprmsk,'cosine');

[maxscores, maxidxs] = maxk(scores_cosine,12);
visualize_align_natimage(CStat, ReprStat.stim.imgfps(maxidxs), maxscores,...
                         'cosine',alpha_msk_flip,msktype,figmtg3,figmtg3+1)
[minscores, minidxs] = mink(scores_cosine,12);
visualize_align_natimage(CStat, ReprStat.stim.imgfps(minidxs), minscores,...
                         'cosine_neg',alpha_msk_flip,msktype,figmtg3,figmtg3+1)


figure(figh5);clf;hold on;set(figh5,'pos',[680   557   860   420])
plot(targAct_sel(Reprmsk),'k','linewidth',2.5)
plot(normReprMat(Reprmsk,maxidxs),'linewidth',1)
xticks(1:sum(Reprmsk));xlim([.5,sum(Reprmsk)+.5])
xticklabels(ReprStat.units.unit_name_arr(Reprmsk))
legend(["Target Pattern";string(ReprStat.stim.imgname_uniq(maxidxs))],"location",'eastoutside')
ylabel("zscore (evk - bsl activation)")
title(CStat.meta.explabel,'interpreter','none','FontSize',12)
saveallform(CStat.meta.figdir,compose("best_align_natimg_%s_pop_pattern","cosine"),figh5,["jpg","fig","pdf"])

scores_dot = scorePopulationVec_direction(ReprStat.resp.meanMat - ReprStat.resp.bslmean,...
    targAct_sel, meanVec, stdVec, Reprmsk,'dot');

[maxscores, maxidxs] = maxk(scores_dot,12);
visualize_align_natimage(CStat, ReprStat.stim.imgfps(maxidxs), maxscores,...
                         'dot',alpha_msk_flip,msktype,figmtg3,figmtg3+1)
[minscores, minidxs] = mink(scores_dot,12);
visualize_align_natimage(CStat, ReprStat.stim.imgfps(minidxs), minscores,...
                         'dot_neg',alpha_msk_flip,msktype,figmtg3,figmtg3+1)
% figure(figmtg3);
% montage(ReprStat.stim.imgfps(maxidxs))
% disp(maxscores)
% %%
% figure(figmtg3);
% montage(ReprStat.stim.imgfps(maxidxs))
% disp(maxscores)

%%
%%
% (targAct_sel' * normReprMat)./vecnorm(normReprMat,2,1)./vecnorm(targAct_sel,2,1)
% %%
% scores_tmp = sum(targAct_sel(Reprmsk) .* normReprMat(Reprmsk,:),1)./vecnorm(normReprMat(Reprmsk,:),2,1)./vecnorm(targAct_sel(Reprmsk),2,1);
% scores_tmp_dot = sum(targAct_sel(Reprmsk) .* normReprMat(Reprmsk,:),1);
end
%%
figure;plot(ReprStat.resp.meanMat./std(ReprStat.resp.meanMat,1,2))
%%
figure;
plot((ReprStat.resp.meanMat - ReprStat.resp.bslmean)./std(ReprStat.resp.meanMat,1,2))
%%
figure;
plot((ReprStat.resp.meanMat - mean(ReprStat.resp.meanMat,2))./std(ReprStat.resp.meanMat,1,2))
%%

function visualize_align_natimage(CStat, imgfps, imscores, score_mode, alpha_msk_flip,msktype,figh,figh2)
    if nargin==1,  
        figh = figure('pos',[511    72   900   760]);
        figh2 = figure('pos',[511    72   970   760]);
    else, 
        figh = figure(figh);set(figh,'pos',[511    72   900   760]);
        figh2 = figure(figh2);set(figh2,'pos',[511    72   970   760]);
    end 
    imgalphamsk = alpha_msk_flip.(msktype);
    best_img_col = {};
    mask_img_col = {};
    % best_imnam_col = strings();
    for imgi = 1:numel(imgfps)
        img = imread(imgfps{imgi});
        best_img_col{imgi} = imresize(img, [256,256]); % note natural images are not necessaruly 256
        mask_img_col{imgi} = single(best_img_col{imgi}) / 255.0 .* imgalphamsk;
    end
    imscores = reshape(imscores,size(best_img_col));
    % TOBE refactored ! 
    figh = figure(figh);%axs = subplot(111);
    montage(best_img_col,'Parent',gca)
    % set(figh,'pos',[511    72   900   760])
    title(compose("%s\nScore mode %s",CStat.meta.explabel,score_mode),'interpreter','none','FontSize',12)
    saveallform(CStat.meta.figdir,compose("best_align_natimg_%s",score_mode),figh,["jpg","pdf"])

    figh = figure(figh);%axs = subplot(111);
    montage(mask_img_col,'Parent',gca)
    % set(figh,'pos',[511    72   900   760])
    title(compose("%s\nScore mode %s Mask Type %s",CStat.meta.explabel,score_mode,msktype),'interpreter','none','FontSize',12)
    saveallform(CStat.meta.figdir,compose("best_align_natimg_%s_wRFmsk_",score_mode)+msktype,figh,["jpg","pdf"])
    
    figh2 = figure(figh2);%axs = subplot(111);
    [scoreframe_imgs, Clim] = score_frame_image_arr(best_img_col, imscores);
    montage(scoreframe_imgs,'Parent',gca)
    caxis([Clim]); cb = colorbar(); cb.Label.Interpreter = 'none';
    cb.Label.String = score_mode; cb.Label.FontSize = 12;
    % set(figh2,'pos',[511    72   970   760])
    title(compose("%s\nScore mode %s",CStat.meta.explabel,score_mode),'interpreter','none','FontSize',12)
    saveallform(CStat.meta.figdir,compose("best_align_natimg_%s_score_framed",score_mode),figh2,["jpg","pdf"])

    figh2 = figure(figh2);%axs = subplot(111);
    [scoreframe_imgs, Clim] = score_frame_image_arr(mask_img_col, imscores);
    montage(scoreframe_imgs,'Parent',gca)
    caxis([Clim]); cb = colorbar(); cb.Label.Interpreter = 'none';
    cb.Label.String = score_mode; cb.Label.FontSize = 12;
    % set(figh2,'pos',[511    72   970   760])
    title(compose("%s\nScore mode %s Mask Type %s",CStat.meta.explabel,score_mode,msktype),'interpreter','none','FontSize',12)
    saveallform(CStat.meta.figdir,compose("best_align_natimg_%s_score_framed_wRFmsk_",score_mode)+msktype,figh2,["jpg","pdf"])

end
