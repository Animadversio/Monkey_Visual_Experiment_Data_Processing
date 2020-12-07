%% Frame Locator
%%
Animal = "Alfa"; Set_Path;
ftr = find(contains(ExpRecord.ephysFN,"Alfa-27102020"));
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr([1,3,5,6,7,8]),Animal);
%%
D = torchImDist();
%%
Trials_mov = Trials_new{3};
rasters_mov = rasters_new{3};
meta_mov = meta_new{3};
% Window
wdw = meta_mov.rasterWindow;
% Get View Time and the Frame Ticks for marking
viewTime = Trials_mov.TrialRecord.User.viewTime;
vid = VideoReader(fullfile(meta_mov.stimuli, Trials_mov.imageName{1}+".avi"));
vidLenth = vid.Duration * 1000;
fps = vid.FrameRate; 
frameN = int32(viewTime / (1000/fps) + 1);
frameTick = (1000/fps)*double(0:frameN); 
% Get movie names and sort the trials into movies
movnm = string(unique(Trials_mov.imageName));
[movnm_sorted, sortedProp, sortIdx] = sortMovieNames(movnm);
mov_idx_arr = arrayfun(@(mv)find(contains(Trials_mov.imageName, mv)),movnm_sorted,'Uni',0);
%%
psthmov_mean = cellfun(@(idx) mean(rasters_mov(:, :, idx),3), mov_idx_arr, 'Uni', 0);
errmov_mean = cellfun(@(idx) std(rasters_mov(:, :, idx),1,3) / sqrt(numel(idx)), mov_idx_arr, 'Uni', 0);
%%
vid = VideoReader(fullfile(meta_mov.stimuli, movnm_sorted{1}+".avi"));
frames = {};
for i = 1:100
frames{end+1} = vid.readFrame();
end
%%
frames_tsr = cell2mat(reshape(frames,1,1,1,[]));

%%
meta_cmb = meta_new{4};
imageName_cmb = cat(1, Trials_new{4}.imageName, Trials_new{5}.imageName, Trials_new{6}.imageName);
rasters_cmb = cat(3, rasters_new{4:6});
stimparts = split(meta_new{4}.stimuli,'\');
spikeID = meta_new{4}.spikeID;
unit_name_arr = generate_unit_labels_new(spikeID, meta_new{4}.unitID);
unit_num_arr = meta_new{4}.unitID;
% unit_name_arr = generate_unit_labels(spikeID);
% [activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, spikeID, rasters_cmb);
%% Predict 
uniq_imgnm = unique(imageName_cmb);
[idx_arr_nos, imgnm_arr_nos, idx_arr_cls, imgnm_arr_cls, ...
    eig_id_arr_nos, dist_arr_nos, eig_id_arr_cls, dist_arr_cls] = parse_image_idx_arr_hess(imageName_cmb);
%%
idx_arr = [idx_arr_cls; idx_arr_nos];
imgnm_arr = [imgnm_arr_cls; imgnm_arr_nos];
psth_mean = cellfun(@(idx)mean(rasters_cmb(:,:,idx),3), idx_arr,'Uni',0);
psth_sem = cellfun(@(idx)std(rasters_cmb(:,:,idx),1,3)/sqrt(numel(idx)), idx_arr,'Uni',0);
%%
keyimgs = {};
for nm = imgnm_arr_cls(1,5:9)
    keyimgs{end+1} = imread(fullfile(meta_cmb.stimuli, nm+".jpg"));
end
keyimgs_tsr = cell2mat(reshape(keyimgs,1,1,1,[]));
%% Show the frames vs the static images. 
figure(3);T=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile(1)
montage(frames_tsr(:,:,:,1:50))
nexttile(2)
montage(keyimgs_tsr(:,:,:,:))
title(T,'Movie Frame Sequence Static Images Comparison 1')
% saveas(3,fullfile(figdir,"frame_img_cmp_class_eig0.png"))
%%  Movie frame - Static image distance matrix
imfrdistmat = [];
for i = 1:numel(keyimgs)
imfrdistmat(i,:) = D.distance(keyimgs_tsr(:,:,:,i), frames_tsr);
end
% Closest frame id, now assume that's the same for each movie
[mindist, minfrid] = mink(imfrdistmat,2,2); % 2 closest frames in the movie. 
% Translate frame id into onset and offset timing of that frame
matchTON = (minfrid - 1) * 100 / 3;
matchTOFF = (minfrid) * 100 / 3;
%% Frame, Static Image Response Comparison.
%  This is super interesting 
MvrspDelayWdw = [80:200];
ImgrspDelayWdw = [80:200];
iCh = 17; iU = 1;
chid = find(meta_cmb.spikeID==iCh & meta_cmb.unitID==iU);%10; % 
rspmat_static = cellfun(@(psth) mean(psth(chid, ImgrspDelayWdw),[1,2]), psth_mean);
chid = find(meta_mov.spikeID==iCh & meta_mov.unitID==iU); % find(meta_mov.spikeID==9)
wdw = meta_mov.rasterWindow;
matchONidx = int32(matchTON) - wdw(1) + 1;
movie_subRsp = [];
for iMv = 1:numel(psthmov_mean) % movies
    for iTic = 1:size(matchONidx, 1) % static frames in this movie
        rsps = [];
        iTons = double(matchONidx(iTic,:));
        for iTon = iTons % if there is multiple matches we average it or pick one.
           rsps(end+1) = mean(psthmov_mean{iMv}(chid, iTon+rspDelayWdw));
        end
        movie_subRsp(iMv, iTic) = mean(rsps);%(rsp1 + rsp2) / 2;
    end
end
% discard the center column for now as repitition suppression can happen. 
rspmat_key = rspmat_static(:,6:9);
movie_subRsp_key = movie_subRsp(:,2:5);
[cc,pp] = corr(rspmat_key(:),movie_subRsp_key(:)) % examine corelation for this cell pair.

%% See how the response calculation window can be tuned!
MvrspDelayWdw = [60:200];
ImgrspDelayWdw = [60:200];
MovImgCorrStats = repmat(struct(),1,numel(meta_cmb.spikeID)); % collect results here. 
% iCh = 17; iU = 1;
% chid = find(meta_cmb.spikeID==iCh & meta_cmb.unitID==iU);%10; % 
for chid = 1:numel(meta_cmb.spikeID) % loop through static image units
iCh = meta_cmb.spikeID(chid);
iU = meta_cmb.unitID(chid);
rspmat_static = cellfun(@(psth) mean(psth(chid, ImgrspDelayWdw),[1,2]), psth_mean);

chid2 = find(meta_mov.spikeID==iCh & meta_mov.unitID==iU); % find(meta_mov.spikeID==9)
wdw = meta_mov.rasterWindow;
matchONidx = int32(matchTON) - wdw(1) + 1;
movie_subRsp = [];
for iMv = 1:numel(psthmov_mean)
    for iTic = 1:size(matchONidx, 1)
        rsps = [];
        iTons = double(matchONidx(iTic,:));
        for iTon = iTons
           rsps(end+1) = mean(psthmov_mean{iMv}(chid2, iTon+MvrspDelayWdw));
        end
        movie_subRsp(iMv, iTic) = mean(rsps);%(rsp1 + rsp2) / 2;
    end
end
% Correlate the non-center image responses
rspmat_key = rspmat_static(:,6:9);
movie_subRsp_key = movie_subRsp(:,2:5);
[cc,pp] = corr(rspmat_key(:),movie_subRsp_key(:));
MovImgCorrStats(chid).corr = cc;
MovImgCorrStats(chid).corr_P = pp;
MovImgCorrStats(chid).rspmat_movie = movie_subRsp;
MovImgCorrStats(chid).rspmat_static = rspmat_static;
MovImgCorrStats(chid).iCh = iCh;
MovImgCorrStats(chid).iU = iU;
end
corr_P_arr = arrayfun(@(S)S.corr_P,MovImgCorrStats);
corr_arr = arrayfun(@(S)S.corr,MovImgCorrStats);
spikeID = meta_cmb.spikeID;
% print results 
fprintf("Response Delay window movie [%d,%d] image [%d,%d] ms\n",MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end))
fprintf("%d / %d channels has significant(0.01) correlation between Movie and Image Response\n",sum(corr_P_arr<0.01),numel(spikeID))
fprintf("%d / %d channels has significant(0.001) correlation\n",sum(corr_P_arr<0.001),numel(spikeID))
fprintf("%d / %d IT channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & spikeID'<=32),sum(spikeID<=32))
fprintf("%d / %d V1 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & spikeID'>=33 & spikeID'<=48),sum(spikeID'>=33 & spikeID'<=48))
fprintf("%d / %d V4 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & spikeID'>=49),sum(spikeID'>=49))

%%
figdir = "E:\OneDrive - Washington University in St. Louis\MovieDynamics\2020-10-27-Alfa-Chan09-1";
save(fullfile(figdir,"MovieImageCorrStat.mat"),'MovImgCorrStats','MvrspDelayWdw','ImgrspDelayWdw')
%% Summary histogram of correlation in all arrays
Pthresh1 = min(corr_arr(corr_P_arr<0.01));
Pthresh2 = min(corr_arr(corr_P_arr<0.05));
figure(2);clf
histogram(corr_arr,10)
vline([Pthresh1,Pthresh2],'r-.',{"P 0.01","P 0.05"})
xlabel("Correlation of Response");box off
title(compose("Correlation of Movie-Image Firing Rate Rsp %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s_2.png",Animal)))
savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s_2.fig",Animal)))
%% Summary histogram for correlation area separated
figure(2);clf;hold on
histogram(corr_arr(spikeID'>=33 & spikeID'<=48),10,'FaceAlpha',0.4)
histogram(corr_arr(spikeID'>=49),10,'FaceAlpha',0.4)
histogram(corr_arr(spikeID'<=32),15,'FaceAlpha',0.4)
vline([Pthresh1,Pthresh2],'r-.',{"P 0.01","P 0.05"})
xlabel("Correlation of Response")
legend(["V1","V4","IT"]) 
title(compose("Correlation of Movie-Image Firing Rate Rsp %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas_2.png",Animal)))
savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas_2.fig",Animal)))
%% Plot correlation for each unit
figure(1);
for chid = 1:numel(MovImgCorrStats)
    iCh = MovImgCorrStats(chid).iCh;
    iU = MovImgCorrStats(chid).iU;
    static_rspmat_key = MovImgCorrStats(chid).rspmat_static(:,6:9);
    movie_rspmat_key = MovImgCorrStats(chid).rspmat_movie(:,2:5);
    scatter(static_rspmat_key(:), movie_rspmat_key(:))
    title(compose("Correlation of Movie-Image Firing Rate Rsp %s Chan%d Unit%s\nCorr Coef %.3f(%.1e)\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, iCh, char(64+iU), MovImgCorrStats(chid).corr, MovImgCorrStats(chid).corr_P, ...
            MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
    xlabel("Static Firing Rate");ylabel("Movie Firing Rate")
    saveas(1,fullfile(figdir,compose("movie_image_corr_%s_%d%s.png",Animal,iCh,char(64+iU))))
end
%%
figdir = "E:\OneDrive - Washington University in St. Louis\MovieDynamics\2020-10-27-Alfa-Chan09-1";
MovPSTHWdw = [1:200];
ImgPSTHWdw = [1:200];
wdw = meta_mov.rasterWindow;
matchONidx = int32(matchTON) - wdw(1) + 1; % Timing for each frame
% MovImgCorrStats = repmat(struct(),1,numel(meta_cmb.spikeID)); % collect results here. 
figure(4);clf;set(4,'pos',[410    62   800   920])
nrow = numel(psthmov_mean); ncol = size(matchONidx, 1); % matched frame number in each movie
T = tiledlayout(nrow, ncol, 'Padding', 'compact', 'TileSpacing', 'compact');
xlabel(T, "Distance from center")
ylabel(T, "Eigen Axes")
% Set up label for first col and last row
ax_col = {};
for iMv = 1:nrow
    for iTic = 1:ncol
    ax_col{iMv, iTic} = nexttile((iMv-1)*5+iTic); 
    end
end
for c = 1:ncol
    ci = strfind(imgnm_arr{nrow,c},'lin');
    xlabel(ax_col{nrow,c}, imgnm_arr{nrow,c}(ci:end),'Interp','none')
end
for r = 1:nrow
    ci = strfind(imgnm_arr{r,1},'lin');
    ylabel(ax_col{r,1}, imgnm_arr{r,1}(1:ci-2),'Interp','none')
end
% Pooling in data 
for chid = 1:numel(meta_cmb.spikeID) % loop through static image units
iCh = meta_cmb.spikeID(chid);
iU = meta_cmb.unitID(chid);
static_PSTH = cellfun(@(psth) psth(chid, ImgPSTHWdw), psth_mean(:,5:9), 'Uni', 0);
chid2 = find(meta_mov.spikeID==iCh & meta_mov.unitID==iU); % find(meta_mov.spikeID==9)
movie_frPSTH = cell(numel(psthmov_mean), size(matchONidx, 1));
for iMv = 1:numel(psthmov_mean)
    for iTic = 1:size(matchONidx, 1)
        iTons = double(matchONidx(iTic,:));
        for iTon = iTons
           movie_frPSTH{iMv, iTic}(end+1,:) = psthmov_mean{iMv}(chid2, iTon+MovPSTHWdw);
        end
    end
end
rspmat_key = cellfun(@mean, static_PSTH);%rspmat_static(:,2:5);
rspmat_key = rspmat_key(:,2:5);
movie_subRsp_key = cellfun(@(P)mean(P,'all'), movie_frPSTH);% movie_subRsp(:,2:5);
movie_subRsp_key = movie_subRsp_key(:,2:5);
[cc,pp] = corr(rspmat_key(:),movie_subRsp_key(:));

YMAX = max(cellfun(@(PS,PM)max(prctile(PS,99,'all'),prctile(PM,99,'all')),static_PSTH,movie_frPSTH),[],'all');
set(0,'CurrentFigure',4)
title(T, compose("%s Chan %d Unit %s Corr %.3f (%.1e)\n Image PSTH Window Image ON + [%d, %d]ms\n Movie PSTH Window Frame ON + [%d, %d]ms", ...
    Animal, iCh, char(64+iU), cc, pp, ...
    ImgPSTHWdw(1), ImgPSTHWdw(end), MovPSTHWdw(1), MovPSTHWdw(end)))
for iMv = 1:numel(psthmov_mean)
    for iTic = 1:size(matchONidx, 1)
    set(gcf,'CurrentAxes',ax_col{iMv, iTic}); cla; hold on
    plot(static_PSTH{iMv, iTic}(:))
    for iMat = 1:size(movie_frPSTH{iMv, iTic}, 1)
       plot(movie_frPSTH{iMv, iTic}(iMat,:))
    end
    ylim([0,YMAX])
    if iMv ~= numel(psthmov_mean), xticklabels([]); end
    if iTic ~= 1, yticklabels([]); end
    end
end
% MovImgCorrStats(chid).corr = cc;
% MovImgCorrStats(chid).corr_P = pp;
% MovImgCorrStats(chid).rspmat_movie = movie_subRsp;
% MovImgCorrStats(chid).rspmat_static = rspmat_static;
% MovImgCorrStats(chid).iCh = iCh;
% MovImgCorrStats(chid).iU = iU;
saveas(4,fullfile(figdir,compose("movie_image_PSTH_corr_%s_%d%s_2.png",Animal,iCh,char(64+iU))))
% savefig(4,fullfile(figdir,compose("movie_image_PSTH_corr_%s_%d%s.fig",Animal,iCh,char(64+iU))))
% pause
end
%% Data computation PSTH correlation
MovPSTHWdw = [1:200];
ImgPSTHWdw = [1:200];
wdw = meta_mov.rasterWindow;
matchONidx = int32(matchTON) - wdw(1) + 1; % Timing for each frame

figure(5); T = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for chid = 1:numel(meta_cmb.spikeID) % loop through static image units
iCh = meta_cmb.spikeID(chid);
iU = meta_cmb.unitID(chid);
static_PSTH = cellfun(@(psth) psth(chid, ImgPSTHWdw), psth_mean(:,5:9), 'Uni', 0);
chid2 = find(meta_mov.spikeID==iCh & meta_mov.unitID==iU); % Find correponding channel in Movie
movie_frPSTH = cell(numel(psthmov_mean), size(matchONidx, 1));
for iMv = 1:numel(psthmov_mean)
    for iTic = 1:size(matchONidx, 1)
        iTons = double(matchONidx(iTic,:));
        for iTon = iTons
           movie_frPSTH{iMv, iTic}(end+1,:) = psthmov_mean{iMv}(chid2, iTon+MovPSTHWdw);
        end
    end
end
PSTH_corrmat2 = cellfun(@(psthM,psthS)corr(psthM(1,:)', psthS'),movie_frPSTH,static_PSTH,'uni',1);
PSTH_corrmat1 = cellfun(@(psthM,psthS)corr(psthM(2,:)', psthS'),movie_frPSTH,static_PSTH,'uni',1);
PSTH_corrmat_m = cellfun(@(psthM,psthS)corr(mean(psthM',2), psthS'),movie_frPSTH,static_PSTH,'uni',1);
rspmat_key = cellfun(@mean, static_PSTH);%rspmat_static(:,2:5);
rspmat_key = rspmat_key(:,2:5);
movie_subRsp_key = cellfun(@(P)mean(P,'all'), movie_frPSTH);% movie_subRsp(:,2:5);
movie_subRsp_key = movie_subRsp_key(:,2:5);
[cc,pp] = corr(rspmat_key(:),movie_subRsp_key(:)); % Correlation of firing rate

CLIM = prctile([PSTH_corrmat2;PSTH_corrmat1;PSTH_corrmat_m],[2,98],'all')';
nexttile(1);
imagesc(PSTH_corrmat1);axis image;caxis(CLIM);colorbar();xticklabels();yticklabels();
title("First Half Occurance")
nexttile(2);
imagesc(PSTH_corrmat2);axis image;caxis(CLIM);colorbar();xticklabels();
title("2nd Half Occurance")
nexttile(3);
imagesc(PSTH_corrmat_m);axis image;caxis(CLIM);colorbar();xticklabels();
title("Mean Occurance")
title(T,compose("%s Chan %d Unit %s Firing Rate Corr %.3f (%.1e)\n Image PSTH Window Image ON + [%d, %d]ms\n Movie PSTH Window Frame ON + [%d, %d]ms", ...
    Animal, iCh, char(64+iU), cc, pp, ...
    ImgPSTHWdw(1), ImgPSTHWdw(end), MovPSTHWdw(1), MovPSTHWdw(end)))
saveas(5,fullfile(figdir, compose("movie_image_PSTH_corrmat_%s_%d%s.png", Animal, iCh, char(64+iU))))

end

function [sortedMovnm, sortedRows, sortIdx] = sortMovieNames(movnm)
% Sort the movie names in **lexicoGraphical order**, by space and by eigen idx. 
% This assumes the names to have the structure like "class_eig17_shortshort"
eigi_cell = cellfun(@(eigi){str2double(eigi{1})},regexp(movnm,"_eig(\d*)_",'tokens'));
space_cell = cellfun(@(eigi){eigi{1}{1}},regexp(movnm,"(.*)_eig(\d*)_",'tokens'));
[sortedRows, sortIdx] = sortrows([space_cell,eigi_cell]);
sortedMovnm = movnm(sortIdx);
end