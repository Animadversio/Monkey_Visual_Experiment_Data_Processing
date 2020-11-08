
Trials_mov = Trials_new{3};
rasters_mov = rasters_new{3};
meta_mov = meta_new{3};

Trials_img = Trials_new{4};
rasters_img = rasters_new{4};
meta_img = meta_new{4};
% figdir will be an input argument.
figroot = "E:\OneDrive - Washington University in St. Louis\MovieDynamics";
figdir = fullfile(figroot, "2020-11-04-Alfa-Chan28-1");
mkdir(figdir)

%% Meta information of movies & spike Id
spikeID_mv = meta_mov.spikeID;
unitID_mv = meta_mov.unitID;
unit_str_mv = generate_unit_labels_new(meta_mov.spikeID, meta_mov.unitID);
wdw = meta_mov.rasterWindow;
% Get View Time and the Frame Ticks for marking
viewTime = Trials_mov.TrialRecord.User.viewTime;
vid = VideoReader(fullfile(meta_mov.stimuli, Trials_mov.imageName{1}+".avi"));
vidLenth = vid.Duration * 1000;
fps = vid.FrameRate; 
frameN = int32(viewTime / (1000/fps) + 1);
frameTick = (1000/fps)*double(0:frameN); 
% Get movie names and Sort the trials into movies
movnm = string(unique(Trials_mov.imageName));
[movnm_sorted, sortedProp, sortIdx] = sortMovieNames(movnm); % Sort again because the naming convention of Hessian Movies
mov_idx_arr = arrayfun(@(mv)find(contains(Trials_mov.imageName, mv)),movnm_sorted,'Uni',0);

psthmov_mean = cellfun(@(idx) mean(rasters_mov(:, :, idx),3), mov_idx_arr, 'Uni',0);
psthmov_sem = cellfun(@(idx) std(rasters_mov(:, :, idx),1,3) / sqrt(numel(idx)), mov_idx_arr, 'Uni',0);

%% Meta information of spike Id
spikeID_im = meta_img.spikeID;
unitID_im = meta_img.unitID;
unit_str_im = generate_unit_labels_new(meta_img.spikeID, meta_img.unitID);
% Sort Static images and trials
[idx_arr_nos, imgnm_arr_nos, idx_arr_cls, imgnm_arr_cls, ...
    eig_id_arr_nos, dist_arr_nos, eig_id_arr_cls, dist_arr_cls] = parse_image_idx_arr_hess(Trials_img.imageName);
idx_arr = [idx_arr_cls; idx_arr_nos];
imgnm_arr = [imgnm_arr_cls; imgnm_arr_nos];
psthimg_mean = cellfun(@(idx)mean(rasters_img(:,:,idx),3), idx_arr,'Uni',0);
psthimg_sem = cellfun(@(idx)std(rasters_img(:,:,idx),1,3)/sqrt(numel(idx)), idx_arr,'Uni',0);
centcol = find(dist_arr_nos==0); % cent col index 
imgnm_per_mov = mat2cell(imgnm_arr(:, centcol:end),ones(12,1),5); % Cell array of the images in 
nImgPerMv = max(cellfun(@numel, imgnm_per_mov)); % may be better ways to do so. 
%% Matching process
% Matchin one movie, image pairs
% [matchfrid, matchTON, matchTOFF] = closest_Kframe_locate(movnm_sorted(1), {imgnm_arr(1, centcol:end)}, meta_mov.stimuli, meta_img.stimuli);
% Matching all movie image pairs
[matchfrids, matchTONs, matchTOFFs] = closest_Kframe_locate(movnm_sorted(:), imgnm_per_mov, ...
                                        meta_mov.stimuli, meta_img.stimuli); % halfSep = false,
%% Save the meta information for functions to use.
Stats.Animal = Animal; 
Stats.movnm = movnm_sorted;
Stats.imgnm_arr = imgnm_arr;
Stats.img_idx_arr = idx_arr;
Stats.imgnm_per_mov = imgnm_per_mov;
Stats.nImgPerMv = nImgPerMv;
Stats.centcol = centcol;

Stats.matchfrids = matchfrids;
Stats.matchTONs = matchTONs;
Stats.matchTOFFs = matchTOFFs;

Stats.unit_str_im = unit_str_im;
Stats.unit_str_mv = unit_str_mv;
Stats.spikeID_im = spikeID_im; 
Stats.unitID_im = unitID_im; 
Stats.spikeID_mv = spikeID_mv; 
Stats.unitID_mv = unitID_mv; 
Stats.MvRstrWdw = wdw; % Movie Raster window. Need this to interpret the raster's timeline.
%%
ImgrspDelayWdw = [81:200]; MvrspDelayWdw = [81:200];
[MovImgCorrStats, corr_arr, corr_P_arr, corr_sep_arr, corr_sep_P_arr] = HessMovImgMatchCorr(ImgrspDelayWdw, MvrspDelayWdw, psthimg_mean, psthmov_mean, Stats);
% [MovImgCorrStats, corr_arr, corr_P_arr] = HessMovImgMatchCorr([81:200], [-59:60], psthimg_mean, psthmov_mean, Stats);
% [MovImgCorrStats, corr_arr, corr_P_arr] = HessMovImgMatchCorr([81:200], [120:240], psthimg_mean, psthmov_mean, Stats);
%%
plot_channel_corr(Stats, MovImgCorrStats, ImgrspDelayWdw, MvrspDelayWdw);
%% Compute Tuning Statistics for Static images for each channel for filtering purpose.
Stats.imgTuneStats = {};
for chid = 1:numel(spikeID_im)
    Stats.imgTuneStats{chid} = calc_tune_stats(cellfun(@(idx) rasters_img(chid, :, idx), idx_arr(:, centcol:end), 'Uni', 0));
end
Stats.imgTuneStats = cell2mat(Stats.imgTuneStats);
%% Verbal summary for the table. 
Tab = struct2table(Stats.imgTuneStats);
ITmsk = spikeID_im<=32; V1msk = spikeID_im>=33 & spikeID_im<=48; V4msk = spikeID_im>=49;
fprintf("Response Delay window image [%d,%d] ms\n", 51, 200)
fprintf("%d / %d channels are visually responsive (T test p < 0.01). %d / %d for p < 0.001\n",sum(Tab.t_P<0.01),numel(Tab.t_P), sum(Tab.t_P<0.001),numel(Tab.t_P))
fprintf("%d / %d IT,  %d / %d V1, %d / %d V4 channels are visually responsive (T test p < 0.01)\n", sum(Tab.t_P<0.01 & ITmsk), sum(ITmsk), ...
											sum(Tab.t_P<0.01 & V1msk),sum(V1msk), sum(Tab.t_P<0.01 & V4msk), sum(V4msk))
fprintf("%d / %d channels have response modulated by image identity (ANOVA p < 0.01). %d / %d for p < 0.001\n",sum(Tab.F_P<0.01),numel(Tab.F_P),sum(Tab.F_P<0.001),numel(Tab.F_P))
fprintf("%d / %d IT,  %d / %d V1, %d / %d V4 channels are modulated (ANOVA p < 0.01)\n", sum(Tab.F_P<0.01 & ITmsk), sum(ITmsk), ...
											sum(Tab.F_P<0.01 & V1msk),sum(V1msk), sum(Tab.F_P<0.01 & V4msk), sum(V4msk))
fprintf("CONTROL: %d / %d channels' baseline activity modulated by image identity (ANOVA p < 0.01). %d / %d for p < 0.001\n",sum(Tab.F_P_bsl<0.01),numel(Tab.F_P_bsl),sum(Tab.F_P_bsl<0.001),numel(Tab.F_P_bsl))
fprintf("%d / %d IT,  %d / %d V1, %d / %d V4 channels are modulated (ANOVA p < 0.01)\n", sum(Tab.F_P_bsl<0.01 & ITmsk), sum(ITmsk), ...
											sum(Tab.F_P_bsl<0.01 & V1msk),sum(V1msk), sum(Tab.F_P_bsl<0.01 & V4msk), sum(V4msk))
%
FsignfMsk = Tab.t_P<0.01 & Tab.F_P<0.01; 
FsignfCorrMsk = Tab.t_P<0.01 & Tab.F_P<0.01 & corr_P_arr' < 0.01; 
fprintf("Response Delay window image [%d,%d] ms\n", 51, 200)
fprintf("From %d channels in total\n%d / %d selective channels (F, p<0.01) have a correlated selectivity in Movie vs Static image. \n", numel(spikeID_im),sum(FsignfCorrMsk), sum(FsignfMsk)) % %d / %d for p < 0.001
fprintf("%d / %d IT,  %d / %d V1, %d / %d V4 selective channels (F, p<0.01) have a significantly correlated (p<0.01) selectivity.\n", ...
	sum(FsignfCorrMsk & ITmsk), sum(FsignfMsk & ITmsk), sum(FsignfCorrMsk & V1msk), sum(FsignfMsk & V1msk), sum(FsignfCorrMsk & V4msk), sum(FsignfMsk & V4msk))
%%
writetable(Tab, fullfile(figdir, "TuneStats.csv"))
save(fullfile(figdir, "ExpStats.mat"), 'Stats', 'MovImgCorrStats','ImgrspDelayWdw','MvrspDelayWdw')
%% Graphic summary on the Distribution of Correlations
Pthresh1 = min(corr_arr(corr_P_arr < 0.01));
figure(2);clf;hold on
histogram(corr_arr(V1msk & Tab.F_P<0.01),10,'FaceAlpha',0.4)
histogram(corr_arr(V4msk & Tab.F_P<0.01),10,'FaceAlpha',0.4)
histogram(corr_arr(ITmsk & Tab.F_P<0.01),15,'FaceAlpha',0.4)
xlim([-0.3,1])
vline([Pthresh1],'r-.',{"P 0.01"})
xlabel("Correlation of Response")
legend(["V1","V4","IT"]) 
title(compose("Correlation of Movie-Image Firing Rate Rsp %s Selective Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas_sel.png",Animal)))
savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas_sel.fig",Animal)))

figure(2);clf;hold on
histogram(corr_arr(V1msk),10,'FaceAlpha',0.4)
histogram(corr_arr(V4msk),10,'FaceAlpha',0.4)
histogram(corr_arr(ITmsk),15,'FaceAlpha',0.4)
vline([Pthresh1],'r-.',{"P 0.01"})
xlabel("Correlation of Response")
legend(["V1","V4","IT"]) 
title(compose("Correlation of Movie-Image Firing Rate Rsp %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas.png",Animal)))
savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas.fig",Animal)))

figure(2);clf;hold on
histogram(corr_arr(Tab.F_P<0.01),10,'FaceAlpha',0.4)
histogram(corr_arr,15,'FaceAlpha',0.4)
vline([Pthresh1],'r-.',{"P 0.01"})
xlabel("Correlation of Response")
legend(["selective", "all"]) 
title(compose("Correlation of Movie-Image Firing Rate Rsp %s Selective Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s_sel.png",Animal)))
savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s_sel.fig",Animal)))

figure(2);clf;hold on
histogram(corr_arr,15,'FaceAlpha',0.4)
vline([Pthresh1],'r-.',{"P 0.01"})
xlabel("Correlation of Response")
legend(["all"]) 
title(compose("Correlation of Movie-Image Firing Rate Rsp %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s.png",Animal)))
savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s.fig",Animal)))

%%
% wvfrm_distmat = pdist2(meta_img.wvfms, meta_mov.wvfms,'correlation');
% figure;imagesc(wvfrm_distmat)
function [sortedMovnm, sortedRows, sortIdx] = sortMovieNames(movnm)
% Sort the movie names in **lexicoGraphical order**, by space and by eigen idx. 
% This assumes the names to have the structure like "class_eig17_shortshort"
eigi_cell = cellfun(@(eigi){str2double(eigi{1})},regexp(movnm,"_eig(\d*)_",'tokens'));
space_cell = cellfun(@(eigi){eigi{1}{1}},regexp(movnm,"(.*)_eig(\d*)_",'tokens'));
[sortedRows, sortIdx] = sortrows([space_cell,eigi_cell]);
sortedMovnm = movnm(sortIdx);
end

function [MovImgCorrStats, corr_arr, corr_P_arr, corr_sep_arr, corr_sep_P_arr] = HessMovImgMatchCorr(ImgrspDelayWdw, MvrspDelayWdw, psthimg_mean, psthmov_mean, S)
% Note: This is limited to Hessian Movies. The correlation part need to be
%   rewrite for non-Hessian
% 
% Input parameters:
%   ImgrspDelayWdw / MvrspDelayWdw: index array of the window to compute
%       e.g. [61:200], [-59:60]
%   S: Stats formed in beforehand. 

MovImgCorrStats = repmat(struct(),1,numel(S.spikeID_im)); % collect results here. 
for chid = 1:numel(S.spikeID_im) % loop through static image units (Loop through movie is also fine). 
iCh = S.spikeID_im(chid); 
iU = S.unitID_im(chid); 

rspmat_static = cellfun(@(psth) mean(psth(chid, ImgrspDelayWdw),[1,2]), psthimg_mean); 
psth_static = cellfun(@(psth) psth(chid, ImgrspDelayWdw), psthimg_mean, 'Uni', 0); 

chid2 = find(S.spikeID_mv==iCh & S.unitID_mv==iU); % Corresponding unit in movie exp. 
rspmat_movie = nan(numel(S.movnm), S.nImgPerMv); 
psth_movie = cell(numel(S.movnm), S.nImgPerMv); 
for iMv = 1:numel(S.movnm) % movies
	matchONidx = round(S.matchTONs{iMv}) - S.MvRstrWdw(1) + 1; % find the index in the raster matrix. Shifting the window by wdw(1) - 1;
    % rspmat_movie(iMv, :, :) = arrayfun(@(iTon) mean(psthmov_mean{iMv}(chid2, iTon+MvrspDelayWdw)), matchONidx); 
    for iSta = 1:size(matchONidx, 1) % Static images Embedded in the movie
        for iOccr = 1:size(matchONidx,2) % Occurance of that image 
           iTon = matchONidx(iSta, iOccr);
           rspmat_movie(iMv, iSta, iOccr) = mean(psthmov_mean{iMv}(chid2, iTon+MvrspDelayWdw));
           psth_movie{iMv, iSta}(end+1,:) = psthmov_mean{iMv}(chid2, iTon+MvrspDelayWdw);
        end
    end
end
% rspmat_movie is # Movie - by - # Static in the Movie - by - # Occurence (12, 5, 2)
% Correlate the non-center image responses
rspmat_key = rspmat_static(:,S.centcol+1:end);
rspmat_movie_key = mean(rspmat_movie(:,2:5,:), 3);
[cc,pp] = corr(rspmat_key(:),rspmat_movie_key(:));
MovImgCorrStats(chid).corr = cc;
MovImgCorrStats(chid).corr_P = pp;

rspmat_movie_vecs = reshape(rspmat_movie(:,2:5,:),[],2);
[cc,pp] = corr(rspmat_key(:), rspmat_movie_vecs);
MovImgCorrStats(chid).corr_sep = cc;
MovImgCorrStats(chid).corr_sep_P = pp;
MovImgCorrStats(chid).rspmat_movie = rspmat_movie;
MovImgCorrStats(chid).rspmat_static = rspmat_static(:,S.centcol+1:end);
MovImgCorrStats(chid).psth_movie = psth_movie; 
MovImgCorrStats(chid).psth_static = psth_static(:,S.centcol+1:end); 
MovImgCorrStats(chid).iCh = iCh;
MovImgCorrStats(chid).iU = iU;
end
corr_arr = arrayfun(@(S)S.corr, MovImgCorrStats); 
corr_P_arr = arrayfun(@(S)S.corr_P, MovImgCorrStats); 
corr_sep_arr = cell2mat(arrayfun(@(S)S.corr_sep, MovImgCorrStats','Uni',0)); % each column is an occurence
corr_sep_P_arr = cell2mat(arrayfun(@(S)S.corr_sep_P, MovImgCorrStats','Uni',0));

fprintf("Response Delay window movie [%d,%d] image [%d,%d] ms\n",MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end))
fprintf("%d / %d channels has significant(0.01) correlation between Movie and Image Response\n",sum(corr_P_arr<0.01),numel(S.spikeID_im))
fprintf("%d / %d channels has significant(0.001) correlation\n",sum(corr_P_arr<0.001),numel(S.spikeID_im))
fprintf("%d / %d IT channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & S.spikeID_im'<=32),sum(S.spikeID_im<=32))
fprintf("%d / %d V1 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & S.spikeID_im'>=33 & S.spikeID_im'<=48),sum(S.spikeID_im'>=33 & S.spikeID_im'<=48))
fprintf("%d / %d V4 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & S.spikeID_im'>=49),sum(S.spikeID_im'>=49))
fprintf("%s ",S.unit_str_im(corr_P_arr<0.01))
fprintf("\n")
end

function h = plot_channel_corr(Stats, MovImgCorrStats, ImgrspDelayWdw, MvrspDelayWdw)
h = figure;set(h,'pos',[1000         462         560         520])
for chid = 1:numel(MovImgCorrStats)
	set(0, 'CurrentFigure', h);hold off;
    iCh = MovImgCorrStats(chid).iCh;
    iU = MovImgCorrStats(chid).iU;
    static_rspmat_key = MovImgCorrStats(chid).rspmat_static;
    movie_rspmat_key = MovImgCorrStats(chid).rspmat_movie;
    scatter(static_rspmat_key(:), movie_rspmat_key(:));
    axis equal; addDiagonal();
    title(compose("Correlation of Movie-Image Firing Rate Rsp %s Chan%d Unit%s\nCorr Coef %.3f(%.1e)\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Stats.Animal, iCh, char(64+iU), MovImgCorrStats(chid).corr, MovImgCorrStats(chid).corr_P, ...
            MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
    xlabel("Static Firing Rate");ylabel("Movie Firing Rate")
    saveas(h,fullfile(Stats.figdir,compose("movie_image_corr_%s_%d%s.png",Stats.Animal,iCh,char(64+iU)))) % Stats.unit_str_im(chid)
end
end

function addDiagonal()
ax = gca;
XLIM = xlim(ax);YLIM = ylim(ax);
% MAX = min(YLIM(2),XLIM(2));
% MIN = max(YLIM(1),XLIM(1));
LIM = [min(XLIM(1),YLIM(1)), max(XLIM(2),YLIM(2))];
hold on; line(LIM,LIM,'Color','r')
xlim(max(0, LIM));ylim(max(0, LIM));
end

% %% Get the Response Mat and PSTH in that period, clipping machine for psth
% MvrspDelayWdw = [60:200];
% ImgrspDelayWdw = [60:200];
% MovImgCorrStats = repmat(struct(),1,numel(spikeID_im)); % collect results here. 
% for chid = 1:numel(spikeID_im) % loop through static image units (Loop through movie is also fine). 
% iCh = spikeID_im(chid); 
% iU = unitID_im(chid); 
% 
% rspmat_static = cellfun(@(psth) mean(psth(chid, ImgrspDelayWdw),[1,2]), psthimg_mean); 
% psth_static = cellfun(@(psth) psth(chid, ImgrspDelayWdw), psthimg_mean, 'Uni', 0); 
% 
% chid2 = find(spikeID_mv==iCh & unitID_mv==iU); % Corresponding unit in movie exp. 
% rspmat_movie = nan(numel(movnm_sorted), nImgPerMv, 2); % rspmat_movie is # Movie - by - # Static in the Movie - by - # Occurence (12, 5, 2)
% psth_movie = cell(numel(movnm_sorted), nImgPerMv);
% for iMv = 1:numel(movnm_sorted) % movies
% 	matchONidx = round(matchTONs{iMv}) - wdw(1) + 1; % find the index in the raster matrix. Shifting the window by wdw(1) - 1;
%     % rspmat_movie(iMv, :, :) = arrayfun(@(iTon)mean(psthmov_mean{iMv}(chid2, iTon+MvrspDelayWdw)), matchONidx); 
%     for iSta = 1:size(matchONidx, 1) % Static images Embedded in the movie
%         for iOccr = size(matchONidx,2) % Occurance of that image 
%            iTon = matchONidx(iSta, iOccr);
%            rspmat_movie(iMv, iSta, iOccr) = mean(psthmov_mean{iMv}(chid2, iTon+MvrspDelayWdw));
%            psth_movie{iMv, iSta}(end+1,:) = psthmov_mean{iMv}(chid2, iTon+MvrspDelayWdw);
%         end
%     end
% end
% % Correlate the non-center image responses
% rspmat_key = rspmat_static(:,6:9);
% rspmat_movie_key = mean(rspmat_movie(:,2:5,:), 3);
% [cc,pp] = corr(rspmat_key(:),rspmat_movie_key(:));
% MovImgCorrStats(chid).corr = cc;
% MovImgCorrStats(chid).corr_P = pp;
% MovImgCorrStats(chid).rspmat_movie = rspmat_movie;
% MovImgCorrStats(chid).rspmat_static = rspmat_static;
% MovImgCorrStats(chid).iCh = iCh;
% MovImgCorrStats(chid).iU = iU;
% end
% corr_P_arr = arrayfun(@(S)S.corr_P,MovImgCorrStats);
% corr_arr = arrayfun(@(S)S.corr,MovImgCorrStats);