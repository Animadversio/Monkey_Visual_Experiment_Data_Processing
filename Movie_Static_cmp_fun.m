%% Very similar to Movie Hess Static cmp fun. But For natural images. 
%  The major difference is the naming convention of frames and movies. 
% Note this code is written extremely concise! 
% 
clearvars -except Trials* meta* rasters*
%%
Animal = "Alfa"; Set_Path; setMatlabTitle("MovieStaticDynamics")
ftr = find(contains(ExpRecord.ephysFN,"Alfa-16122020"));%"Alfa-06112020"
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr, Animal);
%%
figroot = "E:\OneDrive - Washington University in St. Louis\MovieDynamics";
Trials_mov = Trials_new{2};
rasters_mov = rasters_new{2};
meta_mov = meta_new{2};
Trials_img = Trials_new{5};
rasters_img = rasters_new{5};
meta_img = meta_new{5};
% figdir will be an input argument.
prefchan = ExpRecord.pref_chan(contains(ExpRecord.ephysFN,meta_mov.ephysFN));
stimparts = split(meta_mov.stimuli,'\');
fdrname = compose("%s-Chan%02d",stimparts{end},prefchan);
figdir = fullfile(figroot, fdrname);
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
movnm_sorted = movnm;
mov_idx_arr = arrayfun(@(mv)find(contains(Trials_mov.imageName, mv)),movnm_sorted,'Uni',0);

psthmov_all = cellfun(@(idx) rasters_mov(:, :, idx), mov_idx_arr, 'Uni',0);
psthmov_mean = cellfun(@(idx) mean(rasters_mov(:, :, idx),3), mov_idx_arr, 'Uni',0);
psthmov_sem = cellfun(@(idx) std(rasters_mov(:, :, idx),1,3) / sqrt(numel(idx)), mov_idx_arr, 'Uni',0);
%% Meta information of Images & spike Id
spikeID_im = meta_img.spikeID;
unitID_im = meta_img.unitID;
unit_str_im = generate_unit_labels_new(meta_img.spikeID, meta_img.unitID);
% Sort Static images and trials
[imgnm_per_mov, idx_arr, imgnm_arr, mvnms_uniq] = parse_frame_idx_arr(Trials_img.imageName);
[img_mvnm_arr, img_frnum_arr, uniq_imgnms] = parse_frame_name(Trials_img.imageName);
nImgPerMv = max(cellfun(@numel, imgnm_per_mov)); 

psthimg_all = cellfun(@(idx) rasters_img(:, :, idx), idx_arr, 'Uni',0);
psthimg_mean = cellfun(@(idx) mean(rasters_img(:,:,idx),3), idx_arr,'Uni',0);
psthimg_sem = cellfun(@(idx) std(rasters_img(:,:,idx),1,3)/sqrt(numel(idx)), idx_arr,'Uni',0);
% imgnm_per_mov = mat2cell(imgnm_arr, ones(size(imgnm_arr,1),1), size(imgnm_arr,2)); % Cell array of the images in 
%% Matching images into frames process (takes a few seconds)
% Matching all movie image pairs. Currently requires LPIPS and torch to be
% available. 
[matchfrids, matchTONs, matchTOFFs] = closest_Kframe_locate(movnm_sorted(:), imgnm_per_mov, ...
                                        meta_mov.stimuli, meta_img.stimuli, false, 1); % halfSep = false,
%% Records the computed info into a structure
Stats.figdir = figdir;
Stats.Animal = Animal;
Stats.movnm = movnm_sorted;
Stats.mov_idx_arr = mov_idx_arr;
Stats.imgnm_arr = imgnm_arr;
Stats.img_idx_arr = idx_arr;
Stats.imgnm_per_mov = imgnm_per_mov;
Stats.nImgPerMv = nImgPerMv;
% Static-Movie Matching
Stats.matchfrids = matchfrids;
Stats.matchTONs = matchTONs;
Stats.matchTOFFs = matchTOFFs;
% units meta.
Stats.meta_im = meta_img; % copy the meta for img and mov
Stats.meta_mv = meta_mov;
Stats.unit_str_im = unit_str_im;
Stats.unit_str_mv = unit_str_mv;
Stats.spikeID_im = spikeID_im; 
Stats.unitID_im = unitID_im; 
Stats.spikeID_mv = spikeID_mv; 
Stats.unitID_mv = unitID_mv; 
Stats.MvRstrWdw = wdw; % Movie Raster window. Need this to interpret the raster's timeline.
%% Compute basic tuning property
Stats.imgTuneStats = {};
for chid = 1:numel(spikeID_im)
    Stats.imgTuneStats{chid} = calc_tune_stats(cellfun(@(idx) rasters_img(chid, :, idx), idx_arr, 'Uni', 0));
end
Stats.imgTuneStats = cell2mat(Stats.imgTuneStats);
%%
Tab = struct2table(Stats.imgTuneStats);
writetable(Tab, fullfile(figdir, "TuneStats.csv"))
ITmsk = spikeID_im<=32; V1msk = spikeID_im>=33 & spikeID_im<=48; V4msk = spikeID_im>=49;
diary(fullfile(figdir, "SummmaryOutput.log"))
fprintf("Static Image: Response Delay window [%d,%d] ms\n", 51, 200)
fprintf("%d / %d channels are visually responsive (T test p < 0.01). %d / %d for p < 0.001\n",sum(Tab.t_P<0.01),numel(Tab.t_P), sum(Tab.t_P<0.001),numel(Tab.t_P))
fprintf("%d / %d IT,  %d / %d V1, %d / %d V4 channels are visually responsive (T test p < 0.01)\n", sum(Tab.t_P<0.01 & ITmsk), sum(ITmsk), ...
											sum(Tab.t_P<0.01 & V1msk),sum(V1msk), sum(Tab.t_P<0.01 & V4msk), sum(V4msk))
fprintf("%d / %d channels have response modulated by image identity (ANOVA p < 0.01). %d / %d for p < 0.001\n",sum(Tab.F_P<0.01),numel(Tab.F_P),sum(Tab.F_P<0.001),numel(Tab.F_P))
fprintf("%d / %d IT,  %d / %d V1, %d / %d V4 channels are modulated (ANOVA p < 0.01)\n", sum(Tab.F_P<0.01 & ITmsk), sum(ITmsk), ...
											sum(Tab.F_P<0.01 & V1msk),sum(V1msk), sum(Tab.F_P<0.01 & V4msk), sum(V4msk))
fprintf("CONTROL: %d / %d channels' baseline activity modulated by image identity (ANOVA p < 0.01). %d / %d for p < 0.001\n",sum(Tab.F_P_bsl<0.01),numel(Tab.F_P_bsl),sum(Tab.F_P_bsl<0.001),numel(Tab.F_P_bsl))
fprintf("%d / %d IT,  %d / %d V1, %d / %d V4 channels are modulated (ANOVA p < 0.01)\n", sum(Tab.F_P_bsl<0.01 & ITmsk), sum(ITmsk), ...
											sum(Tab.F_P_bsl<0.01 & V1msk),sum(V1msk), sum(Tab.F_P_bsl<0.01 & V4msk), sum(V4msk))
%%
% ImgrspDelayWdw = [80:200]; MvrspDelayWdw = [80:200];%[-20:100];%
% [MovImgCorrStats, corr_arr, corr_P_arr, corr_sep_arr, corr_sep_P_arr] = MovImgMatchCorr(ImgrspDelayWdw, MvrspDelayWdw, psthimg_mean, psthmov_mean, Stats);
ImgrspDelayWdw = [80:200]; MvrspDelayWdw = [80:200];%[-20:100];%
[MovImgCorrStats] = MovImgCorrVariability(ImgrspDelayWdw, MvrspDelayWdw, psthimg_all, psthmov_all, Stats);
visualSummary(MovImgCorrStats, Stats, ImgrspDelayWdw, MvrspDelayWdw)
ImgrspDelayWdw = [80:200]; MvrspDelayWdw = [-20:100];%[80:200];%;%
[MovImgCorrStats] = MovImgCorrVariability(ImgrspDelayWdw, MvrspDelayWdw, psthimg_all, psthmov_all, Stats);
visualSummary(MovImgCorrStats, Stats, ImgrspDelayWdw, MvrspDelayWdw)
diary off
%%
FsignfMsk = Tab.t_P<0.01 & Tab.F_P<0.01; 
FsignfCorrMsk = Tab.t_P<0.01 & Tab.F_P<0.01 & corr_P_arr' < 0.01; 
fprintf("Response Delay window image [%d,%d] ms\n", 51, 200)
fprintf("From %d channels in total\n%d / %d selective channels (F, p<0.01) have a correlated selectivity in Movie vs Static image. \n", numel(spikeID_im),sum(FsignfCorrMsk), sum(FsignfMsk)) % %d / %d for p < 0.001
fprintf("%d / %d IT,  %d / %d V1, %d / %d V4 selective channels (F, p<0.01) have a significantly correlated (p<0.01) selectivity.\n", ...
	sum(FsignfCorrMsk & ITmsk), sum(FsignfMsk & ITmsk), sum(FsignfCorrMsk & V1msk), sum(FsignfMsk & V1msk), sum(FsignfCorrMsk & V4msk), sum(FsignfMsk & V4msk))

%% Save the computation result! 
save(fullfile(figdir, "ExpStats.mat"), 'Stats', 'MovImgCorrStats', 'ImgrspDelayWdw', 'MvrspDelayWdw') % this statistic is really comprehensive! 
writetable(Tab, fullfile(figdir, "TuneStats.csv"))
%% Visualize individual channel for all  
plot_channel_corr(Stats, MovImgCorrStats, ImgrspDelayWdw, MvrspDelayWdw, true); % with error bar
plot_channel_corr(Stats, MovImgCorrStats, ImgrspDelayWdw, MvrspDelayWdw, false); % without error bar

%%
function [imgnm_per_movie, idx_arr, imgnm_arr, mvnms_uniq] = parse_frame_idx_arr(imageName)
% Equivalent to the `parse_image_idx_arr_hess` for non-Hessian experiments
% The key outputs are an idx_arr, 
% Outputs: 
%  idx_arr: the trial indexing for images organized in a 2d cell array. 
%       each row is one movie. 
uniq_nms = string(unique(imageName));
matchtoks = regexp(uniq_nms,"(.*)_(\d*)",'tokens');
mvnm_arr = cellfun(@(tok)tok{1}{1}, matchtoks, 'Uni', 0); 
frnum_arr = cellfun(@(tok)str2double(tok{1}{2}), matchtoks, 'Uni', 1); 
mvnms_uniq = unique(mvnm_arr);
imgnm_per_movie = cellfun(@(nm) uniq_nms(strcmp(mvnm_arr, nm))', mvnms_uniq, 'uni',0);
ncol = max(cellfun(@numel, imgnm_per_movie));
imgnm_arr = cell(numel(mvnms_uniq), ncol); % nMovie by nFrames
idx_arr = cell(numel(mvnms_uniq), ncol);
for iMv = 1:numel(mvnms_uniq)
    for iMat = 1:numel(imgnm_per_movie{iMv})
        imgnm = imgnm_per_movie{iMv}(iMat);
        imgnm_arr{iMv, iMat} = imgnm;
        idx_arr{iMv, iMat} = find(strcmp(imageName, imgnm));
    end
end
end

function [mvnm_arr, frnum_arr, uniq_nms] = parse_frame_name(imageName)
uniq_nms = string(unique(imageName));
matchtoks = regexp(uniq_nms,"(.*)_(\d*)",'tokens');
mvnm_arr = string(cellfun(@(tok)tok{1}{1}, matchtoks, 'Uni', 0)); 
frnum_arr = cellfun(@(tok)str2double(tok{1}{2}), matchtoks, 'Uni', 1); 
end

function [MovImgCorrStats, corr_arr, corr_P_arr, corr_sep_arr, corr_sep_P_arr] = MovImgMatchCorr(ImgrspDelayWdw, MvrspDelayWdw, psthimg_mean, psthmov_mean, S)
% Note: This is limited to Hessian Movies. The correlation part need to be
%   rewrite for non-Hessian.
% Note: This function is now obsolete use MovImgCorrVariability instead.
% 
% Input parameters:
%   ImgrspDelayWdw / MvrspDelayWdw: index array of the window to compute
%       e.g. [61:200], [-59:60]
%   psthimg_mean: Trial averaged PSTH of each image. An cell array of shape
%       (nMovie, nImginMovie) eahc array in it is of shape (nChannel, nTime)
%   psthmov_mean: Trial averaged PSTH of each movie. An cell array of shape
%       (nMovie, ) each array in it is of shape (nChannel, nTime)
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
rspmat_key = rspmat_static(:,:);
rspmat_movie_key = mean(rspmat_movie(:,:,:), 3);
[cc,pp] = corr(rspmat_key(:),rspmat_movie_key(:), 'rows', 'complete');
MovImgCorrStats(chid).corr = cc;
MovImgCorrStats(chid).corr_P = pp;
% Paired ttest of image - frame response. 2 sided.
[H,P,~,STAT] = ttest(rspmat_key(:),rspmat_movie_key(:)); 
MovImgCorrStats(chid).T = STAT.tstat; % t > 0 show the response to image > movie
MovImgCorrStats(chid).T_P = P;
% Correlation of image response to response of each occurance of the image
rspmat_movie_vecs = reshape(rspmat_movie(:,:,:),[], size(rspmat_movie,3));
[cc,pp] = corr(rspmat_key(:), rspmat_movie_vecs, 'rows', 'complete');
MovImgCorrStats(chid).corr_sep = cc;
MovImgCorrStats(chid).corr_sep_P = pp;
MovImgCorrStats(chid).rspmat_movie = rspmat_movie;
MovImgCorrStats(chid).rspmat_static = rspmat_static;
MovImgCorrStats(chid).psth_movie = psth_movie; 
MovImgCorrStats(chid).psth_static = psth_static; 
MovImgCorrStats(chid).iCh = iCh;
MovImgCorrStats(chid).iU = iU;
end
corr_arr = arrayfun(@(S)S.corr, MovImgCorrStats); 
corr_P_arr = arrayfun(@(S)S.corr_P, MovImgCorrStats); 
corr_sep_arr = cell2mat(arrayfun(@(S)S.corr_sep, MovImgCorrStats','Uni',0)); % each column is an occurence
corr_sep_P_arr = cell2mat(arrayfun(@(S)S.corr_sep_P, MovImgCorrStats','Uni',0));
T_arr = arrayfun(@(S)S.T, MovImgCorrStats); 
T_P_arr = arrayfun(@(S)S.T_P, MovImgCorrStats); 
% Print summary string
ITmsk = S.spikeID_im'<=32; V1msk = S.spikeID_im'>=33 & S.spikeID_im'<=48; V4msk = S.spikeID_im'>=49;
fprintf("Response Delay window movie [%d,%d] image [%d,%d] ms\n",MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end))
fprintf("%d / %d channels has significant(0.01) correlation between Movie and Image Response\n",sum(corr_P_arr<0.01),numel(S.spikeID_im))
fprintf("%d / %d channels has significant(0.001) correlation\n",sum(corr_P_arr<0.001),numel(S.spikeID_im))
fprintf("%d / %d IT channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & ITmsk),sum(ITmsk))
fprintf("%d / %d V1 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & V1msk),sum(V1msk))
fprintf("%d / %d V4 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & V4msk),sum(V4msk))
fprintf("%s ",S.unit_str_im(corr_P_arr<0.01)) % list of channels that are strongly correlated 
fprintf("\n")
fprintf("%d / %d channels has significant (0.01) difference between firing rate of paired Movie and Image Response\n",sum(T_P_arr<0.01),numel(S.spikeID_im))
fprintf("In IT, %d / %d Static > Movie, %d / %d Movie > Static  significantly\n",sum(T_P_arr<0.01 & T_arr >0 & ITmsk),sum(ITmsk),sum(T_P_arr<0.01 & T_arr <0 & ITmsk),sum(ITmsk))
fprintf("In V1, %d / %d Static > Movie, %d / %d Movie > Static  significantly\n",sum(T_P_arr<0.01 & T_arr >0 & V1msk),sum(V1msk),sum(T_P_arr<0.01 & T_arr <0 & V1msk),sum(V1msk))
fprintf("In V4, %d / %d Static > Movie, %d / %d Movie > Static  significantly\n",sum(T_P_arr<0.01 & T_arr >0 & V4msk),sum(V4msk),sum(T_P_arr<0.01 & T_arr <0 & V4msk),sum(V4msk))
end

function [MovImgCorrStats] = MovImgCorrVariability(ImgrspDelayWdw, MvrspDelayWdw, psthimg_all, psthmov_all, S)
% Note: This is not limited to Hessian Movies. The correlation part has
%   been rewriten for non-Hessian
% 
% Input parameters:
%   ImgrspDelayWdw / MvrspDelayWdw: index array of the window to compute
%       e.g. [61:200], [-59:60]. Relative to image onset and frame onset.
%   psthimg_all: Single Trial PSTH of each image. An cell array of shape
%       (nMovie, nImginMovie) each array in it is of shape (nChannel, nTime, nTrials)
%   psthmov_all: Single Trial PSTH of each movie. An cell array of shape
%       (nMovie, ) each array in it is of shape (nChannel, nTime, nTrials)
%   S: Stats formed in beforehand. 
MovImgCorrStats = repmat(struct(),1,numel(S.spikeID_im)); % collect results here. 
for chid = 1:numel(S.spikeID_im) % loop through static image units (Loop through movie is also fine). 
iCh = S.spikeID_im(chid); 
iU = S.unitID_im(chid); 

tr_rspmat_static = cellfun(@(psth) squeeze(mean(psth(chid, ImgrspDelayWdw, :),[1,2])), psthimg_all, 'Uni', 0); % cell (nMovie, nImginMovie) with (nTrials,) shape elements
tr_psth_static = cellfun(@(psth) psth(chid, ImgrspDelayWdw, :), psthimg_all, 'Uni', 0); % cell (nMovie, nImginMovie) with (1, TimeWindow, nTrials,) shape elements
emptyImgEntry = find(cellfun(@isempty, psthimg_all));
for idx = emptyImgEntry', tr_psth_static{idx} = []; end
chid2 = find(S.spikeID_mv==iCh & S.unitID_mv==iU); % Corresponding unit in movie exp. 
tr_rspmat_movie = cell(numel(S.movnm), S.nImgPerMv); 
tr_psth_movie = cell(numel(S.movnm), S.nImgPerMv); 
for iMv = 1:numel(S.movnm) % movies
	matchONidx = round(S.matchTONs{iMv}) - S.MvRstrWdw(1) + 1; % find the index in the raster matrix. Shifting the window by wdw(1) - 1;
    % rspmat_movie(iMv, :, :) = arrayfun(@(iTon) mean(psthmov_mean{iMv}(chid2, iTon+MvrspDelayWdw)), matchONidx); 
    for iSta = 1:size(matchONidx, 1) % Static images Embedded in the movie
        for iOccr = 1:size(matchONidx,2) % Occurance of that image 
           iTon = matchONidx(iSta, iOccr);
           tr_rspmat_movie{iMv, iSta, iOccr} = squeeze(mean(psthmov_all{iMv}(chid2, iTon+MvrspDelayWdw, :),[1,2])); % (nTrials, ) shaped
           tr_psth_movie{iMv, iSta, iOccr} = psthmov_all{iMv}(chid2, iTon+MvrspDelayWdw, :); % (1, TimeWindow, nTrials, ) shaped elements
        end
    end
end
% ANOVA of individual trials
anova_mov = anova_cells(tr_rspmat_movie);
anova_img = anova_cells(tr_rspmat_static);
MovImgCorrStats(chid).F_im = anova_img.F;
MovImgCorrStats(chid).F_P_im = anova_img.F_P;
MovImgCorrStats(chid).F_mv = anova_mov.F;
MovImgCorrStats(chid).F_P_mv = anova_mov.F_P;
% Collapse the trial dimension to get std and sem of each group
% rspmat_movie is # Movie - by - # Static in the Movie - by - # Occurence (12, 5, 2)
rspmat_movie = cellfun(@mean, tr_rspmat_movie);
rspstd_movie = cellfun(@std, tr_rspmat_movie);
rspsem_movie = cellfun(@(rsps) std(rsps)/sqrt(numel(rsps)), tr_rspmat_movie);
rspmat_static = cellfun(@mean, tr_rspmat_static);
rspstd_static = cellfun(@std, tr_rspmat_static);
rspsem_static = cellfun(@(rsps) std(rsps)/sqrt(numel(rsps)), tr_rspmat_static);
psth_movie = cellfun(@(psth)mean(psth, [3]), tr_psth_movie, 'Uni', 0);
psth_static = cellfun(@(psth)mean(psth, [3]), tr_psth_static, 'Uni', 0);
% Correlate the non-center image responses
rspmat_static_key = rspmat_static(:,:);
rspmat_movie_key = mean(rspmat_movie(:,:,:), 3);
[cc,pp] = corr(rspmat_static_key(:),rspmat_movie_key(:), 'rows', 'complete');
MovImgCorrStats(chid).corr = cc;
MovImgCorrStats(chid).corr_P = pp;
% Paired ttest of image - frame response. 2 sided.
[H,P,~,STAT] = ttest(rspmat_static_key(:),rspmat_movie_key(:)); 
MovImgCorrStats(chid).T = STAT.tstat; % t > 0 show the response to image > movie
MovImgCorrStats(chid).T_P = P;
% Correlation of image response to response of each occurance of the image
rspmat_movie_vecs = reshape(rspmat_movie(:,:,:),[], size(rspmat_movie,3));
[cc,pp] = corr(rspmat_static_key(:), rspmat_movie_vecs, 'rows', 'complete');
% Save some of response info for plotting.
MovImgCorrStats(chid).corr_sep = cc;
MovImgCorrStats(chid).corr_sep_P = pp;
MovImgCorrStats(chid).rspmat_movie = rspmat_movie; 
MovImgCorrStats(chid).rspstd_movie = rspstd_movie; 
MovImgCorrStats(chid).rspsem_movie = rspsem_movie; 
MovImgCorrStats(chid).rspmat_static = rspmat_static; 
MovImgCorrStats(chid).rspstd_static = rspstd_static; 
MovImgCorrStats(chid).rspsem_static = rspsem_static; 
MovImgCorrStats(chid).psth_movie = psth_movie; 
MovImgCorrStats(chid).psth_static = psth_static; 
MovImgCorrStats(chid).iCh = iCh;
MovImgCorrStats(chid).iU = iU;
end
% summarize channel stats into population
corr_arr = arrayfun(@(S)S.corr, MovImgCorrStats); 
corr_P_arr = arrayfun(@(S)S.corr_P, MovImgCorrStats); 
corr_sep_arr = cell2mat(arrayfun(@(S)S.corr_sep, MovImgCorrStats','Uni',0)); % each column is an occurence
corr_sep_P_arr = cell2mat(arrayfun(@(S)S.corr_sep_P, MovImgCorrStats','Uni',0));
T_arr = arrayfun(@(S)S.T, MovImgCorrStats); 
T_P_arr = arrayfun(@(S)S.T_P, MovImgCorrStats); 
F_im_arr = arrayfun(@(S)S.F_im, MovImgCorrStats); 
F_P_im_arr = arrayfun(@(S)S.F_P_im, MovImgCorrStats); 
F_mv_arr = arrayfun(@(S)S.F_mv, MovImgCorrStats); 
F_P_mv_arr = arrayfun(@(S)S.F_P_mv, MovImgCorrStats); 
% Print summary string
ITmsk = S.spikeID_im'<=32; V1msk = S.spikeID_im'>=33 & S.spikeID_im'<=48; V4msk = S.spikeID_im'>=49;
fprintf("Response Delay window movie [%d,%d] image [%d,%d] ms\n",MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end))
fprintf("%d / %d channels has significant(0.01) correlation between Movie and Image Response\n",sum(corr_P_arr<0.01),numel(S.spikeID_im))
fprintf("%d / %d channels has significant(0.001) correlation\n",sum(corr_P_arr<0.001),numel(S.spikeID_im))
fprintf("%d / %d IT channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & ITmsk),sum(ITmsk))
fprintf("%d / %d V1 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & V1msk),sum(V1msk))
fprintf("%d / %d V4 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & V4msk),sum(V4msk))
fprintf("%s ",S.unit_str_im(corr_P_arr<0.01)) % list of channels that are strongly correlated 
fprintf("\n")
fprintf("%d / %d channels has significant (0.01) difference between firing rate of paired Movie and Image Response\n",sum(T_P_arr<0.01),numel(S.spikeID_im))
fprintf("In IT, %d / %d Static > Movie, %d / %d Movie > Static  significantly\n",sum(T_P_arr<0.01 & T_arr >0 & ITmsk),sum(ITmsk),sum(T_P_arr<0.01 & T_arr <0 & ITmsk),sum(ITmsk))
fprintf("In V1, %d / %d Static > Movie, %d / %d Movie > Static  significantly\n",sum(T_P_arr<0.01 & T_arr >0 & V1msk),sum(V1msk),sum(T_P_arr<0.01 & T_arr <0 & V1msk),sum(V1msk))
fprintf("In V4, %d / %d Static > Movie, %d / %d Movie > Static  significantly\n",sum(T_P_arr<0.01 & T_arr >0 & V4msk),sum(V4msk),sum(T_P_arr<0.01 & T_arr <0 & V4msk),sum(V4msk))
fprintf("\n")
fprintf("%d / %d channels are significantly (0.01) modulated by this set of static images\n",sum(F_P_im_arr<0.01),numel(S.spikeID_im))
fprintf("IT: %d / %d ,V1: %d / %d ,V4: %d / %d \n", sum(F_P_im_arr<0.01 & ITmsk),sum(ITmsk),...
		   sum(F_P_im_arr<0.01 & V1msk),sum(V1msk), sum(F_P_im_arr<0.01 & V4msk),sum(V4msk))
fprintf("%d / %d channels are significantly (0.01) modulated by this set of images embedded in movies\n",sum(F_P_mv_arr<0.01),numel(S.spikeID_im))
fprintf("IT: %d / %d ,V1: %d / %d ,V4: %d / %d \n", sum(F_P_mv_arr<0.01 & ITmsk),sum(ITmsk),...
		   sum(F_P_mv_arr<0.01 & V1msk),sum(V1msk), sum(F_P_mv_arr<0.01 & V4msk),sum(V4msk))
wdwstr = compose("im_%d_%d_mv_%d_%d",ImgrspDelayWdw(1), ImgrspDelayWdw(end), MvrspDelayWdw(1), MvrspDelayWdw(end));
end

% Visualization and summary
function visualSummary(MovImgCorrStats, Stats, ImgrspDelayWdw, MvrspDelayWdw,figdir)
% Visually summarize the `MovImgCorrStats` computed beforehand. Majorly histograms ploted with different separation. 
if nargin <5, figdir = Stats.figdir; end
ITmsk = Stats.spikeID_im<=32; V1msk = Stats.spikeID_im>=33 & Stats.spikeID_im<=48; V4msk = Stats.spikeID_im>=49;
Tab = struct2table(Stats.imgTuneStats); Animal = Stats.Animal;
wdwstr = compose("im_%d_%d_mv_%d_%d",ImgrspDelayWdw(1), ImgrspDelayWdw(end), MvrspDelayWdw(1), MvrspDelayWdw(end));
save(fullfile(figdir, compose("CorrStats_%s.mat",wdwstr)), 'Stats', 'MovImgCorrStats', 'ImgrspDelayWdw', 'MvrspDelayWdw')
% Correlation of activation
corr_arr = arrayfun(@(S)S.corr, MovImgCorrStats);
corr_P_arr = arrayfun(@(S)S.corr_P, MovImgCorrStats);
ccPthrsh = max(corr_arr(corr_P_arr > 0.01));
figure(2);clf;hold on
plotmsk = Tab.F_P<0.01;
histogram(corr_arr(V1msk & plotmsk),10,'FaceAlpha',0.4)
histogram(corr_arr(V4msk & plotmsk),10,'FaceAlpha',0.4)
histogram(corr_arr(ITmsk & plotmsk),15,'FaceAlpha',0.4)
xlim([-0.3,1])
vline([ccPthrsh],'r-.',{"P 0.01"})
xlabel("Correlation of Response")
legend(["V1","V4","IT"]) 
title(compose("Correlation of Movie-Image Firing Rate Rsp %s Selective Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movImgCorrHist_%s_areas_sel_%s.png",Animal,wdwstr)))
savefig(2,fullfile(figdir,compose("movImgCorrHist_%s_areas_sel_%s.fig",Animal,wdwstr)))

figure(2);clf;hold on
histogram(corr_arr(V1msk),10,'FaceAlpha',0.4)
histogram(corr_arr(V4msk),10,'FaceAlpha',0.4)
histogram(corr_arr(ITmsk),15,'FaceAlpha',0.4)
vline([ccPthrsh],'r-.',{"P 0.01"})
xlabel("Correlation of Response")
legend(["V1","V4","IT"]) 
title(compose("Correlation of Movie-Image Firing Rate Rsp %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movImgCorrHist_%s_areas_%s.png",Animal,wdwstr)))
savefig(2,fullfile(figdir,compose("movImgCorrHist_%s_areas_%s.fig",Animal,wdwstr)))

figure(2);clf;hold on
histogram(corr_arr(Tab.F_P<0.01),10,'FaceAlpha',0.4)
histogram(corr_arr,15,'FaceAlpha',0.4)
vline([ccPthrsh],'r-.',{"P 0.01"})
xlabel("Correlation of Response")
legend(["selective", "all"]) 
title(compose("Correlation of Movie-Image Firing Rate Rsp %s Selective Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movImgCorrHist_%s_sel_%s.png",Animal,wdwstr)))
savefig(2,fullfile(figdir,compose("movImgCorrHist_%s_sel_%s.fig",Animal,wdwstr)))

figure(2);clf;hold on
histogram(corr_arr,15,'FaceAlpha',0.4)
vline([ccPthrsh],'r-.',{"P 0.01"})
xlabel("Correlation of Response")
legend(["all"]) 
title(compose("Correlation of Movie-Image Firing Rate Rsp %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(2,fullfile(figdir,compose("movImgCorrHist_%s_%s.png",Animal,wdwstr)))
savefig(2,fullfile(figdir,compose("movImgCorrHist_%s_%s.fig",Animal,wdwstr)))
%% Movie-Static T comparison on Population Level
T_arr = arrayfun(@(S)S.T, MovImgCorrStats); 
T_P_arr = arrayfun(@(S)S.T_P, MovImgCorrStats); 
T_Pthrsh = max(abs(T_arr(T_P_arr>0.01)));
figure(3);clf;hold on
histogram(T_arr(V1msk),20,'FaceAlpha',0.4)
histogram(T_arr(V4msk),20,'FaceAlpha',0.6)
histogram(T_arr(ITmsk),25,'FaceAlpha',0.4)
vline([T_Pthrsh, -T_Pthrsh],'r-.',{"P 0.01"})
xlabel("T statistics: Response to Static > Movie")
legend(["V1","V4","IT"]) 
title(compose("Paired T comparison Movie-Image Firing Rate Rsp\n Movie %s Image %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Stats.meta_mv.ephysFN, Stats.meta_im.ephysFN, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(3,fullfile(figdir,compose("movImgTstatsHist_%s_areas_%s.png",Animal,wdwstr)))
savefig(3,fullfile(figdir,compose("movImgTstatsHist_%s_areas_%s.fig",Animal,wdwstr)))

figure(3);clf; hold on
histogram(T_arr,25,'FaceAlpha',0.4)
histogram(T_arr(Tab.F_P<0.01),25,'FaceAlpha',0.5)
vline([T_Pthrsh, -T_Pthrsh],'r-.',{"P 0.01"})
xlabel("T statistics: Response to Static > Movie")
legend(["all","sel"]) 
title(compose("Paired T comparison Movie-Image Firing Rate Rsp\n Movie %s Image %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
            Stats.meta_mv.ephysFN, Stats.meta_im.ephysFN, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
saveas(3,fullfile(figdir,compose("movImgTstatsHist_%s_%s.png",Animal,wdwstr)))
savefig(3,fullfile(figdir,compose("movImgTstatsHist_%s_%s.fig",Animal,wdwstr)))
end

function h = plot_channel_corr(Stats, MovImgCorrStats, ImgrspDelayWdw, MvrspDelayWdw, errorbarflag)
if nargin == 4, errorbarflag = true; end
h = figure;set(h,'pos',[1000         462         560         520])
for chid = 1:numel(MovImgCorrStats)
	set(0, 'CurrentFigure', h);hold off;
    iCh = MovImgCorrStats(chid).iCh;
    iU = MovImgCorrStats(chid).iU;
    unitstr = char(64+iU); if iU==0, unitstr='U'; end %Stats.unit_str_im(chid);
    static_rspmat_key = MovImgCorrStats(chid).rspmat_static;
    movie_rspmat_key = MovImgCorrStats(chid).rspmat_movie;
    if errorbarflag
    static_rspsem_key = MovImgCorrStats(chid).rspsem_static;
    movie_rspsem_key = MovImgCorrStats(chid).rspsem_movie;
    xerr = static_rspsem_key(:);
    yerr = movie_rspsem_key(:);
    errorbar(static_rspmat_key(:), movie_rspmat_key(:), yerr, yerr, ... % note yneg, ypos error first.
    													xerr, xerr,'o');  % then xneg, xpos error . 
    box off
    else
    scatter(static_rspmat_key(:), movie_rspmat_key(:));
    end
    axis equal; addDiagonal();
    title(compose("Correlation of Movie-Image Firing Rate Rsp\n Movie: %s Image: %s\n Chan%d Unit%s Img~Mov Corr %.3f(%.1e)\n Img - Mov paired T %.2f(%.1e)\nANOVA F: Img %.3f(%.1e) Mov %.3f(%.1e)\nDelay Window: Img %d, %d Mov %d, %d",...
            Stats.meta_mv.ephysFN, Stats.meta_im.ephysFN, iCh, unitstr, MovImgCorrStats(chid).corr, MovImgCorrStats(chid).corr_P, ...
            MovImgCorrStats(chid).T, MovImgCorrStats(chid).T_P, ...
            MovImgCorrStats(chid).F_im, MovImgCorrStats(chid).F_P_im, MovImgCorrStats(chid).F_mv, MovImgCorrStats(chid).F_P_mv, ...
            ImgrspDelayWdw(1), ImgrspDelayWdw(end), MvrspDelayWdw(1), MvrspDelayWdw(end)))
    xlabel("Static Firing Rate");ylabel("Movie Firing Rate")
    if errorbarflag
    saveas(h,fullfile(Stats.figdir,compose("movie_image_corr_%s_%d%s_err.png",Stats.Animal,iCh,unitstr))) % char(64+iU) Stats.unit_str_im(chid) 
    else
    saveas(h,fullfile(Stats.figdir,compose("movie_image_corr_%s_%d%s.png",Stats.Animal,iCh,unitstr))) % char(64+iU) Stats.unit_str_im(chid)
    end
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

% %% Visualize the Correlation Statistic on Population Level
% corr_arr = arrayfun(@(S)S.corr, MovImgCorrStats);
% corr_P_arr = arrayfun(@(S)S.corr_P, MovImgCorrStats);
% Pthresh1 = max(corr_arr(corr_P_arr > 0.01));
% figure(2);clf;hold on
% histogram(corr_arr(V1msk & Tab.F_P<0.01),10,'FaceAlpha',0.4)
% histogram(corr_arr(V4msk & Tab.F_P<0.01),10,'FaceAlpha',0.4)
% histogram(corr_arr(ITmsk & Tab.F_P<0.01),15,'FaceAlpha',0.4)
% xlim([-0.3,1])
% vline([Pthresh1],'r-.',{"P 0.01"})
% xlabel("Correlation of Response")
% legend(["V1","V4","IT"]) 
% title(compose("Correlation of Movie-Image Firing Rate Rsp %s Selective Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
%             Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
% saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas_sel.png",Animal)))
% savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas_sel.fig",Animal)))
% 
% figure(2);clf;hold on
% histogram(corr_arr(V1msk),10,'FaceAlpha',0.4)
% histogram(corr_arr(V4msk),10,'FaceAlpha',0.4)
% histogram(corr_arr(ITmsk),15,'FaceAlpha',0.4)
% vline([Pthresh1],'r-.',{"P 0.01"})
% xlabel("Correlation of Response")
% legend(["V1","V4","IT"]) 
% title(compose("Correlation of Movie-Image Firing Rate Rsp %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
%             Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
% saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas.png",Animal)))
% savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s_areas.fig",Animal)))
% 
% figure(2);clf;hold on
% histogram(corr_arr(Tab.F_P<0.01),10,'FaceAlpha',0.4)
% histogram(corr_arr,15,'FaceAlpha',0.4)
% vline([Pthresh1],'r-.',{"P 0.01"})
% xlabel("Correlation of Response")
% legend(["selective", "all"]) 
% title(compose("Correlation of Movie-Image Firing Rate Rsp %s Selective Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
%             Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
% saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s_sel.png",Animal)))
% savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s_sel.fig",Animal)))
% 
% figure(2);clf;hold on
% histogram(corr_arr,15,'FaceAlpha',0.4)
% vline([Pthresh1],'r-.',{"P 0.01"})
% xlabel("Correlation of Response")
% legend(["all"]) 
% title(compose("Correlation of Movie-Image Firing Rate Rsp %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
%             Animal, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
% saveas(2,fullfile(figdir,compose("movieImageCorrHist_%s.png",Animal)))
% savefig(2,fullfile(figdir,compose("movieImageCorrHist_%s.fig",Animal)))
% %% Visualize the Movie-Static T comparison on Population Level
% T_arr = arrayfun(@(S)S.T, MovImgCorrStats); 
% T_P_arr = arrayfun(@(S)S.T_P, MovImgCorrStats); 
% 
% Pthresh = max(abs(T_arr(T_P_arr>0.01)));
% figure(3);clf;hold on
% histogram(T_arr(V1msk),20,'FaceAlpha',0.4)
% histogram(T_arr(V4msk),20,'FaceAlpha',0.6)
% histogram(T_arr(ITmsk),25,'FaceAlpha',0.4)
% vline([Pthresh, -Pthresh1],'r-.',{"P 0.01"})
% xlabel("T statistics: Response to Static > Movie")
% legend(["V1","V4","IT"]) 
% title(compose("Paired T comparison Movie-Image Firing Rate Rsp\n Movie %s Image %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
%             Stats.meta_mv.ephysFN, Stats.meta_im.ephysFN, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
% saveas(3,fullfile(figdir,compose("movieImageTstatsHist_%s_areas.png",Animal)))
% savefig(3,fullfile(figdir,compose("movieImageTstatsHist_%s_areas.fig",Animal)))
% 
% figure(3);clf; hold on
% histogram(T_arr,25,'FaceAlpha',0.4)
% histogram(T_arr(Tab.F_P<0.01),25,'FaceAlpha',0.5)
% vline([Pthresh, -Pthresh1],'r-.',{"P 0.01"})
% xlabel("T statistics: Response to Static > Movie")
% legend(["all","sel"]) 
% title(compose("Paired T comparison Movie-Image Firing Rate Rsp\n Movie %s Image %s All Channels\nMovie Delay Window %d, %d\n Image Delay Window %d, %d",...
%             Stats.meta_mv.ephysFN, Stats.meta_im.ephysFN, MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end)))
% saveas(3,fullfile(figdir,compose("movieImageTstatsHist_%s.png",Animal)))
% savefig(3,fullfile(figdir,compose("movieImageTstatsHist_%s.fig",Animal)))