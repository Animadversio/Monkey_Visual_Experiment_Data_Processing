Animal = "Alfa"; Set_Path; 
ftr = find(contains(ExpRecord.ephysFN,"Alfa-04112020"));
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr, Animal);
%%
Trials_mov = Trials_new{5};
rasters_mov = rasters_new{5};
meta_mov = meta_new{5};
Trials_img = Trials_new{6};
rasters_img = rasters_new{6};
meta_img = meta_new{6};
% figdir will be an input argument.
figroot = "E:\OneDrive - Washington University in St. Louis\MovieDynamics";
figdir = fullfile(figroot, "2020-11-04-Alfa-Chan28-2");
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

psthmov_mean = cellfun(@(idx) mean(rasters_mov(:, :, idx),3), mov_idx_arr, 'Uni',0);
psthmov_sem = cellfun(@(idx) std(rasters_mov(:, :, idx),1,3) / sqrt(numel(idx)), mov_idx_arr, 'Uni',0);
%%
spikeID_im = meta_img.spikeID;
unitID_im = meta_img.unitID;
unit_str_im = generate_unit_labels_new(meta_img.spikeID, meta_img.unitID);
% Sort Static images and trials
[imgnm_per_mov, idx_arr, imgnm_arr, mvnms_uniq] = parse_frame_idx_arr(Trials_img.imageName);
[img_mvnm_arr, img_frnum_arr, uniq_imgnms] = parse_frame_name(Trials_img.imageName);
nImgPerMv = max(cellfun(@numel, imgnm_per_mov)); 

psthimg_mean = cellfun(@(idx)mean(rasters_img(:,:,idx),3), idx_arr,'Uni',0);
psthimg_sem = cellfun(@(idx)std(rasters_img(:,:,idx),1,3)/sqrt(numel(idx)), idx_arr,'Uni',0);
% imgnm_per_mov = mat2cell(imgnm_arr, ones(size(imgnm_arr,1),1), size(imgnm_arr,2)); % Cell array of the images in 
%% Matching process
% Matching all movie image pairs
[matchfrids, matchTONs, matchTOFFs] = closest_Kframe_locate(movnm_sorted(:), imgnm_per_mov, ...
                                        meta_mov.stimuli, meta_img.stimuli, false, 1); % halfSep = false,
%%
Stats.movnm = movnm_sorted;
Stats.imgnm_arr = imgnm_arr;
Stats.img_idx_arr = idx_arr;
Stats.imgnm_per_mov = imgnm_per_mov;
Stats.nImgPerMv = nImgPerMv;

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
[MovImgCorrStats, corr_arr, corr_P_arr, corr_sep_arr, corr_sep_P_arr] = MovImgMatchCorr(ImgrspDelayWdw, MvrspDelayWdw, psthimg_mean, psthmov_mean, Stats);
%%
Stats.imgTuneStats = {};
for chid = 1:numel(spikeID_im)
    Stats.imgTuneStats{chid} = calc_tune_stats(cellfun(@(idx) rasters_img(chid, :, idx), idx_arr, 'Uni', 0));
end
Stats.imgTuneStats = cell2mat(Stats.imgTuneStats);
%%
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


function [imgnm_per_movie, idx_arr, imgnm_arr, mvnms_uniq] = parse_frame_idx_arr(imageName)
uniq_nms = string(unique(imageName));
matchtoks = regexp(uniq_nms,"(.*)_(\d*)",'tokens');
mvnm_arr = cellfun(@(tok)tok{1}{1}, matchtoks, 'Uni', 0); 
frnum_arr = cellfun(@(tok)str2double(tok{1}{2}), matchtoks, 'Uni', 1); 
mvnms_uniq = unique(mvnm_arr);
imgnm_per_movie = cellfun(@(nm) uniq_nms(contains(mvnm_arr, nm))', mvnms_uniq, 'uni',0);
ncol = max(cellfun(@numel, imgnm_per_movie));
imgnm_arr = cell(numel(mvnms_uniq), ncol); % nMovie by nFrames
idx_arr = cell(numel(mvnms_uniq), ncol);
for iMv = 1:numel(mvnms_uniq)
    for iMat = 1:numel(mvnms_uniq{iMv})
        imgnm = imgnm_per_movie{iMv}(iMat);
        imgnm_arr{iMv, iMat} = imgnm;
        idx_arr{iMv, iMat} = find(contains(imageName, imgnm));
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
rspmat_key = rspmat_static(:,:);
rspmat_movie_key = mean(rspmat_movie(:,:,:), 3);
[cc,pp] = corr(rspmat_key(:),rspmat_movie_key(:));
MovImgCorrStats(chid).corr = cc;
MovImgCorrStats(chid).corr_P = pp;

rspmat_movie_vecs = reshape(rspmat_movie(:,:,:),[], size(rspmat_movie,3));
[cc,pp] = corr(rspmat_key(:), rspmat_movie_vecs);
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

fprintf("Response Delay window movie [%d,%d] image [%d,%d] ms\n",MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end))
fprintf("%d / %d channels has significant(0.01) correlation between Movie and Image Response\n",sum(corr_P_arr<0.01),numel(S.spikeID_im))
fprintf("%d / %d channels has significant(0.001) correlation\n",sum(corr_P_arr<0.001),numel(S.spikeID_im))
fprintf("%d / %d IT channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & S.spikeID_im'<=32),sum(S.spikeID_im<=32))
fprintf("%d / %d V1 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & S.spikeID_im'>=33 & S.spikeID_im'<=48),sum(S.spikeID_im'>=33 & S.spikeID_im'<=48))
fprintf("%d / %d V4 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & S.spikeID_im'>=49),sum(S.spikeID_im'>=49))
fprintf("%s ",S.unit_str_im(corr_P_arr<0.01))
fprintf("\n")
end