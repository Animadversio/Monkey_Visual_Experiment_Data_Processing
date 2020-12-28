
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
assert(exist(figdir,'dir')>0, "Exp folder not exist, no where to load the files.")
% mkdir(figdir)
%%  
D = load(fullfile(figdir, "ExpStats.mat"));%, 'Stats', 'MovImgCorrStats','ImgrspDelayWdw','MvrspDelayWdw')
S = D.Stats;
psthmov_all = cellfun(@(idx) rasters_mov(:,:,idx), S.mov_idx_arr,'Uni',0);
psthimg_all = cellfun(@(idx) rasters_img(:,:,idx), S.img_idx_arr,'Uni',0);
%% Is there a best time window (length, delay) for each channel to be correlated? 
% loop
MovDelays = -120:10:120;nMvDel = numel(MovDelays);
ImgDelays = 40:10:100;nImDel = numel(ImgDelays);
WdwLens = 10:10:150;nWdw = numel(WdwLens);
% kerWids = [5, 11, 21, 31, 51, 81, 121];
CorrStat_col = cell(numel(MovDelays),numel(WdwLens),numel(ImgDelays)); 
wdw_col = cell(numel(MovDelays),numel(WdwLens),numel(ImgDelays));
% workers = parpool(6);
tic
parfor i = 1:nMvDel%numel(MovDelays)
for j = 1:nWdw%numel(WdwLens)
for k = 1:nImDel%numel(ImgDelays)
ImgrspDelayWdw = ImgDelays(k):ImgDelays(k)+WdwLens(j);
MvrspDelayWdw = MovDelays(i):MovDelays(i)+WdwLens(j);
if ImgrspDelayWdw(end) > 200 || ImgrspDelayWdw(1) < 1, continue; end
[MovImgCorrStats] = MovImgCorrVariability_fast(ImgrspDelayWdw, MvrspDelayWdw, psthimg_all, psthmov_all, Stats); % channels
CorrStat_col{i,j,k} = MovImgCorrStats;
wdw_col{i,j,k} = [ImgrspDelayWdw(1),ImgrspDelayWdw(end);
                  MvrspDelayWdw(1),MvrspDelayWdw(end)];
toc
end
end
end
%% Use the template to fill the blanks for CorrStat 
% tmpl = MovImgCorrStats;
tmpl = arrayfun(@(M)struct('corr', nan, 'corr_P', nan, 'T', nan, 'T_P', nan,'iCh',M.iCh,'iU',M.iU),MovImgCorrStats);
CorrStat_col_fill = CorrStat_col;
CorrStat_col_fill(cellfun(@isempty, CorrStat_col)) = {tmpl};
CorrStat_tsr = cell2mat(cellfun(@(C)reshape(C,1,1,1,[]),CorrStat_col_fill(:,:,:),'Uni',0));
cc_tsr = arrayfun(@(C)C.corr,CorrStat_tsr);
%%
imdelay = 60;
k = find(ImgDelays == imdelay);
% cell2mat(CorrStat_col(:,1:end-1,k))
CorrStat_arr = cell2mat(cellfun(@(C)reshape(C,1,1,[]),CorrStat_col(:,1:end-1,k),'Uni',0));
cc_tsr = arrayfun(@(C)C.corr,CorrStat_arr);
%%
save(fullfile(figdir, "corr_wdw_optim.mat"),'CorrStat_col','wdw_col','S')
%%
h=figure;set(h,'pos',[16         353        2538         508])
T=tiledlayout(1,nImDel,'pad','compact','TileSpacing','compact');
for iCh = 1:size(cc_tsr,4)
    for k = 1:nImDel
    nexttile(k)%1,nImDel,
    ccmat = cc_tsr(:,:,k,iCh);
    imagesc(ccmat);hold on 
    [cc_max, idx_max] = max(ccmat,[],'all','linear');
    [ri,cj] = ind2sub(size(ccmat), idx_max);
    % MovDelays(ri),WdwLens(cj)
    bestImWdw = wdw_col{ri,cj,k}(1,:);
    bestMvWdw = wdw_col{ri,cj,k}(2,:);
    title(compose("%s Max cc %.3f\n Window onset delay from image onset %d\nBest window Im:[%d,%d]ms Mv:[%d,%d]ms",...
        S.unit_str_im(iCh),cc_max,imdelay,bestImWdw(1),bestImWdw(2),bestMvWdw(1),bestMvWdw(2)));
    xlabel("Window Length");xticks(1:nWdw);xticklabels(WdwLens)
    ylabel("Window onset delay from movie frame onset");yticks(1:nMvDel);yticklabels(MovDelays)
    axis image
    colorbar()
    hold off
    end
    saveas(h,fullfile(figdir,compose("RateCorr_optim_wdw_Ch%s.png",S.unit_str_im(iCh))))
end
%%
bestCorrWdw = repmat(struct(),1,numel(MovImgCorrStats));
for iCh = 1:size(cc_tsr,4)
    [cc_max, idx_max] = max(cc_tsr(:,:,:,iCh),[],'all','linear');
    [ri,cj,lk] = ind2sub(size(cc_tsr), idx_max);
    bestCorrWdw(iCh).iCh = S.spikeID_im(iCh);
    bestCorrWdw(iCh).iU = S.unitID_im(iCh);
    bestCorrWdw(iCh).cc_max = cc_max;
    bestCorrWdw(iCh).ImWdw = wdw_col{ri,cj,lk}(1,:);
    bestCorrWdw(iCh).MvWdw = wdw_col{ri,cj,lk}(2,:);
    bestCorrWdw(iCh).WdwL = wdw_col{ri,cj,lk}(1,2) - wdw_col{ri,cj,lk}(1,1);
end
CTab = struct2table(bestCorrWdw);
%%
figure;
msk = (CTab.iU > 0) & (CTab.cc_max > 0.5);
scatter(CTab.iCh(msk), CTab.MvWdw(msk,1))
ylabel("Onset of Movie Window");xlabel("Channel")

%%
cmap = jet;CLIM = [0.3,1];
h2=figure;hold on
msk = (CTab.iU > 0) & (CTab.cc_max > 0.5);
% plot([CTab.iCh(msk),CTab.iCh(msk)]', [CTab.MvWdw(msk,:)]', 'color',[0.5,0.5,0.5])
cvalue = (CTab.cc_max - CLIM(1)) / (CLIM(2)-CLIM(1));
cseq = cmap(max(1,round(cvalue * 256)),:);
scatter(CTab.iCh(msk), CTab.MvWdw(msk,1), 36*CTab.iU(msk).^2)
for ri = find(msk)'
    plot([CTab.iCh(ri),CTab.iCh(ri)]',[CTab.MvWdw(ri,:)]', 'color',[cseq(ri,:),0.6],'LineWidth',2.5)
end
vline([32.5,48.5])
ylabel("Onset of Movie Window");xlabel("Channel");hold off
title("Optimal Window Length vs Channel Number")
saveas(h2,fullfile(figdir,"Optim_Window_sumall.png"))
function [MovImgCorrS] = MovImgCorrVariability_fast(ImgrspDelayWdw, MvrspDelayWdw, psthimg_all, psthmov_all, S, chanels)
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
% 
% Output parameter: 
%   MovImgCorrS: An abreviated version of MovImgCorrStats, just to compare
%       the correlation and discard much information.
if nargin<6, chanels = 1:numel(S.spikeID_im); end
MovImgCorrS = repmat(struct(),1,numel(S.spikeID_im)); % collect results here. 
for chid = chanels%1:numel(S.spikeID_im) % loop through static image units (Loop through movie is also fine). 
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
% % ANOVA of individual trials
% anova_mov = anova_cells(tr_rspmat_movie);
% anova_img = anova_cells(tr_rspmat_static);
% MovImgCorrStats(chid).F_im = anova_img.F;
% MovImgCorrStats(chid).F_P_im = anova_img.F_P;
% MovImgCorrStats(chid).F_mv = anova_mov.F;
% MovImgCorrStats(chid).F_P_mv = anova_mov.F_P;
% Collapse the trial dimension to get std and sem of each group
% rspmat_movie is # Movie - by - # Static in the Movie - by - # Occurence.
% e.g. (12, 5, 2) in Hessian image case.
rspmat_movie = cellfun(@mean, tr_rspmat_movie);
% rspstd_movie = cellfun(@std, tr_rspmat_movie);
% rspsem_movie = cellfun(@(rsps) std(rsps)/sqrt(numel(rsps)), tr_rspmat_movie);
rspmat_static = cellfun(@mean, tr_rspmat_static);
% rspstd_static = cellfun(@std, tr_rspmat_static);
% rspsem_static = cellfun(@(rsps) std(rsps)/sqrt(numel(rsps)), tr_rspmat_static);
% psth_movie = cellfun(@(psth)mean(psth, [3]), tr_psth_movie, 'Uni', 0);
% psth_static = cellfun(@(psth)mean(psth, [3]), tr_psth_static, 'Uni', 0);
% Correlate the non-center image responses
rspmat_static_key = rspmat_static(:,:);
rspmat_movie_key = mean(rspmat_movie(:,:,:), 3);
[cc,pp] = corr(rspmat_static_key(:),rspmat_movie_key(:), 'rows', 'complete');
MovImgCorrS(chid).corr = cc;
MovImgCorrS(chid).corr_P = pp;
% Paired ttest of image - frame response. 2 sided.
[~,P,~,STAT] = ttest(rspmat_static_key(:),rspmat_movie_key(:)); 
MovImgCorrS(chid).T = STAT.tstat; % t > 0 show the response to image > movie
MovImgCorrS(chid).T_P = P;
% % Correlation of image response to response of each occurance of the image
% rspmat_movie_vecs = reshape(rspmat_movie(:,:,:),[], size(rspmat_movie,3));
% [cc,pp] = corr(rspmat_static_key(:), rspmat_movie_vecs, 'rows', 'complete');
% % Save some of response info for plotting.
% MovImgCorrStats(chid).corr_sep = cc;
% MovImgCorrStats(chid).corr_sep_P = pp;
% MovImgCorrStats(chid).rspmat_movie = rspmat_movie; 
% MovImgCorrStats(chid).rspstd_movie = rspstd_movie; 
% MovImgCorrStats(chid).rspsem_movie = rspsem_movie; 
% MovImgCorrStats(chid).rspmat_static = rspmat_static; 
% MovImgCorrStats(chid).rspstd_static = rspstd_static; 
% MovImgCorrStats(chid).rspsem_static = rspsem_static; 
% MovImgCorrStats(chid).psth_movie = psth_movie; 
% MovImgCorrStats(chid).psth_static = psth_static; 
MovImgCorrS(chid).iCh = iCh;
MovImgCorrS(chid).iU = iU;
end
% % summarize channel stats into population
% corr_arr = arrayfun(@(S)S.corr, MovImgCorrStats); 
% corr_P_arr = arrayfun(@(S)S.corr_P, MovImgCorrStats); 
% corr_sep_arr = cell2mat(arrayfun(@(S)S.corr_sep, MovImgCorrStats','Uni',0)); % each column is an occurence
% corr_sep_P_arr = cell2mat(arrayfun(@(S)S.corr_sep_P, MovImgCorrStats','Uni',0));
% T_arr = arrayfun(@(S)S.T, MovImgCorrStats); 
% T_P_arr = arrayfun(@(S)S.T_P, MovImgCorrStats); 
% F_im_arr = arrayfun(@(S)S.F_im, MovImgCorrStats); 
% F_P_im_arr = arrayfun(@(S)S.F_P_im, MovImgCorrStats); 
% F_mv_arr = arrayfun(@(S)S.F_mv, MovImgCorrStats); 
% F_P_mv_arr = arrayfun(@(S)S.F_P_mv, MovImgCorrStats); 
% % Print summary string
% ITmsk = S.spikeID_im'<=32; V1msk = S.spikeID_im'>=33 & S.spikeID_im'<=48; V4msk = S.spikeID_im'>=49;
% fprintf("Response Delay window movie [%d,%d] image [%d,%d] ms\n",MvrspDelayWdw(1), MvrspDelayWdw(end), ImgrspDelayWdw(1), ImgrspDelayWdw(end))
% fprintf("%d / %d channels has significant(0.01) correlation between Movie and Image Response\n",sum(corr_P_arr<0.01),numel(S.spikeID_im))
% fprintf("%d / %d channels has significant(0.001) correlation\n",sum(corr_P_arr<0.001),numel(S.spikeID_im))
% fprintf("%d / %d IT channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & ITmsk),sum(ITmsk))
% fprintf("%d / %d V1 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & V1msk),sum(V1msk))
% fprintf("%d / %d V4 channels has significant (0.01) correlation\n",sum(corr_P_arr<0.01 & V4msk),sum(V4msk))
% fprintf("%s ",S.unit_str_im(corr_P_arr<0.01)) % list of channels that are strongly correlated 
% fprintf("\n")
% fprintf("%d / %d channels has significant (0.01) difference between firing rate of paired Movie and Image Response\n",sum(T_P_arr<0.01),numel(S.spikeID_im))
% fprintf("In IT, %d / %d Static > Movie, %d / %d Movie > Static  significantly\n",sum(T_P_arr<0.01 & T_arr >0 & ITmsk),sum(ITmsk),sum(T_P_arr<0.01 & T_arr <0 & ITmsk),sum(ITmsk))
% fprintf("In V1, %d / %d Static > Movie, %d / %d Movie > Static  significantly\n",sum(T_P_arr<0.01 & T_arr >0 & V1msk),sum(V1msk),sum(T_P_arr<0.01 & T_arr <0 & V1msk),sum(V1msk))
% fprintf("In V4, %d / %d Static > Movie, %d / %d Movie > Static  significantly\n",sum(T_P_arr<0.01 & T_arr >0 & V4msk),sum(V4msk),sum(T_P_arr<0.01 & T_arr <0 & V4msk),sum(V4msk))
% fprintf("\n")
% fprintf("%d / %d channels are significantly (0.01) modulated by this set of static images\n",sum(F_P_im_arr<0.01),numel(S.spikeID_im))
% fprintf("IT: %d / %d ,V1: %d / %d ,V4: %d / %d \n", sum(F_P_im_arr<0.01 & ITmsk),sum(ITmsk),...
% 		   sum(F_P_im_arr<0.01 & V1msk),sum(V1msk), sum(F_P_im_arr<0.01 & V4msk),sum(V4msk))
% fprintf("%d / %d channels are significantly (0.01) modulated by this set of images embedded in movies\n",sum(F_P_mv_arr<0.01),numel(S.spikeID_im))
% fprintf("IT: %d / %d ,V1: %d / %d ,V4: %d / %d \n", sum(F_P_mv_arr<0.01 & ITmsk),sum(ITmsk),...
% 		   sum(F_P_mv_arr<0.01 & V1msk),sum(V1msk), sum(F_P_mv_arr<0.01 & V4msk),sum(V4msk))
% wdwstr = compose("im_%d_%d_mv_%d_%d",ImgrspDelayWdw(1), ImgrspDelayWdw(end), MvrspDelayWdw(1), MvrspDelayWdw(end));
end