%%
Animal = "Alfa"; Set_Path;
ftr = find(contains(ExpRecord.ephysFN,"Alfa-27102020"));
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr([1,3,5,6,7,8]),Animal);
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
figroot = "E:\OneDrive - Washington University in St. Louis\MovieDynamics";
figdir = fullfile(figroot, "2020-10-27-Alfa-Chan09-1");
mkdir(figdir)
%%
figure(1);
iChs = 41:76;
for iMv = 1:numel(mov_idx_arr)
imagesc(wdw(1)+1:wdw(2), iChs, mean(rasters_mov(iChs, :, mov_idx_arr{iMv}),3))
title(movnm_sorted(iMv),'Interp','none')
xlim([-50,2000])
pause
end
%% Combine the 3 adjacent experiments into one. 
% 'Alfa-27102020-007', 'Alfa-27102020-008', 'Alfa-27102020-009'
imageName_cmb = cat(1, Trials_new{4}.imageName, Trials_new{5}.imageName, Trials_new{6}.imageName);
rasters_cmb = cat(3, rasters_new{4:6});
stimparts = split(meta_new{4}.stimuli,'\');
spikeID = meta_new{4}.spikeID;
unit_name_arr = generate_unit_labels(spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, spikeID, rasters_cmb);
%% Predict 
uniq_imgnm = unique(imageName_cmb);
% cellfun(@(eigi){str2double(eigi{1})},regexp(uniq_imgnm,"(.*)_eig(\d*)_lin([\d.]*)",'tokens'));
[idx_arr_nos, imgnm_arr_nos, idx_arr_cls, imgnm_arr_cls, ...
    eig_id_arr_nos, dist_arr_nos, eig_id_arr_cls, dist_arr_cls] = parse_image_idx_arr_hess(imageName_cmb);
%%
idx_arr = [idx_arr_cls; idx_arr_nos];
imgnm_arr = [imgnm_arr_cls; imgnm_arr_nos];
psth_mean = cellfun(@(idx)mean(rasters_cmb(:,:,idx),3), idx_arr,'Uni',0);
psth_sem = cellfun(@(idx)std(rasters_cmb(:,:,idx),1,3)/sqrt(numel(idx)), idx_arr,'Uni',0);
%% Set up the figure layout for faster update of static image rsp later. 
h=figure(2);clf
nrow = size(idx_arr,1); ncol = size(idx_arr,2);
T=tiledlayout(nrow,ncol,'TileSpacing','compact','Padding','compact');
for r = 1:nrow
    for c = 1:ncol
    ax_col{r,c} = nexttile(c+(r-1)*ncol);
    lin_col{r,c} = plot(1:200);
    pat_col{r,c} = patch([1:200,fliplr(1:200)], [1:200,fliplr(1:200)], 'k','FaceAlpha',0.15,'EdgeColor','none');
    box off; xticklabels([]); yticklabels([])
    end
end
xticklabels(ax_col{nrow,1}, 'auto') % Add the ticks to somewhere. 
yticklabels(ax_col{nrow,1}, 'auto')
xlabel(T, "Distance from center")
ylabel(T, "Eigen Axes")
TIT = title(T, "");
% Set up label for first col and last row
for c = 1:ncol
    ci = strfind(imgnm_arr{nrow,c},'lin');
    xlabel(ax_col{nrow,c}, imgnm_arr{nrow,c}(ci:end),'Interp','none')
end
for r = 1:nrow
    ci = strfind(imgnm_arr{r,1},'lin');
    ylabel(ax_col{r,1}, strrep(imgnm_arr{r,1}(1:ci-2),'_',' '),'Interp','none')
end
%%
figdir = "E:\OneDrive - Washington University in St. Louis\MovieDynamics\2020-10-27-Alfa-Chan09-1";
smth = 3;
for iCh = 1:76
MaxAmp = max(cellfun(@(psth,sem)max(psth(iCh,:)+sem(iCh)),psth_mean,psth_sem),[],'all');
for r = 1:nrow
    for c = 1:ncol
    psth = movmean(psth_mean{r,c}(iCh,:),smth); err = movmean(psth_sem{r,c}(iCh,:),smth);
    lin_col{r,c}.YData = psth;
    pat_col{r,c}.YData = [psth - err, fliplr(psth + err)];
    ylim(ax_col{r,c}, [0,MaxAmp])
    end
end
chan_str = unit_name_arr(iCh);
TIT.String = compose("%s %s Ch %s",stimparts{end},Animal,chan_str);
saveas(h,fullfile(figdir,compose("%s_Ch%s_frame_psth_avg.png",Animal,chan_str)))
end
%% 
h2=figure(3);clf
ncol = numel(movnm_sorted);
T2=tiledlayout(ncol, 1,'TileSpacing','compact','Padding','compact');

for iCh = 1:76
MaxAmp = max(cellfun(@(idx)max(mean(rasters_mov(iCh, :, idx),3) + ...
        std(rasters_mov(iCh, :, idx),1,3) / sqrt(numel(idx))), mov_idx_arr),[],'all');
for iMv = 1:ncol
    nexttile(iMv)
    psthmov = mean(rasters_mov(iCh, :, mov_idx_arr{iMv}),3);
    errmov = std(rasters_mov(iCh, :, mov_idx_arr{iMv}),1,3) / sqrt(numel(mov_idx_arr{iMv}));
    plot(wdw(1)+1:wdw(2), psthmov)
    patch([wdw(1)+1:wdw(2),fliplr(wdw(1)+1:wdw(2))], [psthmov-errmov,fliplr(psthmov+errmov)], 'k','FaceAlpha',0.15,'EdgeColor','none');
    ylabel(movnm_sorted{iMv}(1:end-6),'Interp','none')
    xlim([-100,2000]);ylim([0,MaxAmp]);vline([0.0,viewTime],'-.r')
    if iMv ~= ncol, box off; xticklabels([]); end
end
chan_str = unit_name_arr(iCh);
title(T2, compose("%s %s Ch %s",stimparts{end},Animal,chan_str));
saveas(h2,fullfile(figdir,compose("%s_Ch%s_movie_psth_avg.png",Animal,chan_str)))
end
%%
mov_act_mean = cell2mat(arrayfun(@(iMv)mean(rasters_mov(:, 450:1800, mov_idx_arr{iMv}),[2,3]),[1:12],'Uni',0));
mov_bsl_mean = cell2mat(arrayfun(@(iMv)mean(rasters_mov(:, [1:300,2000:2300], mov_idx_arr{iMv}),[2,3]),[1:12],'Uni',0));
%%
img_act_mean = cell2mat(cellfun(@(psth)reshape(mean(psth(:,50:200),2),1,1,[]),psth_mean,'Uni',0));
img_bsl_mean = cell2mat(cellfun(@(psth)reshape(mean(psth(:,1:45),2),1,1,[]),psth_mean,'Uni',0));
img_act_mean_half = squeeze(mean(img_act_mean(:,6:9,:),2))';
%%
% img_act_mean_half(spikeID==11,:)
% mov_act_mean(meta_mov.spikeID==11,:)
chid = 9;
[cval, pval] = corr(img_act_mean_half(spikeID==chid,:)', mov_act_mean(meta_mov.spikeID==chid,:)','Type','Spearman')
%%
figure(5);
subplot(121)
scatter(img_act_mean_half(10,:), mov_act_mean(11,:))
[cval, pval] = corr(img_act_mean_half(10,:)', mov_act_mean(11,:)','Type','Spearman');
xlabel("Image (Key Frames) Mean Resp");ylabel("Movie Mean Resp");title(compose("Alfa 9A corr %.3f (p=%.1e)",cval,pval))
subplot(122)
scatter(img_act_mean_half(11,:), mov_act_mean(12,:))
[cval, pval] = corr(img_act_mean_half(11,:)', mov_act_mean(12,:)','Type','Spearman');
xlabel("Image (Key Frames) Mean Resp");ylabel("Movie Mean Resp");title(compose("Alfa 9B corr %.3f (p=%.1e)",cval,pval))
suptitle("Correlation of Average Response to Movie and Image")
savefig(5, fullfile(figdir, "Movie_Img_Tuning_Corr_Chan9.fig"))
saveas(5, fullfile(figdir, "Movie_Img_Tuning_Corr_Chan9.png"))
%%
figure(6);
subplot(121)
scatter(img_act_mean_half(10,:), mov_act_mean(12,:))
[cval, pval] = corr(img_act_mean_half(10,:)', mov_act_mean(12,:)','Type','Spearman');
xlabel("Image (Key Frames) Mean Resp");ylabel("Movie Mean Resp");title(compose("Alfa 9A(Img) ~ Alfa 9B(Mov) corr %.3f (p=%.1e)",cval,pval))
subplot(122)
scatter(img_act_mean_half(11,:), mov_act_mean(11,:))
[cval, pval] = corr(img_act_mean_half(11,:)', mov_act_mean(11,:)','Type','Spearman');
xlabel("Image (Key Frames) Mean Resp");ylabel("Movie Mean Resp");title(compose("Alfa 9B(Img) ~ Alfa 9A(Mov) corr %.3f (p=%.1e)",cval,pval))
suptitle("Control Correlation of Average Response to Movie and Image")
savefig(6, fullfile(figdir, "Movie_Img_Tuning_Corr_Chan9_crosCtrl.fig"))
saveas(6, fullfile(figdir, "Movie_Img_Tuning_Corr_Chan9_crosCtrl.png"))
%%
mov_act_mean_extr = cell2mat(cellfun(@(idx)mean(rasters_mov(:, 1183:1483, idx),[2,3]),mov_idx_arr','Uni',0));
img_act_mean_extr = cell2mat(cellfun(@(psth)reshape(mean(psth(:,50:200),2),1,[]),psth_mean(:,end),'Uni',0))';
%%
chid = 9;
[cval, pval] = corr(img_act_mean_extr(spikeID==chid,:)', mov_act_mean_extr(meta_mov.spikeID==chid,:)','Type','Spearman')
%%

function [sortedMovnm, sortedRows, sortIdx] = sortMovieNames(movnm)
% Sort the movie names in **lexicoGraphical order**, by space and by eigen idx. 
% This assumes the names to have the structure like "class_eig17_shortshort"
eigi_cell = cellfun(@(eigi){str2double(eigi{1})},regexp(movnm,"_eig(\d*)_",'tokens'));
space_cell = cellfun(@(eigi){eigi{1}{1}},regexp(movnm,"(.*)_eig(\d*)_",'tokens'));
[sortedRows, sortIdx] = sortrows([space_cell,eigi_cell]);
sortedMovnm = movnm(sortIdx);
end