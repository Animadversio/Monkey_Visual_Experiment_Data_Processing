%% Movie Selectivity Analysis, Trial by Trial
Animal="Alfa";Set_Path;setMatlabTitle("MovieStaticDynamics")
ftr = find(contains(ExpRecord.ephysFN,"21102020"));
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr(4:7),Animal);

%% Evolved Image
rasters_evo = rasters_new{1};
Trials_evo = Trials_new{1};
meta_evo = meta_new{1};
%%
lastBGidx = find(contains(Trials_evo.imageName,"block031_thread001_gen")); 
%% Static Image
rasters = rasters_new{4};
Trials = Trials_new{4};
meta = meta_new{4};
%% Find the indices for the target static image. 
% centidx = find((contains(Trials.imageName,"0.00") | contains(Trials.imageName,"0.16")) &...
%                 ~contains(Trials.imageName,"eig0"));
centidx = find(contains(Trials.imageName,"0.00"));
iCh = find(meta.spikeID==7);
meanpsth = mean(rasters(iCh, :, centidx),3);
figure;clf;hold on
plot(1:200, meanpsth);
plot(1:200, squeeze(rasters(iCh, :, centidx)), 'Color', [0,0,0,0.1], 'LineWidth',0.2);
%% Movie Part
Trials_mov = Trials_new{3};
rasters_mov = rasters_new{3};
meta_mov = meta_new{3};
unit_name_arr = generate_unit_labels(meta_mov.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta_mov.spikeID, rasters);

vid = VideoReader(fullfile(meta_mov.stimuli, Trials_mov.imageName{1}+".avi"));
vidLenth = vid.Duration * 1000;
fps = vid.FrameRate; 
frameN = int32(viewTime / (1000/fps) + 1);
frameTick = (1000/fps)*double(0:frameN); 
%%
movnm = string(unique(Trials_mov.imageName));
[movnm_sorted, sortedProp, sortIdx] = sortMovieNames(movnm);
idx_arr = arrayfun(@(mv)find(contains(Trials_mov.imageName, mv)),movnm_sorted,'Uni',0);
%%
figdir = "E:\OneDrive - Washington University in St. Louis\MovieDynamics\2020-10-21-Alfa-Chan07-2";
for channum = 49:64 % 1:64
iChs = find(meta_mov.spikeID==channum);
for iMv = 1:numel(movnm_sorted)
figure(4); clf; 
ax = subplot('Position', [0.065, 0.09, 0.9, 0.84]);hold on; xlabel("Time (ms)"); ylabel("Firing Rate")
for iCh = iChs'
meanpsth_mov = mean(rasters_mov(iCh, :, idx_arr{iMv}),3);
sempsth_mov = std(rasters_mov(iCh, :, idx_arr{iMv}),1,3) / sqrt(numel(idx_arr{iMv}));
% meanpsth_evo = mean(rasters_evo(iCh, :, lastBGidx),3);
% sempsth_evo = std(rasters_evo(iCh, :, lastBGidx),1,3) / sqrt(numel(lastBGidx));
meanpsth = mean(rasters(iCh, :, centidx),3);
sempsth = std(rasters(iCh, :, centidx),1,3) / sqrt(numel(centidx));
%%
plot([-250+1:2500],meanpsth_mov, 'r')
plot([-250+1:2500],squeeze(rasters_mov(iCh, :, idx_arr{iMv})), 'Color', [1,0,0,0.1], 'LineWidth',0.2);
plot(1:200, meanpsth, 'k');
plot(1:200, squeeze(rasters(iCh, :, centidx)), 'Color', [0,0,0,0.1], 'LineWidth',0.2);
% plot(1:200, meanpsth_evo, 'blue');
% plot(1:200, squeeze(rasters_evo(iCh, :, lastBGidx)), 'Color', [0,0,1,0.1], 'LineWidth',0.2);
end
xlim([-50,350])
title(compose("%s\n%s corr=%.3f",movnm(iMv),"Chan"+num2str(channum),corr(meanpsth_mov(251:450)', meanpsth')),'Interpreter','none')
vline([0,100],'-.k',"Static Image")
vline([0:100/3:350],':r',"Movie")
saveas(4,fullfile(figdir,compose("%s_%d_%s_static_cmp.png", Animal, channum, movnm(iMv))))
% savefig(4,fullfile(figdir,compose("%s_%d_%s_static_cmp.png", Animal, channum, movnm(iMv))))
% pause;
end
end
%% 
preT = meta.rasterWindow(1) - meta_mov.rasterWindow(1); % the raster window of movie experiment extends before that of image exp. 
movieImgCorr = nan(numel(meta.spikeID),numel(movnm_sorted));
for iCh = 1:numel(meta.spikeID)
    for iMv = 1:numel(movnm_sorted)
    meanpsth_mov = mean(rasters_mov(iCh, :, idx_arr{iMv}),3);
    meanpsth = mean(rasters(iCh, :, centidx),3);
    movieImgCorr(iCh, iMv) = corr(meanpsth_mov(preT+[1:200])', meanpsth');
    end
end
corrtab = array2table(movieImgCorr,'VariableNames',movnm_sorted,'RowNames',unit_name_arr);
writetable(corrtab, fullfile(figdir, "MovieImageRspCorr.csv"), 'WriteRowNames', true)
% save(fullfile(figdir, "MovieImageRspCorr.mat"), 'movieImgCorr')
%%
function [sortedMovnm, sortedRows, sortIdx] = sortMovieNames(movnm)
% Sort the movie names in lexicoGraphical order, by space and eigen idx. 
% This assumes the names to have the structure like "class_eig17_shortshort"
eigi_cell = cellfun(@(eigi){str2double(eigi{1})},regexp(movnm,"_eig(\d*)_",'tokens'));
space_cell = cellfun(@(eigi){eigi{1}{1}},regexp(movnm,"(.*)_eig(\d*)_",'tokens'));
[sortedRows, sortIdx] = sortrows([space_cell,eigi_cell]);
sortedMovnm = movnm(sortIdx);
end