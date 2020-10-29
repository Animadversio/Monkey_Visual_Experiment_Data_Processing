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
%% Predict 
uniq_imgnm = unique(imageName_cmb);
% cellfun(@(eigi){str2double(eigi{1})},regexp(uniq_imgnm,"(.*)_eig(\d*)_lin([\d.]*)",'tokens'));
[idx_arr_nos, imgnm_arr_nos, idx_arr_cls, imgnm_arr_cls] = parse_image_idx_arr_hess(imageName_cmb);
%%
idx_arr = [idx_arr_cls; idx_arr_nos];
imgnm_arr = [imgnm_arr_cls; imgnm_arr_nos];
psth_mean = cellfun(@(idx)mean(rasters_cmb(:,:,idx),3), idx_arr,'Uni',0);
psth_sem = cellfun(@(idx)std(rasters_cmb(:,:,idx),1,3)/sqrt(numel(idx)), idx_arr,'Uni',0);
%%
figure;
tiledlayout(size(idx_arr,1),size(idx_arr,2),'TileSpacing','compact','Padding','compact');


function [sortedMovnm, sortedRows, sortIdx] = sortMovieNames(movnm)
% Sort the movie names in **lexicoGraphical order**, by space and by eigen idx. 
% This assumes the names to have the structure like "class_eig17_shortshort"
eigi_cell = cellfun(@(eigi){str2double(eigi{1})},regexp(movnm,"_eig(\d*)_",'tokens'));
space_cell = cellfun(@(eigi){eigi{1}{1}},regexp(movnm,"(.*)_eig(\d*)_",'tokens'));
[sortedRows, sortIdx] = sortrows([space_cell,eigi_cell]);
sortedMovnm = movnm(sortIdx);
end