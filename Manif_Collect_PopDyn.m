%% Collect the population neural dynamics on manifold stats into a struct array.
%  it keeps every unit response psth to every manifold image. only mean
%  over trials. So the resulting file will be large (300Mb)
ManifDyn = repmat(struct(),1,46);
for Expi = 1:46
%load(fullfile("S:\Data-Ephys-MAT",Stats(Expi).meta.ephysFN+"_formatted.mat"))
rasters = rasters_new{Expi};
%
fprintf("%d\n",Expi)
si = 1;
% manif_psth_all = cellfun(@(idx)rasters(:, :, idx), Stats(Expi).manif.idx_grid{si}, "UniformOutput", false);
manif_psth_avg = cellfun(@(idx)mean(rasters(:, :, idx),3), Stats(Expi).manif.idx_grid{si}, "UniformOutput", false);
manif_psth_std = cellfun(@(idx)std(rasters(:, :, idx),0,3), Stats(Expi).manif.idx_grid{si}, "UniformOutput", false);
manif_psth_avg_tsr = cell2mat(reshape(manif_psth_avg,1,1,11,11)); % reshape and transform, so size 86, 200, 11, 11, key data! 
manif_psth_std_tsr = cell2mat(reshape(manif_psth_std,1,1,11,11)); % reshape and transform, so size 86, 200, 11, 11, key data! 
ManifDyn(Expi).psth_tsr = manif_psth_avg_tsr;
ManifDyn(Expi).psth_std_tsr = manif_psth_std_tsr;
ManifDyn(Expi).rep_num = cellfun(@(idx)length(idx), Stats(Expi).manif.idx_grid{si});
ManifDyn(Expi).mean_rep = mean(ManifDyn(Expi).rep_num,"all");
end
MatStats_path = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Mat_Statistics";
save(fullfile(MatStats_path, compose("%s_ManifPopDynamics.mat", Animal)), 'ManifDyn')
