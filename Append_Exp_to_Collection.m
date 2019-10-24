%% 
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw();
% 
load("D:\\Manifold_Exps.mat")
%%
lfps{end+1}=lfps_new{1};
meta{end+1}=meta_new{1};
Trials{end+1}=Trials_new{1};
rasters{end+1}=rasters_new{1}; 
savefast("D:\\Manifold_Exps.mat", 'lfps', 'meta', 'rasters', 'Trials')