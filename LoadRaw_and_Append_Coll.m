function [meta_new,rasters_new,lfps_new,Trials_new] = LoadRaw_and_Append_Coll(rowi, Animal, meta_new,rasters_new,lfps_new,Trials_new)
% Keep the content of meta_new and append new trials to it. 
% Functional version of Append_Exp_to_Collection.m script
% Obsolete @ 2022
[meta_new2,rasters_new2,lfps_new2,Trials_new2] = Project_Manifold_Beto_loadRaw(rowi, Animal);

meta_new = [meta_new, meta_new2];
rasters_new = [rasters_new, rasters_new2];
lfps_new = [lfps_new, lfps_new2];
Trials_new = [Trials_new, Trials_new2];
clear meta_new2 rasters_new2 lfps_new2 Trials_new2

end