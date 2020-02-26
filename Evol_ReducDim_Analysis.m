%%
Set_Path;
expftr = contains(ExpSpecTable_Aug.expControlFN,"200226");
% expftr = contains(ExpSpecTable_Aug.expControlFN,"generate") & ...
%      contains(ExpSpecTable_Aug.Exp_collection, "Optimizer_cmp");
%  ExpSpecTable_Aug.Expi<=5 & ExpSpecTable_Aug.Expi>=4 & ...
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr)); 
%%
