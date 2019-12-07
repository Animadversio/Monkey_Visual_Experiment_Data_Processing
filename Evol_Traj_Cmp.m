%% Evolution Compare 
% Comparing Evolution Experiments of different size 
% Final Score and Evolution time course. 

expid = find(ExpSpecTable_Aug.Expi > 27 & contains(ExpSpecTable_Aug.expControlFN,'generate'));
expid(1) = 67; expid(2) = 64;
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(expid);
%% 
global  Trials rasters channel sphere_norm ang_step Reps

Triali = 1;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};

exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi_tab = ExpSpecTable_Aug.Expi(exp_rowi);
assert(Expi_tab == Expi, "Data Expi doesn't match that in exp record! Check your indexing or record.")
Expi = Expi_tab; 
fprintf("Processing Exp %d, %s\n", Expi, meta.comments)
%%
