
Animal = "Caos"; Set_Path;
currows = find(contains(ExpRecord.expControlFN,["231106"])); 
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(currows(1:end), Animal, false);
bhvfns = ExpRecord.expControlFN(currows);
% saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_Cosine";
saveroot = "E:\OneDrive - Harvard University\Evol_Cosine";
%%
global Animal
%% Extract and Visualize Cosine Experiments. 
evoidx = contains(bhvfns,"reconstruction_cosine") & ~cellfun(@isempty,meta_new');
CStats = Evol_Cosine_Collect_Stats_fun(meta_new(evoidx), rasters_new(evoidx), Trials_new(evoidx));