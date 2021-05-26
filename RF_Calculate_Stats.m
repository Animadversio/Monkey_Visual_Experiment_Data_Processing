Animal = "Alfa"; Set_Path;
expftr = ...%(contains(ExpRecord.Exp_collection,"Manifold") &...
            contains(ExpRecord.expControlFN, "rfMapper");%&...ExpRecord.Expi > 0); 
            %
rowis = find(expftr); 
rowis = rowis(rowis>=find(strcmp(ExpRecord.ephysFN,"Alfa-05012021-004")));
%%
loadExperiments(rowis,Animal,true,true); % sort the spikes
%% 
expftr = contains(ExpRecord.expControlFN, "rfMapper");%&...ExpRecord.Expi > 0); 
rowis = find(expftr); 
rowis = rowis(rowis>=find(strcmp(ExpRecord.ephysFN,"Alfa-05012021-004")));
[meta_new,rasters_new,lfps_new,Trials_new] = loadExperiments(rowis,Animal,false,true);%Project_Manifold_Beto_loadRaw(rowis,Animal,false,true);
%%
"Alfa-05012021-004"