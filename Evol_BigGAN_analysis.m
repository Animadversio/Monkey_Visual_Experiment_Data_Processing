%% Evol BigGAN Analysis
Animal = "Beto";Set_Path;
expftr = contains(ExpRecord.expControlFN,"200728");% & contains(ExpRecord.Exp_collection,"BigGAN");
fllist = find(expftr);
[meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(fllist(1:end),Animal,false);%find(expftr)%find(expftr)
%%
