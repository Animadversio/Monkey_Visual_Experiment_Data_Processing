%% 0
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw();
%%
load("D:\\Manifold_Evolv_Exps.mat")
%
lfps{end+1}=lfps_new{2};
meta{end+1}=meta_new{2};
Trials{end+1}=Trials_new{2};
rasters{end+1}=rasters_new{2}; 
savefast("D:\\Manifold_Evolv_Exps.mat", 'lfps', 'meta', 'rasters', 'Trials')
%%
load("D:\\Manifold_Exps.mat")
%
lfps{end+1}=lfps_new{1};
meta{end+1}=meta_new{1};
Trials{end+1}=Trials_new{1};
rasters{end+1}=rasters_new{1}; 

savefast("D:\\Manifold_Exps.mat", 'lfps', 'meta', 'rasters', 'Trials')
%%
storedStruct.meta{4}.stimuli = '\\storage1.ris.wustl.edu\crponce\ActiveStimuli\2019-Manifold\beto-191008a\backup_10_08_2019_12_14_29\PC_imgs';
storedStruct.meta{5}.stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-09a-beto' ;
storedStruct.meta{6}.stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191010a\backup_10_10_2019_13_10_15\PC_imgs' ;

%%
load("D:\\Manifold_Exps.mat")
meta{4}.stimuli = '\\storage1.ris.wustl.edu\crponce\ActiveStimuli\2019-Manifold\beto-191008a\backup_10_08_2019_12_14_29\PC_imgs';
meta{5}.stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-09a-beto' ;
meta{6}.stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191010a\backup_10_10_2019_13_10_15\PC_imgs' ;

for Expi = 1: numel(meta_new)/2
lfps{end+1}=lfps_new{2 * Expi - 1};
meta{end+1}=meta_new{2 * Expi - 1};
Trials{end+1}=Trials_new{2 * Expi - 1};
rasters{end+1}=rasters_new{2 * Expi - 1}; 
disp(length(lfps))
end
savefast("D:\\Manifold_Exps.mat", 'lfps', 'meta', 'rasters', 'Trials')
%%
load("D:\\Manifold_Evolv_Exps.mat")
for Expi = 1: numel(meta_new)/2
lfps{end+1}=lfps_new{2 * Expi};
meta{end+1}=meta_new{2 * Expi};
Trials{end+1}=Trials_new{2 * Expi};
rasters{end+1}=rasters_new{2 * Expi}; 
end
savefast("D:\\Manifold_Evolv_Exps.mat", 'lfps', 'meta', 'rasters', 'Trials')
%%
storedStruct.meta = meta;
storedStruct.Trials = Trials;
storedStruct.rasters = rasters;
storedStruct.lfps = lfps;