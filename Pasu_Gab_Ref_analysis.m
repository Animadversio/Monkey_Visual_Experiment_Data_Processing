%% Pasu_reference Analysis
%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats')
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats')
load(fullfile(mat_dir, Animal+"_Manif_PasuStats.mat"), 'PasuStats') % pasupathy data for the whole population
%% 
ui = 1;
Stats(Expi).ref.pasu_psths