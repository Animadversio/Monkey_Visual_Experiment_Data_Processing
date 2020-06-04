%% Comprehensive View
addpath utils
addpath DNN
addpath VisGUI
%% Prepare tools
global G 
G = FC6Generator("matlabGANfc6.mat");
%%
pe = pyenv('Version','C:\Users\binxu\.conda\envs\caffe36\python.exe'); % Note the python env could not be changed in a matlab session
py.importlib.import_module('numpy');
%% Animal specific data.
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
% MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
load(fullfile(MatStats_path, compose("%s_Manif_RFstats.mat", Animal)), 'RFStats')
%% 
Expi = 33;
ExpType = "Evol";
option = struct();
corrFeatTsr_Anim_fun(EStats,ExpType,Animal, Expi)
Evol_Stats_Anim_fun(Stats, EStats,  Animal, Expi)
pause;
%%
ExpType = "Manif";
corrFeatTsr_Anim_fun(EStats,ExpType,Animal,Expi,struct("save",true,"sum_method","L1"))
Manif_Stats_Anim_fun(Stats, EStats, Animal,Expi)
% pause;
%%
for Expi=15:length(EStats)
[x, estimRF] = RF_View_fun(Stats, EStats, RFStats, Animal, Expi);
pause
end

