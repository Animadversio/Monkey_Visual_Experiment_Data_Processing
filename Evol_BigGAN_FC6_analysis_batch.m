Animal = "Alfa";
Set_Path;
ftrrows = find(contains(ExpRecord.ephysFN,["Alfa-25112020-003"]));%"Alfa-27102020-003", "Alfa-27102020-004"
% find(contains(ExpRecord.expControlFN,"generate_") & (ExpRecord.Exp_collection=="BigGAN_fc6" |...
%                ExpRecord.Exp_collection=="BigGAN_FC6"));
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows, Animal, false);
BFEStats = Evol_BigGAN_FC6_Collect_Stats_fun(meta_new, rasters_new, Trials_new);
Evol_BigGAN_FC6_Animation_fun(BFEStats)