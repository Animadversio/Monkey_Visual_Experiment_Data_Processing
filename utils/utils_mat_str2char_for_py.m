MatStats_path = "O:\Mat_Statistics";
Animal = "Beto";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
for Expi = 1:numel(EStats)
    EStats(Expi).ref.impaths = EStats(Expi).ref.impaths';
    EStats(Expi).ref.impaths_chr = cellstr(EStats(Expi).ref.impaths);
end
%%
save(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')