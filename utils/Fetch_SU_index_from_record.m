%%  Re-assign meta comments to files. 
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
Animal = "Beto";Set_Path;
% Put the manifold experiments' comments in the Stats 
expftr = (contains(ExpRecord.Exp_collection,"Manifold") &...
            contains(ExpRecord.expControlFN, "selectivity")&...
            ExpRecord.Expi > 0);
rowis = find(expftr);
for Expi = 1:length(Stats)
    assert(Expi == ExpRecord.Expi(rowis(Expi)))
    Stats(Expi).meta.comments = ExpRecord.comments{rowis(Expi)};
end
% Put the evolution experiments' comments in the EStats.
expftr = (contains(ExpRecord.Exp_collection,"Manifold") &...
            contains(ExpRecord.expControlFN, "generate")&...
            ExpRecord.Expi > 0);
rowis = find(expftr);
for Expi = 1:length(EStats)
    assert(Expi == ExpRecord.Expi(rowis(Expi)))
    EStats(Expi).meta.comments = ExpRecord.comments{rowis(Expi)};
end
%% Channel Quality by record in the comments of exp notes
SU_index = zeros(length(Stats),1);
for Expi = 1:length(Stats)
    fprintf("======Exp %d ======\n",Expi)
    disp(EStats(Expi).meta.comments)
    disp(Stats(Expi).meta.comments)
    Q = input("Channel Quality number\n");
    SU_index(Expi) = Q;
end
%% Try to fill in the missing value from the original notebook
fprintf("The following Exps are missing the SU index\n")
missing_idx = find(isnan(SU_index));
arrayfun(@(Expi)string(Stats(Expi).meta.expControlFN),missing_idx)
%%
for Expi = 1:length(Stats)
    Stats(Expi).meta.SUidx = SU_index(Expi);
    EStats(Expi).meta.SUidx = SU_index(Expi);
end
%%
save(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
save(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')