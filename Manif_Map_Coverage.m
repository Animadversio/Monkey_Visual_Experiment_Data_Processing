%% Manif_Map_Coverage
Set_Path;
mat_dir = "O:\Mat_Statistics";
tabdir = "O:\Manif_MapSmooth\popstats";
sumdir = "O:\Manif_MapSmooth\summary";
Animal="Alfa"; 
alfa_StatsTab_sum = readtable(fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));
Animal="Beto"; 
beto_StatsTab_sum = readtable(fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));
StatsTab_sum = [alfa_StatsTab_sum;beto_StatsTab_sum];
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir,Animal+"_Evol_stats.mat"),'EStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
EStats_all.(Animal) = EStats;
Stats_all.(Animal) = Stats;
MapVarStats_all.(Animal) = MapVarStats;
end
%%
Expi = 3; si = 1;
actcol = MapVarStats_all.Alfa(Expi).manif.act_col{si};
bslcol = MapVarStats_all.Alfa(Expi).manif.bsl_col{si};
bslmat = cell2mat(bslcol(:)');
bslmean = mean(bslmat,2);
bslstd = std(bslmat,1,2);
bslsem = sem(bslmat,2);
%%
actmap_all = arrayfun(@(iCh)cellfun(@(A)...
    mean(A(iCh,:),2),actcol,'uni',1),1:size(bslmat,1),'uni',0);
%%
tab = [];
for iCh = 1:size(bslmat,1)
actmap_col = cellfun(@(A) squeeze(A(iCh,:)),actcol,'uni',0);
anovaStats = anova_cells(actmap_col);
tab(Expi).F = anovaStats.F;
tab(Expi).F_P = anovaStats.F_P;
tab(Expi).F_df = anovaStats.STATS.df;
tab(Expi).FSTATS = anovaStats.STATS;
end

