%% Pasu Manif Plot
%%
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir, Animal+"_Manif_PasuStats.mat"))
%%
for Expi = 11:45
t_p_col = arrayfun(@(c)c.t_p,PasuStats(Expi).ref.pasu_stats);
F_p_col = arrayfun(@(c)c.anova_p,PasuStats(Expi).ref.pasu_stats);
unit_name_arr = PasuStats(Expi).units.unit_name_arr;
fprintf("\nExp%d\n",Expi)
fprintf("%s ", unit_name_arr(t_p_col<0.001 & F_p_col<0.001)')
end