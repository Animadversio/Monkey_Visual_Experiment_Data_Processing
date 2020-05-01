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

%%
t_table = zeros(64,45);
F_table = zeros(64,45);
for Expi = 11:45
t_col = arrayfun(@(c)c.t,PasuStats(Expi).ref.pasu_stats);
F_col = arrayfun(@(c)c.anova_F,PasuStats(Expi).ref.pasu_stats);
spikeID = PasuStats(Expi).units.spikeID;
fprintf("\nExp%d\n",Expi)
t_table(spikeID, Expi) = t_col;
F_table(spikeID, Expi) = F_col;
end
%%
figure
subplot(121)
imagesc(t_table)
xlabel("Exp num")
colorbar()
subplot(122)
imagesc(F_table)
xlabel("Exp num")
colorbar()