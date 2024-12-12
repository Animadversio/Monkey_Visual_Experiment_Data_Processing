%% Compare Evolved image in BigGAN and FC6 in batch
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
for Animal = ["Both"] %, "Beto""Beto"
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
end
%%
S_col = {};
for Expi = 1:numel(BFEStats)
if isempty(BFEStats(Expi).evol)
S_part = struct("Expi",Expi,"space_names1","","space_names2","","optim_names1","","optim_names2","");
else
S_part = struct("Expi",Expi,"space_names1",BFEStats(Expi).evol.space_names(1),...
                            "space_names2",BFEStats(Expi).evol.space_names(2),...
                            "optim_names1",BFEStats(Expi).evol.optim_names(1),...
                            "optim_names2",BFEStats(Expi).evol.optim_names(2));
end
S_col{end+1} = S_part;
end
%%
tab = struct2table(cat(1,S_col{:}));
writetable(tab, "E:\OneDrive - Harvard University\Manuscript_BigGAN\Stats_tables\meta_optimizer_info.csv")
