Set_Path;
% Collect Manifold Experiments specs 
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Evolution_Exp";

expftr = contains(ExpSpecTable_Aug.expControlFN,"generate") & ...
     contains(ExpSpecTable_Aug.Exp_collection, "Manifold");
%  ExpSpecTable_Aug.Expi<=5 & ExpSpecTable_Aug.Expi>=4 & ...
Expi_arr = ExpSpecTable_Aug(expftr,"Expi").Expi;
Pasu_msk = Expi_arr > 10;
pref_chan_arr = ExpSpecTable_Aug(expftr,"pref_chan").pref_chan; % indexed by Expi number
ExpSpecTable_Aug(expftr,"stim_size").stim_size
%%
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr)); 
%% Break down Experiment numbers 
ITftr = 1<=pref_chan_arr & pref_chan_arr<=32;
V4ftr = 49<=pref_chan_arr & pref_chan_arr<=64;
V1ftr = 33<=pref_chan_arr & pref_chan_arr<=48;
fprintf("Manifold experiment series breakdown:\n")
fprintf("IT exps: %d\t(%d with pasu, %d with subspaces)\n", sum(ITftr), sum(ITftr & Pasu_msk), sum(ITftr & ~Pasu_msk))
disp(pref_chan_arr(ITftr)')
fprintf("V4 exps: %d\t(%d with pasu, %d with subspaces)\n", sum(V4ftr), sum(V4ftr & Pasu_msk), sum(V4ftr & ~Pasu_msk))
disp(pref_chan_arr(V4ftr)')
fprintf("V1 exps: %d\t(%d with pasu, %d with subspaces)\n", sum(V1ftr), sum(V1ftr & Pasu_msk), sum(V1ftr & ~Pasu_msk))
disp(pref_chan_arr(V1ftr)')
%% Print exp specs from loaded mat files. 
for Triali = 1:length(meta_new)
    Expi = ExpSpecTable_Aug.Expi(find(contains(ExpSpecTable_Aug.ephysFN, meta_new{Triali}.ephysFN)));
    cfg = Trials_new{Triali}.TrialRecord.User.evoConfiguration;
    fprintf("Exp %d\t evol_chan %d (unit %d)\t position [%.1f, %.1f]\t size %.1f \n",Expi, cfg{1}, cfg{4}, cfg{2}, cfg{3})
end
