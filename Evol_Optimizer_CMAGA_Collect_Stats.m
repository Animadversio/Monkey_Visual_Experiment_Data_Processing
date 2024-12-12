%% Process experiment with CMA vs GA optimizers. 
%  loading => Collect Stats => Generate summary figure 
!ExpRecordBackup.bat 
%%
Animal = "Both";Set_Path;
expftr = contains(ExpRecord.expControlFN,"generate") & ...
        contains(ExpRecord.Exp_collection, "CMAGA_cmp");
row_idx = find(expftr);
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(row_idx, Animal, false);
%%
