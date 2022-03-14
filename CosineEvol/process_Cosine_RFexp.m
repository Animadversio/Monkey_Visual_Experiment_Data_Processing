%% Process Cosine RF Experiments 
Animal = "Alfa";Set_Path;
saveroot = "O:\Evol_Cosine";
strt_row = find(strcmp(ExpRecord.ephysFN,'Alfa-28012021-001'));
end_row  = find(strcmp(ExpRecord.ephysFN,'Alfa-16042021-005'));
%% Analyze the RF stats experiments 
rfrows = find(contains(ExpRecord.expControlFN,'rf'));
rfrows = rfrows((rfrows <= end_row) & (rfrows >= strt_row));
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(rfrows(3:end), Animal, false);
%%
RFS_col = RF_Calc_Stats_fun(meta_new(:), rasters_new(:), Trials_new(:));
%Alfa-02022021-001
%%
saveroot = "O:\Evol_Cosine";
for iRF = 1:numel(RFS_col)
    RFStat = RFS_col(iRF);
    maskS = RFStats_indiv_chan_gen_mask(RFStat);
    expdir = fullfile(saveroot, compose("%s-%s-RF",datestr(RFStat.meta.datetime,"yyyy-mm-dd"),RFStat.Animal));
    mkdir(expdir)
    save(fullfile(expdir,'RFStat.mat'),'RFStat')
    save(fullfile(expdir,'maskStat.mat'),'maskS')
end
%% Analyze the Response distribution in Selectivity Experiments
selrows = find(contains(ExpRecord.expControlFN,'sel') & contains(ExpRecord.Exp_collection,'PopulEvol'));
selrows = selrows((selrows <= end_row) & (selrows >= strt_row));
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(selrows(2:end), Animal, false);
%%
SelS_col = selectivity_Collect_Stats_fun(meta_new(1:end), rasters_new(1:end), Trials_new(1:end)); 
%%
for i = 1:numel(SelS_col)
    seldir = fullfile(saveroot,SelS_col(i).meta.fdrnm);
    SelS_col(i).meta.figdir = seldir; ReprStat = SelS_col(i);
    mkdir(seldir);save(fullfile(seldir,"ReprStat.mat"),'ReprStat')
    fprintf("Selectivity Exp stats saved to %s\n",fullfile(seldir,"ReprStat.mat"))
end
%% Visualize response distribution of all channels 
visusalize_resp_distri_allchan(SelS_col,16);
%%
ReprStats = SelS_col;
save(fullfile(saveroot,"summary","Alfa_ReprStats.mat"),'ReprStats')
%% Analyze the Response distribution in Selectivity Experiments
evolrows = find(contains(ExpRecord.expControlFN,'generate') & contains(ExpRecord.Exp_collection,'PopulEvol'));
evolrows = evolrows((evolrows <= end_row) & (evolrows >= strt_row));
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(evolrows(33:end), Animal, false);%1:32
%%
CosStats = Evol_Cosine_Collect_Stats_fun(meta_new, rasters_new, Trials_new);
%%
