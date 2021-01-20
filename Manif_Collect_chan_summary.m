mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa"; 
load(fullfile(mat_dir, Animal+'_Evol_stats.mat')) 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat')) 
Manif_dir = "E:\OneDrive - Washington University in St. Louis\PC_space_tuning";
%%
SummaryStats = repmat(struct(),1,numel(EStats)); 
for Expi = 1:45
pref_chan = EStats(Expi).evol.pref_chan;
if Animal == "Alfa"
    expoutdir = compose("Alfa_Exp%d_chan%02d", Expi, pref_chan);
elseif Animal == "Beto"
    expoutdir = compose("Exp%d_chan%02d", Expi, pref_chan);
end
%%
load(fullfile(Manif_dir,expoutdir,'Basic_Stats.mat'),'Stat_summary');
Stat_arr = cell2mat(Stat_summary); 
fprintf("%s Exp %d Stats shape (%d,%d)\n",Animal,Expi,size(Stat_arr))
if Stats(Expi).manif.subsp_n == 3
SummaryStats(Expi).manif = Stat_arr;
SummaryStats(Expi).ref = [];
elseif Stats(Expi).manif.subsp_n == 3
SummaryStats(Expi).manif = Stat_arr(:,1);
SummaryStats(Expi).ref = Stat_arr(:,2);
end
end
%% 
save(fullfile(mat_dir, Animal+'_Manif_SummaryStats.mat'),'SummaryStats')
%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat')) 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat')) 
%%
wdw_vect = [[1, 10]+[0:5:190]'; [1, 20] + 10*[0:18]'; [1, 50] + [0:50:150]'; [51,200]];
PrecSummaryStats = repmat(struct(),1,numel(Stats)); 
%%
tic
for Triali = 1:numel(meta_new) %  10 %
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN) & ExpRecord.Exp_collection=="Manifold");
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);Expi=Expi(end); % hack this for beto exp 35
if isnan(Expi) || ~all(contains(ExpRecord.expControlFN(exp_rowi),'selectivity')) ...
        || ~all(contains(ExpRecord.Exp_collection(exp_rowi),'Manifold'))
    % add this filter to process a sequence of Trials_new 
    keyboard
    continue
end
fprintf("Processing  Exp %d:\n",Expi)
fprintf([ExpRecord.comments{exp_rowi},'\n'])
nCh = size(rasters, 1);
for si = 1:Stats(Expi).manif.subsp_n
psth_all_arr = cellfun(@(idx) rasters(:,:,idx), Stats(Expi).manif.idx_grid{si},'Uni',false);
PrecSummaryStats(Expi).manif{si} = arrayfun(@(iCh) calc_tune_stats( ...
    cellfun(@(psth) psth(iCh,:,:),psth_all_arr,'Uni',0), wdw_vect ), 1:nCh);
end
PrecSummaryStats(Expi).gabor = [];
PrecSummaryStats(Expi).pasu  = [];
if Stats(Expi).ref.didPasu
psth_all_arr = cellfun(@(idx) rasters(:,:,idx), Stats(Expi).ref.pasu_idx_grid, 'Uni',false);
PrecSummaryStats(Expi).pasu = arrayfun(@(iCh) calc_tune_stats( ...
    cellfun(@(psth) psth(iCh,:,:),psth_all_arr,'Uni',0), wdw_vect ), 1:nCh);
end
if Stats(Expi).ref.didGabor
psth_all_arr = cellfun(@(idx) rasters(:,:,idx), Stats(Expi).ref.gab_idx_grid, 'Uni',false);
PrecSummaryStats(Expi).gabor = arrayfun(@(iCh) calc_tune_stats( ...
    cellfun(@(psth) psth(iCh,:,:),psth_all_arr,'Uni',0), wdw_vect ), 1:nCh);
end
toc
end
%%
save(fullfile(mat_dir, Animal+'_Manif_PrecSummaryStats.mat'), 'PrecSummaryStats', 'wdw_vect')
%%
% Data = load(fullfile(mat_dir, "Beto"+'_Manif_PrecSummaryStats.mat'),'PrecSummaryStats')
%%
load(fullfile(mat_dir, Animal+'_Manif_PrecSummaryStats.mat'))