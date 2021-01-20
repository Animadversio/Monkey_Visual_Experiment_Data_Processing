load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
%%
for Triali = 1:length(meta_new)
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
pref_chan = Trials.TrialRecord.User.prefChan;
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% note if the evolved channel is marked as null 0 then this can fail! 
pref_chan_id = find(meta.spikeID==pref_chan & ... % the id in the raster and lfps matrix 
                unit_num_arr > 0); % match for unit number

if sum(contains(Trials.imageName, "gab")) > 12
Stats(Expi).ref.didGabor = true;
[gab_idx_grid,~,~,~] = build_Gabor_idx_grid(Trials.imageName);
gab_psths_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), gab_idx_grid, ...
                'Uni', 0);
keyboard
Stats(Expi).ref.gab_psths = gab_psths_col;
Stats(Expi).ref.gab_idx_grid = gab_idx_grid;
end
end
%%
save(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
save(fullfile("D:\",Animal+"_Manif_stats.mat"),'Stats')

function [idx_grid, id_mat, p1_mat, p2_mat] = build_Gabor_idx_grid(imgnms)
ori = 0:30:150; sf = [0.5, 1];
idx_grid = cell(2,6);
id_mat = nan(2, 6);
p1_mat = nan(2, 6);
p2_mat = nan(2, 6);
for i = 1:numel(sf) 
    for j =  1:numel(ori) 
        cur_fn = sprintf('gab_ori_%.1f_%.1f', ori(j), sf(i));
        img_idx = find(contains(imgnms, cur_fn));
        id = 6 * (i - 1) + j;
        id_mat(i, j) = id; 
        idx_grid{i, j} = img_idx; 
    end
end
end