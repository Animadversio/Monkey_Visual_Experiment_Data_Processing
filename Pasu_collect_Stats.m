%% Collect Response of Pasupathy images and Manifold 
% (time averaged to compress, both in baseline and response period)
%%
Animal = "Alfa";Set_Path;
%expftr = (contains(ExpRecord.expControlFN,"200319"));
expftr = (contains(ExpRecord.Exp_collection,"Manifold") &...
            contains(ExpRecord.expControlFN, "selectivity"));%&...
            %ExpRecord.Expi > 10);
rowis = find(expftr);
% Expi_col = [1,2,3,6,7,10,11,12,19,27];
% Expi_col = [1,3,4,5,8,9,10,11,12,13,15,16,17,18,19,20,21,22];
% assert(all(ExpRecord.Expi(rowis(Expi_col))==Expi_col')) % assert you are getting what you want. 
[meta_new,rasters_new,~,Trials_new] = loadExperiments(rowis,Animal);
%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
PasuStats = repmat(struct(), 1, length(meta_new));
%load(compose("D:\\%s_Manif_PasuStats.mat", Animal), 'PasuStats')
load(fullfile(mat_dir, compose("%s_Manif_PasuStats.mat", Animal)), 'PasuStats')
%%
for Triali = 1:length(meta_new)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN) & ExpRecord.Exp_collection=="Manifold");
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);Expi=Expi(end); % hack this for beto exp 35
if isnan(Expi) || ~contains(ExpRecord.expControlFN{exp_rowi},'selectivity') ...
        || ~contains(ExpRecord.Exp_collection{exp_rowi},'Manifold')
    % add this filter to process a sequence of Trials_new 
    keyboard
    continue
end
fprintf("Processing  Exp %d:\n",Expi)
fprintf([ExpRecord.comments{exp_rowi},'\n'])
% savepath = fullfile(result_dir, sprintf("%s_Exp%02d", Animal, Expi));
% mkdir(savepath)
PasuStats(Expi).Animal = Animal;
PasuStats(Expi).Expi = Expi;
PasuStats(Expi).imageName = Trials.imageName;
PasuStats(Expi).meta = meta;
%%
pref_chan = Trials.TrialRecord.User.prefChan;
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% note if the evolved channel is marked as null 0 then this can fail! 
% unit_in_pref_chan = 0;
pref_chan_id = find(meta.spikeID==pref_chan & ... % the id in the raster and lfps matrix 
                unit_num_arr > 0); % match for unit number

PasuStats(Expi).units.pref_chan = pref_chan;
PasuStats(Expi).units.unit_name_arr = unit_name_arr;
PasuStats(Expi).units.unit_num_arr = unit_num_arr;
PasuStats(Expi).units.activ_msk = activ_msk;
PasuStats(Expi).units.spikeID = meta.spikeID;
PasuStats(Expi).units.pref_chan_id = pref_chan_id;

if Animal == "Beto"
    if Expi <= 10, subsp_n = 3; else, subsp_n = 1;end
elseif Animal == "Alfa"
    subsp_n = 1;
end
subsp_templ = {'norm_%d_PC2_%d_PC3_%d', 'norm_%d_PC49_%d_PC50_%d','norm_%d_RND1_%d_RND2_%d'};
sphere_norm = infer_norm_from_imgnm(Trials.imageName);
PasuStats(Expi).manif.sphere_norm = sphere_norm;
PasuStats(Expi).manif.subsp_n = subsp_n;
for subsp_i = 1:subsp_n
name_pattern = subsp_templ{subsp_i};
[idx_grid, ~,~,~] = build_idx_grid(Trials.imageName, name_pattern, sphere_norm);
%%
psths_col = cellfun(@(idx) rasters(:, :, idx), idx_grid, 'UniformOutput', false);
subsp_act = cellfun(@(psth) mean(psth(:, 51:200, :),2), psths_col,'UniformOutput',false);
subsp_bsl = cellfun(@(psth) mean(psth(:, 1:40, :),2), psths_col,'UniformOutput',false);
% activ_map = cellfun(@(c) mean(c(1,50:200,:), [2,3]), psths_col, 'UniformOutput', true);
PasuStats(Expi).manif.idx_grid{subsp_i} = idx_grid;
PasuStats(Expi).manif.act{subsp_i} = subsp_act;
PasuStats(Expi).manif.bsl{subsp_i} = subsp_bsl;
end
%%
PasuStats(Expi).ref.didGabor = false;
PasuStats(Expi).ref.didPasu = false;
if sum(contains(Trials.imageName, "pasu")) > 186 
PasuStats(Expi).ref.didPasu = true;
[pasu_idx_grid,~,~,~] = build_Pasu_idx_grid(Trials.imageName);
pasu_psths_col = cellfun(@(idx) rasters(:, :, idx), pasu_idx_grid, ...
                'UniformOutput', false); % record Pasu response in all channels
PasuStats(Expi).ref.pasu_act = cellfun(@(psth) mean(psth(:, 51:200, :),2), pasu_psths_col,'UniformOutput',false);
PasuStats(Expi).ref.pasu_bsl = cellfun(@(psth) mean(psth(:, 1:40, :),2), pasu_psths_col,'UniformOutput',false);
end
if sum(contains(Trials.imageName, "gab")) > 12
PasuStats(Expi).ref.didGabor = true;
[gab_idx_grid,~,~,~] = build_Gabor_idx_grid(Trials.imageName);
gab_psths_col = cellfun(@(idx) rasters(:, :, idx), gab_idx_grid, ...
                'Uni', false); % record Gabor response in all channels
PasuStats(Expi).ref.gab_act = cellfun(@(psth) mean(psth(:, 51:200, :),2), gab_psths_col,'Uni',false);
PasuStats(Expi).ref.gab_bsl = cellfun(@(psth) mean(psth(:, 1:40, :),2), gab_psths_col,'Uni',false);
end
end
%%
% For Alfa From 34 to 46 the gabors have only 4 images, so very small power
%% Compute summary statistics for these. 
summary = struct('anova_F',nan,'anova_p',nan,'anova_F_bsl',nan,'anova_p_bsl',nan,'t',nan,'t_p',nan,'t_CI',nan(2,1));
for Expi = 1:numel(meta_new)
fprintf("Processing  Exp %d:\n",Expi)
PasuStats(Expi).ref.pasu_stats = repmat(summary,0,0);
PasuStats(Expi).ref.pasu_strs =  repmat("",0,0);
PasuStats(Expi).ref.gab_stats = repmat(summary,0,0);
PasuStats(Expi).ref.gab_strs =  repmat("",0,0);
for chan_j = 1:length(PasuStats(Expi).units.spikeID)
    if PasuStats(Expi).ref.didPasu
    act_cell = cellfun(@(c)squeeze(c(chan_j,:,:)),PasuStats(Expi).ref.pasu_act,'Uni',false);
    bsl_cell = cellfun(@(c)squeeze(c(chan_j,:,:)),PasuStats(Expi).ref.pasu_bsl,'Uni',false);
    [summary, stat_str] = calc_tuning_stats(act_cell, bsl_cell);
    PasuStats(Expi).ref.pasu_stats(chan_j) = summary;
    PasuStats(Expi).ref.pasu_strs(chan_j) = stat_str;
    end
    if PasuStats(Expi).ref.didGabor
    act_cell = cellfun(@(c)squeeze(c(chan_j,:,:)),PasuStats(Expi).ref.gab_act,'Uni',false);
    bsl_cell = cellfun(@(c)squeeze(c(chan_j,:,:)),PasuStats(Expi).ref.gab_bsl,'Uni',false);
    [summary, stat_str] = calc_tuning_stats(act_cell, bsl_cell);
    PasuStats(Expi).ref.gab_stats(chan_j) = summary;
    PasuStats(Expi).ref.gab_strs(chan_j) = stat_str;
    end
end
end
%%
save(compose("D:\\%s_Manif_PasuStats.mat", Animal), 'PasuStats')
save(fullfile(mat_dir, compose("%s_Manif_PasuStats.mat", Animal)), 'PasuStats')
%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
save(compose("D:\\%s_Manif_stats.mat", Animal), 'Stats')
save(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')

%%
[pasu_idx_grid,id_grid,p1_grid,p2_grid] = build_Pasu_idx_grid(Trials.imageName);
pasu_psths_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), pasu_idx_grid, ...
                'UniformOutput', false);
uniti=1;
act_mat = cellfun(@(psth) squeeze(mean(psth(uniti, 51:200, :), [2,3])), pasu_psths_col, ...
                'UniformOutput', true);
act_col = cellfun(@(psth) squeeze(mean(psth(uniti, 51:200, :), [2])), pasu_psths_col, ...
                'UniformOutput', false);
anovan()
%%
% psths_mean = cellfun(@(c) mean(c, 3), psths_col, 'UniformOutput', false);
% psths_1 = cellfun(@(c) squeeze(c(1,:)), psths_mean, 'UniformOutput', false);
% %%
% activ_map = cellfun(@(c) mean(c(1,50:200,:), [2,3]), psths_col, 'UniformOutput', true);
% figure
% imagesc(activ_map)
% axis equal tight; 
% %%
% figure
% for i = 1:size(idx_grid, 1)
% for j = 1:size(idx_grid, 2)
% subplottight(11,11,j+(i-1)*11);
% plot(psths_1{i,j})
% end
% end
%%
function [idx_grid, id_mat, p1_mat, p2_mat] = build_idx_grid(imgnms, name_pattern, sphere_norm)
ang_step = 18;
idx_grid = cell(11,11);
id_mat = nan(11, 11);
p1_mat = nan(11, 11);
p2_mat = nan(11, 11);
for i =  -5:5 % i code for PC2 
    for j =  -5:5 % j code for PC3
        cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
        img_idx = find(contains(imgnms, cur_fn));
        if j ~= 5 && j ~= -5
            id = 11 * i + j;
        elseif j == 5 
            id = 5;
        elseif j == -5
            id = -5;
        end
        id_mat(i+6, j+6) = id; 
        p1_mat(i+6, j+6) = i*ang_step;
        p2_mat(i+6, j+6) = j*ang_step;
        idx_grid{i+6, j+6} = img_idx; 
    end
end
end

function [idx_grid, id_mat, p1_mat, p2_mat] = build_Pasu_idx_grid(imgnms)
idx_grid = cell(51, 4);
id_mat = nan(51, 4);
p1_mat = nan(51, 4);
p2_mat = nan(51, 4);
for i = 1:51
    for j = 1:4
        cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
        img_idx = find(contains(imgnms, cur_fn));
        if i==1 || i==2 || i==3
            id = 4 * (i-1) + 1;
        else 
            id = 4 * (i-1) + j;
        end
        id_mat(i, j) = id;
        p1_mat(i, j) = i;
        p2_mat(i, j) = 2*j-1;
        idx_grid{i, j} = img_idx;
    end
end
end

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


function [summary, stat_str] = calc_tuning_stats(act_cell, bsl_cell) 
    idx_mat = reshape(1:numel(act_cell),size(act_cell));
    idx_mat = arrayfun(@(idx){idx},idx_mat);   
    idx_cell = cellfun(@(idx,act)repmat(idx,length(act),1),idx_mat,act_cell,'Uni',false); % create the idx cell array of same size
    act_vec = cat(1, act_cell{:});
    bsl_vec = cat(1, bsl_cell{:});
    idx_vec = cat(1, idx_cell{:});
    % Do statistics
    [p,tbl,stats] = anova1(act_vec,idx_vec,'off');
    stats.F = tbl{2,5};
    stats.p = p;
    summary.anova_F = stats.F;
    summary.anova_p = stats.p; 
    
    [p,tbl,stats_bsl] = anova1(bsl_vec,idx_vec,'off');
    stats_bsl.F = tbl{2,5};
    stats_bsl.p = p;
    summary.anova_F_bsl = tbl{2,5};
    summary.anova_p_bsl = p; 
    %
    [~,P,CI,STATS] = ttest(act_vec, bsl_vec);
    summary.t = STATS.tstat;
    summary.t_p = P;
    summary.t_CI = CI;
    % visualize
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.1f,%.1f] t=%.2f(%.2e)\n' ...
            'Evoked Modulation: All image, F=%.2f(%.2e)\n' ...
            'Baseline Modulation: All image, F=%.2f(%.2e)'],CI(1), CI(2), STATS.tstat, P, stats.F, stats.p,stats_bsl.F,stats_bsl.p);
end



function [score_mat, bsl_mat, summary, stat_str] = get_stats_from_result(name_pattern, sphere_norm)
    global  Trials rasters channel ang_step Reps
    score_mat = nan(11,11,Reps); 
    bsl_mat = nan(11,11,Reps);  % make sure the 3 dimension has larger size than repitition number! 
    cnt_mat = zeros(11,11); 
    id_mat = zeros(11,11); % record the id correspond to i,j
    for i =  -5:5 % i code for PC2 
        for j =  -5:5 % j code for PC3
            cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
            img_idx = find(contains(Trials.imageName, cur_fn));
            cnt_mat(i+6, j+6) = length(img_idx);
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 50:150, :))- mean(psths(1, 1:50, :)) );
            baseline = squeeze(mean(psths(1, 1:50, :)));
            score_mat(i+6, j+6, 1:length(img_idx)) = scores; % first idx is PC2. 
            bsl_mat(i+6, j+6, 1:length(img_idx)) = baseline;
            if j ~= 5 && j ~= -5
                id = 11 * i + j;
            elseif j == 5 
                id = 5;
            elseif j == -5
                id = -5;
            end
            id_mat(i+6, j+6) = id; 
        end
    end
    [theta_mat, phi_mat] = meshgrid(ang_step*(-5:5), ang_step*(-5:5));
    mean_fr_mat = bsl_mat + score_mat;
    id_vec_nan = reshape(repmat(id_mat, 1,1, Reps), 1, []);
    score_vec_nan = reshape(score_mat, 1, []);
    bsl_vec_nan = reshape(bsl_mat, 1, []);
    mean_fr_vec_nan = bsl_vec_nan + score_vec_nan;
    % Do statistics
    [p,tbl,stats] = anova1(score_vec_nan, id_vec_nan, 'off');
    stats.F = tbl{2,5};
    stats.p = p;
    summary.anova_F = stats.F;
    summary.anova_p = stats.p; 
    %
    [p2,tbl2,stats2] = anovan(score_vec_nan, {reshape(repmat(theta_mat, 1,1,Reps),1,[]), ...
                              reshape(repmat(phi_mat, 1,1,Reps),1,[])}, 'model', 'interaction','display' ,'off');
    stats2.p = p2;
    stats2.F = [tbl2{2:4,6}];
    summary.anova2_p = p2; % p for theta, phi and interaction 
    summary.anova2_F = [tbl2{2:4,6}]; % F for theta, phi and interaction
    %
    [~,P,CI] = ttest(mean_fr_vec_nan, bsl_vec_nan);
    summary.t_p = P;
    summary.t_CI = CI;
    % visualize
    stat_str = sprintf(['Evoked vs bsl(50ms) CI[%.f,%.f](%.3f)\n' ...
            'Modulation: All image, F=%.2f(%.3f)\n'...
            'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],CI(1), CI(2), P, stats.F, stats.p, ...
            stats2.F(1),stats2.p(1), stats2.F(2),stats2.p(2),stats2.F(3),stats2.p(3));

end