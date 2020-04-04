clearvars -except meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug ExpRecord
%%
Animal = "Alfa";Set_Path;
%expftr = (contains(ExpRecord.expControlFN,"200319"));
expftr = (contains(ExpRecord.Exp_collection,"Manifold") &...
            contains(ExpRecord.expControlFN, "selectivity"));
rowis = find(expftr);
% Expi_col = [1,2,3,6,7,10,11,12,19,27];
Expi_col = [1,3,4,5,8,9,10,11,12,13,15,16,17,18,19,20,21,22];
assert(all(ExpRecord.Expi(rowis(Expi_col))==Expi_col')) % assert you are getting what you want. 
[meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(rowis(Expi_col(11:end)),Animal);
%lfps_new
%%
global  Trials rasters channel ang_step Reps
result_dir = "E:\Manif_SUHash";
for Triali = 2:length(Expi_col)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);
assert(sum(Expi==Expi_col) > 0)
fprintf("Processing  Exp %d:\n",Expi)
fprintf([ExpRecord.comments{exp_rowi},'\n'])
savepath = fullfile(result_dir, sprintf("%s_Exp%02d", Animal, Expi));
mkdir(savepath)
pref_chan = Trials.TrialRecord.User.prefChan;
%unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% note if the evolved channel is marked as null 0 then this can fail! 
% unit_in_pref_chan = 0;
pref_chan_id = find(meta.spikeID==pref_chan & ... % the id in the raster and lfps matrix 
                unit_num_arr > 0); % match for unit number
assert(length(pref_chan_id) > 1) % We are testing for those that has a hash unit attached to the evolved unit.
%%
ang_step = 18;
Reps = 15;
sphere_norm = infer_norm_from_imgnm(Trials.imageName);
% collect the processed data here
if Animal == "Beto"
    if Expi <= 10, subsp_n = 3; else, subsp_n = 1;end
elseif Animal == "Alfa"
    subsp_n = 1;
end
subsp_templ = {'norm_%d_PC2_%d_PC3_%d', 'norm_%d_PC49_%d_PC50_%d','norm_%d_RND1_%d_RND2_%d'};
subsp_suff = {'PC23',"PC4950",'RND12'};
subsp_axis = ["PC2","PC3";"PC49","PC50";"RND1","RND2"];
for subsp_i = 1:subsp_n
name_templ = subsp_templ{subsp_i};
suffix = subsp_suff{subsp_i};
score_col = {};
bsl_col = {};
ch_label_col = {};
stats_col = {};
for i = 1:numel(pref_chan_id)
    channel = pref_chan_id(i);
    chan_label_str = sprintf("%s Exp%d Channel %s", Animal, Expi, unit_name_arr{channel});
    ch_label_col{i} = chan_label_str;
    [score_mat, bsl_mat, summary, stat_str] = get_stats_from_result(name_templ, sphere_norm);
    score_col{i} = score_mat;
    bsl_col{i} = bsl_mat;
    stats_col{i} = stat_str;
end
%%
figure(2);clf;set(2,'position', [295         168        1559         793])
for uniti = 1:numel(pref_chan_id)
    ax1 = subplottight(1,numel(pref_chan_id),uniti,0.05,0.13);
    channel = pref_chan_id(uniti);
    chan_label_str = ch_label_col{uniti};
    imagesc(-90:18:90, -90:18:90, nanmean(score_col{uniti}, 3)) % sum(score_mat,3)./cnt_mat
    ylabel(subsp_axis(subsp_i,1) + " degree");xlabel(subsp_axis(subsp_i,2) + " degree")
    title([chan_label_str, compose("Tuning map on %s subspace", suffix), stats_col{uniti}])
    shading flat;axis image;colorbar
end
%%
figure(6);set(6,'position',[275         162        1083         474])
for uniti = 1:numel(pref_chan_id)
mean_act = nanmean(score_col{uniti}, 3 );
[max_score, max_idx ] = max(mean_act,[],'all','linear');
[r_idx,c_idx] = ind2sub(size(mean_act), max_idx);
ang_PC2 = ang_step * (-5 + r_idx -1);
ang_PC3 = ang_step * (-5 + c_idx -1);
[ang_dist, cos_dist] = dist_grid_on_sphere(ang_PC2, ang_PC3);
subplot(1,numel(pref_chan_id),uniti)
scatter(cos_dist(:), mean_act(:))
[b,bint,~,~,stats] = regress(mean_act(:), [cos_dist(:), ones(numel(cos_dist), 1)]);
xlabel("1 - cos dist to the maxima")
ylabel("mean activity of image")
title(sprintf("%s %s\nslope=%.2f (%.2f,%.2f)\n intercp=%.3f (%.2f,%.2f)\nR2=%.3f", ...
    ch_label_col{uniti}, suffix, b(1), bint(1,1), bint(1,2), b(2), bint(2,1), bint(2,2), stats(1)))
end
%%
figure(7);set(7,'position',[275         162        1083         474])
for uniti = 1:numel(pref_chan_id)
mean_act = nanmean(score_col{uniti}, 3 );
[max_score, max_idx ] = max(mean_act,[],'all','linear');
[r_idx,c_idx] = ind2sub(size(mean_act), max_idx);
ang_PC2 = ang_step * (-5 + r_idx -1);
ang_PC3 = ang_step * (-5 + c_idx -1);
[ang_dist, cos_dist] = dist_grid_on_sphere(ang_PC2, ang_PC3);
subplot(1,numel(pref_chan_id),uniti)
scatter(ang_dist(:), mean_act(:))
[b,bint,~,~,stats] = regress(mean_act(:), [ang_dist(:), ones(numel(cos_dist), 1)]);
xlabel("angular dist (rad) to the maxima")
ylabel("mean activity of image")
title(sprintf("%s %s\nslope=%.2f (%.2f,%.2f)\n intercp=%.3f (%.2f,%.2f)\nR2=%.3f", ...
    ch_label_col{uniti}, suffix, b(1), bint(1,1), bint(1,2), b(2), bint(2,1), bint(2,2), stats(1)))
end
%%
figure(8);set(8,'position',[275         162        1083         474])
for uniti = 1:numel(pref_chan_id)
mean_act = nanmean(score_col{uniti}, 3 );
[max_score, max_idx ] = max(mean_act,[],'all','linear');
[r_idx,c_idx] = ind2sub(size(mean_act), max_idx);
ang_PC2 = ang_step * (-5 + r_idx -1);
ang_PC3 = ang_step * (-5 + c_idx -1);
[ang_dist, cos_dist] = dist_grid_on_sphere(ang_PC2, ang_PC3);
subplot(1,numel(pref_chan_id),uniti)
scatter(ang_dist(:), mean_act(:)/max_score)
[b,bint,~,~,stats] = regress(mean_act(:)/max_score, [ang_dist(:), ones(numel(cos_dist), 1)]);
xlabel("angular dist (rad) to the maxima")
ylabel("mean activity of image (normalized by max)")
title(sprintf("%s %s\nslope=%.2f (%.2f,%.2f)\n intercp=%.3f (%.2f,%.2f)\nR2=%.3f", ...
    ch_label_col{uniti}, suffix, b(1), bint(1,1), bint(1,2), b(2), bint(2,1), bint(2,2), stats(1)))
end
%%
figure(9);set(9,'position',[275         162        1083         474])
for uniti = 1:numel(pref_chan_id)
mean_act = nanmean(score_col{uniti}, 3 );
[max_score, max_idx ] = max(mean_act,[],'all','linear');
[r_idx,c_idx] = ind2sub(size(mean_act), max_idx);
ang_PC2 = ang_step * (-5 + r_idx -1);
ang_PC3 = ang_step * (-5 + c_idx -1);
[ang_dist, cos_dist] = dist_grid_on_sphere(ang_PC2, ang_PC3);
subplot(1,numel(pref_chan_id),uniti)
scatter(cos_dist(:), mean_act(:)/max_score)
[b,bint,~,~,stats] = regress(mean_act(:)/max_score, [cos_dist(:), ones(numel(cos_dist), 1)]);
xlabel("1 - cos dist to the maxima")
ylabel("mean activity of image (normalized by max)")
title(sprintf("%s %s\nslope=%.2f (%.2f,%.2f)\n intercp=%.3f (%.2f,%.2f)\nR2=%.3f", ...
    ch_label_col{uniti}, suffix, b(1), bint(1,1), bint(1,2), b(2), bint(2,1), bint(2,2), stats(1)))
end
%%
figure(10);set(10,'position',[275         162        1083         474])
for uniti = 1:numel(pref_chan_id)
mean_act = nanmean(score_col{uniti}, 3 );
[max_score, max_idx ] = max(mean_act,[],'all','linear');
subplot(1,numel(pref_chan_id),uniti)
histogram(mean_act(:), 15) % /max_score
xlabel("mean activity of image")
ylabel("entry number")
title(sprintf("%s %s", ch_label_col{uniti}, suffix))
end
saveas(2, fullfile(savepath, compose('%s_Exp%02d_ManifRspCh%02d_%s.jpg', Animal, Expi, pref_chan, suffix)))
saveas(6, fullfile(savepath, compose('%s_Exp%02d_scor-cos-corr%02d_%s.jpg', Animal, Expi, pref_chan, suffix)))
saveas(7, fullfile(savepath, compose('%s_Exp%02d_scor-ang-corr%02d_%s.jpg', Animal, Expi, pref_chan, suffix)))
saveas(8, fullfile(savepath, compose('%s_Exp%02d_normscor-ang-corr%02d_%s.jpg', Animal, Expi, pref_chan, suffix)))
saveas(9, fullfile(savepath, compose('%s_Exp%02d_normscor-cos-corr%02d_%s.jpg', Animal, Expi, pref_chan, suffix)))
saveas(10, fullfile(savepath, compose('%s_Exp%02d_scor-hist%02d_%s.jpg', Animal, Expi, pref_chan, suffix)))
%%
end
end
% scatter(cos_dist(:), mean_act(:))
%
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