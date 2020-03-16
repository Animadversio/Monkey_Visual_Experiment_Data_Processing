clearvars -except meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug 
%% Estimate Noise level with Manifold Experiment Data
Animal="Beto";Set_Path;
% expftr = contains(ExpSpecTable_Aug.expControlFN,"200303");
expftr = ExpSpecTable_Aug.Expi<=5 & ExpSpecTable_Aug.Expi>=1 & ...
    contains(ExpSpecTable_Aug.expControlFN,"selectivity") & ...
     contains(ExpSpecTable_Aug.Exp_collection, "Manifold");
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr),Animal); 
%%
global  Trials rasters sphere_norm ang_step Reps
ang_step = 18;
Reps = 11;
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Manifold_Noise_Estim";
for Triali = 1:5
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpSpecTable_Aug.Expi(exp_rowi);
fprintf("Processing  Exp %d:\n",Expi)
disp(ExpSpecTable_Aug.comments(exp_rowi))
pref_chan = Trials.TrialRecord.User.prefChan;
unit_name_arr = generate_unit_labels(meta.spikeID);
unit_in_pref_chan = 1;
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
pref_chan_id = find(meta.spikeID==pref_chan & ... % the id in the raster and lfps matrix 
                    unit_num_arr==unit_in_pref_chan); % match for unit number

Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
savepath = fullfile(result_dir, compose("Manifold%02d_chan%02d", Expi, pref_chan(1)));
mkdir(savepath);
% Parse out the sphere norm from image names
tmp = cellfun(@(c) regexp(c,"norm_(?<norm>\d*)_",'names'), Trials.imageName, 'UniformOutput', false);      
extnorms = cellfun(@(c) str2num(c.norm), tmp(~cellfun('isempty',tmp)));
sphere_norm = mode(extnorms);
%% Channel loop
%channel_j = pref_chan_id + 1;
% collect statistics across channel see the distribution
cc_col = nan(size(rasters,1));
slope_col = nan(size(rasters,1));
evk_cc_col = nan(size(rasters,1));
evk_slope_col = nan(size(rasters,1));
fout = fopen(fullfile(savepath, "rate_std_corr_stats.txt"),'w');
for channel_j = 1:size(rasters,1)
% for each channel collect these statistics for each image identity, and form a vector. 
var_vect_col = [];
std_vect_col = [];
mean_vect_col = [];
resid_col = [];
evk_var_vect_col = [];
evk_std_vect_col = [];
evk_mean_vect_col = [];
evk_resid_col = [];
for pattern = string({'norm_%d_PC2_%d_PC3_%d', 'norm_%d_PC49_%d_PC50_%d', 'norm_%d_RND1_%d_RND2_%d'})
% [score_mat, bsl_mat, summary, stat_str] = get_stats_from_result('norm_%d_PC2_%d_PC3_%d', pref_chan_id);
[score_mat, bsl_mat, summary, stat_str] = get_stats_from_result(pattern, channel_j);
%%
% use the score i.e. difference of rates. 
score_mean = nanmean(score_mat, 3);
score_std = nanstd(score_mat, 0, 3);
score_var = nanvar(score_mat, 0, 3);
score_mean(:, 1) = mean(score_mean(:, 1));
score_mean(:, end) = mean(score_mean(:, end));
var_vect = [nanvar(score_mat(:, 1, :), 0, "all"),  reshape(score_var(:, 2:end-1), 1, []),  nanvar(score_mat(:, end, :), 0, "all")];
std_vect = [nanstd(score_mat(:, 1, :), 0, "all"), reshape(score_std(:, 2:end-1), 1, []), nanstd(score_mat(:, end, :), 0, "all")];
mean_vect = [nanmean(score_mat(:, 1, :), "all"), reshape(score_mean(:, 2:end-1), 1, []), nanmean(score_mat(:, end, :), "all")];
res_vect = rmmissing(reshape(score_mat-score_mean,1,[]));

var_vect_col = [var_vect_col, var_vect];
std_vect_col = [std_vect_col, std_vect];
mean_vect_col = [mean_vect_col, mean_vect];
resid_col = [resid_col, res_vect];
% use the evoked rate to compute the same thing
evk_mat = score_mat + bsl_mat;
evk_mean = nanmean(evk_mat, 3);
evk_std = nanstd(evk_mat, 0, 3);
evk_var = nanvar(evk_mat, 0, 3);
evk_mean(:, 1) = mean(evk_mat(:, 1));
evk_mean(:, end) = mean(evk_mat(:, end));
evk_var_vect = [nanvar(evk_mat(:, 1, :), 0, "all"),  reshape(evk_var(:, 2:end-1), 1, []),  nanvar(evk_mat(:, end, :), 0, "all")];
evk_std_vect = [nanstd(evk_mat(:, 1, :), 0, "all"), reshape(evk_std(:, 2:end-1), 1, []), nanstd(evk_mat(:, end, :), 0, "all")];
evk_mean_vect = [nanmean(evk_mat(:, 1, :), "all"), reshape(evk_mean(:, 2:end-1), 1, []), nanmean(evk_mat(:, end, :), "all")];
evk_res_vect = rmmissing(reshape(evk_mat-evk_mean,1,[]));

evk_var_vect_col = [evk_var_vect_col, evk_var_vect];
evk_std_vect_col = [evk_std_vect_col, evk_std_vect];
evk_mean_vect_col = [evk_mean_vect_col, evk_mean_vect];
evk_resid_col = [evk_resid_col, evk_res_vect];
end
%
ccoef = correlation(mean_vect_col, std_vect_col);
[b,bint]=regress(std_vect_col',[mean_vect_col',ones(length(mean_vect_col),1)]);
ccoef2 = correlation(evk_mean_vect_col, evk_std_vect_col);
[b2,bint2]=regress(evk_std_vect_col',[evk_mean_vect_col',ones(length(evk_mean_vect_col),1)]);
% Plot the residue and scatter for each channel
% figure(13);clf;set(13,'position',[503         404        1156         500])
% subplot(121)
% scatter(mean_vect_col, std_vect_col)
% title([sprintf("corr coef %.3f",ccoef),sprintf("slope %.2f intercept %.1f", b(1), b(2))])
% xlabel("mean response")
% ylabel("response std")
% subplot(122)
% histogram(resid_col)
% xlabel("residue response")
% suptitle([sprintf("Manifold Exp%d chan%02d", Expi, pref_chan),sprintf("Unit %s, Statisitcs of Score", unit_name_arr(channel_j))])
%% Plot the residue and scatter for each channel
% figure(14);clf;set(14,'position',[503         404        1156         500])
% subplot(121)
% scatter(evk_mean_vect_col, evk_std_vect_col)
% title([sprintf("corr coef %.3f",ccoef2),sprintf("slope %.2f intercept %.1f", b2(1), b2(2))])
% xlabel("mean response")
% ylabel("response std")
% subplot(122)
% histogram(evk_resid_col)
% xlabel("residue response rate")
% suptitle([sprintf("Manifold Exp%d chan%02d", Expi, pref_chan),sprintf("Unit %s, Statisitcs of Evoked Response", unit_name_arr(channel_j))])
%%
% saveas(14, fullfile(savepath, compose("evk_noise_estim_%s.png", unit_name_arr(channel_j))))
cc_col(channel_j) = ccoef; 
slope_col(channel_j) = b(1);
evk_cc_col(channel_j) = ccoef2; 
evk_slope_col(channel_j) = b2(1);
% write readable statistics
fprintf("Chan %s\t Evoked: cc %.3f, slp %.3f(%.3f,%.3f) interc %.1f(%.1f,%.1f)  Score: cc %.3f, slp %.3f(%.3f,%.3f) interc %.1f(%.1f,%.1f)\n",unit_name_arr(channel_j),...
    ccoef2,b2(1),bint2(1,1),bint2(1,2),b2(2),bint2(2,1),bint2(2,2),...
    ccoef,b(1),bint(1,1),bint(1,2),b(2),bint(2,1),bint(2,2))
fprintf(fout, "Chan %s\t Evoked: cc %.3f, slp %.3f(%.3f,%.3f) interc %.1f(%.1f,%.1f)  Score: cc %.3f, slp %.3f(%.3f,%.3f) interc %.1f(%.1f,%.1f)\n",unit_name_arr(channel_j),...
    ccoef2,b2(1),bint2(1,1),bint2(1,2),b2(2),bint2(2,1),bint2(2,2),...
    ccoef,b(1),bint(1,1),bint(1,2),b(2),bint(2,1),bint(2,2)) ;
end
fclose(fout);
save(fullfile(savepath,"noise_stats.mat"),"cc_col","slope_col","evk_cc_col","evk_slope_col")
% Plot the pooled statistics across channel
figure(15);clf;set(15,'position',[ 214         122        1587         856]);
subplot(221)
histogram(cc_col,30)
xlabel("correlation coefficient");title(["Score (evoked rate - baseline rate)"]);
subplot(222)
histogram(slope_col,30)
xlabel("slope");title(["Score (evoked rate - baseline rate)"]);
subplot(223)
histogram(evk_cc_col,30)
xlabel("correlation coefficient");title(["Firing rate in 50-200 Evoked window"]);
subplot(224)
histogram(evk_slope_col,30)
xlabel("slope");title(["Firing rate in 50-200 Evoked window"]);
suptitle([sprintf("Manifold %s", Exp_label_str), "correlation and slope between std and mean"])
saveas(15, fullfile(savepath, compose("noise_estim_summary.png")))
end
%%
figure
subplot(231)
imagesc(score_mean)
axis image
subplot(232)
imagesc(score_std)
axis image
subplot(233)
imagesc(score_var)
axis image
subplot(234)
scatter(mean_vect, std_vect)
xlabel("mean response")
ylabel("score std")
subplot(235)
histogram(score_mat - score_mean)
%%
figure, imagesc(nanstd(bsl_mat,0,3))

% Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
% savepath = fullfile(result_dir, compose("Evol%02d_chan%02d", Expi, pref_chan(1)));
% mkdir(savepath);
% unit_in_pref_chan = cell2mat(Trials_new{1}.TrialRecord.User.evoConfiguration(:,4))';
%%



function [score_mat, bsl_mat, summary, stat_str] = get_stats_from_result(name_pattern, channel)
    global  Trials rasters sphere_norm ang_step Reps
    score_mat = nan(11,11,Reps); 
    bsl_mat = nan(11,11,Reps);  % make sure the 3 dimension has larger size than repitition number! 
    cnt_mat = zeros(11,11); 
    id_mat = zeros(11,11); % record the id correspond to i,j
    for i =  -5:5
        for j =  -5:5
            cur_fn = sprintf(name_pattern, sphere_norm, i*ang_step, j*ang_step);
            img_idx = find(contains(Trials.imageName, cur_fn));
            cnt_mat(i+6, j+6) = length(img_idx);
            psths = rasters(channel, :, img_idx);
            scores = squeeze(mean(psths(1, 50:150, :))- mean(psths(1, 1:50, :)) );
            baseline = squeeze(mean(psths(1, 1:50, :)));
            score_mat(i+6, j+6, 1:length(img_idx)) = scores;
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