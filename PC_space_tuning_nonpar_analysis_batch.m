% Batch processing code  calculating the statistics (t, 1way ANOVA, 2way ANOVA) 
% and generate annotated figures for each and every experiment 
%% 
storedStruct = load("D:\\Manifold_Exps.mat");
%%
%load("D:\\Beto64chan-02102019-003_formatted");
norm_arr = [328, 326, 269, 329, 401, 386];
pref_chan_arr = [29, 6, 5, 20, 19, 13]; 
global sphere_norm Trials channel rasters ang_step Reps meta
%%
for Expi=1:6
    rasters = storedStruct.rasters{Expi};
    Trials = storedStruct.Trials{Expi};
    meta = storedStruct.meta{Expi};
    unit_name_arr = generate_unit_labels();
Stat_summary = {};
pref_chan = pref_chan_arr(Expi);
sphere_norm = norm_arr(Expi);%328; % 269 Day3 % 326 Daye % 328 Day 1 [328, 326, 269, 329, 401]
ang_step = 18;
Reps = 10; % can be any number LARGER than the largest repitition. or there may be problems caused by NAN and 0 filling
% savepath = "C:\Users\ponce\OneDrive\Desktop\OneDrive_Binxu\OneDrive\PC_space_tuning\Exp1_chan29";
savepath = sprintf("C:\\Users\\ponce\\OneDrive\\Desktop\\OneDrive_Binxu\\OneDrive\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
mkdir(savepath);
fid = fopen(fullfile(savepath, "Unit_Label.txt"),'w');
fprintf(fid,'%s\n', unit_name_arr{:});
fclose(fid);
for channel = 1:size(rasters,1)
figure(1);clf;set(1, 'position', [304    12   560   577]);
figure(2);clf;set(2, 'position', [304    12   560   577]);
figure(3);clf;set(3, 'position', [304    12   560   577]);
figure(4);clf;set(4, 'position', [73  -40  2418  705]);
chan_label_str = sprintf("Exp%d Channel %s", Expi, unit_name_arr{channel});
%% PC12
[score_mat, ~, summary, stat_str] = get_stats_from_result('norm_%d_PC2_%d_PC3_%d');
Stat_summary{channel, 1} = summary;
disp(summary)
% visualize
figure(1);
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("PC 2 degree")
xlabel("PC 3 degree")
title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
shading flat
axis image
colorbar

figure(4)
subplot(131)
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("PC 2 degree")
xlabel("PC 3 degree")
title([chan_label_str, "Tuning map on PC2 3 subspace", stat_str])
shading flat
axis image
colorbar

%% PC4950
[score_mat, ~, summary, stat_str] = get_stats_from_result('norm_%d_PC49_%d_PC50_%d');
Stat_summary{channel, 2} = summary;
disp(summary)
    
figure(2);
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("PC 49 degree")
xlabel("PC 50 degree")
title([chan_label_str, "Tuning map on PC49 50 subspace", stat_str])
shading flat
axis image
colorbar

figure(4)
subplot(132)
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("PC 49 degree")
xlabel("PC 50 degree")
title([chan_label_str, "Tuning map on PC49 50 subspace", stat_str])
shading flat
axis image
colorbar
%% RND12
[score_mat, ~, summary, stat_str] = get_stats_from_result('norm_%d_RND1_%d_RND2_%d');
Stat_summary{channel, 3} = summary;
disp(summary)

figure(3);clf
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("Rand vector 1 degree")
xlabel("Rand vector 2 degree")
title([chan_label_str, "Tuning map on Random Vector (outside first 50 PCs) subspace", stat_str])
shading flat
axis image
colorbar

figure(4)
subplot(133)
imagesc(-90:18:90, -90:18:90, nanmean(score_mat,3))
ylabel("Rand vector 1 degree")
xlabel("Rand vector 2 degree")
title([chan_label_str, "Tuning map on Random Vector (outside first 50 PCs) subspace", stat_str])
shading flat
axis image
colorbar
%%
saveas(4, fullfile(savepath, sprintf("chan%02d_PC_tune_cmp_stat.png", channel)))
saveas(1, fullfile(savepath, sprintf("chan%02d_PC23_tune_stat.png", channel)))
saveas(2, fullfile(savepath, sprintf("chan%02d_PC4950_tune_stat.png", channel)))
saveas(3, fullfile(savepath, sprintf("chan%02d_RND_tune_stat.png", channel)))
end
save(fullfile(savepath, "Basic_Stats.mat"), 'Stat_summary')
end
%%
function [score_mat, bsl_mat, summary, stat_str] = get_stats_from_result(name_pattern)
    global  Trials rasters channel sphere_norm ang_step Reps
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
            baseline = squeeze(mean(psths(1, 50:150, :)));
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

function unit_name_arr = generate_unit_labels()
% Generate the unit labels 17B from the spikeID variable
global meta
Unit_id = meta.spikeID;
unit_name_arr = {}; % name tag for each unit 
for i = 1:length(Unit_id)
    cur_chan = Unit_id(i);
    if sum(Unit_id == cur_chan) == 1
        unit_name_arr{i} = num2str(cur_chan);
    else
        cur_chan = Unit_id(i);
        rel_idx = find(find(Unit_id == cur_chan) == i);
        unit_name_arr{i} = [num2str(cur_chan), char(64+rel_idx)];
    end
end
end

