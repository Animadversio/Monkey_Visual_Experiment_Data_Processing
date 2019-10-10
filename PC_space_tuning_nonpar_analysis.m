% load(fullfile("\\storage1.ris.wustl.edu\crponce\Active\Data-Ephys-MAT","Beto64chan-03102019-002_formatted.mat"))
%%
load("D:\\Beto64chan-02102019-003_formatted");
%%
img_dir = "\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-02a-beto";
sphere_norm = 328;
ang_step = 18;
%%
channel = 62;
score_mat = nan(11,11,8); 
bsl_mat = nan(11,11,8); 
mean_fr_mat = nan(11,11,8);  % make sure the 3 dimension has larger size than repitition number! 
cnt_mat = zeros(11,11); 
% score_vec = []; % concatenate the scores together 
% bsl_vec = [];
% id_vec = []; % code the i,j as an id for one way anova
% theta_vec = []; 
% phi_vec = []; 

id_mat = zeros(11,11); % record the id correspond to i,j
for i =  -5:5
    for j =  -5:5
        cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', sphere_norm, i*ang_step, j*ang_step);
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
        %id_vec = [id_vec, id*ones(1, length(img_idx))]; % for 1 way anova 
        %score_vec = [score_vec, squeeze(scores)']; % for ttest 
        %bsl_vec = [bsl_vec, baseline'];
        %theta_vec = [theta_vec, i*ang_step*ones(1, length(img_idx))]; % for 2 way anova 
        %phi_vec = [phi_vec, j*ang_step*ones(1, length(img_idx))];
    end
end
[theta_mat, phi_mat] = meshgrid(ang_step*(-5:5), ang_step*(-5:5));
mean_fr_mat = bsl_mat + score_mat;
id_vec_nan = reshape(repmat(id_mat, 1,1,8), 1, []);
score_vec_nan = reshape(score_mat, 1, []);
bsl_vec_nan = reshape(bsl_mat, 1, []);
mean_fr_vec_nan = bsl_vec_nan + score_vec_nan;
%%
[p,tbl,stats] = anova1(score_vec_nan, id_vec_nan);%,'off'
stats.F = tbl{2,5};
stats.p = p;
multcompare(stats);
summary.anova_F = stats.F;
summary.anova_p = stats.p; 
%%
[p,tbl,stats] = anovan(score_vec_nan, {reshape(repmat(theta_mat, 1,1,8),1,[]), reshape(repmat(phi_mat, 1,1,8),1,[])}, 'model', 'interaction');
multcompare(stats);
summary.anova2_p = p; % p for theta, phi and interaction 
summary.anova2_F = [tbl{2:4,6}]; % F for theta, phi and interaction
%%
[~,P,CI] = ttest(mean_fr_vec_nan, bsl_vec_nan);
summary.t_p = P;
summary.t_CI = CI;


