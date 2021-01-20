% Calculate the tuning statistics t,F, p from a matrix of response
% Majorly designed for the original Manifold experiments.
function [summary, stat_str] = calc_tuning_stats(score_mat, bsl_mat, theta_arr, phi_arr)
Reps = size(score_mat, 3);
num_theta = size(score_mat, 1);
num_phi = size(score_mat, 2);
assert(num_theta==length(theta_arr))
assert(num_phi==length(phi_arr))
id_mat = cell(num_theta, num_phi);
for i = 1:num_theta
    theta = theta_arr(i);
    for j = 1:num_phi
        phi = phi_arr(j);
        if phi ~= 90 && phi ~= -90
            id = sprintf("%d_%d", theta, phi);
        elseif phi == 90
            id = sprintf("%d_%d", 0, phi);
        elseif phi ==-90
            id = sprintf("%d_%d", 0, phi);
        end
        id_mat{i, j} = id;
    end
end
id_mat = string(id_mat);
[phi_mat, theta_mat] = meshgrid(phi_arr, theta_arr);
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
        'Theta, F=%.2f(%.3f), Phi, F=%.2f(%.3f), Interact, F=%.2f(%.3f)'],...
        CI(1), CI(2), P, stats.F, stats.p, ...
        stats2.F(1),stats2.p(1), stats2.F(2),stats2.p(2),stats2.F(3),stats2.p(3));
end