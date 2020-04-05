addpath D:\Github\Fit_Spherical_Tuning

ft = fittype( @(theta, phi, psi, kappa, beta, A, x, y) KentFunc(theta, phi, psi, kappa, beta, A, x, y), ...
    'independent', {'x', 'y'},'dependent',{'z'});
%% Fit the Kent Distribution. Collect statistics
tic
Kent_stats = repmat(struct(), 1, length(meta_new));
for i = 1:length(Stats)
% there can be multiple units, or multiple subspace, so we have this cell
% array
stats_grid = cell(numel(Stats(i).units.pref_chan_id), Stats(i).manif.subsp_n);
scr_stats_grid = cell(numel(Stats(i).units.pref_chan_id), Stats(i).manif.subsp_n);
subsp_suff = {'PC23',"PC4950",'RND12'};
subsp_axis = ["PC2","PC3";"PC49","PC50";"RND1","RND2"];
for subsp_i = 1:Stats(i).manif.subsp_n
imgnm_grid = cellfun(@(idx) unique([Stats(i).imageName(idx)]), Stats(i).manif.idx_grid{subsp_i});
psths = Stats(i).manif.psth{subsp_i};
fprintf(subsp_suff{subsp_i}+"\n")
for uniti = 1:1:numel(Stats(i).units.pref_chan_id)
channel_j = Stats(i).units.pref_chan_id(uniti);
chan_label_str = sprintf("%s Exp%d Channel %s", Stats(i).Animal, Stats(i).Expi, ...
            Stats(i).units.unit_name_arr{channel_j});

mean_score_map = cellfun(@(psth) mean(psth(uniti, 51:200, :),[2,3])-mean(psth(uniti, 1:40, :),[2,3]), psths);
[Parameter, gof] = fit_Kent(mean_score_map);
scr_stats_grid{uniti, subsp_i} = gof;
% V = coeffvalues(Parameter);
% CI = confint(Parameter);
% param_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
%                 V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));

mean_act_map = cellfun(@(psth) mean(psth(uniti, 51:200, :),[2,3]), psths);
[Parameter, gof] = fit_Kent(mean_act_map);
stats_grid{uniti, subsp_i} = gof;
V = coeffvalues(Parameter);
CI = confint(Parameter);
param_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
                V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
fprintf(chan_label_str+"\n "+param_str+"\n")
end
end
Kent_stats(i).act_fit = stats_grid;
Kent_stats(i).scr_fit = scr_stats_grid;
end
toc % took 6.88s to run through! not intensive
%%
save('D:\Alfa_Manif_Kent_Fit.mat','Kent_stats')
save(fullfile(savepath, 'Alfa_Manif_Kent_Fit.mat'),'Kent_stats')
%% Plot the summary Scatter and line 
Expi_col = [1,3,4,5,8,9,10,11,12,13,15,16,17,18,19,20,21,22];
KSU_vec = [];
KSU_CI = [];
KHash_vec = [];
KHash_CI = [];
for i = 1:numel(Expi_col)
Expi = Expi_col(i);
SUi = 1; Hai = 2; % index for SU and Ha within each struct
if Expi == 10 && Animal=='Alfa', SUi = 2; Hai =1; end % At Exp 10 SU label swapped....
KSU_vec = [KSU_vec, Kent_stats(Expi).act_fit{SUi}.coef(4)];
KSU_CI = [KSU_CI, Kent_stats(Expi).act_fit{SUi}.confint(:,4)];

KHash_vec = [KHash_vec, Kent_stats(Expi).act_fit{Hai}.coef(4)];
KHash_CI = [KHash_CI, Kent_stats(Expi).act_fit{Hai}.confint(:,4)];
end
%%
i_list = 1:numel(Expi_col);
figure(6);clf;hold on
errorbar(i_list, KSU_vec, KSU_vec-KSU_CI(1,:), KSU_CI(2,:)-KSU_vec,'o')
errorbar(i_list, KHash_vec, KHash_vec-KHash_CI(1,:), KHash_CI(2,:)-KHash_vec,'o')
axis equal tight
xlabel("Exp num");ylabel("kappa (Kentfun)")
xticks(i_list)
xticklabels(Expi_col)

figure(7);clf;hold on
errorbar(KSU_vec, KHash_vec, ...
        KSU_vec-KSU_CI(1,:), KSU_CI(2,:)-KSU_vec,...
        KHash_vec-KHash_CI(1,:), KHash_CI(2,:)-KHash_vec,'o')
xlabel("kappa (SU)");ylabel("kappa (Hash)")
line([0,2],[0,2])
axis equal



%%
errorbar(kap_SU, kapHash)
%% Plot it with spherical plot 
ang_step = 18;
theta_arr = -90:ang_step:90;
phi_arr   = -90:ang_step:90;
[phi_grid, theta_grid] = meshgrid(theta_arr, phi_arr);
figure(3)
ax1 = subplot(121);
sphere_plot(ax1, theta_grid, phi_grid, mean_act_map);
ylabel("RND1(theta)");zlabel("RND2(phi)");cl = caxis;
title([chan_label_str, "Tuning map on RND 1 2 subspace"])%, stat_str])
ax2 = subplot(122);
Kent_sphere_plot(ax2, coeffvalues(Parameter));
ylabel("RND1(theta)");zlabel("RND2(phi)");caxis(cl);
title([chan_label_str, "Kent Fit Tuning map on RND 1 2 subspace", param_str])
%%



%%
mean_act = nanmean(score_col{uniti}, 3 );
[max_score, max_idx ] = max(mean_act,[],'all','linear');
[r_idx,c_idx] = ind2sub(size(mean_act), max_idx);
ang_PC2 = ang_step * (-5 + r_idx -1);
ang_PC3 = ang_step * (-5 + c_idx -1);
[ang_dist, cos_dist] = dist_grid_on_sphere(ang_PC2, ang_PC3);
subplot(1,numel(pref_chan_id),uniti)
scatter(ang_dist(:), mean_act(:))
[b,bint,~,~,stats] = regress(mean_act(:), [ang_dist(:), ones(numel(cos_dist), 1)]);

