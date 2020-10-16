%% BigGAN load codes
Animal="Both"; Set_Path;
% ftrrows = find(contains(ExpRecord.expControlFN,"generate_") & (ExpRecord.Exp_collection=="BigGAN_fc6" |...
%                ExpRecord.Exp_collection=="BigGAN_FC6")) ;
ftrrows = find(contains(ExpRecord.expControlFN,"selectiv") &...
               ExpRecord.Exp_collection=="BigGAN_Hessian") ;
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftrrows, Animal, true, true); % 
%% 
ftrrows = find(contains(ExpRecord.expControlFN,"generate_") & (ExpRecord.Exp_collection=="BigGAN_fc6" |...
               ExpRecord.Exp_collection=="BigGAN_FC6")) ;
%%
[codes_all, img_ids, generations] = load_codes_all(ExpRecord.stimuli{ftrrows(1)}, 2);
%% Load up Hessian matrix
Hdata = py.numpy.load("E:\OneDrive - Washington University in St. Louis\Hessian_summary\BigGAN\H_avg_1000cls.npz");
eva_all = Hdata.get('eigvals_avg').double;
evc_all = Hdata.get('eigvects_avg').double;
evc_all = evc_all(:,end:-1:1);
eva_all = eva_all(end:-1:1);
eva_cls = Hdata.get('eigvals_clas_avg').double;
evc_cls = Hdata.get('eigvects_clas_avg').double;
evc_cls = evc_cls(:,end:-1:1);
eva_cls = eva_cls(end:-1:1);
eva_nos = Hdata.get('eigvals_nois_avg').double;
evc_nos = Hdata.get('eigvects_nois_avg').double;
evc_nos = evc_nos(:,end:-1:1);
eva_nos = eva_nos(end:-1:1);
evc_nos_ag = [evc_nos;zeros(128)];
evc_cls_ag = [zeros(128);evc_cls];
%% Go through all the evolutions using BigGAN.
tic
gen_num_col = [];
last_gen_col = {};
last_gen_mean_col = [];
for Expi = 1:numel(ftrrows)
try
[codes_all, img_ids, generations] = load_codes_all(ExpRecord.stimuli{ftrrows(Expi)}, 2);
catch 
    fprintf("Exp %d folder %s loading Error. Comments:\n", Expi, ExpRecord.stimuli{ftrrows(Expi)})
    fprintf(ExpRecord.comments{ftrrows(Expi)})
   continue
end
gen_num = max(generations);
last_gen_code = codes_all(generations == gen_num, :);
last_gen_col{end+1} = last_gen_code;
last_gen_mean_col(end+1, :) = mean(last_gen_code,1);
gen_num_col(end+1) = gen_num;
end
toc
%%
figure;hold on
for i = 1:size(last_gen_mean_col,1)
   scatter(1:128, last_gen_mean_col(i,:)) 
end

%%
perm_last_gen_mean = [];
for i = 1:size(last_gen_mean_col,1)
perm_last_gen_mean = [perm_last_gen_mean; last_gen_mean_proj(i,randperm(128))];
end
%%
last_gen_mean_proj = last_gen_mean_col * evc_cls;
perm_last_gen_proj = perm_last_gen_mean * evc_cls;
%%
figure;hold on
for i = 1:size(last_gen_mean_col,1)
   scatter(1:128, last_gen_mean_proj(i,:)) 
end
%%
figure;
errorbar(1:128, mean(abs(last_gen_mean_proj(:,:)),1), std(abs(last_gen_mean_proj(:,:)),1))
corr([1:128]', mean(abs(last_gen_mean_proj(:,:)),1)');
%%
cc_idx = corr([1:128]', mean(abs(last_gen_mean_proj(:,:)),1)')
cc_eva = corr(eva_cls', mean(abs(last_gen_mean_proj(:,:)),1)')
cc_logeva = corr(log10(eva_cls)', mean(abs(last_gen_mean_proj(:,:)),1)')

cc_idx_shfl = corr([1:128]', mean(abs(perm_last_gen_proj(:,:)),1)')
cc_eva_shfl = corr(eva_cls', mean(abs(perm_last_gen_proj(:,:)),1)')
cc_logeva_shfl = corr(log10(eva_cls)', mean(abs(perm_last_gen_proj(:,:)),1)')
%% Confidence Interval Computation
cc_idx_dist = []; cc_eva_dist = []; cc_logeva_dist = [];
for triali = 1:1000
perm_last_gen_mean = [];
for i = 1:size(last_gen_mean_col,1)
perm_last_gen_mean = [perm_last_gen_mean; last_gen_mean_proj(i,randperm(128))];
end
perm_last_gen_proj = perm_last_gen_mean * evc_cls;
cc_idx_shfl = corr([1:128]', mean(abs(perm_last_gen_proj(:,:)),1)');
cc_eva_shfl = corr(eva_cls', mean(abs(perm_last_gen_proj(:,:)),1)');
cc_logeva_shfl = corr(log10(eva_cls)', mean(abs(perm_last_gen_proj(:,:)),1)');
cc_idx_dist = [cc_idx_dist, cc_idx_shfl];
cc_eva_dist = [cc_eva_dist, cc_eva_shfl];
cc_logeva_dist = [cc_logeva_dist, cc_logeva_shfl];
end
%% Verbal statement
fprintf("For %d BigGAN Evolution Experiment, we collect the mean last generation gene, and project\n"+...
    " it using the Hessian EigenBasis of Class subspace. The mean absolute value of projection coefficient\n"+...
    " shows a correlation with the eigenIdx; eigenvalue; and log(eigenvalue), \n"+...
    " correlation with eigenIdx (%.3f) eigenvalue (%.3f) and log(eigenvalue) (%.3f), \n"+...
    "while the Correlation with the projection of shuffled code is less pronounced \n"+...
    "[i.e. correlation with the eigenIdx (%.3f, [%.3f, %.3f]) eigenvalue (%.3f, [%.3f, %.3f]) and log(eigenvalue) (%.3f, [%.3f, %.3f])\n"+...
    "1-99 percentile denoted in [ , ]\n",...
    size(last_gen_mean_proj,1), cc_idx, cc_eva, cc_logeva, ...
    mean(cc_idx_dist), prctile(cc_idx_dist,1), prctile(cc_idx_dist,99), ...
    mean(cc_eva_dist), prctile(cc_eva_dist,1), prctile(cc_eva_dist,99), ...
    mean(cc_logeva_dist), prctile(cc_logeva_dist,1), prctile(cc_logeva_dist,99))
%% Plotting shuffling result! 
figdir = "E:\OneDrive - Washington University in St. Louis\HessEvolStruct_BigGAN";
H = figure('position',[300,500, 900, 300]);
T = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile(1)%subplot(131)
histogram(cc_idx_dist, 30)
line([cc_idx, cc_idx], ylim(), 'Color', 'red')
legend(["Permutation Control","Real Value"])
xlabel("Corr with eigen rank")
ylabel("Density")
title(compose("cc=%.3f (p=%.4f)", cc_idx, 1-sum(cc_idx_dist>cc_idx)/numel(cc_idx_dist)))
box off;
nexttile(2)%subplot(132)
histogram(cc_eva_dist, 30)
line([cc_eva, cc_eva], ylim(), 'Color', 'red')
xlabel("Corr with eigval")
title(compose("cc=%.3f (p=%.4f)", cc_eva, sum(cc_eva_dist>cc_eva)/numel(cc_eva_dist)))
box off;
nexttile(3)%subplot(133)
histogram(cc_logeva_dist, 30)
line([cc_logeva, cc_logeva], ylim(), 'Color', 'red')
xlabel("Corr with Log(eigval)")
title(compose("cc=%.3f (p=%.4f)", cc_logeva, sum(cc_logeva_dist>cc_logeva)/numel(cc_logeva_dist)))
box off;
title(T,"Correlation of Projection Amplitude of Evolved Code and Hessian Eigen Structures (BigGAN)")
savefig(H, fullfile(figdir,"BigGAN_proj_eig_corr_permtest.fig"))
saveas(H, fullfile(figdir,"BigGAN_proj_eig_corr_permtest.png"))
%%
figure;hold on 
errorbar(1:128, mean(abs(last_gen_mean_proj(:,:)),1), std(abs(last_gen_mean_proj(:,:)),1)/sqrt( size(last_gen_mean_proj,1)))
ylabel("abs(proj coef)");xlabel("eig rank (BigGAN class space)");
% errorbar(1:128, mean(abs(perm_last_gen_proj(:,:)),1), std(abs(perm_last_gen_proj(:,:)),1)/sqrt( size(last_gen_mean_proj,1)))
% scatter(1:128, mean(abs(perm_last_gen_proj(:,:)),1))