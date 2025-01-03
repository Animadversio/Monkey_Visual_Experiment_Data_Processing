%%
Animal = "Alfa";Set_Path;
ftrrows = find(...
               contains(ExpRecord.expControlFN,["select"])&...
               contains(ExpRecord.Exp_collection,["FC6Evol_Decomp"])&...
               ~isnan(ExpRecord.Expi)...
               );
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows, Animal, false);
%%

sumdir = "O:\corrFeatVis_FactorPredict\summary";
%%
tval_arr = [];
mean_mat = [];
sem_mat = [];
for  iTr = 1:numel(S_col)
    if iTr == 4, continue;end
    fulltsr_gi = find(strcmp(S_col(iTr).stim.grouplabs,'full_tsr'));
    evobest_gi = find(strcmp(S_col(iTr).stim.grouplabs,'evol_best'));
%     [~,P,CI,TST] = ttest2(S_col(iTr).prefresp.group_score_col{fulltsr_gi},S_col(iTr).prefresp.group_score_col{evobest_gi});
%     fprintf("fulltsr - evobest: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
    [tval,pval,sumstr,mean_arr,sem_arr] = ttest2_print(S_col(iTr).prefresp.group_score_col{fulltsr_gi},S_col(iTr).prefresp.group_score_col{evobest_gi},"fulltsr","evolbest");
    tval_arr = [tval_arr,tval];
    mean_mat = [mean_mat;mean_arr];
    sem_mat = [sem_mat;sem_arr];
end
%%
examp_col = [];
for  iTr = 1:numel(S_col)
    if iTr == 4, continue;end
    evobest_gi = find(strcmp(S_col(iTr).stim.grouplabs,'evol_best'));
    fulltsr_gi = find(strcmp(S_col(iTr).stim.grouplabs,'full_tsr'));
    evobest_imgfn = S_col(iTr).stim.imgfn_mapper(S_col(iTr).prefresp.bestimgnms{evobest_gi});
    fulltsr_imgfn = S_col(iTr).stim.imgfn_mapper(S_col(iTr).prefresp.bestimgnms{fulltsr_gi});
    examp_col{iTr,1} = imread(evobest_imgfn);
    examp_col{iTr,2} = imread(fulltsr_imgfn);
    % [~,P,CI,TST] = ttest2(S_col(iTr).prefresp.group_score_col{fulltsr_gi},S_col(iTr).prefresp.group_score_col{evobest_gi});
    % fprintf("fulltsr - evobest: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",P,TST.tstat,TST.df,CI(1),CI(2));
    % [tval,pval,sumstr,mean_arr,sem_arr] = ttest2_print(S_col(iTr).prefresp.group_score_col{fulltsr_gi},S_col(iTr).prefresp.group_score_col{evobest_gi},"fulltsr","evolbest");
    % tval_arr = [tval_arr,tval];
    % mean_mat = [mean_mat;mean_arr];
    % sem_mat = [sem_mat;sem_arr];
end
%%
figure(6); set(6, 'pos', [40   280   2520   420])
montage(examp_col(:),'Size',[2,18])%,'Parent',gca
title("Comparing Evolved Best Image and Re-evolved from the Linear Model")
saveallform(sumdir, "FullTsrOptim_BestEvol_image_cmp")
%%
figure(7); clf; set(7, 'pos', [1000  489  500  500])
%scatter(mean_mat(:,1),mean_mat(:,2));axis image
errorbar(mean_mat(:,1),mean_mat(:,2),sem_mat(:,1),sem_mat(:,1),sem_mat(:,2),sem_mat(:,2),'o');
axis image;box off;
xlabel("Full Tensor Projection")
ylabel("Best in Evolution")
title(compose("Comparing Neeural Responses of Projection of Full \nTensor Model and Best Evolved Img(1 Exp per pnt)"))
add_diagonal(gca,'k:');
saveallform(sumdir, "FullTsrOptim_BestEvol_score_cmp")
%%
h=stripe_simple_plot(tval_arr,"tval","FullTsr-BestEvol",sumdir,"");
%%

EvolFactStats = S_col;
save(fullfile(sumdir, Animal+"_EvolFactStats.mat"), "EvolFactStats")
save(fullfile(matdir, Animal+"_EvolFactStats.mat"), "EvolFactStats")