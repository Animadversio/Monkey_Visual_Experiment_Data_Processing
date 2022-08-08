%% Manifold paper response to reviewers 
% Compare image diversity for Manifold Experiments driven by units from different
% Visual areas. 
% IT driven manifold images are more diverse than V1 driven manifold images. 
% Consistent with CNN results. 
%%
Animal="Beto";Set_Path;
%%
S_col = [];
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_Manif_ImDist.mat"),'ManifImDistStat');
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"))
%%
for Expi = 1:numel(ManifImDistStat)
    S = struct();
    S.Animal = Animal;
    S.Expi = Expi;
    S.prefchan = Stats(Expi).units.pref_chan;
    for dist_metric = ["squ","L2","FC6","FC6_corr","SSIM"]
    distmat = ManifImDistStat(Expi).(dist_metric);
%     distmat = ManifImDistStat(Expi).squ;
%     distmat = ManifImDistStat(Expi).L2;
    distmat_cent = distmat(11:111,11:111);
    msk = triu(nan(101));
    S.(dist_metric+"std") = nanstd(distmat_cent+msk, 1,'all');
    S.(dist_metric+"mean") = nanmean(distmat_cent+msk, 'all');
    end
    fprintf("Exp%d %.3f+-%.3f\n",Expi,S.squmean,S.squstd)
    S_col = [S_col;S];
end
end
%%
imdivers_tab = struct2table(S_col);
%%
ITvec = imdivers_tab.LPIPSmean((imdivers_tab.prefchan<33));
V1vec = imdivers_tab.LPIPSmean((imdivers_tab.prefchan<49) & (imdivers_tab.prefchan>32));
V4vec = imdivers_tab.LPIPSmean((imdivers_tab.prefchan>48));
%%
ttest2_print(ITvec, V1vec, "IT", "V1");
ttest2_print(ITvec, V4vec, "IT", "V4");
ttest2_print(V4vec, V1vec, "V4", "V1");
%%
global figdir 
figdir = "E:\OneDrive - Harvard University\Manuscript_Manifold\Response\Manif_imdiversity";
%%
ITmsk = (imdivers_tab.prefchan<33);
V1msk = (imdivers_tab.prefchan<49) & (imdivers_tab.prefchan>32);
V4msk = (imdivers_tab.prefchan>48);
for dist_metric = ["squ","L2","FC6","FC6_corr","SSIM"]
stripe_plot(imdivers_tab, dist_metric+"mean",{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    "img diversity", "ManifImg_diversity", {[1,2],[2,3],[1,3]})
end
%%
