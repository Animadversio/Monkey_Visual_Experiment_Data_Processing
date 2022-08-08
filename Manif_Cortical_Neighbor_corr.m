% 
Set_Path;
%%
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_CortiDisCorr.mat"))
for Expi = 1:numel(CortiDisCorr)
M = struct();
areavec = area_map(CortiDisCorr(Expi).units.spikeID,"Alfa");
M.ITmsk = areavec=="IT"; M.V4msk = areavec=="V4"; M.V1msk = areavec=="V1";
Fmsk = ([CortiDisCorr(Expi).FStats.F_P]'<0.01) ;
valmsk   = CortiDisCorr(Expi).units.activ_msk;
prefchan = CortiDisCorr(Expi).units.pref_chan;
Cortdist = CortiDisCorr(Expi).cortexDismat;
corrmat = CortiDisCorr(Expi).avgsph_corrmat;
%%
samemsk = abs(Cortdist - 0) < 1E-5;
neighbormsk = abs(Cortdist - 400) < 1E-5;
nonneighbormsk = (~neighbormsk)&(~samemsk)&~isnan(Cortdist);
% diagmsk = nan(size(Cortdist));
%%
% nonneighcorrvec = CortiDisCorr(Expi).avgsph_corrmat(nonneighbormsk);
% corrvec  = CortiDisCorr(Expi).avgsph_corrmat(neighbormsk);
%%
fprintf("%s Exp%d chan %d\n",Animal,Expi, prefchan)
prefchan_id = CortiDisCorr(Expi).units.pref_chan_id(1);
prefneighmsk = (Cortdist(:,prefchan_id) == 400) & valmsk & Fmsk;
pref_neigh_corr = corrmat(prefchan_id, prefneighmsk);
pref_nonneigh_msk = ones(size(corrmat),'logical');
pref_nonneigh_msk(prefchan_id, prefneighmsk) = 0;
pref_nonneigh_msk(prefneighmsk, prefchan_id) = 0;

disp(pref_neigh_corr)
uppermsk = triu(ones(size(Cortdist)),1);
for area = ["V1","V4","IT"]
areamsk = M.(area+"msk");
unitmsk = areamsk & Fmsk;
% ITcorrvec = corrmat((unitmsk*unitmsk') & neighbormsk & uppermsk);
% ITnoncorrvec = corrmat((unitmsk*unitmsk') & nonneighbormsk & uppermsk);
% ttest2_print(ITcorrvec, ITnoncorrvec,area+" neighbor",area+" non-neighbor");
corrmat(prefchan_id, prefneighmsk) = nan;
corrmat(prefneighmsk, prefchan_id) = nan;
ITcorrvec = corrmat((unitmsk*unitmsk') & pref_nonneigh_msk & uppermsk);
% ITnoncorrvec = corrmat((unitmsk*unitmsk') & nonneighbormsk & uppermsk);
% ttest2_print(ITcorrvec, ITnoncorrvec,area+" neighbor",area+" non-neighbor");
disp(compose("%s %.3f+-%.3f (N=%d)",area, mean(ITcorrvec),sem(ITcorrvec),numel(ITcorrvec)))
end
end
end
%%
Animvec_pool = [];
Expivec_pool = [];
corrvec_pool = [];
distvec_pool = [];
areavec_pool = [];
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_CortiDisCorr.mat"))
for Expi = 1:numel(CortiDisCorr)
M = struct();
areavec = area_map(CortiDisCorr(Expi).units.spikeID,"Alfa");
M.ITmsk = areavec=="IT"; M.V4msk = areavec=="V4"; M.V1msk = areavec=="V1";
Fmsk = ([CortiDisCorr(Expi).FStats.F_P]'<0.01) ;
valmsk   = CortiDisCorr(Expi).units.activ_msk;
prefchan = CortiDisCorr(Expi).units.pref_chan;
Cortdist = CortiDisCorr(Expi).cortexDismat;
corrmat  = CortiDisCorr(Expi).avgsph_corrmat;
%%
uppermsk = triu(ones(size(Cortdist)),1);
for area = ["V1","V4","IT"]
unitmsk = M.(area+"msk") & Fmsk; %& valmsk;
corrvec = corrmat((unitmsk*unitmsk') & uppermsk);
distvec = Cortdist((unitmsk*unitmsk') & uppermsk);
corrvec_pool = [corrvec_pool;corrvec];
distvec_pool = [distvec_pool;distvec];
Animvec_pool = [Animvec_pool;repmat(Animal,numel(distvec),1)];
Expivec_pool = [Expivec_pool;Expi*ones(numel(distvec),1)];
areavec_pool = [areavec_pool;repmat(area,numel(distvec),1)];
end
end
end
%%
tab = array2table(Animvec_pool,'VariableNames',["Animal"]);
tab.Expi = Expivec_pool;
tab.area = areavec_pool;
tab.dist = distvec_pool;
tab.mapcorr = corrvec_pool;
writetable(tab, fullfile(savedir, "Both_all_pair_df.csv"));
%%
[cval, pval] = corr(distvec_pool, corrvec_pool,'type','pearson')
[cval, pval] = corr(distvec_pool, corrvec_pool,'type','spearman');
[cval, pval] = corr(distvec_pool(distvec_pool>0), corrvec_pool(distvec_pool>0),'type','pearson')
[cval_excl0, pval_excl0] = corr(distvec_pool(distvec_pool>0), corrvec_pool(distvec_pool>0),'type','spearman');

%%
savedir = "E:\OneDrive - Harvard University\Manuscript_Manifold\Response\Cortical_smoothness";
% [cval, pval] = corr(distvec_pool, corrvec_pool,'type','pearson')
[cval, pval] = corr(distvec_pool, corrvec_pool,'type','spearman');
% [cval, pval] = corr(distvec_pool(distvec_pool>0), corrvec_pool(distvec_pool>0),'type','pearson')
[cval_excl0, pval_excl0] = corr(distvec_pool(distvec_pool>0), corrvec_pool(distvec_pool>0),'type','spearman');
figure('position',[200,200,500,500]);
scatter(distvec_pool, corrvec_pool,...
    'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
xlabel("cortical distance (\mum)")
ylabel("Tuning map correlation")
title(["Relation of Cortical distance and Tuning map similarity",...
       compose("Spearman Corr %.3f(P=%.1e)\n Excluding 0 Spearman Corr %.3f(P=%.1e)",...
       cval, pval,cval_excl0, pval_excl0)])
saveallform(savedir,"overall_corrplot")
%%
% [cval, pval] = corr(distvec_pool, corrvec_pool,'type','spearman');
% [cval, pval] = corr(distvec_pool(distvec_pool>0), corrvec_pool(distvec_pool>0),'type','pearson')
[cval_alfa, pval_alfa] = corr(distvec_pool(Animvec_pool=="Alfa"), corrvec_pool(Animvec_pool=="Alfa"),'type','spearman')
[cval_beto, pval_beto] = corr(distvec_pool(Animvec_pool=="Beto"), corrvec_pool(Animvec_pool=="Beto"),'type','spearman')
%%
[cval_alfa, pval_alfa] = corr(distvec_pool((Animvec_pool=="Alfa")&(distvec_pool>0)), corrvec_pool((Animvec_pool=="Alfa")&(distvec_pool>0)),'type','spearman')
[cval_beto, pval_beto] = corr(distvec_pool((Animvec_pool=="Beto")&(distvec_pool>0)), corrvec_pool((Animvec_pool=="Beto")&(distvec_pool>0)),'type','spearman')

