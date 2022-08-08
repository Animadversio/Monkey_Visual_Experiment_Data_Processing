% reshape(ReprStat.stim.imgname_uniq,6,10)
function calc_sparseness_fun(ReprStats)
for i = 1:numel(ReprStats)
ReprStat = ReprStats(i);
%%
area_vec = area_map(ReprStat.units.spikeID,"Beto_new");
areaidx_vec = area2index(area_vec);
for area = ["V1","V4","IT"]
    areamsk.(area) = area_vec == area;
end
Fmask = [ReprStat.stat.Fstats.F_P]' < 0.005;
stimstr = compose("Image pos [%.1f %.1f] size %.1f deg",...
    ReprStat.stim.impos(1),ReprStat.stim.impos(2),ReprStat.stim.imsize_deg);
%%
% for population correlation, we need to normalize cells properly! 
% respMat_bslnorm = ReprStat.resp.meanMat ./ ReprStat.resp.bslmean;
% evkmat = ReprStat.resp.meanMat - ReprStat.resp.bslmean;
% evkmat = ReprStat.resp.meanMat ./ ReprStat.resp.bslmean;
evkmat = ReprStat.resp.meanMat;
sparseness_vec = (1 - mean(evkmat,2).^2 ./ mean(evkmat.^2,2))/(1 - 1/size(evkmat,2));
%%
[cval,pval] = corr(sparseness_vec(Fmask),areaidx_vec(Fmask),'type','Spearman');
corrstat = struct("corr",cval,"pval",pval);
%%
figure('pos',[680   251   400   540])
% boxplot(sparseness_vec,area_vec,'GroupOrder',["V1","V4","IT"])
boxplot(sparseness_vec(Fmask),area_vec(Fmask),'GroupOrder',["V1","V4","IT"])
ylabel("Sparseness coef")
title([ReprStat.meta.fdrnm, "Single Neuron Sparseness ~ Area", ...
    compose("Spearman cc:%.3f(P=%.1e) N=%d",cval,pval,sum(Fmask)),stimstr])
saveallform(ReprStat.meta.figdir, "Sparseness_summary_Fsig_bar")
%%
[p,anova_tab,Fstats] = anova1(sparseness_vec(Fmask),area_vec(Fmask),'off');
disp(anova_tab)
%%
save(fullfile(ReprStat.meta.figdir, "SparsenessStats.mat"), ...
    'sparseness_vec', 'area_vec', 'anova_tab','Fstats','corrstat')
end
end

function areaidx = area2index(area)
areaidx = nan(size(area));
areaidx(area=="V1") = 1;
areaidx(area=="V4") = 2;
areaidx(area=="IT") = 3;
end