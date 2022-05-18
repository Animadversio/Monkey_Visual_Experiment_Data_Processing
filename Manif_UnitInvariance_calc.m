% reshape(ReprStat.stim.imgname_uniq,6,10)
function calc_invariance_fun(ReprStat)
%%
respmeantsr = reshape(ReprStat.resp.meanMat,[],6,10);
%%
ccmat_tsr = [];
for iCh = 1:size(respmeantsr,1)
ccmat = corr(squeeze(respmeantsr(iCh,:,:))');
ccmat_tsr(iCh,:,:) = ccmat;
end
%%
mask = reshape(diag(nan(1,6)),1,6,6);
%%
ccmat_tsrmask = ccmat_tsr+mask;
unit_invariance = nanmean(ccmat_tsrmask,[2,3]);
area_vec = area_map(ReprStat.units.spikeID,"Beto_new");
areaidx_vec = area2index(area_vec);
[cc,pval] = corr(areaidx_vec, unit_invariance,'type','Spearman');
%%
figure('pos',[680   251   400   540])
boxplot(unit_invariance,area_vec,'GroupOrder',["V1","V4","IT"])
ylabel("Unit Invariance (corr)")
title([ReprStat.meta.fdrnm, "Single Neuron Invariance ~ Area", compose("Spearman cc:%.3f(P=%.1e)",cc,pval)])
saveallform(ReprStat.meta.figdir, "UnitInvariance_summary_bar")
%%
[p,tbl,Fstats] = anova1(unit_invariance,area_vec,'off');
%%
save(fullfile(ReprStat.meta.figdir, "InvarianceStats.mat"), 'unit_invariance', 'area_vec', 'Fstats')
end

function areaidx = area2index(area)
areaidx = nan(size(area));
areaidx(area=="V1") = 1;
areaidx(area=="V4") = 2;
areaidx(area=="IT") = 3;
end