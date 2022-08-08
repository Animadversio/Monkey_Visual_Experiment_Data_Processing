% reshape(ReprStat.stim.imgname_uniq,6,10)
function calc_invariance_fun(ReprStats)
for i = 1:numel(ReprStats)
ReprStat = ReprStats(i);
%%
area_vec = area_map(ReprStat.units.spikeID,"Beto_new");
areaidx_vec = area2index(area_vec);
for area = ["V1","V4","IT"]
    areamsk.(area) = area_vec == area;
end
%% Reshape the mean response to be a tensor 
NObj = 10;
NTfm = 6;
respmeantsr = reshape(ReprStat.resp.meanMat,[],NTfm,NObj);
%%
% for population correlation, we need to normalize cells properly! 
% respMat_bslnorm = ReprStat.resp.meanMat ./ ReprStat.resp.bslmean;
Fmask = [ReprStat.stat.Fstats.F_P]' < 0.005;
respMat_bslnorm = zscore(ReprStat.resp.meanMat,1,2);
respTsr_bslnorm = reshape(respMat_bslnorm,[],NTfm,NObj);
ccmat_obj = struct();
ccmat_obj_mean = struct();
for area = ["V1","V4","IT"]
ccmat_obj.(area) = [];
for iObj = 1:NObj
    ccmat_obj.(area)(iObj,:,:) = corr(squeeze(respTsr_bslnorm(areamsk.(area)&Fmask,:,iObj)));
end
objmask = reshape(diag(nan(1,NTfm)),1,NTfm,NTfm);
maskedccmat = ccmat_obj.(area) + objmask;
cc_per_obj = nanmean(maskedccmat,[2,3]);
ccmat_obj_mean.(area) = cc_per_obj;
fprintf("Area %s %.3f+-%.3f  N=%d\n ",...
    area,mean(cc_per_obj),sem(cc_per_obj),sum(areamsk.(area)&Fmask))
disp(compose("%.2f",cc_per_obj)')
end
%%
ccmat_tsr = [];
for iCh = 1:size(respmeantsr,1)
ccmat = corr(squeeze(respmeantsr(iCh,:,:))');
ccmat_tsr(iCh,:,:) = ccmat;
end
mask = reshape(diag(nan(1,6)),1,6,6);
ccmat_tsrmask = ccmat_tsr+mask;
unit_invariance = nanmean(ccmat_tsrmask,[2,3]);
%%
[cval,pval] = corr(areaidx_vec, unit_invariance,'type','Spearman','Rows','complete');
corrstat = struct("corr",cval,"pval",pval);
%%
stimstr = compose("Image pos [%.1f %.1f] size %.1f deg",...
    ReprStat.stim.impos(1),ReprStat.stim.impos(2),ReprStat.stim.imsize_deg);
figure('pos',[680   251   400   540])
boxplot(unit_invariance,area_vec,'GroupOrder',["V1","V4","IT"])
ylabel("Unit Invariance (corr)")
title([ReprStat.meta.fdrnm, "Single Neuron Invariance ~ Area", ...
    compose("Spearman cc:%.3f(P=%.1e)",cval,pval),stimstr])
saveallform(ReprStat.meta.figdir, "UnitInvariance_summary_bar")
%%
[p,anova_tab,Fstats] = anova1(unit_invariance,area_vec,'off');
%%
save(fullfile(ReprStat.meta.figdir, "InvarianceStats.mat"), ...
    'unit_invariance', 'area_vec', 'anova_tab','Fstats','corrstat','ccmat_obj','ccmat_obj_mean')
end
end

function areaidx = area2index(area)
areaidx = nan(size(area));
areaidx(area=="V1") = 1;
areaidx(area=="V4") = 2;
areaidx(area=="IT") = 3;
end