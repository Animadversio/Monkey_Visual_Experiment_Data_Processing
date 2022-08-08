NObj = 10;
NTfm = 6;
respmeantsr = reshape(ReprStat.resp.meanMat,[],NTfm,NObj);
%%
F_P_arr = [ReprStat.stat.Fstats.F_P]';
%%
cval_arr = nan(numel(ReprStat.meta.spikeID),1);
pval_arr = nan(numel(ReprStat.meta.spikeID),1);
for iCh = 1:numel(ReprStat.meta.spikeID)
if F_P_arr(iCh) >0.001, continue;
end
% unit = find(ReprStat.meta.spikeID==54);
tunemap = squeeze(respmeantsr(iCh,:,:));
%%
bsl = ReprStat.resp.bslmean(iCh);
tunemap_bsl = tunemap - bsl;
objinv = (sum(tunemap_bsl ./ max(tunemap_bsl,[],1),1)-1)/(6-1);
%%
obj_avgact = mean(tunemap_bsl,1);
obj_centact = tunemap_bsl(4,:);
% figure;scatter(obj_avgact, objinv)
[cval,pval] = corr(obj_avgact', objinv');%,'type','Spearman'
cval_arr(iCh) = cval;
pval_arr(iCh) = pval;
fprintf("%02d corr %.3f (%.1e)\n",iCh,cval,pval)
end
sum(pval_arr < 0.05)/sum(~isnan(pval_arr))
%%