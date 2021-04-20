%% Tuning map Spatial Relationship
%  characterize the relationship of tuning landscapes between the neurons
%  in and across the array. 
Animal = "Alfa";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir, Animal+"_ManifMapVarStats.mat"),'MapVarStats');
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
%%
%% Integration Weight Matrix
%  Integrate the area of the patch that each node governs
[phi_grid, theta_grid] = meshgrid(-90:18:90, -90:18:90); 
global Wgrid
phi1_grid = max(phi_grid - 9, -90) /180 *pi;
phi2_grid = min(phi_grid + 9,  90) /180 *pi;
theta1_grid = max(theta_grid - 9,  -90) /180 *pi;
theta2_grid = min(theta_grid + 9,   90) /180 *pi;
Wgrid = abs(sin(phi2_grid) - sin(phi1_grid)).*(theta2_grid - theta1_grid);
%%
CortiDisCorr = repmat(struct(),1,numel(MapVarStats));
%%
figdir = "E:\OneDrive - Washington University in St. Louis\CortiDistCorr";
mkdir(figdir)
% diary(fullfile(figdir,Animal+"_cortdistcorr.log"))
F.plot_corrmat = false;
F.plot_corrcorr = false;
for Expi = 1:numel(EStats)
tic
preflab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
fprintf("Processing %s Exp%d prefchan %s\n",Animal,Expi,preflab)
CortiDisCorr(Expi).units = MapVarStats(Expi).units;
CortiDisCorr(Expi).meta = MapVarStats(Expi).meta;
CortiDisCorr(Expi).Expi = MapVarStats(Expi).Expi;
CortiDisCorr(Expi).Animal = MapVarStats(Expi).Animal;
si = 1;
% Get the FStatistics
FStats = [];
for iCh = 1:numel(MapVarStats(Expi).units.spikeID)
actmap_col = cellfun(@(A)...
    squeeze(A(iCh,:)),MapVarStats(Expi).manif.act_col{si},'uni',0);
FStat = anova_cells(actmap_col);
FStats = [FStats, FStat];
end
CortiDisCorr(Expi).FStats = FStats;
FSigMsk = struct2table(FStats).F_P<1E-2;
% Single trial activity correlation
act_col_vec = reshape(MapVarStats(Expi).manif.act_col{si},1,[]);
actmat = cell2mat(act_col_vec);
sgtr_corrmat = corrcoef(actmat');
% Residue, noise correlation
res_col = cellfun(@(A)A-mean(A), MapVarStats(Expi).manif.act_col{si}, 'uni',0);
resmat = cell2mat(reshape(res_col, 1, []));
res_corrmat = corrcoef(resmat');
% Avg activity correlation
avgact_col = cellfun(@(A)mean(A,2),MapVarStats(Expi).manif.act_col{si},'uni',0);
acgact_tsr = cell2mat(reshape(avgact_col,1,11,11));
corrmat1d = map_corr(acgact_tsr, "corr1d");
corrmat_sph = map_corr(acgact_tsr, "corr_sph");
% corrmat2d = map_corr(acgact_tsr);
% [cval,pval] = corr(corrmat2d(:),cortexDistmat(:),'row','complete');
% fprintf("Tr Avg firing rate (signal) correlation, corr %.3f (%.1e)\n",cval,pval)
%%
[cortexDistmat, V1msk, V4msk, ITmsk] = spikeID2cortexDist(MapVarStats(Expi).units.spikeID, ...
    MapVarStats(Expi).units.unit_num_arr);
CortiDisCorr(Expi).sgtr_corrmat = sgtr_corrmat;
CortiDisCorr(Expi).res_corrmat = res_corrmat;
CortiDisCorr(Expi).avg_corrmat = corrmat1d;
CortiDisCorr(Expi).avgsph_corrmat = corrmat_sph;
CortiDisCorr(Expi).cortexDismat = cortexDistmat;
if F.plot_corrmat
figure(1);clf;
subtightplot(1,5,1);imagesc(sgtr_corrmat);axis image;colorbar;title("Single Trial Correlation")
subtightplot(1,5,2);imagesc(res_corrmat);axis image;colorbar;title("Residue Correlation")
subtightplot(1,5,3);imagesc(corrmat1d);axis image;colorbar;title("Average Correlation")
subtightplot(1,5,4);imagesc(corrmat_sph);axis image;colorbar;title("Average Correlation(Sph)")
subtightplot(1,5,5);imagesc(cortexDistmat);axis image;colorbar;title("Cortical Distance")
suptitle(compose("%s Exp%02d manif PrefChan %s",Animal,Expi,preflab))
fignm = compose("%s_Exp%02d_manif_corrmats",Animal,Expi);
saveas(1,fullfile(figdir,fignm+".png"))
saveas(1,fullfile(figdir,fignm+".pdf"))
savefig(1,fullfile(figdir,fignm+".fig"))
end

[cval,pval] = corr(sgtr_corrmat(:),cortexDistmat(:),'row','complete');
fprintf("Single trial firing rate correlation ~ DistMat, corr %.3f (%.1e)\n",cval,pval)
[cval,pval] = corr(res_corrmat(:),cortexDistmat(:),'row','complete');
fprintf("Single trial Residue (noise) correlation ~ DistMat, corr %.3f (%.1e)\n",cval,pval)
[cval,pval] = corr(corrmat1d(:),cortexDistmat(:),'row','complete');
fprintf("Tr Avg firing rate (signal) correlation ~ DistMat, corr %.3f (%.1e)\n",cval,pval)
[cval,pval] = corr(corrmat_sph(:),cortexDistmat(:),'row','complete');
fprintf("Tr Avg firing rate (signal) correlation ~ DistMat, corr %.3f (%.1e)\n",cval,pval)
%%
fprintf("Correlation of Noise Correlation and cortical distance\n")
R = summary_mat_corr(res_corrmat, cortexDistmat, [], {[],V1msk, V4msk, ITmsk}, ["all","V1area", "V4area", "ITarea"], ...
    ["all","V1", "V4", "IT"]);
CortiDisCorr(Expi).res = R;
%%
fprintf("Correlation of Single Trial Correlation and cortical distance\n")
S = summary_mat_corr(sgtr_corrmat, cortexDistmat, [], {[],V1msk, V4msk, ITmsk}, ["all","V1area", "V4area", "ITarea"], ...
    ["all","V1", "V4", "IT"]);
CortiDisCorr(Expi).sgtr = S;
%%
C = [];
fprintf("Correlation of tuning landscape and cortical distance\n")
C = summary_mat_corr(corrmat1d, cortexDistmat, [], {[], V1msk, V4msk, ITmsk}, ["All", "V1area", "V4area", "ITarea"], ...
    ["all","V1", "V4", "IT"], C);
C = summary_mat_corr(corrmat1d, cortexDistmat, (cortexDistmat>0), {[], V1msk, V4msk, ITmsk}, ...
    ["DiffChan", "V1DiffChan", "V4DiffChan", "ITDiffChan"], ...
    ["allDfCh","V1DfCh", "V4DfCh", "ITDfCh"], C);
C = summary_mat_corr(corrmat1d, cortexDistmat, FSigMsk, {[], V1msk, V4msk, ITmsk}, ...
    ["Fsignif", "V1area Fsignif", "V4area Fsignif", "ITarea Fsignif"], ...
    ["all_F","V1_F", "V4_F", "IT_F"], C);
C = summary_mat_corr(corrmat1d, cortexDistmat, FSigMsk&(cortexDistmat>0)&FSigMsk', ...
    {[], V1msk, V4msk, ITmsk}, ...
    ["DiffChan Fsignif","V1DiffChan Fsignif", "V4DiffChan Fsignif", "ITDiffChan Fsignif"], ...
    ["allDfCh_F","V1DfCh_F", "V4DfCh_F", "ITDfCh_F"], C);
CortiDisCorr(Expi).avg = C;
%%
Sp = [];
fprintf("Spherical Correlation of tuning landscape and cortical distance\n")
Sp = summary_mat_corr(corrmat_sph, cortexDistmat, [], {[], V1msk, V4msk, ITmsk}, ["All", "V1area", "V4area", "ITarea"], ...
    ["all","V1", "V4", "IT"], Sp);
Sp = summary_mat_corr(corrmat_sph, cortexDistmat, (cortexDistmat>0), {[], V1msk, V4msk, ITmsk}, ...
    ["DiffChan", "V1DiffChan", "V4DiffChan", "ITDiffChan"], ...
    ["allDfCh","V1DfCh", "V4DfCh", "ITDfCh"], Sp);
Sp = summary_mat_corr(corrmat_sph, cortexDistmat, FSigMsk, {[], V1msk, V4msk, ITmsk}, ...
    ["Fsignif", "V1area Fsignif", "V4area Fsignif", "ITarea Fsignif"], ...
    ["all_F","V1_F", "V4_F", "IT_F"], Sp);
Sp = summary_mat_corr(corrmat_sph, cortexDistmat, FSigMsk&(cortexDistmat>0)&FSigMsk', ...
    {[], V1msk, V4msk, ITmsk}, ...
    ["DiffChan Fsignif","V1DiffChan Fsignif", "V4DiffChan Fsignif", "ITDiffChan Fsignif"], ...
    ["allDfCh_F","V1DfCh_F", "V4DfCh_F", "ITDfCh_F"], Sp);
CortiDisCorr(Expi).avgsph = Sp;

%%
CC = [];
CC = summary_corrmat(corrmat1d, [], {[], V1msk, V4msk, ITmsk}, ["All", "V1area", "V4area", "ITarea"], ...
    ["cc_all","cc_V1", "cc_V4", "cc_IT"], CC);
CC = summary_corrmat(corrmat1d, FSigMsk, {[], V1msk, V4msk, ITmsk}, ["All", "V1area", "V4area", "ITarea"], ...
    ["cc_all_F","cc_V1_F", "cc_V4_F", "cc_IT_F"], CC);
CC = summary_corrmat(sgtr_corrmat, [], {[], V1msk, V4msk, ITmsk}, ["All", "V1area", "V4area", "ITarea"], ...
    ["sgtr_cc_all","sgtr_cc_V1", "sgtr_cc_V4", "sgtr_cc_IT"], CC);
CC = summary_corrmat(sgtr_corrmat, FSigMsk, {[], V1msk, V4msk, ITmsk}, ["All", "V1area", "V4area", "ITarea"], ...
    ["sgtr_cc_all_F","sgtr_cc_V1_F", "sgtr_cc_V4_F", "sgtr_cc_IT_F"], CC);
CC = summary_corrmat(res_corrmat, [], {[], V1msk, V4msk, ITmsk}, ["All", "V1area", "V4area", "ITarea"], ...
    ["res_cc_all","res_cc_V1", "res_cc_V4", "res_cc_IT"], CC);
CC = summary_corrmat(res_corrmat, FSigMsk, {[], V1msk, V4msk, ITmsk}, ["All", "V1area", "V4area", "ITarea"], ...
    ["res_cc_all_F","res_cc_V1_F", "res_cc_V4_F", "res_cc_IT_F"], CC);
CortiDisCorr(Expi).cc = CC;

if F.plot_corrcorr
h=figure(2);clf;
ax1 = subtightplot(1,3,1);
addScatterWmask(ax1, corrmat1d, cortexDistmat, [], {[], V1msk, V4msk, ITmsk}, ["all","V1", "V4", "IT"])
title(ax1, "All pairs")
ax2 = subtightplot(1,3,2);
addScatterWmask(ax2, corrmat1d, cortexDistmat, FSigMsk, {[], V1msk, V4msk, ITmsk}, ...
    ["Fsignif", "V1area Fsignif", "V4area Fsignif", "ITarea Fsignif"])
title(ax2, "All F significant p<0.01")
ax3 = subtightplot(1,3,3);
addScatterWmask(ax3, corrmat1d, cortexDistmat, FSigMsk&(cortexDistmat>0)&FSigMsk', ...
    {[], V1msk, V4msk, ITmsk}, ["DiffChan Fsignif","V1DiffChan Fsignif", "V4DiffChan Fsignif", "ITDiffChan Fsignif"])
title(ax3, "All Pairs with different channel and F signif")
suptitle(compose("%s Exp%02d manif PrefChan %s",Animal,Expi,preflab))
fignm = compose("%s_Exp%02d_manif_corrDistScatter",Animal,Expi);
saveas(2,fullfile(figdir,fignm+".png"))
saveas(2,fullfile(figdir,fignm+".pdf"))
savefig(2,fullfile(figdir,fignm+".fig"))
end
toc
end
diary off
%%
save(fullfile(mat_dir, Animal+'_CortiDisCorr.mat'), 'CortiDisCorr') 
%%

function S = summary_corrmat(corrmat, commonmsk, masks, labels, entries, S)
if nargin < 6, S = struct(); end
if nargin < 5, entries = labels; end
if size(commonmsk,2)==1 || size(commonmsk,1)==1
    commonmsk = reshape(commonmsk,[],1) & reshape(commonmsk,1,[]) & tril(ones(size(corrmat),'logical'),-1);
elseif isempty(commonmsk)
    commonmsk = tril(ones(size(corrmat),'logical'),-1);%ones(size(corrmat),'logical');
else,
    commonmsk = commonmsk & tril(ones(size(corrmat),'logical'),-1);
end
for mi = 1:numel(masks)
msk = masks{mi};
if size(msk,2)==1 || size(msk,1)==1
finalmsk = reshape(msk,[],1) & commonmsk & reshape(msk,1,[]);
corrvec = reshape(corrmat(finalmsk),[],1);
elseif isempty(msk)
corrvec = reshape(corrmat(commonmsk),[],1);
else
corrvec = reshape(corrmat(msk&commonmsk),[],1);
end
num = sum(~isnan(corrvec));
S.(entries(mi)) = nanmean(corrvec);
S.(entries(mi)+"_sem") = sem(corrvec);
S.(entries(mi)+"_N") = num;
end
end

function S = summary_mat_corr(corrmat, distmat, commonmsk, masks, labels, entries, S)
if nargin < 7, S = struct(); end
if nargin < 6, entries = labels; end
if size(commonmsk,2)==1 || size(commonmsk,1)==1
    commonmsk = reshape(commonmsk,[],1) & reshape(commonmsk,1,[]) & tril(ones(size(corrmat),'logical'),-1);
elseif isempty(commonmsk)
    commonmsk = tril(ones(size(corrmat),'logical'),-1);%ones(size(corrmat),'logical');
else,
    commonmsk = commonmsk & tril(ones(size(corrmat),'logical'),-1);
end
for mi = 1:numel(masks)
msk = masks{mi};
if size(msk,2)==1 || size(msk,1)==1
finalmsk = reshape(msk,[],1) & commonmsk & reshape(msk,1,[]);
corrvec = reshape(corrmat(finalmsk),[],1);
distvec = reshape(distmat(finalmsk),[],1);
elseif isempty(msk)
corrvec = reshape(corrmat(commonmsk),[],1);
distvec = reshape(distmat(commonmsk),[],1);
else
corrvec = reshape(corrmat(msk&commonmsk),[],1);
distvec = reshape(distmat(msk&commonmsk),[],1);
end
num = sum(~isnan(corrvec)&~isnan(distvec));
if num > 0
[cval,pval] = corr(corrvec, distvec,'row','complete','type','spearman');
S.(entries(mi)) = cval;
S.(entries(mi)+"_P") = pval;
else
S.(entries(mi)) = nan;
S.(entries(mi)+"_P") = nan;
end
S.(entries(mi)+"_df") = num;
fprintf("For %s, CorrMat ~ DistMat, corr %.3f (%.1e)\n",labels{mi},cval,pval)
end
end

function addScatterWmask(ax, corrmat, distmat, commonmsk, masks, labels)
legend_strs = [];
hold on;
if size(commonmsk,2)==1 || size(commonmsk,1)==1
    commonmsk = reshape(commonmsk,[],1) & reshape(commonmsk,1,[]) & tril(ones(size(corrmat),'logical'),-1);
elseif isempty(commonmsk)
    commonmsk = tril(ones(size(corrmat),'logical'),-1);%ones(size(corrmat),'logical');
else,
    commonmsk = commonmsk & tril(ones(size(corrmat),'logical'),-1);
end
for mi = 1:numel(masks)
msk = masks{mi};
if size(msk,2)==1 || size(msk,1)==1
finalmsk = reshape(msk,[],1) & commonmsk & reshape(msk,1,[]);
corrvec = reshape(corrmat(finalmsk),[],1);
distvec = reshape(distmat(finalmsk),[],1);
elseif isempty(msk)
corrvec = reshape(corrmat(commonmsk),[],1);
distvec = reshape(distmat(commonmsk),[],1);
else
corrvec = reshape(corrmat(msk&commonmsk),[],1);
distvec = reshape(distmat(msk&commonmsk),[],1);
end
num = sum(~isnan(corrvec)&~isnan(distvec));
if num > 0
[cval,pval] = corr(corrvec, distvec,'row','complete','type','spearman');
else
cval=nan;pval=nan;
end
legend_str = compose("%s %.2f(%.1e) n=%d",labels(mi),cval,pval,num);
scatter(ax, distvec, corrvec,'DisplayName', legend_str)
legend_strs(mi) = legend_str;
end
xlabel(ax,"cortical distance");
ylabel(ax,"Activation corr");
legend(ax,'Location','best')
end

function distmat = map_corr(maps_tsr,type)
if nargin == 1, type="corr2d"; end
chN = size(maps_tsr, 1);
distmat = zeros(chN,chN);
if type == "corr2d"
for i = 1:chN
    for j = 1:chN
    distmat(i,j) = corr2(squeeze(maps_tsr(i,:,:)),squeeze(maps_tsr(j,:,:)));
    end
end
elseif type == "corr_sph"
for i = 1:chN
    for j = i:chN
    distmat(i,j) = sphere_corr(squeeze(maps_tsr(i,:,:)),squeeze(maps_tsr(j,:,:)));
    distmat(j,i) = sphere_corr(squeeze(maps_tsr(i,:,:)),squeeze(maps_tsr(j,:,:)));
    end
end
elseif type == "corr1d"
actmat = reshape(maps_tsr,[],11*11);
distmat = corrcoef(actmat');
end
end

function corr_sph = sphere_corr(map1, map2)
global Wgrid
Wsum = sum(Wgrid,'all');
M1 = sum(map1.*Wgrid,'all')/Wsum;
M2 = sum(map2.*Wgrid,'all')/Wsum;
innprod = nansum((map1 - M1).*(map2 - M2).*Wgrid,'all');
var1 = nansum((map1 - M1).^2.*Wgrid,'all');
var2 = nansum((map2 - M2).^2.*Wgrid,'all');
corr_sph = innprod / sqrt(var1*var2);
end

function mapIOU

end