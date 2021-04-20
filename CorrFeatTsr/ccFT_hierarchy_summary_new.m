%% cc_hierarchy_summary, newer functional version
%  Produce the Higher order summary plot for correlation across hierarchy (WiP)
global sumdir
moviedir = "E:\OneDrive - Washington University in St. Louis\Evol_Manif_Movies"; 
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
sumdir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy\summary";
%% Collect the stats data into struct
Animal = "Beto";
ccHierStats = repmat(struct("Evol",[],"Manif",[]),1,45);
tic
for Expi = 1:45
ExpType = "Manif";
hierdata = load(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t",...
    "poscorr_vox_num", "poscorr_vox_prct", "negcorr_vox_num", "negcorr_vox_prct");
% hierdata.corr_vox_prct = hierdata.corr_vox_num./hierdata.totl_vox_num;
ccHierStats(Expi).Manif = hierdata;
ExpType = "Evol";
hierdata = load(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t",...
    "poscorr_vox_num", "poscorr_vox_prct", "negcorr_vox_num", "negcorr_vox_prct");
% hierdata.corr_vox_prct = hierdata.corr_vox_num./hierdata.totl_vox_num;
ccHierStats(Expi).Evol = hierdata;
end
toc
%%
save(fullfile(hier_savedir,Animal+"_ccHierStats.mat"),"ccHierStats")
%% Create the combined ccHierStats
ccHierStats_cmb = [];
for Animal = ["Alfa","Beto"]
load(fullfile(hier_savedir,Animal+"_ccHierStats.mat"),"ccHierStats")
ccHierStats_cmb = [ccHierStats_cmb, ccHierStats];
end
save(fullfile(hier_savedir,"Both"+"_ccHierStats.mat"),"ccHierStats_cmb")
%% Load the hierarchy Stats
Animal = "Alfa";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
sumdir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy\summary";

load(fullfile(hier_savedir,Animal+"_ccHierStats.mat"),"ccHierStats")
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
% Create prefchan array and the time window vector and masks
prefchan_arr = arrayfun(@(S)S.units.pref_chan,EStats);
V1msk = prefchan_arr <=48 & prefchan_arr>=33;
V4msk = prefchan_arr <=64 & prefchan_arr>=49;
ITmsk = prefchan_arr <=32 & prefchan_arr>= 1;
wdw_vect = [1, 20] + 10 * [0:18]';
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; 
%% Plot Correlated number of voxels across layers with the final time window.
fi=24; 
layernames = ccHierStats(1).Evol.layernames;
EvolCorrVoxPrct = arrayfun(@(H)H.Evol.poscorr_vox_prct(fi,:),ccHierStats,'Uni',0);
EvolCorrVoxPrct = cell2mat(EvolCorrVoxPrct');
ManifCorrVoxPrct = arrayfun(@(H)H.Manif.poscorr_vox_prct(fi,:),ccHierStats,'Uni',0);
ManifCorrVoxPrct = cell2mat(ManifCorrVoxPrct');
xjit = randn(numel(ccHierStats),1)*0.2;
figure;hold on
for i=1:6
    scatter(i+xjit, ManifCorrVoxPrct(:,i))
end
plot([1:6]'+xjit',ManifCorrVoxPrct','color',[0,0,0,0.1])
xticks(1:6); xticklabels(layernames)
% ManifCorrVoxPrct = arrayfun(@(H)H.Manif.corr_vox_prct(fi,:),ccHierStats,'Uni',0);
% ManifCorrVoxPrct = cell2mat(ManifCorrVoxPrct');
%% Compact functional version of plot
%% Plot the dynamics of correlated voxel number across time. 
ExpType = "Evol";varnm = "poscorr_vox_prct";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm),true);
%
ExpType = "Manif";varnm = "poscorr_vox_prct";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm),true);

%% Plot the dynamics of Median Correlation across time. 
ExpType = "Evol";varnm = "med_pos_cc";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Median Corr for Pos Correlated Voxel Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Median Corr for Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm),true);
ExpType = "Manif";varnm = "med_pos_cc";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Median Corr for Pos Correlated Voxel Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
    compose("%s %s Median Corr for Pos Correlated Voxel Across Layers",Animal,ExpType),...
    compose("%s_%s_%s",Animal,ExpType,varnm),true);

%% Gross Hierarchy Summary across exps and layers
varnm = "poscorr_vox_prct";
[h,T] = distr_summary(ccHierStats,varnm,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel",Animal),...
    compose("%s_%s",Animal,varnm));
varnm = "med_pos_cc";
[h,T] = distr_summary(ccHierStats,varnm,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    layernames,compose("%s Median Corr for Pos Correlated Voxel",Animal),...
    compose("%s_%s",Animal,varnm));
%%
varnm = "poscorr_vox_prct";
[h,T] = scatter_summary(ccHierStats,varnm,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel",Animal),...
    compose("%s_%s",Animal,varnm));
%%
Animal = "Both";
tabA = readtable(fullfile(mat_dir, "Alfa"+"_EvolTrajStats.csv"));
tabB = readtable(fullfile(mat_dir, "Beto"+"_EvolTrajStats.csv"));
ExpTab_cmb = cat(1, tabA, tabB); % load hte table for Evol Traject successfulness
sucs_msk = (ExpTab_cmb.t_p_initend<1E-3)&(ExpTab_cmb.t_p_initmax<1E-3);
V1msk = ExpTab_cmb.pref_chan <=48 & ExpTab_cmb.pref_chan>=33;
V4msk = ExpTab_cmb.pref_chan <=64 & ExpTab_cmb.pref_chan>=49;
ITmsk = ExpTab_cmb.pref_chan <=32 & ExpTab_cmb.pref_chan>= 1;
Alfamsk = ExpTab_cmb.Animal=="Alfa"; Betomsk = ExpTab_cmb.Animal=="Beto";
%%
varnm = "poscorr_vox_prct";
[h,T] = distr_summary(ccHierStats_cmb,varnm,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel",Animal),...
    compose("%s_%s",Animal,varnm));
varnm = "med_pos_cc";
[h,T] = distr_summary(ccHierStats_cmb,varnm,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    layernames,compose("%s Median Corr for Pos Correlated Voxel",Animal),...
    compose("%s_%s",Animal,varnm));
%%
varnm = "poscorr_vox_prct";
[h,T] = scatter_summary(ccHierStats_cmb,varnm,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel",Animal),...
    compose("%s_%s",Animal,varnm));
varnm = "med_pos_cc";
[h,T] = scatter_summary(ccHierStats_cmb,varnm,{V1msk,V4msk,ITmsk},["V1","V4","IT"],...
    layernames,compose("%s Median Corr for Pos Correlated Voxel",Animal),...
    compose("%s_%s",Animal,varnm));
%%
varnm = "poscorr_vox_prct";
[h,T] = scatter_summary(ccHierStats_cmb,varnm,{sucs_msk&V1msk,sucs_msk&V4msk,sucs_msk&ITmsk},["V1 Success","V4 Success","IT Success"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel, success Evol Exp",Animal),...
    compose("%s_%s_success",Animal,varnm));
[h,T] = distr_summary(ccHierStats_cmb,varnm,{sucs_msk&V1msk,sucs_msk&V4msk,sucs_msk&ITmsk},["V1 Success","V4 Success","IT Success"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel, success Evol Exp",Animal),...
    compose("%s_%s_success",Animal,varnm));
varnm = "med_pos_cc";
[h,T] = scatter_summary(ccHierStats_cmb,varnm,{sucs_msk&V1msk,sucs_msk&V4msk,sucs_msk&ITmsk},["V1 Success","V4 Success","IT Success"],...
    layernames,compose("%s Median Corr for Pos Correlated Voxel, success Evol Exp",Animal),...
    compose("%s_%s_success",Animal,varnm));
%%
varnm = "poscorr_vox_prct";
[h,T] = scatter_summary(ccHierStats_cmb,varnm,{sucs_msk},["All Success"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel, success Evol Exp",Animal),...
    compose("%s_%s_Allsuccess",Animal,varnm));
[h,T] = distr_summary(ccHierStats_cmb,varnm,{sucs_msk},["All Success"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel, success Evol Exp",Animal),...
    compose("%s_%s_Allsuccess",Animal,varnm));
%%
varnm = "poscorr_vox_prct";
[h,T] = scatter_summary(ccHierStats_cmb,varnm,{sucs_msk,~sucs_msk},["All Success","Not Success"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel, success Evol Exp",Animal),...
    compose("%s_%s_Allsuccess_cmp",Animal,varnm));
[h,T] = distr_summary(ccHierStats_cmb,varnm,{sucs_msk,~sucs_msk},["All Success","Not Success"],...
    layernames,compose("%s Prctage of Pos Correlated Voxel, success Evol Exp",Animal),...
    compose("%s_%s_Allsuccess_cmp",Animal,varnm));
%%
varnm = "med_pos_cc";
[h,T] = scatter_summary(ccHierStats_cmb,varnm,{sucs_msk,~sucs_msk},["All Success","Not Success"],...
    layernames,compose("%s Median Corr for Pos Correlated Voxel, success Evol Exp",Animal),...
    compose("%s_%s_Allsuccess_cmp",Animal,varnm));
[h,T] = distr_summary(ccHierStats_cmb,varnm,{sucs_msk,~sucs_msk},["All Success","Not Success"],...
    layernames,compose("%s Median Corr for Pos Correlated Voxel, success Evol Exp",Animal),...
    compose("%s_%s_Allsuccess_cmp",Animal,varnm));
%%
ExpType = "Evol";varnm = "poscorr_vox_prct";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats_cmb,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm),true);
[h,T]=ccDynam_plot(StatPrct,{Alfamsk,Betomsk},["Alfa","Beto"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s_animsep",Animal,ExpType,varnm));
%
ExpType = "Manif";varnm = "poscorr_vox_prct";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats_cmb,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s",Animal,ExpType,varnm),true);
[h,T]=ccDynam_plot(StatPrct,{Alfamsk,Betomsk},["Alfa","Beto"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s_animsep",Animal,ExpType,varnm));
%%
ExpType = "Evol";varnm = "poscorr_vox_prct";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats_cmb,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk&sucs_msk,V4msk&sucs_msk,ITmsk&sucs_msk},["V1 success","V4 success","IT success"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s_succ",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{Alfamsk&sucs_msk,Betomsk&sucs_msk},["Alfa success","Beto success"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s_animsep_succ",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{sucs_msk,~sucs_msk},["success","Not success"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s_succ_cmp",Animal,ExpType,varnm));

ExpType = "Manif";varnm = "poscorr_vox_prct";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats_cmb,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk&sucs_msk,V4msk&sucs_msk,ITmsk&sucs_msk},["V1 success","V4 success","IT success"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s_succ",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{Alfamsk&sucs_msk,Betomsk&sucs_msk},["Alfa success","Beto success"],layernames,...
	compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
	compose("%s_%s_%s_animsep_succ",Animal,ExpType,varnm));
[h,T]=ccDynam_plot(StatPrct,{sucs_msk,~sucs_msk},["success","Not success"],layernames,...
    compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType),...
    compose("%s_%s_%s_succ_cmp",Animal,ExpType,varnm));
%%
function [h,T]=distr_summary(ccHierStats,varnm,msks,labels,layernames,titstr,savestr)
global sumdir
EvolCorrVoxPrct = arrayfun(@(H)H.Evol.(varnm)(end,:),ccHierStats,'Uni',false);
EvolCorrVoxPrct = cell2mat(EvolCorrVoxPrct');
ManifCorrVoxPrct = arrayfun(@(H)H.Manif.(varnm)(end,:),ccHierStats,'Uni',false);
ManifCorrVoxPrct = cell2mat(ManifCorrVoxPrct');
h = figure();clf;set(h,'pos',[1000, 469, 320*numel(msks), 510]);
T = tiledlayout(1,numel(msks),'Padd','none','TileSp','compact');
for i = 1:numel(msks)
nexttile(T,i);hold on %subtightplot(1,3,1);
msk = msks{i};
if isempty(msk),msk = ones(size(EvolCorrVoxPrct,1),1,'logical');end
evol_m = nanmean(EvolCorrVoxPrct(msk,:),1);
evol_s = sem(EvolCorrVoxPrct(msk,:),1);
manif_m = nanmean(ManifCorrVoxPrct(msk,:),1);
manif_s = sem(ManifCorrVoxPrct(msk,:),1);
b = bar([evol_m;manif_m]');
x = [];for j = 1:numel(b), x = [x ; b(j).XEndPoints];end
errorbar(x',[evol_m;manif_m]',[evol_s;manif_s]','k','linestyle','none')'; % Plot the errorbars

xticks(1:numel(layernames));xticklabels(layernames);%ylim(YLIM)
legend(["Evol","Manif"])
title(compose("%s exp(n=%d)",labels(i),sum(msk)))
end
title(T,titstr+" Across VGG Layers")
linkaxes(T.Children,'y')
savefn = compose("%s_alllayer_gross",savestr);
saveas(h,fullfile(sumdir,savefn+".jpg"))
saveas(h,fullfile(sumdir,savefn+".pdf"))
savefig(h,fullfile(sumdir,savefn+".fig"))
end

function [h,T]=scatter_summary(ccHierStats,varnm,msks,labels,layernames,titstr,savestr)
global sumdir
EvolCorrVoxPrct = arrayfun(@(H)H.Evol.(varnm)(end,:),ccHierStats,'Uni',false);
EvolCorrVoxPrct = cell2mat(EvolCorrVoxPrct');
ManifCorrVoxPrct = arrayfun(@(H)H.Manif.(varnm)(end,:),ccHierStats,'Uni',false);
ManifCorrVoxPrct = cell2mat(ManifCorrVoxPrct');
nLayer = size(EvolCorrVoxPrct,2);
nTot = size(EvolCorrVoxPrct,1);
xjit = 0.1*randn(nTot,1);
OFFSET = 0.4;
h = figure();clf;set(h,'pos',[1000, 469, 360*numel(msks), 510]);
T = tiledlayout(1,numel(msks),'Padd','none','TileSp','compact');
Cord = colororder();
for i = 1:numel(msks)
nexttile(T,i);hold on %subtightplot(1,3,1);
msk = msks{i};
nExp = sum(msk);
for li = 1:nLayer
scatter(li+xjit(msk),EvolCorrVoxPrct(msk,li),'MarkerEdgeColor',Cord(1,:));
scatter(li+OFFSET+xjit(msk),ManifCorrVoxPrct(msk,li),'MarkerEdgeColor',Cord(2,:));
end
plot([1:nLayer]'+xjit(msk)',EvolCorrVoxPrct(msk,:)','color',[0,0,0,0.1])
plot([1:nLayer]'+xjit(msk)'+OFFSET,ManifCorrVoxPrct(msk,:)','color',[0,0,0,0.1])
% bar([nanmean(EvolCorrVoxPrct(msk,:),1);nanmean(ManifCorrVoxPrct(msk,:),1)]')
% errorbar(1:6,mean(CorrVoxPrct(V1msk,:),1),std(CorrVoxPrct(V1msk,:),0,1),'o')
xticks(1:numel(layernames));xticklabels(layernames);%ylim(YLIM)
legend(["Evol","Manif"])
title(compose("%s exp(n=%d)",labels(i),nExp))
end
title(T,titstr+" Across VGG Layers")
linkaxes(T.Children,'y')
savefn = compose("%s_alllayer_scatter",savestr);
saveas(h,fullfile(sumdir,savefn+".jpg"))
saveas(h,fullfile(sumdir,savefn+".pdf"))
savefig(h,fullfile(sumdir,savefn+".fig"))
end

function [h,T]=ccDynam_plot(StatPrct,msks,labels,layernames,titstr,savestr,doerror)
if nargin < 7, doerror = false; end
global sumdir
nLayer = size(StatPrct,2);
nWdw = size(StatPrct,1);
% layernames = ccHierStats(1).Evol.layernames;
h=figure; T = tiledlayout(1,numel(msks),'Padd','none','TileSp','compact');
set(h,'pos',[526,   462,   420*numel(msks)+120,   420])
for i = 1:numel(msks)
nexttile(T,i);
msk = msks{i};
VarTrace_mean = nanmean(StatPrct(:,:,msk),3);
VarTrace_sem = sem(StatPrct(:,:,msk),3);
if ~ doerror
plot([1:19,NaN,20:23,NaN,24],[VarTrace_mean(1:19,:);nan(1,nLayer);VarTrace_mean(20:23,:);nan(1,nLayer);VarTrace_mean(24,:)]...
    , "-o","LineWidth",1.5)
xlabel("Time Windows")
box off;
else
hold on
% meantrmat = [VarTrace_mean(1:19,:);nan(1,nLayer);VarTrace_mean(20:23,:);nan(1,nLayer);VarTrace_mean(24,:)];
% semtrmat = [VarTrace_sem(1:19,:);nan(1,nLayer);VarTrace_sem(20:23,:);nan(1,nLayer);VarTrace_sem(24,:)];
% errorbar(repmat([1:19,NaN,20:23,NaN,24]',1,nLayer),meantrmat,semtrmat...
%     , "-o","LineWidth",1.5)
Cord = colororder();
for li = 1:nLayer
shadedErrorBar([1:19],VarTrace_mean(1:19,li),VarTrace_sem(1:19,li), ...
    'LineProp', {"-o","LineWidth",1.5,'Color',Cord(li,:)})
shadedErrorBar([20:23],VarTrace_mean(20:23,li),VarTrace_sem(20:23,li), ...
    'LineProp', {"-o","LineWidth",1.5,'Color',Cord(li,:),'HandleVisibility','off'})
errorbar([24],VarTrace_mean(24,li),VarTrace_sem(24,li), "-o","LineWidth",1.5,'Color',Cord(li,:),'HandleVisibility','off')
end
xlabel("Time Windows")
end
YLIM = ylim();ylim([max(0,YLIM(1)),YLIM(2)]);
xlim([0.5,nWdw+0.5]);
legend(layernames,'Location','Best')
title(compose("%s exp(n=%d)",labels(i),sum(msk)))
end
title(T,titstr)
linkaxes(T.Children,'y')
savefn = compose("%s_cctraj",savestr);
if doerror, savefn = savefn+"err"; end
savefig(h,fullfile(sumdir, savefn+".fig")) 
saveas(h,fullfile(sumdir, savefn+".png"))
saveas(h,fullfile(sumdir, savefn+".pdf"))
end
