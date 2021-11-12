%% cc_hierarchy_summary
%  Produce the Higher order summary plot for correlation across hierarchy (WiP)
moviedir = "E:\OneDrive - Washington University in St. Louis\Evol_Manif_Movies"; 
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
sumdir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy\summary";
%% Collect the stats data into struct
Animal = "Alfa";
ccHierStats = repmat(struct("Evol",[],"Manif",[]),1,46);
tic
for Expi = 1:46
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
ExpType = "Evol";varnm = "poscorr_vox_prct";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames);
title(T,compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType))
%%
ExpType = "Manif";varnm = "poscorr_vox_prct";
StatPrct = arrayfun(@(H)H.(ExpType).(varnm)(:,:),ccHierStats,'Uni',0);
StatPrct = cell2mat(reshape(StatPrct,1,1,[]));
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames);
title(T,compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType))
%%
[h,T]=ccDynam_plot(StatPrct,{V1msk,V4msk,ITmsk},["V1","V4","IT"],layernames,true);
title(T,compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType))

%%
for fi = 1:24
EvolCorrVoxPrct = arrayfun(@(H)H.Evol.corr_vox_prct(fi,:),ccHierStats,'UniformOutput',false);
EvolCorrVoxPrct = cell2mat(EvolCorrVoxPrct');
ManifCorrVoxPrct = arrayfun(@(H)H.Manif.corr_vox_prct(fi,:),ccHierStats,'UniformOutput',false);
ManifCorrVoxPrct = cell2mat(ManifCorrVoxPrct');
YLIM = [min(EvolCorrVoxPrct,[],'all'), max(EvolCorrVoxPrct,[],'all')];
%
figure(7);clf
subtightplot(1,3,1,0.05,0.07,0.12);hold on
plot(EvolCorrVoxPrct(ITmsk,:)',"-.o","color",[0,0,1]);ylim(YLIM)
xticks([1:6]);xticklabels(["fc8","fc7","fc6","conv5-3","conv4-3","conv3-3"])
title("Evol IT");ylabel("Correlated Voxel Percentage");xlabel("VGG Layer")
subtightplot(1,3,2,0.05,0.07,0.12);hold on
plot(EvolCorrVoxPrct(V4msk,:)',"-.*","color",[1,0,1]);ylim(YLIM)
xticks([1:6]);xticklabels(["fc8","fc7","fc6","conv5-3","conv4-3","conv3-3"])
title("Evol V4");ylabel("Correlated Voxel Percentage");xlabel("VGG Layer")
subtightplot(1,3,3,0.05,0.07,0.12);hold on
plot(EvolCorrVoxPrct(V1msk,:)',"-.^","color",[1,0,0]);ylim(YLIM)
xticks([1:6]);xticklabels(["fc8","fc7","fc6","conv5-3","conv4-3","conv3-3"])
title("Evol V1");ylabel("Correlated Voxel Percentage");xlabel("VGG Layer")
% subtightplot(1,2,2,0.02,0.04,0.05);hold on
% plot(ManifCorrVoxPrct(ITmsk,:)',"-.o","color",[0,0,1])
% plot(ManifCorrVoxPrct(V4msk,:)',"-.*","color",[1,0,1])
% plot(ManifCorrVoxPrct(V1msk,:)',"-.^","color",[1,0,0])
% xticks([1:6]);xticklabels(["fc8","fc7","fc6","conv5-3","conv4-3","conv3-3"])
% title("Manif")
suptitle(compose("%s [%d, %d] ms",Animal,wdw_vect(fi,1),wdw_vect(fi,2)))
pause
end
%% Gross Hierarchy Plot
fi =24;
EvolCorrVoxPrct = arrayfun(@(H)H.Evol.corr_vox_prct(fi,:),ccHierStats,'UniformOutput',false);
EvolCorrVoxPrct = cell2mat(EvolCorrVoxPrct');
ManifCorrVoxPrct = arrayfun(@(H)H.Manif.corr_vox_prct(fi,:),ccHierStats,'UniformOutput',false);
ManifCorrVoxPrct = cell2mat(ManifCorrVoxPrct');
%% Evolution Manifold Combined Plot
YLIM = [0,0.7];
figure(2);clf
subtightplot(1,3,1);hold on 
bar([mean(EvolCorrVoxPrct(V1msk,:),1);mean(ManifCorrVoxPrct(V1msk,:),1)]')
% errorbar(1:6,mean(CorrVoxPrct(V1msk,:),1),std(CorrVoxPrct(V1msk,:),0,1),'o')
xticks(1:6);xticklabels(layernames);ylim(YLIM)
legend(["Evol","Manif"])
title(compose("V1 exp(n=%d)",sum(V1msk)))
subtightplot(1,3,2);hold on 
bar([mean(EvolCorrVoxPrct(V4msk,:),1);mean(ManifCorrVoxPrct(V4msk,:),1)]')
% errorbar(1:6,mean(CorrVoxPrct(V4msk,:),1),std(CorrVoxPrct(V4msk,:),0,1),'o')
legend(["Evol","Manif"])
xticks(1:6);xticklabels(layernames);ylim(YLIM)
title(compose("V4 exp(n=%d)",sum(V4msk)))
subtightplot(1,3,3);hold on 
bar([mean(EvolCorrVoxPrct(ITmsk,:),1);mean(ManifCorrVoxPrct(ITmsk,:),1)]')
% errorbar(1:6,mean(CorrVoxPrct(ITmsk,:),1),std(CorrVoxPrct(ITmsk,:),0,1),'o')
legend(["Evol","Manif"])
xticks(1:6);xticklabels(layernames);ylim(YLIM)
title(compose("IT exp(n=%d)",sum(ITmsk)))
suptitle(compose("%s Evol and Manif Correlated Voxels Across VGG Layers",Animal))
saveas(2,fullfile(hier_savedir,compose("%s_alllayer_gross_mean.jpg",Animal)))
savefig(2,fullfile(hier_savedir,compose("%s_alllayer_gross_mean.fig",Animal)))

%% Separate plot
YLIM = [0,0.5];
ExpType = "Evol";
CorrVoxPrct = arrayfun(@(H)H.(ExpType).corr_vox_prct(fi,:),ccHierStats,'UniformOutput',false);
CorrVoxPrct = cell2mat(CorrVoxPrct');
figure(11);clf
subtightplot(1,3,1);hold on 
bar(mean(CorrVoxPrct(V1msk,:),1))
errorbar(1:6,mean(CorrVoxPrct(V1msk,:),1),std(CorrVoxPrct(V1msk,:),0,1),'o')
xticks(1:6);xticklabels(layernames);ylim(YLIM)
title(compose("V1 exp(n=%d)",sum(V1msk)))
subtightplot(1,3,2);hold on 
bar(mean(CorrVoxPrct(V4msk,:),1))
errorbar(1:6,mean(CorrVoxPrct(V4msk,:),1),std(CorrVoxPrct(V4msk,:),0,1),'o')
xticks(1:6);xticklabels(layernames);ylim(YLIM)
title(compose("V4 exp(n=%d)",sum(V4msk)))
subtightplot(1,3,3);hold on 
bar(mean(CorrVoxPrct(ITmsk,:),1))
errorbar(1:6,mean(CorrVoxPrct(ITmsk,:),1),std(CorrVoxPrct(ITmsk,:),0,1),'o')
xticks(1:6);xticklabels(layernames);ylim(YLIM)
title(compose("IT exp(n=%d)",sum(ITmsk)))
suptitle(compose("%s %s Correlated Voxels Across VGG Layers",Animal,ExpType))
saveas(11,fullfile(hier_savedir,compose("%s_%s_gross_mean.jpg",Animal,ExpType)))
savefig(11,fullfile(hier_savedir,compose("%s_%s_gross_mean.fig",Animal,ExpType)))
%% Visualize Individual Experiment and see videos
Animal = "Beto";
for Expi = 14:45
prefchanlab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
ExpType = "Manif";
load(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t");
corr_vox_prct = corr_vox_num./totl_vox_num;
%
figure(4);clf;hold on 
nLayer = length(layernames);
subtightplot(1,3,1,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[corr_vox_prct(1:19,:);nan(1,nLayer);corr_vox_prct(20:23,:);nan(1,nLayer);corr_vox_prct(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best");xlim([.5,24.5])
title("Correlated Voxel Percents")
subtightplot(1,3,2,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[med_pos_cc(1:19,:);nan(1,nLayer);med_pos_cc(20:23,:);nan(1,nLayer);med_pos_cc(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best"); xlim([.5,24.5])
title("Median Positive Correlated Coefficient")
subtightplot(1,3,3,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[med_neg_cc(1:19,:);nan(1,nLayer);med_neg_cc(20:23,:);nan(1,nLayer);med_neg_cc(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best");xlim([.5,24.5])
title("Median Negative Correlated Coefficient")
suptitle(compose("%s %s Exp%d VGG16 Pref Chan %s\nCorrelation Percent and cc Distribution", Animal, ExpType, Expi, prefchanlab))
saveas(4, fullfile(hier_savedir,compose(Animal+compose("_%s_Exp%d_VGG16_cc.jpg", ExpType, Expi))))

ExpType = "Evol";
load(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t");
corr_vox_prct = corr_vox_num./totl_vox_num;

figure(5);clf;hold on 
nLayer = length(layernames);
subtightplot(1,3,1,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[corr_vox_prct(1:19,:);nan(1,nLayer);corr_vox_prct(20:23,:);nan(1,nLayer);corr_vox_prct(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best");xlim([.5,24.5])
title("Correlated Voxel Percents")
subtightplot(1,3,2,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[med_pos_cc(1:19,:);nan(1,nLayer);med_pos_cc(20:23,:);nan(1,nLayer);med_pos_cc(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best"); xlim([.5,24.5])
title("Median Positive Correlated Coefficient")
subtightplot(1,3,3,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[med_neg_cc(1:19,:);nan(1,nLayer);med_neg_cc(20:23,:);nan(1,nLayer);med_neg_cc(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best");xlim([.5,24.5])
title("Median Negative Correlated Coefficient")
suptitle(compose("%s %s Exp%d VGG16 Pref Chan %s\nCorrelation Percent and cc Distribution", Animal, ExpType, Expi, prefchanlab))
saveas(5, fullfile(hier_savedir,compose(Animal+compose("_%s_Exp%d_VGG16_cc.jpg", ExpType, Expi))))

winopen(fullfile(moviedir, compose("%s_Evol_Exp%02d_Best_PSTH.mov.avi",Animal,Expi)))
winopen(fullfile(moviedir, compose("%s_Manif_Exp%02d_Avg_PSTH.mov.avi",Animal,Expi)))
pause;
end

function [h,T]=ccDynam_plot(StatPrct,msks,labels,layernames,doerror)
if nargin < 5, doerror = false; end

nLayer = size(StatPrct,2);
% layernames = ccHierStats(1).Evol.layernames;
h=figure; T = tiledlayout(1,numel(msks),'Padd','none','TileSp','compact');
for i = 1:numel(msks)
nexttile(T,i);
msk = msks{i};
VarTrace_mean = mean(StatPrct(:,:,msk),3);
VarTrace_sem = sem(StatPrct(:,:,msk),3);
if ~ doerror
plot([1:19,NaN,20:23,NaN,24],[VarTrace_mean(1:19,:);nan(1,nLayer);VarTrace_mean(20:23,:);nan(1,nLayer);VarTrace_mean(24,:)]...
    , "-o","LineWidth",1.5)
else
hold on
meantrmat = [VarTrace_mean(1:19,:);nan(1,nLayer);VarTrace_mean(20:23,:);nan(1,nLayer);VarTrace_mean(24,:)];
semtrmat = [VarTrace_sem(1:19,:);nan(1,nLayer);VarTrace_sem(20:23,:);nan(1,nLayer);VarTrace_sem(24,:)];
% errorbar(repmat([1:19,NaN,20:23,NaN,24]',1,nLayer),meantrmat,semtrmat...
%     , "-o","LineWidth",1.5)
Cord = colororder();
for li = 1:nLayer
shadedErrorBar([1:19,NaN,20:23,NaN,24],meantrmat(:,li),semtrmat(:,li), ...
    'LineProp', {"-o","LineWidth",1.5,'Color',Cord(li,:)})
end
end
legend(layernames,'Location','Best')
title(labels(i))
end
end
