%% ccFT_hierarchy_plot
ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";

Animal = "Alfa";
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
ExpType = "Evol";
layernames = ["fc8", "fc7", "fc6", "conv5_3", "conv4_3", "conv3_3"];
wdwn = 24; layerN = numel(layernames);
tic;
for Expi = 1:numel(EStats)
corr_vox_num = nan(wdwn, layerN); 
corr_vox_prct = nan(wdwn, layerN); 
poscorr_vox_num = nan(wdwn, layerN); 
poscorr_vox_prct = nan(wdwn, layerN); 
negcorr_vox_num = nan(wdwn, layerN); 
negcorr_vox_prct = nan(wdwn, layerN); 
med_pos_cc = nan(wdwn, layerN); 
med_neg_cc = nan(wdwn, layerN); 
mean_pos_t = nan(wdwn, layerN); 
mean_neg_t = nan(wdwn, layerN); 
totl_vox_num = nan(1,layerN);
for li = 1:layerN
layername = layernames(li); %["conv5_3"]%"conv5_3";
if layername=="conv3_3" && Expi ==40,continue; end
outfn = fullfile(ccmat_dir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername));
load(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
posmsk = t_signif_tsr > 5;
negmsk = t_signif_tsr <-5;
signif_n = sum(posmsk | negmsk,[1,2,3]);
possig_n = sum(posmsk,[1,2,3]);
negsig_n = sum(negmsk,[1,2,3]);
vox_num = prod(size(cc_tsr,[1,2,3]));%numel(cc_tsr);
totl_vox_num(li) = vox_num; 
corr_vox_num(:,li) = signif_n;
corr_vox_prct(:,li) = signif_n / vox_num;
poscorr_vox_num(:,li) = possig_n;
poscorr_vox_prct(:,li) = possig_n / vox_num;
negcorr_vox_num(:,li) = negsig_n;
negcorr_vox_prct(:,li) = negsig_n / vox_num;
for fi = 1:wdwn
ctmp = cc_tsr(:,:,:,fi);
ttmp = t_signif_tsr(:,:,:,fi);
med_pos_cc(fi,li) = median(ctmp(posmsk(:,:,:,fi)),'all');
med_neg_cc(fi,li) = median(ctmp(negmsk(:,:,:,fi)),'all');
mean_pos_t(fi,li) = mean(ttmp(posmsk(:,:,:,fi)),'all');
mean_neg_t(fi,li) = mean(ttmp(negmsk(:,:,:,fi)),'all');
end
end
save(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t",...
    "poscorr_vox_num", "poscorr_vox_prct", "negcorr_vox_num", "negcorr_vox_prct");


prefchanlab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
load(fullfile(hier_savedir, compose("%s_%s_Exp%d_VGG16.mat",Animal,ExpType,Expi)),...
    "layernames","wdw_vect","totl_vox_num","corr_vox_num","corr_vox_prct","med_pos_cc","med_neg_cc","mean_pos_t","mean_neg_t",...
    "poscorr_vox_num", "poscorr_vox_prct", "negcorr_vox_num", "negcorr_vox_prct");
% corr_vox_prct = corr_vox_num ./ totl_vox_num;

figure(4);clf;hold on 
layerN = length(layernames);
subtightplot(1,5,1,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[corr_vox_prct(1:19,:);nan(1,layerN);corr_vox_prct(20:23,:);nan(1,layerN);corr_vox_prct(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best");xlim([.5,24.5])
title("Correlated Voxel Percents")
subtightplot(1,5,2,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[poscorr_vox_prct(1:19,:);nan(1,layerN);poscorr_vox_prct(20:23,:);nan(1,layerN);poscorr_vox_prct(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best");xlim([.5,24.5])
title("Positiv Correlated Voxel Percents")
subtightplot(1,5,3,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[negcorr_vox_prct(1:19,:);nan(1,layerN);negcorr_vox_prct(20:23,:);nan(1,layerN);negcorr_vox_prct(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best");xlim([.5,24.5])
title("Negativ Correlated Voxel Percents")
subtightplot(1,5,4,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[med_pos_cc(1:19,:);nan(1,layerN);med_pos_cc(20:23,:);nan(1,layerN);med_pos_cc(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best"); xlim([.5,24.5])
title("Median Positive Correlated Coefficient")
subtightplot(1,5,5,0.04,0.07,0.03)
plot([1:19,NaN,20:23,NaN,24],[med_neg_cc(1:19,:);nan(1,layerN);med_neg_cc(20:23,:);nan(1,layerN);med_neg_cc(24,:)], "-o"...
    ,"LineWidth",1.5)
legend(layernames, 'Location',"Best");xlim([.5,24.5])
title("Median Negative Correlated Coefficient")
suptitle(compose("%s %s Exp%d VGG16 Pref Chan %s\nCorrelation Percent and cc Distribution", Animal, ExpType, Expi, prefchanlab))
    
savenm = compose(Animal+compose("_%s_Exp%d_VGG16_cc", ExpType, Expi));
saveas(4, fullfile(hier_savedir,savenm + ".jpg"))
saveas(4, fullfile(hier_savedir,savenm + ".pdf"))
savefig(4, fullfile(hier_savedir,savenm + ".fig"))
toc;
end