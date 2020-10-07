%% Adapted from Hess_Tune_Analysis
Animal = "Alfa";Set_Path;
expftr = contains(ExpRecord.Exp_collection,"Sel_BigSet");%contains(ExpRecord.expControlFN,["200817"]);%& contains(ExpRecord.expControlFN, "selectivity");; %& ,"200812"contains(ExpRecord.Exp_collection,"BigGAN_Hessian");% & contains(ExpRecord.Exp_collection,"BigGAN");
% expftr = contains(ExpRecord.Exp_collection,"BigGAN_Hessian") & contains(ExpRecord.expControlFN, "selectivity");
fllist = find(expftr);no_return=false;
[meta_new,rasters_new,~,Trials_new] = loadExperiments(fllist(1:end),Animal,no_return);
%%
for Triali = 2%2:7

meta = meta_new{Triali};
rasters = rasters_new{Triali};
% lfps = lfps_new{Triali};
Trials = Trials_new{Triali};

unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
imgname_uniq = unique(Trials.imageName); 
idx_col = cellfun(@(imgnm)find(contains(Trials.imageName, imgnm)), imgname_uniq, 'Uni', false);
cnt_arr = cellfun(@length, idx_col);

bsl_mat = squeeze(mean(rasters(:,1:50,:),2));
rsp_mat = squeeze(mean(rasters(:,51:200,:),2));
bsl_img_arr = cell2mat(cellfun(@(idx) mean(bsl_mat(:,idx),2), idx_col', 'Un', false));
act_img_arr = cell2mat(cellfun(@(idx) mean(rsp_mat(:,idx),2), idx_col', 'Un', false));
rsp_img_arr = cell2mat(cellfun(@(idx) mean(rsp_mat(:,idx) - bsl_mat(:,idx),2), idx_col', 'Un', false));
%%
figure(6);clf;
TL = tiledlayout(9, 9);
for iCh = 1:size(rasters,1)
ax = nexttile;hold on 
histogram(act_img_arr(iCh,:));
histogram(bsl_img_arr(iCh,:));
histogram(rsp_img_arr(iCh,:));
% legend(["Evoked Rate","Baseline","Activation"])
title(ax, unit_name_arr(iCh))
end
end
%%
figure;hold on
histogram(act_img_arr(iCh,:));
histogram(bsl_img_arr(iCh,:));
histogram(rsp_img_arr(iCh,:));
legend(["Evoked Rate","Baseline","Activation"])