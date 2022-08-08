%% function to make the figure of Grant, Optimizer Comparison in vivo
% run the data processing code before running this
% Used in GECCO 2022 paper. 

meanscore_nat_share = nan(size(rasters, 1), length(block_list));
stdscore_nat_share = nan(size(rasters, 1), length(block_list));
for blocki = min(block_arr):max(block_arr)
    nat_msk = row_nat & block_arr == blocki;
    meanscore_nat_share(:, blocki) = mean(scores_tsr(:, nat_msk), 2);
    stdscore_nat_share(:, blocki)  = std(scores_tsr(:, nat_msk), 1, 2) / sqrt(sum(nat_msk));
end
%%
h4 = figure(11);clf;
hold on 
shadedErrorBar(block_list(1:end-1), meanscore_syn(channel_j, 1:end-1, 1), stdscore_syn(channel_j, 1:end-1, 1),...
    'lineprops',{'Color',[1.0,0.6,0.2,0.7]},'transparent',1,'patchSaturation',0.175)
shadedErrorBar(block_list(1:end-1), meanscore_syn(channel_j, 1:end-1, 2), stdscore_syn(channel_j, 1:end-1, 2),...
    'lineprops',{'Color',[0.2,0.6,1.0,0.7]},'transparent',1,'patchSaturation',0.175)

shadedErrorBar(block_list(1:end-1), meanscore_nat_share(channel_j, 1:end-1), stdscore_nat_share(channel_j, 1:end-1),...
    'lineprops',{'Color',[0,1,0,0.7]},'transparent',1,'patchSaturation',0.175)
axis tight
legend(["CMA-ES","GA-Classic","Natural img"],'Location',"Best")
xlabel("generations")
ylabel("Firing Rate Response (Hz)")
title([Exp_label_str, compose('Generation averaged score, channel %s', unit_name_arr{channel_j}), compose("Optimizers"))
%%
save_to_pdf(h4, compose("D:\\Optimizer_cmp_Exp%02d_2.pdf", Expi))

