%% Evol_Optimizer_Summary plot 
%% Analysis Code for Comparing Optimizers on a same unit. 
% much adapted from Evol_Traj_Cmp code, inspired Evol_Traj_analysis code.
% (it's kind of a multi-thread) version of Evol Traj analysis

Set_Path;
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Optimizer_Cmp";
ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xls");
expftr = contains(ExpSpecTable_Aug.expControlFN,"generate") & ...
     contains(ExpSpecTable_Aug.Exp_collection, "Optimizer_cmp");
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr)); 
savepath = fullfile(result_dir, "summary");
mkdir(savepath);
%% Prepare figure frames 
h = figure('Visible','on');set(h,'position',[1          41        2560         963]);
axs{1} = subplot(1,2,1);axs{2} = subplot(1,2,2);
h2 = figure('Visible','on');clf; h2.Position = [  19         235        1779         743];
axs2 = {}; axs2{1} = subplot(1,2,1); axs2{2} = subplot(1,2,2);
h3 = figure('Visible','on');h3.Position = [  782          43        1779         743];
axs3 = {}; axs3{1} = subplot(1,2,1); axs3{2} = subplot(1,2,2);
%%
score_traces = cell(length(meta_new), 2);
block_traces = cell(length(meta_new), 2);
ref_traces = cell(length(meta_new), 2);
ref_block_traces = cell(length(meta_new), 2);
for Triali = 1:length(meta_new)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpSpecTable_Aug.Expi(exp_rowi);
fprintf("Exp %d:\n",Expi)
disp(ExpSpecTable_Aug.comments(exp_rowi))
% assert(Expi_tab == Expi, "Data Expi doesn't match that in exp record! Check your indexing or record.")
%% Sort channel id
% finding spike ID, note for multi-thread optimizer, we will have multiple
% pref_chan for different optimizers 
pref_chan = Trials.TrialRecord.User.prefChan;
assert(pref_chan(1) == pref_chan(2))
pref_chan_id = find(meta.spikeID==pref_chan(1)); % the id in the raster and lfps matrix 
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
unit_name_arr = generate_unit_labels(meta.spikeID);
%% Optimizer Names 
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, string(Trials.TrialRecord.User.evoConfiguration{i,end})];
end
%% Sort the images
imgnm = Trials.imageName;
% seperate the thread natural images and generated images 
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the thread 1 and 2 
row_thread0 = contains(imgnm, compose("thread%03d", 0));
row_thread1 = contains(imgnm, compose("thread%03d", 1));
assert(sum(row_thread0)+sum(row_thread1) == length(imgnm))
thread_msks = {row_thread0, row_thread1}; % store masks in a structure for the ease to iterate
% get the generation number 
block_arr = cell2mat(Trials.block);
% if needed, analyze the image names to see the block and thread
% information. see Evol_Traj_Cmp
%% Compute score for evolution 
block_list = min(block_arr):max(block_arr);
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]

meanscore_syn = nan(size(rasters, 1), length(block_list), 2); % [unit_num, gen_nums, threads]
stdscore_syn = nan(size(rasters, 1), length(block_list), 2); 
meanscore_nat = nan(size(rasters, 1), length(block_list), 2);
stdscore_nat = nan(size(rasters, 1), length(block_list), 2);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        meanscore_syn(:, blocki, threadi) = mean(scores_tsr(:, gen_msk), 2);
        meanscore_nat(:, blocki, threadi) = mean(scores_tsr(:, nat_msk), 2);
        stdscore_syn(:, blocki, threadi)  = std(scores_tsr(:, gen_msk), 1, 2) / sqrt(sum(gen_msk));
        stdscore_nat(:, blocki, threadi)  = std(scores_tsr(:, nat_msk), 1, 2) / sqrt(sum(nat_msk));
    end
    score_traces{Triali,threadi} = scores_tsr(pref_chan_id, row_gen & thread_msks{threadi});
    block_traces{Triali,threadi} = block_arr(row_gen & thread_msks{threadi});
    ref_traces{Triali,threadi} = scores_tsr(pref_chan_id, row_nat & thread_msks{threadi});
    ref_block_traces{Triali,threadi} = block_arr(row_nat & thread_msks{threadi});
end
%% Compute Average PSTH and sem for evolved image
evol_stim_fr = nan(size(rasters, 1), size(rasters, 2), length(block_list), thread_num);
evol_stim_sem = nan(size(rasters, 1), size(rasters, 2), length(block_list), thread_num);
for threadi = 1:thread_num
for blocki = 1:length(block_list)
    gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
    evol_stim_fr(:, :, blocki, threadi) = mean(rasters(:,:, gen_msk),3);
    evol_stim_sem(:, :, blocki, threadi) = std(rasters(:,:, gen_msk),1,3) / sqrt(sum(gen_msk));
end
end
% %% Prepare color sequence 
% MAX_BLOCK_NUM = length(block_list); 
% color_seq = brewermap(MAX_BLOCK_NUM, 'spectral');
% for channel_j = 1:size(rasters, 1)%pref_chan_id%
% %% Plot Mean response with shaded error bar compare
% % channel_j = pref_chan_id;
% % h2 = figure('Visible','on');h2.Position = [  782          43        1779         743];
% % axs = {}; axs{1} = subplot(1,2,1);axs{2} = subplot(1,2,2);
% for threadi = 1:thread_num
% set(0,'CurrentFigure',h2); %clf; %
% set(gcf, "CurrentAxes", axs2{threadi}); cla(axs2{threadi},'reset'); % some plot from shadedErrorBar cannot be cleared by cla!
% hold on 
% shadedErrorBar(block_list(1:end-1), meanscore_syn(channel_j, 1:end-1, threadi), stdscore_syn(channel_j, 1:end-1, threadi),...
%     'lineprops',{'Color',[0,0,0,0.7]},'transparent',1,'patchSaturation',0.075)
% shadedErrorBar(block_list(1:end-1), meanscore_nat(channel_j, 1:end-1, threadi), stdscore_nat(channel_j, 1:end-1, threadi),...
%     'lineprops',{'Color',[0,1,0,0.7]},'transparent',1,'patchSaturation',0.075)
% axis tight
% legend(["Generated img","Natural img"],'Location',"Best")
% xlabel("generations")
% title([Exp_label_str, compose('Generation averaged score, channel %s', unit_name_arr{channel_j}), compose("Optimizer %s", Optim_names(threadi))])
% end
% axs2 = AlignAxisLimits(axs2);
% %% Plot PSTH Evolution Plot 
% % h3 = figure('Visible','on');h3.Position = [  782          43        1779         743];
% % axs3 = {}; axs3{1} = subplot(1,2,1);axs3{2} = subplot(1,2,2);
% for threadi = 1:thread_num
% set(0,'CurrentFigure',h3); %clf; %
% set(gcf, "CurrentAxes", axs3{threadi}); cla(axs3{threadi},'reset'); % some plot from shadedErrorBar cannot be cleared by cla!
% hold on 
% for blocki = 1:length(block_list) - 1
%     shadedErrorBar([],evol_stim_fr(channel_j, :, blocki, threadi),evol_stim_sem(channel_j, :, blocki, threadi),...
%     'lineprops',{'Color',[color_seq(blocki, :),0.85]},'transparent',1,'patchSaturation',0.075)
% end
% axis tight
% % legend(["Generated img","Natural img"])
% xlabel("Post Stimulus Time (ms)")
% title([Exp_label_str, compose('Generation averaged PSTH , channel %s', unit_name_arr{channel_j}), compose("Optimizer %s", Optim_names(threadi))])
% end
% axs3 = AlignAxisLimits(axs3);
%%
% end
end
%% Plot the Image Evolution Trajectory 
%% t test the last 5 generations 
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
scores_thread1 = scores_tsr(pref_chan_id, row_gen & (block_arr > max(block_arr)-6) & thread_msks{1});
scores_thread2 = scores_tsr(pref_chan_id, row_gen & (block_arr > max(block_arr)-6) & thread_msks{2});
[~,Pval,CI] = ttest2(scores_thread1, scores_thread2);
%%
CMA_score = zeros(length(meta_new),1);
CMA_err = zeros(length(meta_new),1);
GA_score = zeros(length(meta_new),1);
GA_err = zeros(length(meta_new),1);
CMA_ref_score = zeros(length(meta_new),1);
CMA_ref_err = zeros(length(meta_new),1);
GA_ref_score = zeros(length(meta_new),1);
GA_ref_err = zeros(length(meta_new),1);
for Triali = 1:length(meta_new)
    scores = score_traces{Triali,1}(block_traces{Triali,1} > max(block_traces{Triali,1})-5);
    CMA_score(Triali) = mean(scores);
    CMA_err(Triali) = std(scores)/sqrt(length(scores));
    scores = score_traces{Triali,2}(block_traces{Triali,2} > max(block_traces{Triali,2})-5);
    GA_score(Triali) = mean(scores);
    GA_err(Triali) = std(scores)/sqrt(length(scores));
    
    scores = ref_traces{Triali,1}(ref_block_traces{Triali,1} > max(ref_block_traces{Triali,1})-5);
    CMA_ref_score(Triali) = mean(scores);
    CMA_ref_err(Triali) = std(scores)/sqrt(length(scores));
    scores = ref_traces{Triali,2}(ref_block_traces{Triali,2} > max(ref_block_traces{Triali,2})-5);
    GA_ref_score(Triali) = mean(scores);
    GA_ref_err(Triali) = std(scores)/sqrt(length(scores));
end
%% Error bar comparison plot
% https://www.mathworks.com/matlabcentral/answers/455796-errorbars-on-scatter-plot
figure(1);clf
scatter(GA_ref_score, CMA_ref_score, 50,'g','filled');
hold on
eb2(1) = errorbar(GA_ref_score,CMA_ref_score,GA_ref_err, 'horizontal', 'LineStyle', 'none');
eb2(2) = errorbar(GA_ref_score,CMA_ref_score,CMA_ref_err, 'vertical', 'LineStyle', 'none');
set(eb2, 'color', 'g', 'LineWidth', 2)
scatter(GA_score, CMA_score, 50,'k','filled');
eb(1) = errorbar(GA_score,CMA_score,GA_err, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(GA_score,CMA_score,CMA_err, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 2)
line([0,350],[0,350])
xlim([0,350]);ylim([0,350]);axis equal
% Add annotation! 
for Triali = 1:length(meta_new)
    pref_chan = Trials_new{Triali}.TrialRecord.User.prefChan;
    Expi = ExpSpecTable_Aug.Expi(contains(ExpSpecTable_Aug.ephysFN, meta_new{Triali}.ephysFN));
    text(GA_score(Triali)-25,CMA_score(Triali)+20, compose("Exp%d\npref chan %d", Expi, pref_chan(1)));
end
title("Score Comparison of last 5 Generations")
ylabel("CMA-ES")
xlabel("GA-classic")
legend(["Synthetic","","","Natural","","",""])
%%
h=figure(1);
save_to_pdf(h,fullfile(savepath,"ScoreCmp.pdf"))