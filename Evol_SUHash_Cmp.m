clearvars -except meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug ExpRecord
%%
Animal = "Beto";Set_Path;
%expftr = (contains(ExpRecord.expControlFN,"200319"));
expftr = contains(ExpRecord.expControlFN,"generate") & ...
    (contains(ExpRecord.Exp_collection,"SUHash"));%find(expftr)
Project_Manifold_Beto_loadRaw(find(expftr),Animal,false);
%% Analysis Code for Comparing Optimizers on a same unit. 
% much adapted from Evol_Traj_Cmp code, inspired Evol_Traj_analysis code.
% (it's kind of a multi-thread) version of Evol Traj analysis

% global block_arr gen_list color_seq row_gen row_nat
% global evol_stim_fr evol_stim_sem meanscore_syn stdscore_syn meanscore_nat stdscore_nat
Animal = "Beto"; Set_Path;
% result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Optimizer_Cmp";
% result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Optimizer_Tuning";

expftr = contains(ExpRecord.expControlFN,"generate") & ...
    (contains(ExpRecord.Exp_collection,"SUHash"));
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr),Animal); %find(expftr)
%% Prepare figure frames 
%result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Optimizer_Tuning";
% result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Optimizer_Cmp";
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Evol_SUHash";
for Triali = [1:4,6:length(meta_new)]
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);
fprintf("Processing  Exp %d:\n",Expi)
fprintf([ExpRecord.comments{exp_rowi},'\n'])
%% Sort channel id
% finding spike ID, note for multi-thread optimizer, we will have multiple
% pref_chan for different optimizers 
pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
if thread_num == 2, assert(pref_chan(1) == pref_chan(2)); end
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
if contains(meta.ephysFN, "Beto"), Animal = "Beto"; elseif contains(meta.ephysFN, "Alfa"), Animal = "Alfa"; end
savepath = fullfile(result_dir, compose("%s_Evol%02d_chan%02d", Animal, Expi, pref_chan(1)));
mkdir(savepath);
if thread_num == 2
    h = figure(1);set(h,'position',[1          41        2560         963],'Visible','on');
    axs{1} = subplot(1,2,1);axs{2} = subplot(1,2,2);
    h2 = figure(2);h2.Visible='on';clf; h2.Position = [  19         235        1779         743];
    axs2 = {}; axs2{1} = subplot(1,2,1); axs2{2} = subplot(1,2,2);
    h3 = figure(3);h3.Visible='on';h3.Position = [  782          43        1779         743];
    axs3 = {}; axs3{1} = subplot(1,2,1); axs3{2} = subplot(1,2,2);
elseif thread_num == 1
    h = figure(1);set(h,'position',[134    46   949   904],'Visible','on');
    axs{1} = subplot(1,1,1);
    h2 = figure(2);h2.Visible='on';clf; h2.Position = [  76   110   899   774];
    axs2 = {}; axs2{1} = subplot(1,1,1); 
    h3 = figure(3);h3.Visible='on';h3.Position = [  76   110   899   774];
    axs3 = {}; axs3{1} = subplot(1,1,1); 
end
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
for i = 1:thread_num
    pref_chan_id(i) = find(meta.spikeID==pref_chan(i) & ... % the id in the raster and lfps matrix 
                    unit_num_arr==unit_in_pref_chan(i)); % match for unit number
end
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
labstr = strcat(Exp_label_str,' (');
for i = 1:thread_num
    labstr = strcat(labstr, unit_name_arr{pref_chan_id(i)},', ');
end
Exp_label_str = strcat(labstr, ') ');
clear labstr
%% Optimizer Names 
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, strrep(string(Trials.TrialRecord.User.evoConfiguration{i,end}),"_"," ")];
end
%% Sort the images
imgnm = Trials.imageName;
% seperate the thread natural images and generated images 
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm) & ...% first 2 characters are not digits
        cellfun(@(c) ~contains(c(end-4:end), "_nat"), imgnm) ; % last 4 characters are not `_nat`
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
%% Get image name array
if ~ contains(meta.stimuli,'N:\') && contains(meta.stimuli(1:7),'Stimuli')
    meta.stimuli = fullfile('N:\', meta.stimuli);
end
imgColl = repmat("", length(block_list)-1, thread_num);
scoreColl = zeros(length(block_list)-1, thread_num);
for threadi = 1:thread_num
    channel_j = pref_chan_id(threadi);
    for blocki = min(block_arr):max(block_arr)-1 % exclude last gen if there is no max
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        [maxScore, maxIdx] = max(scores_tsr(channel_j, gen_msk));
        tmpimgs = imgnm(gen_msk);
        imgfullfn = ls(fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]));
        assert(~isempty(imgfullfn), "Image not found %s",fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]))
        imgColl(blocki, threadi) = fullfile(meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % solved @ Jan.29
        scoreColl(blocki, threadi) = maxScore;
    end
end
% Montage the images
set(0,'CurrentFigure',h); %clf; %
for threadi = 1:thread_num
set(gcf, "CurrentAxes", axs{threadi}); cla(axs{threadi},'reset'); 
montage(imgColl(:,threadi))
title([Exp_label_str, compose('Best Image per Generation'), compose("Optimizer %s", Optim_names(threadi))])
end
% set(gcf, "CurrentAxes", axs{2}); cla(axs{2},'reset'); 
% montage(imgColl(:,2))
% title([Exp_label_str, compose('Best Image per Generation'), compose("Optimizer %s", Optim_names(2))])
% saveas(h, fullfile(savepath, "EvolImageSeq_cmp.png"))
saveas(h, fullfile(result_dir, compose("%s_Evol%02d_EvolImageSeq.png", Animal, Expi)))
%% Prepare color sequence 
MAX_BLOCK_NUM = length(block_list); 
color_seq = brewermap(MAX_BLOCK_NUM, 'spectral');
for channel_j = pref_chan_id%1:size(rasters, 1)%pref_chan_id%
%% Plot Mean response compare figure
%channel_j = pref_chan_id;
% h1 = figure(1);clf
% ax1{1} = subplot(1,2,1);hold on
% plot(block_list, meanscore_syn(channel_j, :, 1), 'LineWidth',2,'Color','k')
% plot(block_list, meanscore_nat(channel_j, :, 1),'LineWidth',2,'Color','g')
% ax1{2} = subplot(1,2,2);hold on
% plot(block_list, meanscore_syn(channel_j, :, 2), 'LineWidth',2,'Color','k')
% plot(block_list, meanscore_nat(channel_j, :, 2),'LineWidth',2,'Color','g')
% ax1 = AlignAxisLimits(ax1);
%% Plot Mean response with shaded error bar compare
% channel_j = pref_chan_id;
% h2 = figure('Visible','on');h2.Position = [  782          43        1779         743];
% axs = {}; axs{1} = subplot(1,2,1);axs{2} = subplot(1,2,2);
for threadi = 1:thread_num
set(0,'CurrentFigure',h2); %clf; %
set(gcf, "CurrentAxes", axs2{threadi}); cla(axs2{threadi},'reset'); % some plot from shadedErrorBar cannot be cleared by cla!
hold on 
shadedErrorBar(block_list(1:end-1), meanscore_syn(channel_j, 1:end-1, threadi), stdscore_syn(channel_j, 1:end-1, threadi),...
    'lineprops',{'Color',[0,0,0,0.7]},'transparent',1,'patchSaturation',0.075)
shadedErrorBar(block_list(1:end-1), meanscore_nat(channel_j, 1:end-1, threadi), stdscore_nat(channel_j, 1:end-1, threadi),...
    'lineprops',{'Color',[0,1,0,0.7]},'transparent',1,'patchSaturation',0.075)
axis tight
legend(["Generated img","Natural img"],'Location',"Best")
xlabel("generations")
title([Exp_label_str, compose('Generation averaged score, channel %s', unit_name_arr{channel_j}), compose("Optimizer %s", Optim_names(threadi))])
end
axs2 = AlignAxisLimits(axs2);
%% Plot PSTH Evolution Plot 
% h3 = figure('Visible','on');h3.Position = [  782          43        1779         743];
% axs3 = {}; axs3{1} = subplot(1,2,1);axs3{2} = subplot(1,2,2);
for threadi = 1:thread_num
set(0,'CurrentFigure',h3); %clf; %
set(gcf, "CurrentAxes", axs3{threadi}); cla(axs3{threadi},'reset'); % some plot from shadedErrorBar cannot be cleared by cla!
hold on 
for blocki = 1:length(block_list) - 1
    shadedErrorBar([],evol_stim_fr(channel_j, :, blocki, threadi),evol_stim_sem(channel_j, :, blocki, threadi),...
    'lineprops',{'Color',[color_seq(blocki, :),0.85]},'transparent',1,'patchSaturation',0.075)
end
axis tight
% legend(["Generated img","Natural img"])
xlabel("Post Stimulus Time (ms)")
title([Exp_label_str, compose('Generation averaged PSTH , channel %s', unit_name_arr{channel_j}), compose("Optimizer %s", Optim_names(threadi))])
end
axs3 = AlignAxisLimits(axs3);
%%
% saveas(h2, fullfile(savepath, compose("score_traj_cmp_chan%s.png", unit_name_arr{channel_j})))
% saveas(h3, fullfile(savepath, compose("Evolv_psth_cmp_chan%s.png", unit_name_arr{channel_j})))
saveas(h2, fullfile(result_dir, compose("%s_Evol%02d_score_traj_chan%s.png", Animal, Expi, unit_name_arr{channel_j})))
saveas(h3, fullfile(result_dir, compose("%s_Evol%02d_Evolv_psth_chan%s.png", Animal, Expi, unit_name_arr{channel_j})))
end
end
%% Plot the Image Evolution Trajectory 
%% t test the last 5 generations 
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
scores_thread1 = scores_tsr(pref_chan_id, row_gen & (block_arr > max(block_arr)-6) & thread_msks{1});
scores_thread2 = scores_tsr(pref_chan_id, row_gen & (block_arr > max(block_arr)-6) & thread_msks{2});
[~,Pval,CI] = ttest2(scores_thread1, scores_thread2);


