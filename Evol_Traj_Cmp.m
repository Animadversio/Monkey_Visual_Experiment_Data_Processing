%% Evolution Compare 
% Comparing Evolution Experiments of different size 
% Final Score and Evolution time course. 
% Get all the pairs! 
expid = find(ExpSpecTable_Aug.Expi >= 28 & ...
    contains(ExpSpecTable_Aug.expControlFN,'generate') & ...
    contains(ExpSpecTable_Aug.Exp_collection,'Manifold'));
expid(1:2) = expid(2:-1:1);
expid(3:4) = expid(4:-1:3);% Switch the first 2 entries so that 1deg experiments go first, 3 deg experiments go second. 
%  Note exp 30 it's 3 deg. 31 it's 1 deg, also need to be switched 
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(expid);
%% Check That the experiment follows the pairs
for Triali = 1:length(meta_new)
    fprintf("Position: %d %d, %d Deg\n",Trials_new{Triali}.XY{1}, Trials_new{Triali}.TrialRecord.User.width_perChan)
end
%% 
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Evolution_Exp";
savepath = fullfile(result_dir, compose("Manifold_Evol_Resize_Cmp"));
mkdir(savepath);
global  Trials rasters 
%%
%h1 = figure('Visible','off');
h1 = figure('Visible','on');clf; h1.Position = [1          41        2560         963];
axs1 = {}; axs1{1} = subplot(1,2,1);axs1{2} = subplot(1,2,2);
h12 = figure('Visible','on');clf; h12.Position = [1          41        2560         963];
axs12 = {}; axs12{1} = subplot(1,2,1);axs12{2} = subplot(1,2,2);
h2 = figure('Visible','on');clf; h2.Position = [  19         235        1779         743];
axs = {}; axs{1} = subplot(1,2,1);axs{2} = subplot(1,2,2);
h3 = figure('Visible','on');h3.Position = [  782          43        1779         743];
axs3 = {}; axs3{1} = subplot(1,2,1);axs3{2} = subplot(1,2,2);
%color_seq = brewermap(length(gen_list), 'spectral');
%%
scores_col = {};
for groupi = 1:9
    rel_i = (groupi - 1) * 2;
    %MAX_BLOCK_NUM = max(max(cell2mat(Trials_new{rel_i+1}.block)), max(cell2mat(Trials_new{rel_i+2}.block)));
    MAX_BLOCK_NUM = max(max(cell2mat(Trials_new{rel_i+1}.block)), max(cell2mat(Trials_new{rel_i+2}.block)));
    color_seq = brewermap(MAX_BLOCK_NUM, 'spectral');
    % MAX_BLOCK_NUM = 50
for Triali = 1:2
meta = meta_new{Triali + rel_i};
rasters = rasters_new{Triali + rel_i};
Trials = Trials_new{Triali + rel_i};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi_tab = ExpSpecTable_Aug.Expi(exp_rowi);
%assert(Expi_tab == Expi, "Data Expi doesn't match that in exp record! Check your indexing or record.")
Expi = Expi_tab; 
fprintf("Processing Exp %d, %s\n", Expi, meta.comments)

%% Sort channel id
pref_chan = Trials.TrialRecord.User.prefChan;
pref_chan_id = find(meta.spikeID==pref_chan); % the id in the raster and lfps matrix 
Exp_label_str = sprintf("Exp%d pref chan %d (center [%.1f %.1f] width %.1f)", Expi, pref_chan, Trials.XY{1}, Trials.TrialRecord.User.width_perChan);
unit_name_arr = generate_unit_labels(meta.spikeID);
%% Sort the block id and the identity of image. 
block_arr = cell2mat(Trials.block);
imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
gen_list = min(block_arr):max(block_arr);
%%
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % debug @ Jan.31
meanscore_syn = zeros(size(rasters, 1), length(gen_list)); % []
stdscore_syn = zeros(size(rasters, 1), length(gen_list));
meanscore_nat = zeros(size(rasters, 1), length(gen_list));
stdscore_nat = zeros(size(rasters, 1), length(gen_list));
for blocki = gen_list
    meanscore_syn(:,blocki) = mean(scores_tsr(:, row_gen & block_arr == blocki), 2);%cat(2, meanscore_syn, tmpscore_syn);
    meanscore_nat(:,blocki) = mean(scores_tsr(:, row_nat & block_arr == blocki), 2);%cat(2, meanscore_nat, tmpscore_nat);
    stdscore_syn(:,blocki)  = std(scores_tsr(:, row_gen & block_arr == blocki), 1, 2) / sqrt(sum(row_gen & block_arr == blocki));%cat(2, stdscore_syn, tmpstdscore_syn);
    stdscore_nat(:,blocki)  = std(scores_tsr(:, row_nat & block_arr == blocki), 1, 2) / sqrt(sum(row_nat & block_arr == blocki));%cat(2, stdscore_nat, tmpstdscore_nat);
    % error fixed on Jan.26 2020, the row_nat is mistaken as row_gen
end 
% Average PSTH and sem for evolved image
evol_stim_fr = zeros(size(rasters, 1), size(rasters, 2), length(gen_list));
evol_stim_sem = zeros(size(rasters, 1), size(rasters, 2), length(gen_list));
for blocki = 1:length(gen_list)
    evol_stim_fr(:, :, blocki) = mean(rasters(:,:, row_gen & block_arr == blocki),3);
    evol_stim_sem(:, :, blocki) = std(rasters(:,:, row_gen & block_arr == blocki),1,3)/sqrt(sum(row_gen & block_arr == blocki));
end
%% Collect Final Generation Data
channel_j = pref_chan_id;
msk = blocki & block_arr <= gen_list(end-1) & block_arr >= gen_list(end-5);
scores_col{groupi,Triali} = scores_tsr(channel_j,msk);
%% Plot the Evolution Image Sequence. 
channel_j = pref_chan_id(1);
imgColl = repmat("", length(gen_list), 1);
scoreColl = zeros(length(gen_list), 1);
for blocki = min(block_arr):max(block_arr)
    gen_msk = row_gen & block_arr == blocki; 
    [maxScore, maxIdx] = max(scores_tsr(channel_j, gen_msk));
    tmpimgs = imgnm(gen_msk);
    imgfullfn = ls(fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]))
    imgColl(blocki) = fullfile(meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % solved @ Jan.29
    scoreColl(blocki) = maxScore;
end
% Montage the images
set(0,'CurrentFigure',h1); %clf; %
set(gcf, "CurrentAxes", axs1{Triali}); cla(axs1{Triali},'reset'); 
montage(imgColl(:,1))
title([Exp_label_str, compose('Best Image per Generation')])
% Montage the images resized to the size of image.
set(0,'CurrentFigure',h12); %clf; %
set(gcf, "CurrentAxes", axs12{Triali}); cla(axs12{Triali},'reset'); 
if Trials.TrialRecord.User.width_perChan == 1
    montage(imgColl(:,1), 'ThumbnailSize', [180,180], 'BorderSize', 180, 'BackgroundColor', [0.5, 0.5, 0.5]) % 180 is manually chosen to get the image size. 
elseif Trials.TrialRecord.User.width_perChan == 3
    montage(imgColl(:,1), 'ThumbnailSize', [180,180], 'BackgroundColor', [0.5, 0.5, 0.5])
end
title([Exp_label_str, compose('Best Image per Generation')])
%% Plot the PSTH evolution
channel_j = pref_chan_id;
set(0,'CurrentFigure',h3); %clf; %
set(gcf, "CurrentAxes", axs3{Triali}); cla(axs3{Triali},'reset'); % some plot from shadedErrorBar cannot be cleared by cla!
hold on 
for blocki = 1:length(gen_list) - 1
    shadedErrorBar([],evol_stim_fr(channel_j, :, blocki),evol_stim_sem(channel_j, :, blocki),...
    'lineprops',{'Color',[color_seq(blocki, :),0.85]},'transparent',1,'patchSaturation',0.075)
end
axis tight
% legend(["Generated img","Natural img"])
xlabel("generations")
title([Exp_label_str, compose('PSTH averaged scores, channel %s', unit_name_arr{channel_j})])
%% Plot the score changing with Standard Deviation. 
channel_j = pref_chan_id;
set(0,'CurrentFigure',h2); %clf; %
set(gcf, "CurrentAxes", axs{Triali}); cla(axs{Triali},'reset'); % some plot from shadedErrorBar cannot be cleared by cla!
hold on 
shadedErrorBar(gen_list(1:end-1), meanscore_syn(channel_j, 1:end-1), stdscore_syn(channel_j, 1:end-1),...
    'lineprops',{'Color',[0,0,0,0.7]},'transparent',1,'patchSaturation',0.075)
shadedErrorBar(gen_list(1:end-1), meanscore_nat(channel_j, 1:end-1), stdscore_nat(channel_j, 1:end-1),...
    'lineprops',{'Color',[0,1,0,0.7]},'transparent',1,'patchSaturation',0.075)
axis tight
legend(["Generated img","Natural img"])
xlabel("generations")
title([Exp_label_str, compose('PSTH averaged scores, channel %s', unit_name_arr{channel_j})])
%%
end
axs = AlignAxisLimits(axs);
axs3 = AlignAxisLimits(axs3);
saveas(h1, fullfile(savepath, compose("EvolImageSeq_cmp_group%02d_Exp%02d.png", groupi, Expi)))
saveas(h12, fullfile(savepath, compose("EvolImageSeq_Rszcmp_group%02d_Exp%02d.png", groupi, Expi)))
%saveas(h2, fullfile(savepath, compose("score_traj_cmp_group%02d_Exp%02d.png", groupi, Expi)))
%saveas(h3, fullfile(savepath, compose("Evolv_psth_cmp_group%02d_Exp%02d.png", groupi, Expi)))
end
%% T test of score
for groupi = 1:9
[H,P,CI] = ttest2(scores_col{groupi,1}, scores_col{groupi,2}); % 1 deg (col 1) vs 3 deg (col 2) 
fprintf("H=%d P=%e CI = [%.2f %.2f]\n",H,P,CI)
end

%% Plot the ttest comparison from all 9 pair of experiments. 
h = figure(1);h.Position = [331         185        1658         793];
clf;hold on 
expname_col = [];
for groupi = 1:9
[H,P,CI] = ttest2(scores_col{groupi,1}, scores_col{groupi,2}); % 1 deg (col 1) vs 3 deg (col 2) 
ttest_str = sprintf("H=%d \nP=%.2e \nCI =[%.1f %.1f]",H,P,CI);
errorbar(3 * groupi + [0, 1],  cellfun(@mean, scores_col(groupi,:)), cellfun(@(c) std(c)/sqrt(length(c)), scores_col(groupi,:)),'LineWidth',2)
text(3*groupi, 200, ttest_str,'FontSize',12)
Expi1 = ExpSpecTable_Aug.Expi(contains(ExpSpecTable_Aug.ephysFN, meta_new{2*groupi-1}.ephysFN));
Expi2 = ExpSpecTable_Aug.Expi(contains(ExpSpecTable_Aug.ephysFN, meta_new{2*groupi}.ephysFN));
W1 = Trials_new{2*groupi-1}.TrialRecord.User.width_perChan;
W2 = Trials_new{2*groupi  }.TrialRecord.User.width_perChan;
pref_chan = Trials_new{2*groupi-1}.TrialRecord.User.prefChan;
exp_str = sprintf("Exp%02d (%d deg)\nVS Exp%02d (%d deg)\n Pref chan %02d",Expi1, W1, Expi2, W2, pref_chan);
text(3*groupi, 280, exp_str,'FontSize',12)
expname_str = sprintf("Exp%02d VS Exp%02d",Expi1,Expi2);
expname_col = [expname_col,expname_str];
end
xticks(3 * [1:9] + 0.5)
xticklabels(expname_col)
xlim([2, 30])
suptitle(["Score Comparison of 1 deg evolution VS 3 deg Evolution","Score computed with [50:200] - [0:40]"])
ylabel("Score of Last 5 Generations Compared (st.e plotted )")
saveas(h,fullfile(savepath, "Final_Stage_Score_Comparison.png"))
% function axs = AlignAxisLimits(axs)
% % Given a group of axis, make their YLim and XLim the same as each other
% %     YLIM = [min([axs{1}.YLim, axs{2}.YLim]), max([axs{1}.YLim, axs{2}.YLim])];
% %     XLIM = [min([axs{1}.XLim, axs{2}.XLim]), max([axs{1}.XLim, axs{2}.XLim])];
% % Now a util function in utils subfolder
%     XLIM = axs{1}.XLim; YLIM = axs{1}.YLim;
%     for i = 2:numel(axs)
%         XLIM = [min(XLIM(1), axs{i}.XLim(1)), max(XLIM(2), axs{i}.XLim(2))];
%         YLIM = [min(YLIM(1), axs{i}.YLim(1)), max(YLIM(2), axs{i}.YLim(2))];
%     end
%     for i = 1:numel(axs)
%         axs{i}.XLim = XLIM; axs{i}.YLim = YLIM;
%     end
% end