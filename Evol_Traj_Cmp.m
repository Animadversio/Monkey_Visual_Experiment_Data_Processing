%% Evolution Compare 
% Comparing Evolution Experiments of different size 
% Final Score and Evolution time course. 
% Get all the pairs! 
% expid = find(ExpSpecTable_Aug.Expi > 27 & contains(ExpSpecTable_Aug.expControlFN,'generate'));
% expid(1) = 67; expid(2) = 64; % Switch the first 2 entries so that 1deg experiments go first, 3 deg experiments go second. 
% [meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(expid);
%% 
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Evolution_Exp";
savepath = fullfile(result_dir, compose("Manifold_Evol_Resize_Cmp"));
mkdir(savepath);
global  Trials rasters 
%h1 = figure('Visible','off');
h2 = figure('Visible','on');clf; h2.Position = [  19         235        1779         743];
axs = {}; axs{1} = subplot(1,2,1);axs{2} = subplot(1,2,2);
h3 = figure('Visible','on');h3.Position = [  782          43        1779         743];
axs3 = {}; axs3{1} = subplot(1,2,1);axs3{2} = subplot(1,2,2);
%color_seq = brewermap(length(gen_list), 'spectral');
for groupi = 9%1:8
    rel_i = -2;%(groupi - 1) * 2;
    %MAX_BLOCK_NUM = max(max(cell2mat(Trials_new{rel_i+1}.block)), max(cell2mat(Trials_new{rel_i+2}.block)));
    MAX_BLOCK_NUM = max(max(cell2mat(Trials_new{2}.block)), max(cell2mat(Trials_new{6}.block)));
    color_seq = brewermap(MAX_BLOCK_NUM, 'spectral');
    % MAX_BLOCK_NUM = 50
for Triali = 1:2
meta = meta_new{Triali * 4 + rel_i};
rasters = rasters_new{Triali * 4 + rel_i};
Trials = Trials_new{Triali * 4 + rel_i};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi_tab = ExpSpecTable_Aug.Expi(exp_rowi);
%assert(Expi_tab == Expi, "Data Expi doesn't match that in exp record! Check your indexing or record.")
Expi = Expi_tab; 
fprintf("Processing Exp %d, %s\n", Expi, meta.comments)

%% Sort channel id
pref_chan = Trials.TrialRecord.User.prefChan;
pref_chan_id = find(meta.spikeID==pref_chan); % the id in the raster and lfps matrix 
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan);
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
scores_tsr = squeeze(mean(rasters(:, 151:200, :), 2) - mean(rasters(:, 1:40, :), 2));
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
%%
channel_j = pref_chan_id;
%%
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
%%
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
saveas(h2, fullfile(savepath, compose("score_traj_cmp_group%02d_Exp%02d.png", groupi, Expi)))
saveas(h3, fullfile(savepath, compose("Evolv_psth_cmp_group%02d_Exp%02d.png", groupi, Expi)))
end
function axs = AlignAxisLimits(axs)
% Given a group of axis, make their YLim and XLim the same as each other
%     YLIM = [min([axs{1}.YLim, axs{2}.YLim]), max([axs{1}.YLim, axs{2}.YLim])];
%     XLIM = [min([axs{1}.XLim, axs{2}.XLim]), max([axs{1}.XLim, axs{2}.XLim])];
% Now a util function in utils subfolder
    XLIM = axs{1}.XLim; YLIM = axs{1}.YLim;
    for i = 2:numel(axs)
        XLIM = [min(XLIM(1), axs{i}.XLim(1)), max(XLIM(2), axs{i}.XLim(2))];
        YLIM = [min(YLIM(1), axs{i}.YLim(1)), max(YLIM(2), axs{i}.YLim(2))];
    end
    for i = 1:numel(axs)
        axs{i}.XLim = XLIM; axs{i}.YLim = YLIM;
    end
end