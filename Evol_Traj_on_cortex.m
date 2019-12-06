%% 
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Evolution_Exp";
global block_arr gen_list color_seq row_gen row_nat
global evol_stim_fr evol_stim_sem meanscore_syn stdscore_syn meanscore_nat stdscore_nat
for Expi = 1:37
Triali = Expi;%Expi - 34;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpSpecTable_Aug.Expi(exp_rowi);
fprintf("Processing Exp %d, %s\n", Expi, meta.comments)
% Prepare the relevant folder and info
pref_chan = Trials.TrialRecord.User.prefChan;
pref_chan_id = find(meta.spikeID==pref_chan);
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan);
savepath = fullfile(result_dir, compose("Manifold_Evol%02d_chan%02d", Expi, pref_chan));
mkdir(savepath);
unit_name_arr = generate_unit_labels(meta.spikeID, savepath); % Generate readable labels for each channel 
% Compute the block structure
imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & contains(imgnm, "block") & cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm);
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
block_arr = zeros(numel(imgnm), 1);
generations = zeros(numel(imgnm), 1);
blocki = 0;geni = 0;
for i = 1:numel(imgnm)
    if row_gen(i)
        matchstr = regexp(imgnm{i}, "block(?<blocki>\d\d\d)_thread(?<threadi>\d\d\d)_gen_(?<imgname>.*)",'names');
        if str2num(matchstr.blocki) == blocki
            % still in same block
        else
            blocki = str2num(matchstr.blocki);
            geni = blocki - 1;
        end
    end
    block_arr(i) = blocki;
    generations(i) = geni;
end
%% Generate Gradual Changing Color Labels 
gen_list = min(block_arr):max(block_arr); % list of generations
% fun = @(m)srgb_to_Lab(m);
% color_seq = maxdistcolor(30,fun); % Use this to generate maximally distinguishable color sequence
color_seq = brewermap(length(gen_list), 'spectral'); 
% Use this to generate gradual changing sequence. match the number of generations 
%% Compute scores for each channel and every image. 
scores_tsr = squeeze(mean(rasters(:, 151:200, :), 2) - mean(rasters(:, 1:40, :), 2));
meanscore_syn = [];
stdscore_syn = [];
meanscore_nat = [];
stdscore_nat = [];
for blocki = min(block_arr):max(block_arr)
    tmpscore_syn = mean(scores_tsr(:, row_gen & block_arr == blocki), 2);cloase   
    tmpscore_nat = mean(scores_tsr(:, row_nat & block_arr == blocki), 2);
    tmpstdscore_syn = std(scores_tsr(:, row_gen & block_arr == blocki), 1, 2) / sqrt(sum(row_gen & block_arr == blocki));
    tmpstdscore_nat = std(scores_tsr(:, row_nat & block_arr == blocki), 1, 2) / sqrt(sum(row_gen & block_arr == blocki));
    meanscore_syn = cat(2, meanscore_syn, tmpscore_syn);
    meanscore_nat = cat(2, meanscore_nat, tmpscore_nat);
    stdscore_syn  = cat(2, stdscore_syn, tmpstdscore_syn);
    stdscore_nat  = cat(2, stdscore_nat, tmpstdscore_nat);
end
%% Compute mean and std PSTH within generation for every channel
gen_list = min(block_arr):max(block_arr);
evol_stim_fr = zeros(size(rasters, 1), size(rasters, 2), length(gen_list));
evol_stim_sem = zeros(size(rasters, 1), size(rasters, 2), length(gen_list));
for blocki = 1:length(gen_list)
    evol_stim_fr(:, :, blocki) = mean(rasters(:,:, row_gen & block_arr == blocki),3);
    evol_stim_sem(:, :, blocki) = std(rasters(:,:, row_gen & block_arr == blocki),1,3)/sqrt(sum(row_gen & block_arr == blocki));
end
% nat_stim_fr = zeros(max(natural_stim_i), size(rasters, 2), size(rasters, 3));
% nat_stim_fr_std = zeros(max(natural_stim_i), size(rasters, 2), size(rasters, 3));
% nat_stim_fr_sem = zeros(max(natural_stim_i), size(rasters, 2), size(rasters, 3));
% for i = 1:max(natural_stim_i)
%     nat_stim_fr(i,:,:) =    mean(rasters(sort_idx(natural_stim_i==i), :, :),1);
%     nat_stim_fr_std(i,:,:) = std(rasters(sort_idx(natural_stim_i==i), :, :),1,1);
%     nat_stim_fr_sem(i,:,:) = nat_stim_fr_std(i,:,:) / sqrt(sum(natural_stim_i==i));
% end
%% Prepare the axes for 3 arrays 
[chan_idxA, chan_idxB] = unit_id2_chan_idx(1:64, meta.spikeID); 
% separate the spikeID into unit A group and unit B group 
Exp_label = sprintf("Exp %d Evolv channel %d", Expi, pref_chan);
figIT = figure('Visible','off');clf;hold on
figV1 = figure('Visible','off');clf;hold on
figV4 = figure('Visible','off');clf;hold on
[ax_arrA,tIT,tV1,tV4] = Cortex_Channel_Tile_Layout_All(figIT, figV1, figV4);
if all(chan_idxA==chan_idxB) % each channel has only 1 unit exactly
    plot2unit = false;
else % some channels have more than 1 unit
    figITB = figure('Visible','off');clf;hold on
    figV1B = figure('Visible','off');clf;hold on
    figV4B = figure('Visible','off');clf;hold on
    [ax_arrB,tITB,tV1B,tV4B] = Cortex_Channel_Tile_Layout_All(figITB, figV1B, figV4B);
    plot2unit = true;
end
%% Plot the figure to target axis 
for arr_chan = 1:64 % array channel! not number in the resulting array
    channel = chan_idxA(arr_chan); % map the array channel id to the matlab channel number
    if isnan(channel) % if the channel has no unit. 
        set(ax_arrA{arr_chan},'Visible','off')
        set(ax_arrB{arr_chan},'Visible','off')
        continue  % then hide the axis and continue
    end
    chan_label_str = sprintf("Ch %s", unit_name_arr{channel} );
    plot_score_error_traj(channel, ax_arrA{arr_chan})
    title(chan_label_str)
    if plot2unit % if more than one unit in this channel
        channel = chan_idxB(arr_chan); % the c
        chan_label_str = sprintf("Ch %s", unit_name_arr{channel} );
        plot_score_error_traj(channel, ax_arrB{arr_chan})
        title(chan_label_str)
    end
end
title(tIT,strcat(Exp_label, " IT array"))
title(tV1,strcat(Exp_label, " V1V2 array"))
title(tV4,strcat(Exp_label, " V4 array"))
if plot2unit
    title(tITB,strcat(Exp_label, " IT array"))
    title(tV1B,strcat(Exp_label, " V1V2 array"))
    title(tV4B,strcat(Exp_label, " V4 array"))
end
% Save everything in the exp folder and the overall folder
saveas(figIT, fullfile(savepath, sprintf("IT_array_Exp%02d_score_trajA.jpg", Expi)))
saveas(figIT, fullfile(result_dir, sprintf("IT_array_Exp%02d_score_trajA.jpg", Expi)))
saveas(figV1, fullfile(savepath, sprintf("V1_array_Exp%02d_score_trajA.jpg", Expi)))
saveas(figV1, fullfile(result_dir, sprintf("V1_array_Exp%02d_score_trajA.jpg", Expi)))
saveas(figV4, fullfile(savepath, sprintf("V4_array_Exp%02d_score_trajA.jpg", Expi)))
saveas(figV4, fullfile(result_dir, sprintf("V4_array_Exp%02d_score_trajA.jpg", Expi)))
if plot2unit
saveas(figITB, fullfile(savepath, sprintf("IT_array_Exp%02d_score_trajB.jpg", Expi)))
saveas(figITB, fullfile(result_dir, sprintf("IT_array_Exp%02d_score_trajB.jpg", Expi)))
saveas(figV1B, fullfile(savepath, sprintf("V1_array_Exp%02d_score_trajB.jpg", Expi)))
saveas(figV1B, fullfile(result_dir, sprintf("V1_array_Exp%02d_score_trajB.jpg", Expi)))
saveas(figV4B, fullfile(savepath, sprintf("V4_array_Exp%02d_score_trajB.jpg", Expi)))
saveas(figV4B, fullfile(result_dir, sprintf("V4_array_Exp%02d_score_trajB.jpg", Expi)))
end
%% Plot the figure to target axis 
for arr_chan = 1:64 % array channel! not number in the resulting array
    channel = chan_idxA(arr_chan); % map the array channel id to the matlab channel number
    if isnan(channel) % if the channel has no unit. 
        set(ax_arrA{arr_chan},'Visible','off')
        set(ax_arrB{arr_chan},'Visible','off')
        continue  % then hide the axis and continue
    end
    chan_label_str = sprintf("Ch %s", unit_name_arr{channel} );
    plot_PSTH_change(channel, ax_arrA{arr_chan})
    title(chan_label_str)
    if plot2unit % if more than one unit in this channel
        channel = chan_idxB(arr_chan); % the c
        chan_label_str = sprintf("Ch %s", unit_name_arr{channel} );
        plot_PSTH_change(channel, ax_arrB{arr_chan})
        title(chan_label_str)
    end
end
title(tIT,strcat(Exp_label, " IT array"))
title(tV1,strcat(Exp_label, " V1V2 array"))
title(tV4,strcat(Exp_label, " V4 array"))
if plot2unit
    title(tITB,strcat(Exp_label, " IT array"))
    title(tV1B,strcat(Exp_label, " V1V2 array"))
    title(tV4B,strcat(Exp_label, " V4 array"))
end
% Save everything in the exp folder and the overall folder
saveas(figIT, fullfile(savepath, sprintf("IT_array_Exp%02d_PSTH_trajA.jpg", Expi)))
saveas(figIT, fullfile(result_dir, sprintf("IT_array_Exp%02d_PSTH_trajA.jpg", Expi)))
saveas(figV1, fullfile(savepath, sprintf("V1_array_Exp%02d_PSTH_trajA.jpg", Expi)))
saveas(figV1, fullfile(result_dir, sprintf("V1_array_Exp%02d_PSTH_trajA.jpg", Expi)))
saveas(figV4, fullfile(savepath, sprintf("V4_array_Exp%02d_PSTH_trajA.jpg", Expi)))
saveas(figV4, fullfile(result_dir, sprintf("V4_array_Exp%02d_PSTH_trajA.jpg", Expi)))
if plot2unit
saveas(figITB, fullfile(savepath, sprintf("IT_array_Exp%02d_PSTH_trajB.jpg", Expi)))
saveas(figITB, fullfile(result_dir, sprintf("IT_array_Exp%02d_PSTH_trajB.jpg", Expi)))
saveas(figV1B, fullfile(savepath, sprintf("V1_array_Exp%02d_PSTH_trajB.jpg", Expi)))
saveas(figV1B, fullfile(result_dir, sprintf("V1_array_Exp%02d_PSTH_trajB.jpg", Expi)))
saveas(figV4B, fullfile(savepath, sprintf("V4_array_Exp%02d_PSTH_trajB.jpg", Expi)))
saveas(figV4B, fullfile(result_dir, sprintf("V4_array_Exp%02d_PSTH_trajB.jpg", Expi)))
end



end
function plot_score_scatter_traj(channel, ax)
    % Given a averaged score matrix, plot a contour heatmap to certain
    % axis! 
    global block_arr gen_list meanscore_syn stdscore_syn meanscore_nat stdscore_nat row_gen row_nat
    set(0, "CurrentFigure", ancestor(ax,'figure'))
    set(gcf, "CurrentAxes", ax);cla(ax);hold on % 
    gen_list = min(block_arr):max(block_arr);
    scatter(ax, block_arr(row_gen), scores_tsr(channel, row_gen))
    scatter(ax, block_arr(row_nat), scores_tsr(channel, row_nat))
    plot(ax, gen_list, meanscore_syn(channel, :), 'LineWidth',2,'Color','k')
    plot(ax, gen_list, meanscore_nat(channel, :),'LineWidth',2,'Color','g')
    legend(["Generated img","Natural img","Gen mean","Nat mean"])
    xlabel("generations")
end
function plot_score_error_traj(channel, ax)
    % Given a averaged score matrix, plot a contour heatmap to certain
    % axis! 
    global block_arr gen_list meanscore_syn stdscore_syn meanscore_nat stdscore_nat
    set(0, "CurrentFigure", ancestor(ax,'figure'))
    set(gcf, "CurrentAxes", ax);cla(ax);hold on 
    axis normal;
    gen_list = min(block_arr):max(block_arr);
    errorbar(ax, gen_list, meanscore_syn(channel, :), stdscore_syn(channel, :), 'LineWidth',2,'Color',[0,0,0,0.5])
    errorbar(ax, gen_list, meanscore_nat(channel, :), stdscore_nat(channel, :), 'LineWidth',2,'Color',[0,1,0,0.5])
    legend(["Gen","Nat"])
    xlabel("generations")
    hold off
end
function plot_PSTH_change(channel, ax)
    % plot PSTH over generations
    global block_arr gen_list color_seq evol_stim_fr evol_stim_sem
    set(0, "CurrentFigure", ancestor(ax,'figure'))
    set(gcf, "CurrentAxes", ax);cla(ax);hold on 
    axis normal;
    for i = 1:length(gen_list)
        shadedErrorBar([],evol_stim_fr(channel, :, i),evol_stim_sem(channel, :, i),...
        'lineprops',{'Color',[color_seq(i, :),0.85]},'transparent',1,'patchSaturation',0.075)
    end
    YL=ylim;YL(1)=0;ylim(YL);
    XL=xlim;XL(1)=0;xlim(XL);
    xlabel("time (ms)")
    legend(gca,'off');
end
