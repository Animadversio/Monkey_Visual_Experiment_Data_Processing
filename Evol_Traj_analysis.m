%% Standard code for Evolution Exp Analysis
% well editted @Jan 30 for batch processing. 
clearvars -except meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug % keep only the codes store data
%%
Set_Path;
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Evolution_Exp";
%ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx");
%%
% 15 has someth special
expftr = ExpSpecTable_Aug.Expi<=14 & ExpSpecTable_Aug.Expi>=1 & ...
    contains(ExpSpecTable_Aug.expControlFN,"generate") & ...
    contains(ExpSpecTable_Aug.Exp_collection,"Manifold");
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr)); 
%%
h   = figure('Visible','on');h.Position = [828 42 1026 954]; % Evolution Image Sequence 
h0  = figure('Visible','on');h0.Position = [828 42 1026 954]; % Evolution Image Sequence with score frame
h1 = figure('Visible','off'); % score line plot + scatter plot for each trial 
h2 = figure('Visible','off'); % shaded Errorbar of score for each generation 
h3 = figure('Visible','off');h3.Position = [ 1128         314         899         505]; % shaded Errorbar of PSTH for each generation
%
for Triali = 1:length(meta_new)
% Fetch the trial info
%Triali = Expi - 26;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpSpecTable_Aug.ephysFN, meta.ephysFN));
Expi = ExpSpecTable_Aug.Expi(exp_rowi);
fprintf("Exp %d:\n",Expi)
disp(ExpSpecTable_Aug.comments(exp_rowi))
% Prepare the relevant folder and info
pref_chan = Trials.TrialRecord.User.prefChan;
pref_chan_id = find(meta.spikeID==pref_chan);
pref_chan_unit = Trials.TrialRecord.User.evoConfiguration{Trials.TrialRecord.User.iConfig.unit}; % which unit is it evolving to
Exp_label_str = sprintf("Exp%d pref chan %d (unit %d)", Expi, pref_chan, pref_chan_unit);
savepath = fullfile(result_dir, compose("Manifold_Evol%02d_chan%02d", Expi, pref_chan));
mkdir(savepath);

unit_name_arr = generate_unit_labels(meta.spikeID, savepath); % Generate readable labels for each channel
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters, savepath);
% Check if the units are active
% Compute the block structure from imagenames
imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
block_arr = cell2mat(Trials.block);
% generations = zeros(numel(imgnm), 1);
% blocki = 0;geni = 0;
% for i = 1:numel(imgnm)
%     if row_gen(i)
%         matchstr = regexp(imgnm{i}, "block(?<blocki>\d\d\d)_thread(?<threadi>\d\d\d)_gen_(?<imgname>.*)",'names');
%         if str2num(matchstr.blocki) == blocki
%             
%         else
%             blocki = str2num(matchstr.blocki);
%             geni = blocki - 1;
%         end
%     end
%     block_arr(i) = blocki;
%     generations(i) = geni;
% end
%% Generate Gradual Changing Color Labels 
block_list = min(block_arr):max(block_arr);
% fun = @(m)srgb_to_Lab(m);
% color_seq = maxdistcolor(30,fun); % Use this to generate maximally distinguishable color sequence
color_seq = brewermap(length(block_list), 'spectral'); % Use this to generate gradual changing sequence
%% 
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
meanscore_syn = nan(size(rasters, 1), length(block_list)); % [unit_num, gen_nums, threads]
stdscore_syn = nan(size(rasters, 1), length(block_list)); 
meanscore_nat = nan(size(rasters, 1), length(block_list));
stdscore_nat = nan(size(rasters, 1), length(block_list));
for blocki = min(block_arr):max(block_arr)
    gen_msk = row_gen & block_arr == blocki;% & thread_msks{threadi}; 
    nat_msk = row_nat & block_arr == blocki;% & thread_msks{threadi};
    meanscore_syn(:, blocki) = mean(scores_tsr(:, gen_msk), 2);
    meanscore_nat(:, blocki) = mean(scores_tsr(:, nat_msk), 2);
    stdscore_syn(:, blocki)  = std(scores_tsr(:, gen_msk), 1, 2) / sqrt(sum(gen_msk));
    stdscore_nat(:, blocki)  = std(scores_tsr(:, nat_msk), 1, 2) / sqrt(sum(nat_msk));
end
%%
block_list = min(block_arr):max(block_arr);
evol_stim_fr = zeros(size(rasters, 1), size(rasters, 2), length(block_list));
evol_stim_sem = zeros(size(rasters, 1), size(rasters, 2), length(block_list));
for blocki = block_list
    gen_msk = row_gen & block_arr == blocki;
    evol_stim_fr(:, :, blocki) = mean(rasters(:,:, gen_msk),3);
    evol_stim_sem(:, :, blocki) = std(rasters(:,:, gen_msk),1,3)/sqrt(sum(gen_msk));
end
% nat_stim_fr = zeros(max(natural_stim_i), size(rasters, 2), size(rasters, 3));
% nat_stim_fr_std = zeros(max(natural_stim_i), size(rasters, 2), size(rasters, 3));
% nat_stim_fr_sem = zeros(max(natural_stim_i), size(rasters, 2), size(rasters, 3));
% for i = 1:max(natural_stim_i)
%     nat_stim_fr(i,:,:) =    mean(rasters(sort_idx(natural_stim_i==i), :, :),1);
%     nat_stim_fr_std(i,:,:) = std(rasters(sort_idx(natural_stim_i==i), :, :),1,1);
%     nat_stim_fr_sem(i,:,:) = nat_stim_fr_std(i,:,:) / sqrt(sum(natural_stim_i==i));
% end
%% Get the Image FileName Sequence
%EMPTYTHR= 1000;
emptychan = ~activ_msk(pref_chan_id(1)); %sum(rasters(pref_chan_id(1),:,1:10),[2,3]) < EMPTYTHR;
if emptychan
    fprintf("There is an empty unit %d in the prefered channel %d\n", pref_chan_id(1), pref_chan)
end
channel_j = pref_chan_id(emptychan + pref_chan_unit);
imgColl = repmat("", length(block_list),1);
scoreColl = zeros(length(block_list),1); % have to be same size with imgColl. 
for blocki = block_list
    gen_msk = row_gen & block_arr == blocki; 
    [maxScore, maxIdx] = max(scores_tsr(channel_j, gen_msk));
    tmpimgs = imgnm(gen_msk);
    imgfullfn = ls(fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]))
    imgColl(blocki) = fullfile(meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % solved @ Jan.29
    scoreColl(blocki) = maxScore;
end
% Montage the images
set(0,'CurrentFigure',h); %clf; %
montage(imgColl(:))
title([Exp_label_str, compose('Best Image per Generation')])
saveas(h, fullfile(savepath, "EvolImageSeq.png"))

set(0,'CurrentFigure',h0);
cLim = [min(scoreColl), max(scoreColl)];
frame_img_list = score_frame_image_arr(imgColl, scoreColl, cLim, parula, 20);
montage(frame_img_list)
colorbar();caxis(cLim)
title([Exp_label_str, compose('Best Image per Generation (Best Score frame)')])
saveas(h0, fullfile(savepath, "EvolImageSeq_BestScore.png"))
for channel_j = 1:size(rasters, 1) %pref_chan_id; % []
%channel_j = pref_chan_id;
set(0,'CurrentFigure',h1); clf; hold on 
scatter(block_arr(row_gen), scores_tsr(channel_j, row_gen))
scatter(block_arr(row_nat), scores_tsr(channel_j, row_nat))
plot(block_list, meanscore_syn(channel_j, :), 'LineWidth',2,'Color','k')
plot(block_list, meanscore_nat(channel_j, :),'LineWidth',2,'Color','g')
legend(["Generated img","Natural img","Gen mean","Nat mean"])
xlabel("generations")
title([Exp_label_str, compose('PSTH averaged scores, channel %s', unit_name_arr{channel_j})])
%saveas(h1,fullfile(savepath,compose("score_traj_chan%d.png",channel_j)))
saveas(h1,fullfile(savepath,compose("score_traj_chan%s.png",unit_name_arr{channel_j}))) % update @Feb.5 to use new unit name label on filename
%h1.Visible='on';
%%
set(0,'CurrentFigure',h2); clf; hold on %
shadedErrorBar(block_list(1:end-1), meanscore_syn(channel_j, 1:end-1), stdscore_syn(channel_j, 1:end-1),...
    'lineprops',{'Color',[0,0,0,0.7]},'transparent',1,'patchSaturation',0.075)
shadedErrorBar(block_list(1:end-1), meanscore_nat(channel_j, 1:end-1), stdscore_nat(channel_j, 1:end-1),...
    'lineprops',{'Color',[0,1,0,0.7]},'transparent',1,'patchSaturation',0.075)
axis tight
legend(["Generated img","Natural img"])
xlabel("generations")
title([Exp_label_str, compose('Generation averaged score, channel %s', unit_name_arr{channel_j})])
%saveas(h2,fullfile(savepath,compose("score_traj_std_chan%d.png",channel_j)))
saveas(h2,fullfile(savepath,compose("score_traj_std_chan%s.png",unit_name_arr{channel_j}))) % update @Feb.5 to use new unit name label on filename
%h2.Visible='on';
%%
set(0,'CurrentFigure',h3); clf;hold on;
block_list = min(block_arr):max(block_arr);
for i = block_list(1:end-1)
    shadedErrorBar([],evol_stim_fr(channel_j, :, i),evol_stim_sem(channel_j, :, i),...
    'lineprops',{'Color',[color_seq(i, :),0.85]},'transparent',1,'patchSaturation',0.075)
    % plot(evol_stim_fr(i,:,channel_j))
end
YL=ylim;YL(1)=0;ylim(YL);
XL=xlim;XL(1)=0;xlim(XL);
xlabel("time (ms)")
title([Exp_label_str, compose('Generation averaged PSTH of Evolved Stimuli channel %s', unit_name_arr{channel_j})])%num2str(meta.spikeID(pref_chan_id))
%saveas(h3,fullfile(savepath,compose("Evolv_psth_chan%d.png",channel_j)))
saveas(h3,fullfile(savepath,compose("Evolv_psth_chan%s.png",unit_name_arr{channel_j}))) % update @Feb.5 to use new unit name label on filename
hold off
%h3.Visible='on';
end
end
