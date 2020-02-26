%% Animation of Evolution 
expftr = contains(ExpSpecTable_Aug.expControlFN,"generate") & ...
     ExpSpecTable_Aug.Expi==7 &...
     contains(ExpSpecTable_Aug.Exp_collection, "Manifold");
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr)); 
%%
Triali=1; 
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

imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
block_arr = cell2mat(Trials.block);
block_list = min(block_arr):max(block_arr);
color_seq = brewermap(length(block_list), 'GnBu'); 

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
block_list = min(block_arr):max(block_arr);
evol_stim_fr = zeros(size(rasters, 1), size(rasters, 2), length(block_list));
evol_stim_sem = zeros(size(rasters, 1), size(rasters, 2), length(block_list));
for blocki = block_list
    gen_msk = row_gen & block_arr == blocki;
    evol_stim_fr(:, :, blocki) = mean(rasters(:,:, gen_msk),3);
    evol_stim_sem(:, :, blocki) = std(rasters(:,:, gen_msk),1,3)/sqrt(sum(gen_msk));
end

%%
channel_j = pref_chan_id(pref_chan_unit);
imgColl = repmat("", length(block_list),1);
scoreColl = zeros(length(block_list),1); % have to be same size with imgColl. 
for blocki = block_list 
    gen_msk = row_gen & block_arr == blocki; 
    [maxScore, maxIdx] = max(scores_tsr(channel_j, gen_msk));
    tmpimgs = imgnm(gen_msk);
    imgfullfn = ls(fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(meta.stimuli, [tmpimgs(maxIdx)+"*"]))
    imgColl(blocki) = fullfile(meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % FIXED @ Jan.29
    scoreColl(blocki) = maxScore;
end
%%
%Fs = {};
h=figure(1);clf;
for blocki = block_list 
    imshow(imgColl(blocki)); 
    title(compose("Evoked Rate %.1f", scoreColl(blocki))) ; 
    Fs(blocki) = getframe(h);
end
%%
v = VideoWriter(fullfile(savepath,'Evol_Best_PSTH.avi'));
v.FrameRate = 3;
open(v);
h2=figure(2);clf;
subplot(212);
ylabel("PSTH (Hz)")
xlabel("time(ms)")
title("Evoked PSTH")
ylim([0,250])
for blocki = block_list(1:end-1)
    subplot(211);
    imshow(imgColl(blocki)); 
    title(compose("Evoked Rate %.1f", scoreColl(blocki))) ; 
    subplot(212);hold on
    shadedErrorBar([],evol_stim_fr(channel_j, :, blocki),evol_stim_sem(channel_j, :, blocki),...
    'lineprops',{'Color',[color_seq(blocki, :),0.85]},'transparent',1,'patchSaturation',0.15)
    drawnow;
    Fs(blocki) = getframe(h2);
    writeVideo(v,Fs(blocki));
end
close(v);
% % Montage the images
% set(0,'CurrentFigure',h); %clf; %
% montage(imgColl(:))
% title([Exp_label_str, compose('Best Image per Generation')])
% saveas(h, fullfile(savepath, "EvolImageSeq.png"))