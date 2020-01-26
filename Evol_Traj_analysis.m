%% Standard code for Evolution Exp Analysis
clearvars -except meta_new rasters_new lfps_new Trials_new % keep only the codes store data
%%
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Evolution_Exp";
Expi = 28:33;

for Expi = 27
% Fetch the trial info
Triali = Expi - 26;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
% Prepare the relevant folder and info
pref_chan = Trials.TrialRecord.User.prefChan;
pref_chan_id = find(meta.spikeID==pref_chan);
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan);
savepath = fullfile(result_dir, compose("Manifold_Evol%02d_chan%02d", Expi, pref_chan));
mkdir(savepath);
unit_name_arr = generate_unit_labels(meta.spikeID, savepath); % Generate readable labels for each channel
% Compute the block structure from imagenames
imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
block_arr = zeros(numel(imgnm), 1);
generations = zeros(numel(imgnm), 1);
blocki = 0;geni = 0;
for i = 1:numel(imgnm)
    if row_gen(i)
        matchstr = regexp(imgnm{i}, "block(?<blocki>\d\d\d)_thread(?<threadi>\d\d\d)_gen_(?<imgname>.*)",'names');
        if str2num(matchstr.blocki) == blocki
            
        else
            blocki = str2num(matchstr.blocki);
            geni = blocki - 1;
        end
    end
    block_arr(i) = blocki;
    generations(i) = geni;
end
%% Generate Gradual Changing Color Labels 
gen_list = min(block_arr):max(block_arr);
% fun = @(m)srgb_to_Lab(m);
% color_seq = maxdistcolor(30,fun); % Use this to generate maximally distinguishable color sequence
color_seq = brewermap(length(gen_list), 'spectral'); % Use this to generate gradual changing sequence
%% 
scores_tsr = squeeze(mean(rasters(:, 151:200, :), 2) - mean(rasters(:, 1:40, :), 2));
meanscore_syn = [];
stdscore_syn = [];
meanscore_nat = [];
stdscore_nat = [];
for blocki = min(block_arr):max(block_arr)
    tmpscore_syn = mean(scores_tsr(:, row_gen & block_arr == blocki), 2);
    tmpscore_nat = mean(scores_tsr(:, row_nat & block_arr == blocki), 2);
    tmpstdscore_syn = std(scores_tsr(:, row_gen & block_arr == blocki), 1, 2) / sqrt(sum(row_gen & block_arr == blocki));
    tmpstdscore_nat = std(scores_tsr(:, row_nat & block_arr == blocki), 1, 2) / sqrt(sum(row_gen & block_arr == blocki));
    meanscore_syn = cat(2, meanscore_syn, tmpscore_syn);
    meanscore_nat = cat(2, meanscore_nat, tmpscore_nat);
    stdscore_syn  = cat(2, stdscore_syn, tmpstdscore_syn);
    stdscore_nat  = cat(2, stdscore_nat, tmpstdscore_nat);
end
% clear tmpscore_syn tmpscore_nat tmpstdscore_syn tmpstdscore_nat
%%
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
h1 = figure('Visible','off');
h2 = figure('Visible','off'); 
h3 = figure('Visible','off');h3.Position = [ 794         328        1233         463];
for channel_j = 1:size(rasters, 1) %pref_chan_id;
%channel_j = pref_chan_id;
set(0,'CurrentFigure',h1); clf; hold on 
scatter(block_arr(row_gen), scores_tsr(channel_j, row_gen))
scatter(block_arr(row_nat), scores_tsr(channel_j, row_nat))
plot(min(block_arr):max(block_arr), meanscore_syn(channel_j, :), 'LineWidth',2,'Color','k')
plot(min(block_arr):max(block_arr), meanscore_nat(channel_j, :),'LineWidth',2,'Color','g')
legend(["Generated img","Natural img","Gen mean","Nat mean"])
xlabel("generations")
title([Exp_label_str, compose('PSTH averaged scores, channel %s', unit_name_arr{channel_j})])
saveas(h1,fullfile(savepath,compose("score_traj_chan%d.png",channel_j)))
%h1.Visible='on';
%%
set(0,'CurrentFigure',h2); clf; hold on %
errorbar(min(block_arr):max(block_arr), meanscore_syn(channel_j, :), stdscore_syn(channel_j, :), 'LineWidth',2,'Color',[0,0,0,0.5])
errorbar(min(block_arr):max(block_arr), meanscore_nat(channel_j, :), stdscore_nat(channel_j, :), 'LineWidth',2,'Color',[0,1,0,0.5])
legend(["Generated img","Natural img"])
xlabel("generations")
title([Exp_label_str, compose('PSTH averaged scores, channel %s', unit_name_arr{channel_j})])
saveas(h2,fullfile(savepath,compose("score_traj_std_chan%d.png",channel_j)))
%h2.Visible='on';
%%
set(0,'CurrentFigure',h3); clf;hold on;
gen_list = min(block_arr):max(block_arr);
for i = 1:length(gen_list)
    shadedErrorBar([],evol_stim_fr(channel_j, :, i),evol_stim_sem(channel_j, :, i),...
    'lineprops',{'Color',[color_seq(i, :),0.85]},'transparent',1,'patchSaturation',0.075)
    % plot(evol_stim_fr(i,:,channel_j))
end
YL=ylim;YL(1)=0;ylim(YL);
XL=xlim;XL(1)=0;xlim(XL);
xlabel("time (ms)")
title([Exp_label_str, compose('Generation averaged PSTH of Evolved Stimuli channel %s', unit_name_arr{channel_j})])%num2str(meta.spikeID(pref_chan_id))
saveas(h3,fullfile(savepath,compose("Evolv_psth_chan%d.png",channel_j)))
hold off
%h3.Visible='on';
end
end
