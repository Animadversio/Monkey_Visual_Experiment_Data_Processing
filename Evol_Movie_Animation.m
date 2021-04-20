%% Evol_Movie Animation
movdir = "O:\Evol_Movies_Anim";
mkdir(movdir)
% for Triali = 1:2 % 3:numel(meta_new)
for Triali = 1:2
%%
rasters = rasters_new{Triali};
meta = meta_new{Triali};
Trials = Trials_new{Triali};
wdw = meta.rasterWindow;
%%
pref_chan = Trials.TrialRecord.User.prefChan;
assert(all(pref_chan==pref_chan(1)))
pref_chan=pref_chan(1);
stimparts = split(meta.stimuli,"\");
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
% fdrnm = compose("%s-%s-Chan%02d-1",datestr(expday,'yyyy-mm-dd'), Animal, pref_chan(1));
fdrnm = compose("%s-Chan%02d", stimparts{end-1}, pref_chan(1));
figdir = fullfile(movdir, fdrnm);
fprintf("Going to create folder %s\n",fdrnm)
mkdir(figdir)
%% Collect basic info (adapted from Evol_Collect_Stats)
pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
imgpos = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,2));
imgsize = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,3));
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
unit_num_arr = meta.unitID;
pref_chan_id = zeros(1, thread_num);
for i = 1:thread_num % note here we use the unit_num_arr which exclude the null channels.
    pref_chan_id(i) = find(meta.spikeID==pref_chan(i) & ... % the id in the raster and lfps matrix 
                    unit_num_arr==unit_in_pref_chan(i)); % match for unit number
end
if size(Trials.TrialRecord.User.evoConfiguration,2) == 5
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, strrep(string(Trials.TrialRecord.User.evoConfiguration{i,end}),"_"," ")]; % use this substitution for better printing effect
end
elseif size(Trials.TrialRecord.User.evoConfiguration,2) == 4 % old fashioned config
    Optim_names = repmat("CMAES", 1,thread_num);
end
% separate nat, gen
imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm) & ...% first 2 characters are not digits
        cellfun(@(c) ~contains(c(end-4:end), "_nat"), imgnm) ; % last 4 characters are not `_nat`
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the threads 1 and 2 (and maybe thread 3 4)
thread_msks = cell(1, thread_num);
for threadi = 1:thread_num
    msk = contains(imgnm, compose("thread%03d", threadi - 1));
    thread_msks{threadi} = msk; % store masks in a structure for the ease to iterate
end
assert(sum(cellfun(@sum, thread_msks)) == length(imgnm))
% separate block, nat, gen
block_arr = cell2mat(Trials.block);
block_list = min(block_arr):max(block_arr);
block_num = length(block_list);
gen_idx_seq = cell(thread_num, block_num); % generated image idx cell as a horizontal array. 
nat_idx_seq = cell(thread_num, block_num);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        gen_idx_seq{threadi, blocki} = find(gen_msk);
        nat_idx_seq{threadi, blocki} = find(nat_msk);
    end
end
%%
isMovie = contains(meta.expControlFN,"Movie");
if isMovie
bslwdw = [-50:45] - wdw(1); actwdw = [51:400] - wdw(1);
else
bslwdw = [1:45] - wdw(1); actwdw = [51:200] - wdw(1);
end
meanscore_syn = cellfun(@(idx) mean(rasters(:,actwdw,idx),[2,3]) - mean(rasters(:,bslwdw,idx),[2,3]), gen_idx_seq', 'uni', 0);
semscore_syn = cellfun(@(idx) squeeze(std(mean(rasters(:,actwdw,idx),2),1,3)) / sqrt(numel(idx)), gen_idx_seq', 'uni', 0);
meanscore_nat = cellfun(@(idx) mean(rasters(:,actwdw,idx),[2,3]) - mean(rasters(:,bslwdw,idx),[2,3]), nat_idx_seq', 'uni', 0);
semscore_nat = cellfun(@(idx) squeeze(std(mean(rasters(:,actwdw,idx),2),1,3)) / sqrt(numel(idx)), nat_idx_seq', 'uni', 0);
meanscore_syn = cell2mat(reshape(meanscore_syn,[],size(gen_idx_seq',1),size(gen_idx_seq',2)));
semscore_syn = cell2mat(reshape(semscore_syn,[],size(gen_idx_seq',1),size(gen_idx_seq',2)));
meanscore_nat = cell2mat(reshape(meanscore_nat,[],size(nat_idx_seq',1),size(nat_idx_seq',2)));
semscore_nat = cell2mat(reshape(semscore_nat,[],size(nat_idx_seq',1),size(nat_idx_seq',2)));
%%
imgfullfn = ls(fullfile(meta.stimuli, Trials.imageName(find(row_gen,1))+"*"));
imgparts = split(imgfullfn,".");
suffix = "."+imgparts{2};
%     imgfullfn = ls(fullfile(meta.stimuli, Trials.imageName(imgidx)+"*"));
%     assert(~isempty(imgfullfn), "Image not found %s",fullfile(meta.stimuli, Trials.imageName(imgidx)+"*"))
imgColl = repmat("", block_num-1,1);
scoreColl = zeros(block_num-1,1);
for blocki = 1:block_num-1
    gen_scores = squeeze(mean(rasters(pref_chan_id, actwdw, gen_idx_seq{blocki}),[2]));
    [maxScore, maxIdx] = max(gen_scores); % Idx in local generation
    imgidx = gen_idx_seq{blocki}(maxIdx);
    imgfullfn= fullfile(meta.stimuli, Trials.imageName(imgidx)+suffix);
    imgColl(blocki) = imgfullfn;%fullfile(meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % FIXED @ Jan.29
    scoreColl(blocki) = maxScore;
end
%%
score_avg = meanscore_syn(pref_chan_id, :);
score_sem = semscore_syn(pref_chan_id, :);
ctrls_avg = meanscore_nat(pref_chan_id, :);
ctrls_sem = semscore_nat(pref_chan_id, :);
psth_col = cellfun(@(idx)rasters(pref_chan_id,:,idx), gen_idx_seq, 'uni', false);
evol_stim_fr = cellfun(@(psth)mean(psth,3),psth_col,'Uni',0);
evol_stim_fr = cell2mat(reshape(evol_stim_fr,1,1,[]));
evol_stim_sem = cellfun(@(psth)std(psth,0,3)/sqrt(size(psth,3)),psth_col,'Uni',0);
evol_stim_sem = cell2mat(reshape(evol_stim_sem,1,1,[]));
%%
savepath = movdir; % fullfile(movdir, compose("%s_Evol_"))
clrseq = brewermap(block_num,'Spectral');%color_seq = EStats(Expi).color_seq;
v = VideoWriter(fullfile(movdir,compose("%s_Evol_Best_PSTH",fdrnm)));%'%s_Evol_Exp%02d_Best_PSTH.mov',Animal,Expi)));
if isMovie, v.FrameRate = 15; else, v.FrameRate = 3;  end
open(v);
h2=figure(2);set(2,'position',[1258         237         662         735]);clf;
ax1 = subplot(2,2,1);
set(ax1,"position",[0.07,0.578,0.439,0.3412]);
ax3 = subplot(2,2,2);
shadedErrorBar([],score_avg(1:end-1),score_sem(1:end-1),...
        'lineprops',{'Color','k','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
shadedErrorBar([],ctrls_avg(1:end-1),ctrls_sem(1:end-1),...
        'lineprops',{'Color','g','LineWidth',1},'transparent',1,'patchSaturation',0.05);
axis tight;
xlabel("Generations");ylabel("Response fr (Hz)");title("Evolution Traj")
ax2 = subplot(2,1,2);
bgpsth = plot(wdw(1)+1:wdw(2),squeeze(evol_stim_fr(:,:,1:end-1)),'Color',[0.7,0.7,0.7]);hold on
sEB = shadedErrorBar(wdw(1)+1:wdw(2),evol_stim_fr(1, :, 1),evol_stim_sem(1, :, 1),...
        'lineprops',{'Color','r','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
ylabel("PSTH (Hz)");xlabel("time(ms)");xlim([-50,400])
ylim([0,max(evol_stim_fr(:,:,1:end-1)+evol_stim_sem(:,:,1:end-1),[],'all')])
title("Evoked PSTH")
ST = suptitle(compose("%s Evol pref chan %s", ...
    fdrnm, unit_name_arr(pref_chan_id)));%Animal, 
GRAYSCREEN = 0.5*ones(256,256,3);
for blocki = 1:block_num-1 
    set(h2,"CurrentAxes",ax2)
    ax2.Title.String = compose("Evoked PSTH Gen%d Evoked Rate %.1f", blocki, scoreColl(blocki));
%     shadedErrorBar([],evol_stim_fr(ui, :, blocki),evol_stim_sem(ui, :, blocki),...
%         'lineprops',{'Color',[color_seq(blocki, :),0.85]},'transparent',1,'patchSaturation',0.15)
    psthcur = evol_stim_fr(1, :, blocki);
    psthsemcur = evol_stim_sem(1, :, blocki);
    uE = psthcur + psthsemcur; lE = psthcur - psthsemcur;
    sEB.mainLine.YData = psthcur;
    sEB.edge(1).YData = lE;
    sEB.edge(2).YData = uE;
    sEB.patch.Vertices(:,2) = [lE,fliplr(uE)];
    set(h2,"CurrentAxes",ax1)%subplot(2,2,1);
    title(compose("Gen%d", blocki)) ;
    if isMovie
    vid = VideoReader(imgColl(blocki));
    imgs = vid.read();%close(vid)
    for fi = 1:vid.NumFrames
    imshow(imgs(:,:,:,fi)); 
    drawnow;
    Fs = getframe(h2);
    writeVideo(v,Fs);
    end
    imshow(GRAYSCREEN); 
    else
    imshow(imgColl(blocki));
    end
    drawnow;
    Fs = getframe(h2);
    writeVideo(v,Fs);
end
close(v);
end