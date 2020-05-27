%% Evol_Stats_Animation 
% adapted from Evol_Animation Suited for Stats and EStats
Animal = "Alfa";
MatStats_path = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
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
Expi = 3;
ui=1;
Window = 50:200;
imgColl = repmat("", EStats(3).evol.block_n-1,1);
scoreColl = zeros(EStats(3).evol.block_n-1,1);
for blocki = 1:EStats(3).evol.block_n-1
    gen_scores = squeeze(mean(EStats(3).evol.psth{blocki}(:,Window,:),[1,2]));
    [maxScore, maxIdx] = max(gen_scores); % Idx in local generation
    imgidx = EStats(3).evol.idx_seq{blocki}(maxIdx);
    imgfullfn = ls(fullfile(EStats(3).meta.stimuli, EStats(3).imageName(imgidx)+"*"));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(EStats(3).meta.stimuli, EStats(3).imageName(imgidx)+"*"))
    imgColl(blocki) = fullfile(EStats(3).meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % FIXED @ Jan.29
    scoreColl(blocki) = maxScore;
end
evol_stim_fr = cellfun(@(psth)mean(psth,3),EStats(3).evol.psth,'UniformOutput',false);
evol_stim_fr = cell2mat(reshape(evol_stim_fr,1,1,[]));
evol_stim_sem = cellfun(@(psth)std(psth,0,3)/sqrt(size(psth,3)),EStats(3).evol.psth,'UniformOutput',false);
evol_stim_sem = cell2mat(reshape(evol_stim_sem,1,1,[]));
%%
result_dir = "E:\OneDrive - Washington University in St. Louis\Evol_Manif_Movies";
savepath = result_dir;%fullfile(result_dir, compose("%s_Evol_"))
color_seq = EStats(3).color_seq;
v = VideoWriter(fullfile(savepath,compose('%s_Evol_Exp%02d_Best_PSTH.mov',Animal,Expi)));
v.FrameRate = 3;
open(v);
h2=figure(2);clf;
subplot(212);
ylabel("PSTH (Hz)")
xlabel("time(ms)")
title("Evoked PSTH")
ylim([0,max(evol_stim_fr+evol_stim_sem,[],'all')])
for blocki = 1:EStats(3).evol.block_n-1
    subplot(211);
    imshow(imgColl(blocki)); 
    title(compose("Gen%d Evoked Rate %.1f", blocki, scoreColl(blocki))) ; 
    subplot(212);%hold on
    shadedErrorBar([],evol_stim_fr(ui, :, blocki),evol_stim_sem(ui, :, blocki),...
        'lineprops',{'Color',[color_seq(blocki, :),0.85]},'transparent',1,'patchSaturation',0.15)
    drawnow;
    Fs(blocki) = getframe(h2);
    writeVideo(v,Fs(blocki));
end
close(v);