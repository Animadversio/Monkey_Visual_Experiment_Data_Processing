%% Evol_Stats_Animation 
% adapted from Evol_Animation Suited for generate movies for Stats and EStats
Animal = "Beto";
MatStats_path = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')

%% Evolution Exp 
for Expi = 27:length(EStats)
fprintf("Processing Evolution Exp %d\n",Expi)
ui=1;
Window = 50:200;
imgColl = repmat("", EStats(Expi).evol.block_n-1,1);
scoreColl = zeros(EStats(Expi).evol.block_n-1,1);
for blocki = 1:EStats(Expi).evol.block_n-1
    gen_scores = squeeze(mean(EStats(Expi).evol.psth{blocki}(:,Window,:),[1,2]));
    [maxScore, maxIdx] = max(gen_scores); % Idx in local generation
    imgidx = EStats(Expi).evol.idx_seq{blocki}(maxIdx);
    imgfullfn = ls(fullfile(EStats(Expi).meta.stimuli, EStats(Expi).imageName(imgidx)+"*"));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(EStats(Expi).meta.stimuli, EStats(Expi).imageName(imgidx)+"*"))
    imgColl(blocki) = fullfile(EStats(Expi).meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % FIXED @ Jan.29
    scoreColl(blocki) = maxScore;
end
evol_stim_fr = cellfun(@(psth)mean(psth,3),EStats(Expi).evol.psth,'UniformOutput',false);
evol_stim_fr = cell2mat(reshape(evol_stim_fr,1,1,[]));
evol_stim_sem = cellfun(@(psth)std(psth,0,3)/sqrt(size(psth,3)),EStats(Expi).evol.psth,'UniformOutput',false);
evol_stim_sem = cell2mat(reshape(evol_stim_sem,1,1,[]));
% Scalor score for evolution
score_avg = cellfun(@(psth)mean(psth(:,51:200,:),'all') - mean(psth(:,1:50,:),'all'),EStats(Expi).evol.psth);
score_sem = cellfun(@(psth)std(squeeze(mean(psth(:,51:200,:),[1,2]) ))...
    /sqrt(size(psth,3)),EStats(Expi).evol.psth);%- mean(psth(:,1:50,:),[1,2])
%% Generate Movies
result_dir = "E:\OneDrive - Washington University in St. Louis\Evol_Manif_Movies";
savepath = result_dir; % fullfile(result_dir, compose("%s_Evol_"))
color_seq = EStats(Expi).color_seq;
v = VideoWriter(fullfile(savepath,compose('%s_Evol_Exp%02d_Best_PSTH.mov',Animal,Expi)));
v.FrameRate = 4;
open(v);
h2=figure(2);set(2,'position',[1258         237         662         735]);clf;
ax1 = subplot(2,2,1);
set(ax1,"position",[0.07,0.578,0.439,0.3412]);
ax3 = subplot(2,2,2);
shadedErrorBar([],score_avg(1:end-1),score_sem(1:end-1),...
        'lineprops',{'Color','k','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
axis tight;
xlabel("Generations");ylabel("Response fr (Hz)");title("Evolution Traj")
ax2 = subplot(2,1,2);
bgpsth = plot(squeeze(evol_stim_fr(:,:,1:end-1)),'Color',[0.7,0.7,0.7]);hold on
sEB = shadedErrorBar([],evol_stim_fr(1, :, 1),evol_stim_sem(1, :, 1),...
        'lineprops',{'Color','r','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
ylabel("PSTH (Hz)");xlabel("time(ms)")
ylim([0,max(evol_stim_fr(:,:,1:end-1)+evol_stim_sem(:,:,1:end-1),[],'all')])
title("Evoked PSTH")
ST = suptitle(compose("%s Evol (Manif) Exp %02d pref chan %s", ...
    Animal, Expi, EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id)));
for blocki = 1:EStats(Expi).evol.block_n-1
    set(h2,"CurrentAxes",ax1)%subplot(2,2,1);
    imshow(imgColl(blocki)); 
    title(compose("Gen%d", blocki)) ; 
    set(h2,"CurrentAxes",ax2)
    ax2.Title.String = compose("Evoked PSTH Gen%d Evoked Rate %.1f", blocki, scoreColl(blocki));
%     shadedErrorBar([],evol_stim_fr(ui, :, blocki),evol_stim_sem(ui, :, blocki),...
%         'lineprops',{'Color',[color_seq(blocki, :),0.85]},'transparent',1,'patchSaturation',0.15)
    psthcur = evol_stim_fr(ui, :, blocki);
    psthsemcur = evol_stim_sem(ui, :, blocki);
    uE = psthcur + psthsemcur; lE = psthcur - psthsemcur;
    sEB.mainLine.YData = psthcur;
    sEB.edge(1).YData = lE;
    sEB.edge(2).YData = uE;
    sEB.patch.Vertices(:,2) = [lE,fliplr(uE)];
    drawnow;
    Fs = getframe(h2);
    writeVideo(v,Fs);
end
close(v);
end

