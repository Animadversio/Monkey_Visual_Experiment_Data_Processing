function Evol_BigGAN_FC6_Animation_fun(EStats)
result_dir = "E:\OneDrive - Washington University in St. Louis\BigGAN_Evol_Movies";
ExpType = "BigGAN_FC6";

for Expi = 1:length(EStats)
% fprintf("Processing BigGAN-FC6 Evolution Exp %d\n",Expi)
Animal = EStats(Expi).Animal;
stimparts = split(EStats(Expi).meta.stimuli,"\");
expday = datetime(EStats(Expi).meta.expControlFN(1:6),'InputFormat','yyMMdd');
fprintf("Processing BigGAN-FC6 Evolution Exp %d\n",Expi)
thread_n = EStats(Expi).evol.thread_num;
% ui = EStats(Expi).evol.unit_in_pref_chan(1);
% assert(all(ui==EStats(Expi).evol.unit_in_pref_chan))
% the following snippet is to get rid of 0 unit (null channel)
pref_chan = EStats(Expi).evol.pref_chan(1);
prefchan_id = find((EStats(Expi).units.spikeID == EStats(Expi).evol.pref_chan(1))); % Note unit U will be included here. 
unit_in_pref_chan = EStats(Expi).evol.unit_in_pref_chan; % this numbering corresponds to A,B,C... U is not included. 
% chid = find((EStats(Expi).units.unit_num_arr == unit_in_pref_chan(1)) & (EStats(Expi).units.spikeID == EStats(Expi).evol.pref_chan(1))); 
ui = unit_in_pref_chan(1); %find(prefchan_id==chid);
Window = 50:200;
% Find the image name and score for the best image in each block.
% Compute this before hand for faster video generation 
imgColl = repmat("", EStats(Expi).evol.block_n-1,thread_n);
scoreColl = zeros(EStats(Expi).evol.block_n-1,thread_n);
meanscoreColl = zeros(EStats(Expi).evol.block_n-1,thread_n);
for thr_i = 1:thread_n
for blocki = 1:EStats(Expi).evol.block_n-1
    gen_scores = squeeze(mean(EStats(Expi).evol.psth{thr_i,blocki}(ui,Window,:),[1,2]));
    [maxScore, maxIdx] = max(gen_scores); % Idx in local generation
    imgidx = EStats(Expi).evol.idx_seq{thr_i,blocki}(maxIdx);
    imgfullfn = ls(fullfile(EStats(Expi).meta.stimuli, EStats(Expi).imageName(imgidx)+"*"));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(EStats(Expi).meta.stimuli, EStats(Expi).imageName(imgidx)+"*"))
    imgColl(blocki,thr_i) = fullfile(EStats(Expi).meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % FIXED @ Jan.29
    scoreColl(blocki,thr_i) = maxScore;
    meanscoreColl(blocki,thr_i) = mean(gen_scores);
end
end
%% Generation averaged psth and sem
evol_stim_fr = cellfun(@(psth)mean(psth,3),EStats(Expi).evol.psth,'UniformOutput',false);
evol_stim_fr = cell2mat(reshape(evol_stim_fr',1,1,[],thread_n));
evol_stim_sem = cellfun(@(psth)std(psth,0,3)/sqrt(size(psth,3)),EStats(Expi).evol.psth,'UniformOutput',false);
evol_stim_sem = cell2mat(reshape(evol_stim_sem',1,1,[],thread_n));
% Scalor score for evolution
score_avg = cellfun(@(psth)mean(psth(:,51:200,:),'all') - mean(psth(:,1:50,:),'all'),EStats(Expi).evol.psth);
score_sem = cellfun(@(psth)std(squeeze(mean(psth(:,51:200,:),[1,2])))...
    /sqrt(size(psth,3)),EStats(Expi).evol.psth); % - mean(psth(:,1:50,:),[1,2])
%% Generate Movies
savepath = result_dir; mkdir(savepath) % fullfile(result_dir, compose("%s_Evol_"))
color_seq = EStats(Expi).color_seq;
% v = VideoWriter(fullfile(savepath,compose('%s_BGFC6Evol_Exp%02d_Best_PSTH',Animal,Expi)));
v = VideoWriter(fullfile(savepath,compose('%s_BGFC6Evol_%s_Best_PSTH',Animal,stimparts{end-1})));
v.FrameRate = 4;
open(v);
% Set up figure canvas, get the object for use.
h3=figure(4);set(4,'position',[263         148        1341         735]);clf;
ax1_1 = subtightplot(2,4,1,0.07); % Axes 1 images
% set(ax1_1,"position",[0.07,0.578,0.439,0.3412]);
title(compose("%s",EStats(Expi).evol.optim_names(1)))
scoreYLIM = [0,max(score_avg(:,1:end-1)+score_sem(:,1:end-1),[],'all')]; % Preset and align the YLIM
ax3_1 = subtightplot(2,4,2,0.07); % Axes 2 score trajectory
shadedErrorBar([],score_avg(1,1:end-1),score_sem(1,1:end-1),...
        'lineprops',{'Color','k','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
xlabel("Generations");ylabel("Response fr (Hz)");axis tight;ylim(scoreYLIM)
title(EStats(Expi).evol.optim_names(1)+" Evol Traj");
psthYLIM = [0,max(evol_stim_fr(ui,:,1:end-1,:)+evol_stim_sem(ui,:,1:end-1,:),[],'all')];% Preset and align the YLIM
ax2_1 = subtightplot(2,2,3,0.07,0.07,0.05); % Axes 3 PSTH evolution
bgpsth_1 = plot(squeeze(evol_stim_fr(ui,:,1:end-1,1)),'Color',[0.7,0.7,0.7]);hold on
sEB_1 = shadedErrorBar([],evol_stim_fr(ui, :, 1, 1),evol_stim_sem(ui, :, 1, 1),...
        'lineprops',{'Color','r','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%
ylabel("PSTH (Hz)");xlabel("time(ms)");ylim(psthYLIM)
title("Evoked PSTH")
if thread_n==2
ax1_2 = subtightplot(2,4,3,0.07); % Axes 1 images
% set(ax1_2,"position",[0.07,0.578,0.439,0.3412]);
title(compose("%s",EStats(Expi).evol.optim_names(2)))
ax3_2 = subtightplot(2,4,4,0.07); % Axes 2 score trajectory
shadedErrorBar([],score_avg(2,1:end-1),score_sem(2,1:end-1),...
        'lineprops',{'Color','k','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
xlabel("Generations");ylabel("Response fr (Hz)");axis tight;ylim(scoreYLIM)
title(EStats(Expi).evol.optim_names(2)+" Evol Traj");
ax2_2 = subtightplot(2,2,4,0.07,0.07,0.05); % Axes 3 PSTH evolution
bgpsth_2 = plot(squeeze(evol_stim_fr(ui,:,1:end-1,2)),'Color',[0.7,0.7,0.7]);hold on
sEB_2 = shadedErrorBar([],evol_stim_fr(ui, :, 1, 2),evol_stim_sem(ui, :, 1, 2),...
        'lineprops',{'Color','r','LineWidth',1.5},'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
ylabel("PSTH (Hz)");xlabel("time(ms)");ylim(psthYLIM)
title("Evoked PSTH")
end
stimparts = split(EStats(Expi).meta.stimuli,"\");
expday = datetime(EStats(Expi).meta.expControlFN(1:6),'InputFormat','yyMMdd');
% fdrnm = compose("%s-Chan%02d", stimparts{end-1}, pref_chan(1));
ST = suptitle(compose("%s BigGAN FC6 Evol Exp %02d pref chan %s", ...
    stimparts{end-1}, Expi, EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id(1))));
% Dynamically change the components on figure;
for blocki = 1:EStats(Expi).evol.block_n-1
    set(h3,"CurrentAxes",ax1_1)  % change the image shown
    imshow(imgColl{blocki,1}); 
    title(compose("Gen%d Best Rate %.1f", blocki, scoreColl(blocki,1))) ; 
    if thread_n==2,set(h3,"CurrentAxes",ax1_2)
    imshow(imgColl{blocki,2}); 
    title(compose("Gen%d BestRate %.1f", blocki, scoreColl(blocki,2))) ; 
    end
    % Change the title of PSTH
    ax2_1.Title.String = compose("Evoked PSTH Gen%d Evoked Rate %.1f", blocki, meanscoreColl(blocki,1));
    if thread_n==2
    ax2_2.Title.String = compose("Evoked PSTH Gen%d Evoked Rate %.1f", blocki, meanscoreColl(blocki,2));
    end
    % Change the data of current PSTH shaded Errorbar.
    psthcur = evol_stim_fr(ui, :, blocki, 1);
    psthsemcur = evol_stim_sem(ui, :, blocki, 1);
    uE = psthcur + psthsemcur; lE = psthcur - psthsemcur;
    sEB_1.mainLine.YData = psthcur;
    sEB_1.edge(1).YData = lE;
    sEB_1.edge(2).YData = uE;
    sEB_1.patch.Vertices(:,2) = [lE,fliplr(uE)];
    if thread_n==2
    psthcur = evol_stim_fr(ui, :, blocki, 2);
    psthsemcur = evol_stim_sem(ui, :, blocki, 2);
    uE = psthcur + psthsemcur; lE = psthcur - psthsemcur;
    sEB_2.mainLine.YData = psthcur;
    sEB_2.edge(1).YData = lE;
    sEB_2.edge(2).YData = uE;
    sEB_2.patch.Vertices(:,2) = [lE,fliplr(uE)];
    end
    drawnow;
%     pause(0.1)
    Fs = getframe(h3);
    writeVideo(v,Fs);
end
close(v);
end
end