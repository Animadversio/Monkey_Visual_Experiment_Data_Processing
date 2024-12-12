%% This script demonstrate the efficiency of using functionalized script to do analysis.
%  It uses 
%     Evol_BigGAN_FC6_Animation_fun.m
%     Evol_BigGAN_FC6_Collect_Stats_fun.m
%% StyleGAN 
Animal="Alfa";Set_Path;
ftr = find(contains(ExpRecord.Exp_collection, "StyleGAN_Evol") );
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr,Animal); % 4 till now 2024
%%
SGEStats = Evol_BigGAN_FC6_Collect_Stats_fun(meta_new(1:end), rasters_new(1:end), Trials_new(1:end));
%%
save(fullfile("E:\OneDrive - Harvard University\Evol_StyleGAN_cmp\summary","EvolStatsAll"),'SGEStats')
%%
videopaths = Evol_BigGAN_FC6_Animation_fun(SGEStats);

%%
sumdir = "E:\OneDrive - Harvard University\Evol_StyleGAN_cmp\summary";
Cord = [[0,0,1];
        [1,0,0]];
failures = struct('iter', {}, 'str', {});
P = struct();
P.plot_img = true;
P.plot_traj = true;
for Expi = 1%:numel(SGEStats) % failure_Expis %5,
% try
BFES = SGEStats(Expi);
stimpath = BFES.meta.stimuli;
prefchan = BFES.evol.pref_chan(1);
area  = area_map(prefchan);
% expdir = fullfile(rootdir,compose("%s_Exp%02d_Ch%02d",Animal,Expi,prefchan));
% mkdir(expdir)
expdir = BFES.meta.figdir;
mkdir(expdir)
expstr = compose("%s Exp%02d PrefChan %02d (%s)\nThr1: %s x %s Thr2: %s x %s\n%s",Animal,Expi,prefchan,area,...
                 BFES.evol.space_cfg{1}{1},BFES.evol.optim_names(1),...
                 BFES.evol.space_cfg{2}{1},BFES.evol.optim_names(2), BFES.meta.ephysFN);
expstr = strrep(expstr,"_","-");
fprintf(expstr+"\n")
%% Get the responses
bsl_col = cellfun(@(P)squeeze(mean(P(:,1:45,:),[1,2])), BFES.evol.psth, 'uni',0);
act_col = cellfun(@(P)squeeze(mean(P(:,51:200,:),[1,2])), BFES.evol.psth, 'uni',0);
refbsl_col = cellfun(@(P)squeeze(mean(P(:,1:45,:),[1,2])), BFES.ref.psth, 'uni',0);
refact_col = cellfun(@(P)squeeze(mean(P(:,51:200,:),[1,2])), BFES.ref.psth, 'uni',0);
bsl_vec = cellfun(@mean, bsl_col);
sembsl_vec = cellfun(@sem, bsl_col);

meanact_vec = cellfun(@mean, act_col);
semact_vec = cellfun(@sem, act_col);
stdact_vec = cellfun(@std, act_col);

refact_vec = cellfun(@mean, refact_col);
refsem_vec = cellfun(@sem, refact_col);
refstd_vec = cellfun(@std, refact_col);

nGen = size(meanact_vec,2);
if P.plot_traj
%% Evol Trajectory mean
figure(1);clf;hold off;set(1,'pos',[680   400   460   480])
for iThr = 1:2
    % plot(meanact_vec(iThr,1:end-1),'linestyle','-','color',Cord(iThr,:),'Linewidth',1.5);
    % shadedErrorBar([],refact_vec(iThr,1:end-1),refsem_vec(iThr,1:end-1),...
    %     'lineProps',{'Color',Cord(iThr+2,:),'Linewidth',1.,'linestyle','-.'},'patchSaturation',0.7)
    shadedErrorBar([],meanact_vec(iThr,1:end-1),stdact_vec(iThr,1:end-1),...
        'lineProps',{'Color',Cord(iThr,:),'Linewidth',2},'patchSaturation',0.3)
    hold on 
    plot(refact_vec(iThr,1:end-1),'Color',Cord(iThr,:),'Linewidth',1,'linestyle','-.')
    hold on 
    plot(bsl_vec(iThr,1:end-1),'linestyle',':','color',Cord(iThr,:),'Linewidth',1.)
end
xlabel("Generations")
ylabel("Firing Rate (event/sec)")
title(compose("%s\nEvolution Trajectory",expstr))
% legend(["FC6","FC6 baseline","BigGAN","BigGAN baseline"],'location','best')
legend(["FC6","FC6 ref","FC6 baseline","BigGAN","BigGAN ref","BigGAN baseline"],'location','best')
xlim([0,nGen+1])
%%
figure(3);clf;hold off;set(3,'pos',[680   400   460   480])
for iThr = 1:2
    % plot(meanact_vec(iThr,1:end-1),'linestyle','-','color',Cord(iThr,:),'Linewidth',1.5);
    % shadedErrorBar([],refact_vec(iThr,1:end-1),refsem_vec(iThr,1:end-1),...
    %     'lineProps',{'Color',Cord(iThr+2,:),'Linewidth',1.,'linestyle','-.'},'patchSaturation',0.7)
    shadedErrorBar([],meanact_vec(iThr,1:end-1),semact_vec(iThr,1:end-1),...
        'lineProps',{'Color',Cord(iThr,:),'Linewidth',2},'patchSaturation',0.3)
    hold on 
    shadedErrorBar([],refact_vec(iThr,1:end-1),refsem_vec(iThr,1:end-1),...
        'lineProps',{'Color',Cord(iThr,:),'Linewidth',1,'linestyle','-.'},'patchSaturation',0.1)
%     plot(refact_vec(iThr,1:end-1),'Color',Cord(iThr+2,:),'Linewidth',1,'linestyle','-')
    hold on 
%     plot(bsl_vec(iThr,1:end-1),'linestyle',':','color',Cord(iThr,:),'Linewidth',1.)
    shadedErrorBar([],bsl_vec(iThr,1:end-1),sembsl_vec(iThr,1:end-1),...
        'lineProps',{'Color',Cord(iThr,:),'Linewidth',1,'linestyle',':'},'patchSaturation',0.1)
end
xlabel("Generations")
ylabel("Firing Rate (event/sec)")
title(compose("%s\nEvolution Trajectory",expstr))
% legend(["FC6","FC6 baseline","BigGAN","BigGAN baseline"],'location','best')
legend(["FC6","FC6 ref image","FC6 baseline","BigGAN","BigGAN ref image","BigGAN baseline"],'location','southeast')
xlim([0,nGen+1])
%% Evol Trajectory std
% Cord = colororder;
% figure(3);clf;hold off;set(3,'pos',[680   558   460   420])
% for iThr = 1:2
% plot(stdact_vec(iThr,1:end-1),'linestyle','-','color',Cord(iThr,:),'Linewidth',1.5);
% hold on
% end
% xlabel("Generations")
% ylabel("Firing Rate (event/sec)")
% title(compose("%s\nEvolution Trajectory",expstr))
% legend(["FC6","BigGAN"],'location','best')
%% Single trial trajectory
[genvec1, actvec1] = actcell2traj(act_col(1,1:end-1));
[genvec2, actvec2] = actcell2traj(act_col(2,1:end-1));
figure(2);clf;hold off;set(2,'pos',[1140   400    460   480])
for iThr = 1:2
    [genvec, actvec] = actcell2traj(act_col(iThr,1:end-1));
    scatter(genvec,actvec,'markerfacecolor',Cord(iThr,:),...
        'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
    hold on
    plot(bsl_vec(iThr,1:end-1),'color',Cord(iThr,:),...
        'linestyle',':','Linewidth',2)
end
legend(["FC6","FC6 baseline","BigGAN","BigGAN baseline"],'location','best')
xlabel("Generations")
ylabel("Firing Rate (event/sec)")
title(compose("%s\nEvolution Trajectory (single trial)",expstr))
xlim([0,nGen+1])

saveallform([expdir], compose("%s_Exp%02d_EvolTraj",Animal,Expi), 1)
saveallform([sumdir], compose("%s_Exp%02d_EvolTraj",Animal,Expi), 1, ["png"])
saveallform([expdir], compose("%s_Exp%02d_EvolTraj_sgtr",Animal,Expi), 2)
saveallform([sumdir], compose("%s_Exp%02d_EvolTraj_sgtr",Animal,Expi), 2, ["png"])
saveallform([expdir], compose("%s_Exp%02d_EvolTraj_sem_shaded",Animal,Expi), 3)
saveallform([sumdir], compose("%s_Exp%02d_EvolTraj_sem_shaded",Animal,Expi), 3, ["png"])
end


%% Preferred images across Evolution
if P.plot_img
% BFES.imageName(BFES.evol.idx_seq{2,3})
bestFCimg_col = {};
bestBGimg_col = {};
meanFCimg_col = {};
meanBGimg_col = {};
WtmeanFCimg_col = {};
WtmeanBGimg_col = {};
for igen = 1:size(BFES.evol.idx_seq,2)-1
%     [FCim_gen, FCimnames] = loadEvolImages(stimpath, 0, igen);
%     [BGim_gen, BGimnames] = loadEvolImages(stimpath, 1, igen);
    tic
    FCimnms = BFES.imageName(BFES.evol.idx_seq{1,igen});
    BGimnms = BFES.imageName(BFES.evol.idx_seq{2,igen});
    FCimfps = fullfile(stimpath, strcat(FCimnms,".bmp"));
    BGimfps = fullfile(stimpath, strcat(BGimnms,".bmp"));
    FCim_gen = cellfun(@(fp) imread(fp), FCimfps, "uni", 0);
    BGim_gen = cellfun(@(fp) imread(fp), BGimfps, "uni", 0);
    toc
    %%
    FCact = act_col{1,igen};
    BGact = act_col{2,igen};
    % meanFCimg = mean(cat(4, FCim_gen{:}) / 255.0, 4); % Buggy! this will qunatize pixel to 0, 1
    % meanBGimg = mean(cat(4, BGim_gen{:}) / 255.0, 4); % Buggy! this will qunatize pixel to 0, 1
    meanFCimg = mean(double(cat(4, FCim_gen{:})) / 255.0, 4);
    meanBGimg = mean(double(cat(4, BGim_gen{:})) / 255.0, 4);
    WtmeanFCimg = sum(double(cat(4, FCim_gen{:})) / 255.0 ...
                       .* reshape(FCact,1,1,1,[]), 4) / sum(FCact);
    WtmeanBGimg = sum(double(cat(4, BGim_gen{:})) / 255.0 ...
                       .* reshape(BGact,1,1,1,[]), 4) / sum(BGact);
    % uniformly averaged images
    meanFCimg_col{end+1} = meanFCimg;
    meanBGimg_col{end+1} = meanBGimg;
    % activation averaged images
    WtmeanFCimg_col{end+1} = WtmeanFCimg;
    WtmeanBGimg_col{end+1} = WtmeanBGimg;
    % get the max activation image, very noisy
    [~,maxidx] = max(FCact);
    bestFCimg_col{end+1} = FCim_gen{maxidx};
    [~,maxidx] = max(BGact);
    bestBGimg_col{end+1} = BGim_gen{maxidx};
    % figure(3);
    % subplot(121);montage(FCim_gen)
    % subplot(122);montage(BGim_gen)
    % figure(4);
    % subplot(121);imshow(meanFCimg)
    % subplot(122);imshow(meanBGimg)
    % pause
end
%%
% from the first to last - 1 generation, the mean activation, 2% and 98% to avoid the max min unstable
CLIM = prctile(meanact_vec(:,1:end-1),[2,98],'all')'; 
frame_FCimg_col = score_frame_image_arr(WtmeanFCimg_col, meanact_vec(1,1:end-1), CLIM);%, cmap, LineWidth)
frame_BGimg_col = score_frame_image_arr(WtmeanBGimg_col, meanact_vec(2,1:end-1), CLIM);%, cmap, LineWidth)
frame_FCbestimg_col = score_frame_image_arr(bestFCimg_col, meanact_vec(1,1:end-1), CLIM);%, cmap, LineWidth)
frame_BGbestimg_col = score_frame_image_arr(bestBGimg_col, meanact_vec(2,1:end-1), CLIM);%, cmap, LineWidth)
%%
FCbestmtg = imtile(bestFCimg_col,8);
BGbestmtg = imtile(bestBGimg_col,8);
frame_FCbestmtg = imtile(frame_FCbestimg_col,8);
frame_BGbestmtg = imtile(frame_BGbestimg_col,8);
FCmtg = imtile(WtmeanFCimg_col,8);
BGmtg = imtile(WtmeanBGimg_col,8);
frame_FCmtg = imtile(frame_FCimg_col,8);
frame_BGmtg = imtile(frame_BGimg_col,8);
figure(5);set(5,'pos',[ 81         100       1720         750])
T = tiledlayout(1,2,'tilesp','compact','pad','compact');
nexttile(1);imshow(FCmtg) % montage(meanFCimg_col)
nexttile(2);imshow(BGmtg) % montage(meanBGimg_col)
title(T,compose("Image Evolution Color lim [%.2f,%.2f]\n%s",CLIM(1),CLIM(2),expstr))
figure(6);set(6,'pos',[ 81         100       1720         750])
T = tiledlayout(1,2,'tilesp','compact','pad','compact');
nexttile(T,1);imshow(frame_FCmtg)
nexttile(T,2);imshow(frame_BGmtg)
title(T,compose("Image Evolution Color lim [%.2f,%.2f]\n%s",CLIM(1),CLIM(2),expstr))
figure(7);set(7,'pos',[ 81         100       1720         750])
T = tiledlayout(1,2,'tilesp','compact','pad','compact');
nexttile(T,1);imshow(frame_FCbestmtg)
nexttile(T,2);imshow(frame_BGbestmtg)
title(T,compose("Image Evolution Color lim [%.2f,%.2f]\n%s",CLIM(1),CLIM(2),expstr))
%%
saveallform([sumdir,expdir], compose("%s_Exp%02d_FC6BGImageEvol_cmp",Animal,Expi), 5, ["png","pdf"])
saveallform([sumdir,expdir], compose("%s_Exp%02d_FC6BGImageEvol_framed_cmp",Animal,Expi), 6, ["png","pdf"])
saveallform([sumdir,expdir], compose("%s_Exp%02d_FC6BGImageEvol_bestframed_cmp",Animal,Expi), 7, ["png","pdf"])

imwrite(FCmtg, fullfile(expdir,compose("%s_Exp%02d_FC6ImageEvol.png",Animal,Expi)))
imwrite(BGmtg, fullfile(expdir,compose("%s_Exp%02d_BGImageEvol.png",Animal,Expi)))
imwrite(frame_FCmtg, fullfile(expdir,compose("%s_Exp%02d_FC6ImageEvol_framed.png",Animal,Expi)))
imwrite(frame_BGmtg, fullfile(expdir,compose("%s_Exp%02d_BGImageEvol_framed.png",Animal,Expi)))
imwrite(FCbestmtg, fullfile(expdir,compose("%s_Exp%02d_FC6ImageEvol_best.png",Animal,Expi)))
imwrite(BGbestmtg, fullfile(expdir,compose("%s_Exp%02d_BGImageEvol_best.png",Animal,Expi)))
imwrite(frame_FCbestmtg, fullfile(expdir,compose("%s_Exp%02d_FC6ImageEvol_best_framed.png",Animal,Expi)))
imwrite(frame_BGbestmtg, fullfile(expdir,compose("%s_Exp%02d_BGImageEvol_best_framed.png",Animal,Expi)))
end
%% 
% catch ME
%      failures(end + 1).iter = Expi;
%      failures(end).str  = getReport(ME);
%      disp(ME)
%      save("D:\BigGAN_run_failure.mat",'failures')
%      fileID = fopen('D:\BigGAN_run_failure.txt','a+');
%      nbytes = fprintf(fileID,"Exp %d\n",Expi);
%      nbytes = fprintf(fileID,getReport(ME));
% end
end




%%
result_dir = "E:\OneDrive - Harvard University\StyleGAN_Evol_Movies";
savepath = result_dir; mkdir(savepath) % fullfile(result_dir, compose("%s_Evol_"))
ExpType = "StyleGAN";
EStats = SGEStats;

for Expi = 1:length(EStats)
fprintf("Processing StyleGAN Evolution Exp %d\n",Expi)
thread_n = EStats(Expi).evol.thread_num;
% ui = EStats(Expi).evol.unit_in_pref_chan(1);
% assert(all(ui==EStats(Expi).evol.unit_in_pref_chan))
% the following snippet is to get rid of 0 unit (null channel)
prefchan_id = find((EStats(Expi).units.spikeID == EStats(Expi).evol.pref_chan(1))); % Note unit U will be included here. 
unit_in_pref_chan = EStats(Expi).evol.unit_in_pref_chan; % this numbering corresponds to A,B,C... U is not included. 
% chid = find((EStats(Expi).units.unit_num_arr == unit_in_pref_chan(1)) & (EStats(Expi).units.spikeID == EStats(Expi).evol.pref_chan(1))); 
ui = unit_in_pref_chan(1); %find(prefchan_id==chid);
Window = 50:200;
% Find the image name and score for the best image in each block.
% Compute this before hand for faster video generation 
% Sort a image name list for each gen each thread. 
imgColl = repmat("", EStats(Expi).evol.block_n-1,thread_n);
sortImgColl = cell(EStats(Expi).evol.block_n-1,thread_n);
sortScoreColl = cell(EStats(Expi).evol.block_n-1,thread_n);
maxscoreColl = zeros(EStats(Expi).evol.block_n-1,thread_n);
minscoreColl = zeros(EStats(Expi).evol.block_n-1,thread_n);
meanscoreColl = zeros(EStats(Expi).evol.block_n-1,thread_n);
for thr_i = 1:thread_n
for blocki = 1:EStats(Expi).evol.block_n-1
    gen_scores = squeeze(mean(EStats(Expi).evol.psth{thr_i,blocki}(ui,Window,:),[1,2]));
    [maxScore, maxIdx] = max(gen_scores); % Idx in local generation
    maximgidx = EStats(Expi).evol.idx_seq{thr_i,blocki}(maxIdx); % get the trial index array for this block, get the trial index for block max
    imgfullfn = ls(fullfile(EStats(Expi).meta.stimuli, EStats(Expi).imageName(maximgidx)+"*"));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(EStats(Expi).meta.stimuli, EStats(Expi).imageName(maximgidx)+"*"))
%     imgColl(blocki,thr_i) = fullfile(EStats(Expi).meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % FIXED @ Jan.29
    maxscoreColl(blocki,thr_i) = max(gen_scores);
    minscoreColl(blocki,thr_i) = min(gen_scores);
    meanscoreColl(blocki,thr_i) = mean(gen_scores);
    suffix = imgfullfn(end-3:end);
    [sortScore, sortIdx] = sort(gen_scores,'Descend');
    sortimgidx = EStats(Expi).evol.idx_seq{thr_i,blocki}(sortIdx); % get the trial index array for this block, get the trial index for block max
    sortImgColl{blocki,thr_i} = arrayfun(@(idx)fullfile(EStats(Expi).meta.stimuli, ...
        string(EStats(Expi).imageName{idx})+suffix), sortimgidx); % full image path for sorted image name list for each block each thread
    sortScoreColl{blocki,thr_i} = sortScore;
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
%% Generate Movies with Color coded frame
v = VideoWriter(fullfile(savepath,compose('%s_BGFC6Evol_Exp%02d_imbatch_frame',Animal,Expi)));
v.FrameRate = 2;
open(v);
% Set up figure canvas, get the object for use.
h4=figure(6);set(6,'position',[263          34        1611         849]);clf;
T = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1); ax2 = nexttile(2);
stimparts = split(EStats(Expi).meta.stimuli,"\");
% expday = datetime(EStats.meta.expControlFN(1:6),'InputFormat','yyMMdd');
% fdrnm = compose("%s-Chan%02d", stimparts{end-1}, pref_chan(1));
ST = sgtitle(compose("%s BigGAN FC6 Evol Exp %02d pref chan %s", ...
    stimparts{end-1}, Expi, EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id(1))));
% Plot the 2 batch of images
for blocki = 1:EStats(Expi).evol.block_n-1
    set(h4,"CurrentAxes",ax1)  % change the image shown
    sortScore = sortScoreColl{blocki,1};
    imglist1 = score_frame_image_arr(sortImgColl{blocki,1}, sortScoreColl{blocki,1}, [minscoreColl(blocki,1), maxscoreColl(blocki,1)], parula); 
    montage(imglist1)
    title(compose("Gen%d Rate Max %.1f Mean %.1f Min %.1f", blocki, maxscoreColl(blocki,1), ...
        meanscoreColl(blocki,1), minscoreColl(blocki,1) )) ; 
    if thread_n==2,set(h4,"CurrentAxes",ax2)
    imglist2 = score_frame_image_arr(sortImgColl{blocki,2}, sortScoreColl{blocki,2}, [minscoreColl(blocki,2), maxscoreColl(blocki,2)], parula); 
    montage(imglist2)
    title(compose("Gen%d Rate Max %.1f Mean %.1f Min %.1f", blocki, maxscoreColl(blocki,2), ...
        meanscoreColl(blocki,2), minscoreColl(blocki,2) )) ; 
    end
    drawnow;
%     pause(0.1)
    Fs = getframe(h4);
    writeVideo(v,Fs);
end
close(v);
%% Generate Movies
color_seq = EStats(Expi).color_seq;
v = VideoWriter(fullfile(savepath,compose('%s_BGFC6Evol_Exp%02d_imbatch',Animal,Expi)));
v.FrameRate = 2;
open(v);
% Set up figure canvas, get the object for use.
h4=figure(5);set(5,'position',[263          34        1611         849]);clf;
T = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1); ax2 = nexttile(2);
stimparts = split(EStats(Expi).meta.stimuli,"\");
ST = sgtitle(compose("%s BigGAN FC6 Evol Exp %02d pref chan %s", ...
    stimparts{end-1}, Expi, EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id(1))));
% Plot the 2 batch of images
for blocki = 1:EStats(Expi).evol.block_n-1
    set(h4,"CurrentAxes",ax1)  % change the image shown
    montage(sortImgColl{blocki,1}); 
    title(compose("Gen%d Rate Max %.1f Mean %.1f Min %.1f", blocki, maxscoreColl(blocki,1), ...
        meanscoreColl(blocki,1), minscoreColl(blocki,1) )) ; 
    if thread_n==2,set(h4,"CurrentAxes",ax2)
    montage(sortImgColl{blocki,2}); 
    title(compose("Gen%d Rate Max %.1f Mean %.1f Min %.1f", blocki, maxscoreColl(blocki,2), ...
        meanscoreColl(blocki,2), minscoreColl(blocki,2) )) ; 
    end
    drawnow;
%     pause(0.1)
    Fs = getframe(h4);
    writeVideo(v,Fs);
end
close(v);
end





%%




%%
% StyleGAN 
Animal="Alfa";Set_Path;
ftr = find(contains(ExpRecord.ephysFN,"21102020") & contains(ExpRecord.Exp_collection, "StyleGAN_Evol") );
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr,Animal);
%%
Animal="Alfa";Set_Path;
ftr = find(contains(ExpRecord.ephysFN,"23102020") & contains(ExpRecord.Exp_collection, "StyleGAN_Evol") );
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr,Animal);
%%
saveroot = "E:\OneDrive - Washington University in St. Louis\StyleGAN_evol";
Triali = 1;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
%%
pref_chan = Trials.TrialRecord.User.prefChan;
assert(all(pref_chan==pref_chan(1)))
pref_chan=pref_chan(1);

stimparts = split(meta.stimuli,"\");
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
% fdrnm = compose("%s-%s-Chan%02d-1",datestr(expday,'yyyy-mm-dd'), Animal, pref_chan(1));
fdrnm = compose("%s-Chan%02d", stimparts{end-1}, pref_chan(1));
figdir = fullfile(saveroot, fdrnm);
mkdir(figdir)
%% Collect basic info (adapted from Evol_Collect_Stats)
pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
imgpos = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,2));
imgsize = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,3));
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
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

imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm) & ...% first 2 characters are not digits
        cellfun(@(c) ~contains(c(end-4:end), "_nat"), imgnm) ; % last 4 characters are not `_nat`
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the thread 1 and 2 (and maybe thread 3 4)
thread_msks = cell(1, thread_num);
for threadi = 1:thread_num
    msk = contains(imgnm, compose("thread%03d", threadi - 1));
    thread_msks{threadi} = msk; % store masks in a structure for the ease to iterate
end
assert(sum(cellfun(@sum, thread_msks)) == length(imgnm))

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
%% Activation matrix. The simplest stats. 
act_mat = squeeze(mean(rasters(:,51:200,:),2));
%% Score Trajectory Plot. 
gen_clr = {'red','magenta'}; nat_clr = {'green','cyan'};
figure(6);
iCh=10;
for iCh = 1:numel(meta.spikeID)
cla(gca, 'reset')
gen_act_mean = cellfun(@(idx)mean(act_mat(iCh,idx),2), gen_idx_seq(:,1:end-1));
gen_act_sem = cellfun(@(idx)std(act_mat(iCh,idx),1,2)/sqrt(length(idx)), gen_idx_seq(:,1:end-1));
nat_act_mean = cellfun(@(idx)mean(act_mat(iCh,idx),2), nat_idx_seq(:,1:end-1));
nat_act_sem = cellfun(@(idx)std(act_mat(iCh,idx),1,2)/sqrt(length(idx)), nat_idx_seq(:,1:end-1));
% plot(block_list(1:end-1), gen_act_mean')
for thri = 1:thread_num
   shadedErrorBar(block_list(1:end-1), gen_act_mean(thri,:), gen_act_sem(thri,:),'lineProps',...
       {'Color', gen_clr{thri},'LineWidth',2},'patchSaturation',0.15)
end
for thri = 1:thread_num
   shadedErrorBar(block_list(1:end-1), nat_act_mean(thri,:), nat_act_sem(thri,:),'lineProps',...
       {'Color', nat_clr{thri},'LineWidth',2},'patchSaturation',0.15)
end
legend(["thr1 evo","thr2 evo","thr1 nat","thr2 nat"])
xlabel("Generation");ylabel("Evoked Firing Rate");box off
title(compose("%s %s score traj\n Thr1 Face Thr2 ImageNet",Animal,unit_name_arr(iCh)))
saveas(6,fullfile(figdir,compose("%s_%s_score_traj.png",Animal,unit_name_arr(iCh))))
% pause
end


function [genvec, actvec] = actcell2traj(actcol)
actvec = cat(1,actcol{:});
gencol = arrayfun(@(geni) geni * ones(numel(actcol{geni}),1), 1:numel(actcol),'uni',0);
genvec = cat(1,gencol{:});
end

%%
function ImageEvol_pair_cmp(FCimg_col, BGimg_col, actmat, prefix, suffix, expdir, rootdir)
CLIM = prctile(actmat(:,:),[2,98],'all')'; 
frame_FCimg_col = score_frame_image_arr(FCimg_col, actmat(1,:), CLIM);%, cmap, LineWidth)
frame_BGimg_col = score_frame_image_arr(BGimg_col, actmat(2,:), CLIM);%, cmap, LineWidth)
FCmtg = imtile(FCimg_col,8);
BGmtg = imtile(BGimg_col,8);
frame_FCmtg = imtile(frame_FCimg_col,8);
frame_BGmtg = imtile(frame_BGimg_col,8);
imwrite(FCmtg, fullfile(expdir,compose("%s_FC6ImageEvol%s.png",prefix,suffix)))
imwrite(BGmtg, fullfile(expdir,compose("%s_BGImageEvol%s.png",prefix,suffix)))

end