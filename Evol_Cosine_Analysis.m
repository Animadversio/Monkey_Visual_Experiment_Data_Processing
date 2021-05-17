%%  This is all in one analysis for Cosine Evolution Experiment
%%
!ExpRecordBackup.bat 
%%
Animal = "Alfa";Set_Path;
ftrrows = find(contains(ExpRecord.expControlFN,["generate_BigGAN_cosine"]));%&contains(ExpRecord.expControlFN,["210218","210408","210413","210415"]));%...
    %&~contains(ExpRecord.ephysFN,["Alfa-02032021","Alfa-04032021"]));%&..."Alfa-18022021"
%     contains(ExpRecord.ephysFN,["Alfa-12022021"]));%,"Alfa-09022021""Alfa-27102020-003", "Alfa-27102020-004"
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows, Animal, false);
%%
saveroot = "O:\Evol_Cosine";
mkdir(saveroot)
refcoldir = "N:\Stimuli\2020-CosineEvol\RefCollection";
targmap = get_refimg_map(refcoldir);
%%
global figdir explabel
MAXUNUM = 4;
P.plotCriter = true;
P.plotPopuEvol = true;
P.plotImageSeq = true;
P.plotMovieSum = true;
figh0 = figure(5);
figh = figure(1);
movh = figure(2);
mtgh = figure(3);
h = figure(4);
for iTr = 1:numel(meta_new)%numel(meta_new)-6:
meta = meta_new{iTr};
rasters = rasters_new{iTr};
Trials = Trials_new{iTr};
imageName = string(Trials.imageName);
stdActMat = Trials.TrialRecord.User.stdActMat;
meanActMat = Trials.TrialRecord.User.meanActMat;
targetActMat = Trials.TrialRecord.User.targetActMat;
maskMat = Trials.TrialRecord.User.maskMat;
target_cfg = Trials.TrialRecord.User.target_cfg;
score_mode = Trials.TrialRecord.User.score_mode;
scores_rec = Trials.TrialRecord.User.scores_record; % recorded score from online Exp 
natscores_rec = Trials.TrialRecord.User.natscores_record; % recorded score for natural images
genN_rec = arrayfun(@(b)b*ones(size(scores_rec{b})),[1:numel(scores_rec)]','uni',0); % format gen number as the score
scores_vec = cell2mat(scores_rec);
genN_vec = cell2mat(genN_rec);
% Sort trial idx into different masks for usage below
[thread_msks,thread_num,row_gen,row_nat,gen_idx_seq,nat_idx_seq,block_arr,block_num] =...
    sort_trial_into_masks(Trials); 
%% Create dir for stimuli and title labels
stimparts = split(meta.stimuli,"\");
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
fdrnm = compose("%s-%s", stimparts{end-1}, score_mode{1});
% fdrnm = compose("%s-%s-Chan%02d-1",datestr(expday,'yyyy-mm-dd'), Animal, pref_chan(1));
figdir = fullfile(saveroot, fdrnm);
if exist(figdir,'dir'),warning("%s figure directory exist! Beware",figdir);end
mkdir(figdir)
explabel = stimparts{end-1};
for iThr = 1:thread_num
spacenm = Trials.TrialRecord.User.space_cfg{iThr}{1};
optimnm = Trials.TrialRecord.User.evoConfiguration{iThr,5};
imgpos = Trials.TrialRecord.User.evoConfiguration{iThr,2};
imgsize = Trials.TrialRecord.User.evoConfiguration{iThr,3};
explabel = explabel + compose(" %s (trg: %s)\n%s (%s)",score_mode{iThr},target_cfg{iThr}{2},spacenm,optimnm);
explabel = explabel + compose("   Pos [%.1f,%.1f] Size %.1f deg",imgpos,imgsize);
end
%% Online computed score
figh0 = figure;hold on 
score_avg = cellfun(@mean,scores_rec);
score_sem = cellfun(@sem,scores_rec);
errorbar(1:numel(score_avg),score_avg,score_sem);
natscore_avg = cellfun(@mean,natscores_rec);
natscore_sem = cellfun(@sem,natscores_rec);
errorbar(1:numel(score_avg),natscore_avg,natscore_sem);
ylabel("Online Score");xlabel("Generation")
legend(["Gen","Nat"])
title(compose("Scores Online\n")+explabel,'interpreter','none')
saveas(figh0,fullfile(figdir,"online_scoretraj.png"))
saveas(figh0,fullfile(figdir,"online_scoretraj.pdf"))
% save(fullfile(figdir,"EvolStat.mat"),'EvolStat')
% fprintf("ExpStats saved to %s EvolStat.mat!\n",figdir)
%% Form the activation
act_mat = squeeze(mean(rasters(:,51:200,:),2));
% bsl_mat = squeeze(mean(rasters(:,1:50,:),[2,3])); % baseline matrix for each channel
bsl_mat = squeeze(mean(rasters(:,1:50,:),[2,3])); % baseline matrix for each channel
act_tsr = nan(64,MAXUNUM+1,size(rasters,3));
bslMat = nan(64,MAXUNUM+1);
for iCh = 1:size(rasters,1)
chani = meta.spikeID(iCh);
uniti = meta.unitID(iCh);
act_tsr(chani,uniti+1,:) = act_mat(iCh,:);
bslMat(chani,uniti+1) = bsl_mat(iCh);
end
%%
chan_arr = [1:64]';
unit_arr = 0:MAXUNUM;
chanMat = repmat(chan_arr,1,MAXUNUM+1);
unitMat = repmat(unit_arr,64,1);
for iThr = 1:thread_num
stdActMat = Trials.TrialRecord.User.stdActMat{iThr};
meanActMat = Trials.TrialRecord.User.meanActMat{iThr};
targetActMat = Trials.TrialRecord.User.targetActMat{iThr};
FMat = Trials.TrialRecord.User.maskMat{iThr}; % Mask in the ReprTsr file, F significance
mode = score_mode{iThr};
targetName = target_cfg{iThr}{2};
targimg = imread(targmap(targetName)); % get the target image
objMask = parse_mode2mask(mode); % Parse out the chanXunit mask used to evaluate objective. (area specific)
maskMat = objMask & FMat & ~isnan(targetActMat); % Channel Selection & F significant
% score_recalc = squeeze(nansum(((act_tsr - meanActMat)./stdActMat).*....*targetActMat.*maskMat,[1,2]));
%                       ((targetActMat- meanActMat)./stdActMat).*maskMat,[1,2]));
% Note to baseline subtract! or correlation will be bad.
vec_reprs = reshape(act_tsr(repmat(maskMat,1,1,size(act_tsr,3))),[],size(act_tsr,3)) - bslMat(maskMat); 
vec_targs = targetActMat(maskMat);
vec_norm_mean = meanActMat(maskMat);
vec_norm_std = stdActMat(maskMat);
vec_reprs_norm = (vec_reprs-vec_norm_mean)./vec_norm_std; % Z scored repr 
vec_targs_norm = (vec_targs-vec_norm_mean)./vec_norm_std; % Z scored target
vecchan = chanMat(maskMat);
vecunit = unitMat(maskMat);
%% Visualize the evolution of population pattern. 
if P.plotPopuEvol
figh = popuEvol(vec_reprs,vec_targs,vecchan,vecunit,gen_idx_seq,nat_idx_seq,block_arr,figh);
ylabel("Raw Activation")
saveas(figh,fullfile(figdir,"RawAct_Pattern_Evolution.png"))
saveas(figh,fullfile(figdir,"RawAct_Pattern_Evolution.pdf"))
figh = popuEvol(vec_reprs_norm,vec_targs_norm,vecchan,vecunit,gen_idx_seq,nat_idx_seq,block_arr,figh);
ylabel("Z Activation")
saveas(figh,fullfile(figdir,"NormAct_Pattern_Evolution.png"))
saveas(figh,fullfile(figdir,"NormAct_Pattern_Evolution.pdf"))
end
%% Q: What's the right criterion?
%% Different Matching Criterions to compare with target.
MSE = -nanmean((vec_targs_norm - vec_reprs_norm).^2,1);
innprod = (vec_targs_norm)'*(vec_reprs_norm);
corr_norm2_arr = corr(vec_reprs_norm, vec_targs_norm);
repr_norm_norm = norm_axis(vec_reprs_norm, 1);
if P.plotCriter
% Correlations with various normalization scheme
% corr_arr = corr(vec_reprs,vec_targs);
% h=plotStatsTraj(corr_arr,"Raw Repr (no normalized) Correlation","raw_corr",gen_idx_seq,nat_idx_seq,block_arr);
% corr_norm1_arr = corr((vec_reprs-vec_norm_mean)./vec_norm_std, vec_targs);
% h=plotStatsTraj(corr_norm1_arr,"Z Repr (Evol normalized) Correlation","Z_1norm_corr",gen_idx_seq,nat_idx_seq,block_arr);
corr_norm2_arr = corr(vec_reprs_norm, vec_targs_norm);%(vec_targs-vec_norm_mean)./vec_norm_std);
h=plotStatsTraj(corr_norm2_arr,"Correlation of (Both normalized) Z Repr","Z_2norm_corr",gen_idx_seq,nat_idx_seq,block_arr,h);
repr_norm_norm = norm_axis(vec_reprs_norm, 1);
h=plotStatsTraj(repr_norm_norm,"Norm Evol Z Repr","EvolZ_norm",gen_idx_seq,nat_idx_seq,block_arr,h);
% Inner Product
innprod = (vec_targs_norm)'*(vec_reprs_norm);
h=plotStatsTraj(innprod,"Inner Prod of (Both normalized) Z Repr","Z_Innprod",gen_idx_seq,nat_idx_seq,block_arr,h);
innprod = (vec_targs)'*((vec_reprs-vec_norm_mean)./vec_norm_std);
h=plotStatsTraj(innprod,"Inner Prod of (Evol normalized) Repr","Z_1norm_Innprod",gen_idx_seq,nat_idx_seq,block_arr,h);
innprod = (vec_targs)'*vec_reprs;
h=plotStatsTraj(innprod,"Inner Prod of (no normalized) Repr","raw_Innprod",gen_idx_seq,nat_idx_seq,block_arr,h);
% MSE and L1 distance
MSE = -mean((vec_targs_norm - vec_reprs_norm).^2,1);
h=plotStatsTraj(MSE,"MSE of (Both normalized) Repr","-Z MSE",gen_idx_seq,nat_idx_seq,block_arr,h);
L1 = -mean(abs(vec_targs_norm - vec_reprs_norm),1);
h=plotStatsTraj(L1,"L1 of (Both normalized) Repr","-Z L1",gen_idx_seq,nat_idx_seq,block_arr,h);
end
%%
h=plotStatsTraj(innprod,"Inner Prod of (Evol normalized) Repr","Z_1norm_Innprod",gen_idx_seq,nat_idx_seq,block_arr,h);

%% Visualize Image evolution
if P.plotImageSeq
% ExempN = 3;
% imgscore_vec = MSE; scorestr="MSE";
refnmMap = get_refimg_map(meta.stimuli, true); % mapping from ref image name to the full path. true to find in parent folder
mtgh = imgevol_movie(MSE, "MSE", gen_idx_seq, nat_idx_seq, imageName, meta.stimuli, refnmMap, Animal, 3,mtgh);
mtgh = imgevol_movie(corr_norm2_arr', "corr_norm", gen_idx_seq, nat_idx_seq, imageName, meta.stimuli, refnmMap, Animal, 3, mtgh);
mtgh = imgevol_movie(innprod, "dot_norm", gen_idx_seq, nat_idx_seq, imageName, meta.stimuli, refnmMap, Animal, 3, mtgh);
end
%
if P.plotMovieSum
movh = summary_movie(MSE, "MSE", vec_reprs_norm, vec_targs_norm, vecchan,...
    gen_idx_seq, nat_idx_seq, imageName, meta.stimuli, targimg, Animal, 4, movh);
movh = summary_movie(corr_norm2_arr, "corr", vec_reprs_norm, vec_targs_norm, vecchan,...
    gen_idx_seq, nat_idx_seq, imageName, meta.stimuli, targimg, Animal, 4, movh);
end
end
end

%%
[sortchan,sortidx] = sort(vecchan);
getChanX(sortchan)

%% 
function movh = summary_movie(imgscore_vec, scorestr, reprMat, targVec, chanArr, ...
    gen_idx_seq, nat_idx_seq, imageName, stimpath, targimg, Animal, ExempN, h)
% Generate overall movie summary for each experiment: 
%   Showing the population vector, target image, evolved image and scoring. 
%   
global figdir explabel
if nargin <=11, ExempN=4; end % best 3 vs worst 3.
if nargin <=12, movh=figure; else, movh = figure(h); end
% ExempN = 4;
% movh = figure; 
v = VideoWriter(fullfile(figdir,compose('%s_PopEvol_summary_%s',Animal,scorestr)));
v.FrameRate = 2; open(v);
set(movh,'pos',[781  171  1248  660]);clf(movh)
T=tiledlayout(movh,2,3,'Padding','compact','TileSp','compact');
ax0=nexttile(T,1);ax1=nexttile(T,2);ax2=nexttile(T,3);ax3=nexttile(T,4,[1,3]);
title(T,explabel,'interp','none')
block_num = size(gen_idx_seq,2);
% Panel for Target img
set(movh,'CurrentAxes',ax0)
imshow(targimg)
title(ax0,compose("Target"))
Cord = colororder;
clrseq = brewermap(block_num-1,'Spectral');%from red to blue/purple-ish
% row_gen = cell2mat(reshape(gen_idx_seq(1:end-1),[],1));
% row_nat = cell2mat(reshape(nat_idx_seq(1:end-1),[],1));
% Panel for overall score
mean_score_gen = cellfun(@(idx)mean(imgscore_vec(idx)),gen_idx_seq(1:end-1));
sem_score_gen = cellfun(@(idx)sem(imgscore_vec(idx)),gen_idx_seq(1:end-1));
mean_score_nat = cellfun(@(idx)mean(imgscore_vec(idx)),nat_idx_seq(1:end-1));
sem_score_nat = cellfun(@(idx)sem(imgscore_vec(idx)),nat_idx_seq(1:end-1));
set(movh,'CurrentAxes',ax2)
% scatter(block_arr(row_gen),imgscore_vec(row_gen),16,'MarkerEdgeAlpha',0.6);hold on
% scatter(block_arr(row_nat),imgscore_vec(row_nat),16,'MarkerEdgeAlpha',0.6);
errorbar([],mean_score_gen,sem_score_gen,'color',[Cord(1,:),0.6],'LineWidth',1.5);hold on
errorbar([],mean_score_nat,sem_score_nat,'color',[Cord(2,:),0.6],'LineWidth',1.5)
YLIM=ylim();L = plot([0,0],YLIM,'-.r','LineWidth',1,'HandleVis','off');
xlabel("Generation");ylabel(scorestr);
legend('Gen','Nat','Location','best','AutoUpdate','off');box('off')
set(movh,'CurrentAxes',ax3)
% Population activity pattern
[sortchan, sortidx] = sort(chanArr);
sortchanX = getChanX(sortchan); % X posistion to plot each chanel
reprBlkMat = cell2mat(cellfun(@(idx)mean(reprMat(:,idx),2),gen_idx_seq,'uni',0));
reprBlkMat_sem = cell2mat(cellfun(@(idx)sem(reprMat(:,idx),2),gen_idx_seq,'uni',0));
reprBlkMat_nat = cell2mat(cellfun(@(idx)mean(reprMat(:,idx),2),nat_idx_seq,'uni',0));
TGTL = plot(ax3,sortchanX, targVec(sortidx),'color',[0,0,0,0.6],'LineWidth',2.5); hold on
YLIM3 = ylim(ax3);
plot([32.5,32.5],YLIM3,'-.r','LineWidth',0.5,'HandleVis','off') % IT V1 sep
plot([48.5,48.5],YLIM3,'-.r','LineWidth',0.5,'HandleVis','off') % V1 V4 sep 
xlim([0,65])
xlabel("Channels");ylabel("Activation");
title(compose("PopRepr Evol (Dot line: natural, Solid line: generated, Black: Target Repr)"),'interpreter','none')
for blocki = 1:block_num-1
L.XData = [blocki,blocki];
legend(ax2,'Gen','Nat','Location','southeast');
% evol pattern
plot(ax3,sortchanX,reprBlkMat(sortidx,blocki),'color',[clrseq(blocki,:),0.5],'LineWidth',1.5) % sortidx
plot(ax3,sortchanX,reprBlkMat_nat(sortidx,blocki),':','color',[clrseq(blocki,:),0.3],'LineWidth',1.5) % sortidx
% uistack(TGTL,'top'); % move target to the top layer?
% Curate best few images
[maxScore, maxIdx] = sort(imgscore_vec(gen_idx_seq{blocki}),'descend');
score_mean = mean_score_gen(blocki);%mean(imgscore_vec(gen_idx_seq{blocki}));
score_sem = sem_score_gen(blocki);%sem(imgscore_vec(gen_idx_seq{blocki}));
title(ax2,compose('Block %d %s scores mean %.2f sem %.2f',blocki,scorestr,score_mean,score_sem),'interp','none')    
Exemplar_score = [maxScore(1:ExempN)];%, maxScore(end-ExempN+1:end)
Exemplar = imageName(gen_idx_seq{blocki}([maxIdx(1:ExempN)]));
for i = 1:ExempN    
    imshow(imread(fullfile(stimpath,Exemplar(i)+".bmp")),'Parent',ax1)
    title(ax1,compose('%s scores Cur %.2f',scorestr,Exemplar_score(i)),'interp','none')
    pause(0.1)
    Fs = getframe(movh);
    writeVideo(v,Fs);
end
end
close(v);
end

function [sortchanX] = getChanX(sortchan)
sortchanX = sortchan;
for i = 1:numel(sortchan)
    totN = sum(sortchan==sortchan(i));
    preN = sum(sortchan(1:i)==sortchan(i));
    sortchanX(i) = sortchanX(i) + (preN - 1) / totN;
end
end

function refnmMap = get_refimg_map(stim_dir,parent)
if nargin==1, parent=false; end
if parent,
parent_dir = fileparts(stim_dir);
else
parent_dir = stim_dir;
end
refimgs = dir(parent_dir);
refimgs = refimgs(~arrayfun(@(R)R.isdir,refimgs));
refnmMap = containers.Map();
for i = 1:numel(refimgs)
    [~,fn,ext] = fileparts(refimgs(i).name);
    refnmMap(string(fn)) = string(fullfile(refimgs(i).folder,refimgs(i).name));
end
end

function mtgh = imgevol_movie(imgscore_vec, scorestr, gen_idx_seq, nat_idx_seq, imageName, stimpath, refnmMap, Animal, ExempN, h)
% movie of image popuplation (best and worst `ExempN` examples showed.)
if nargin <=8, ExempN=3; end % best 3 vs worst 3.
if nargin <=9, mtgh=figure; else, mtgh = h; clf(mtgh);end
global figdir explabel
v = VideoWriter(fullfile(figdir,compose('%s_Pop_ImageEvol_%s',Animal,scorestr)));
v.FrameRate = 2; open(v);
T=tiledlayout(mtgh,2,1,'Padding','compact');set(mtgh,'pos',[731  189  1600  660])
block_num = size(gen_idx_seq,2);
for blocki = 1:block_num-1
    [maxScore, maxIdx] = sort(imgscore_vec(gen_idx_seq{blocki}),'descend');
    [maxScore_nat, maxIdx_nat] = sort(imgscore_vec(nat_idx_seq{blocki}),'descend');
    score_mean = mean(imgscore_vec(gen_idx_seq{blocki}));
    score_sem = sem(imgscore_vec(gen_idx_seq{blocki}));
    score_nat_mean = mean(imgscore_vec(nat_idx_seq{blocki}));
    score_nat_sem = sem(imgscore_vec(nat_idx_seq{blocki}));
    Exemplar_score = [maxScore(1:ExempN), maxScore(end-ExempN+1:end)];
    Exemplar = imageName(gen_idx_seq{blocki}([maxIdx(1:ExempN),maxIdx(end-ExempN+1:end)]));
    Exemplar_nat_score = [maxScore_nat(1:ExempN), maxScore_nat(end-ExempN+1:end)];
    Exemplar_nat = imageName(nat_idx_seq{blocki}([maxIdx_nat(1:ExempN),maxIdx_nat(end-ExempN+1:end)]));
    Exemplar_nat_nmprt = split(Exemplar_nat,"_thread00"); % map the ref image name (w\o suffix) to full path
    Exemplar_nat_nmprt = Exemplar_nat_nmprt(:,1);
    reffullpath = arrayfun(@(F)refnmMap(F),Exemplar_nat_nmprt); 
    ax1=nexttile(T,1);
    montage(fullfile(stimpath,Exemplar+".bmp"),'Size',[1,ExempN*2],'ThumbnailSize',[256,256],'BorderSize',2)%,'Parent',ax1
    title(ax1,compose('Generated, %s scores %s mean %.2f sem %.2f',scorestr,num2str(Exemplar_score,"%.2f, "),score_mean,score_sem),'interp','none')
    ax2=nexttile(T,2);
    montage(reffullpath,'Size',[1,ExempN*2],'ThumbnailSize',[256,256],'BorderSize',2)%,'Parent',ax1
    title(ax2,compose('Natural %s scores %s mean %.2f sem %.2f',scorestr,num2str(Exemplar_nat_score,"%.2f, "),score_nat_mean,score_nat_sem),'interp','none')
    title(T,explabel+compose("Block %d",blocki),'interp','none')
    Fs = getframe(mtgh);
    writeVideo(v,Fs);
end
close(v);
end
function figh = popuEvol(reprMat,targVec,chanArr,unitId,gen_idx_seq,nat_idx_seq,block_arr,figh)
% This creates figure of population response pattern evolving.
global figdir explabel
if nargin<8, figh=[]; end
if isempty(figh), figh = figure;else, clf(figh); end
reprBlkMat = cell2mat(cellfun(@(idx)mean(reprMat(:,idx),2),gen_idx_seq,'uni',0));
reprBlkMat_sem = cell2mat(cellfun(@(idx)sem(reprMat(:,idx),2),gen_idx_seq,'uni',0));
reprBlkMat_nat = cell2mat(cellfun(@(idx)mean(reprMat(:,idx),2),nat_idx_seq,'uni',0));
blocklist = min(block_arr):max(block_arr);
blockN = max(block_arr);
clrseq = brewermap(blockN-1,'Spectral');%from red to blue/purple-ish
[sortedchan, sortidx] = sort(chanArr);%sortidx to make sure the channel follows from 1:64
figure(figh);hold on;set(figh,'pos',[472   499   895   445])
for i =1:blockN-1
    plot(reprBlkMat(sortidx,i),'color',[clrseq(i,:),0.5],'LineWidth',1.5) % sortidx
    plot(reprBlkMat_nat(sortidx,i),':','color',[clrseq(i,:),0.3],'LineWidth',1.5) % sortidx
    %shadedErrorBar(-249:500,psthavg_col{threadi,i},psthsem_col{threadi,i},'lineProps',{'color',[clrseq(i,:),0.7],'LineWidth',1.5})
end
plot(targVec(sortidx),'color',[0,0,0,0.75],'LineWidth',2.5);
ITV4sep = sum(chanArr<33)+0.5;
V1V4sep = sum(chanArr<49)+0.5;
vline([ITV4sep,V1V4sep],'-.r')
xlim([0,size(reprBlkMat,1)+1])
xlabel("Channels");ylabel("Activation");
title(compose("PopRepr Evol\n%s",explabel),'interpreter','none')
end

function chanmsk = parse_mode2mask(mode, chan_arr, unit_arr)
% Parse the score_mode string into a channel mask of which channels are
% used in the exp.
% chanmsk = parse_mode2mask("corr_V4IT");
MAXUNUM = 4;
if nargin ==1, chan_arr=[1:64]'; unit_arr = 0:MAXUNUM; end
chanMat = repmat(chan_arr,1,numel(unit_arr));
unitMat = repmat(unit_arr,numel(chan_arr),1);
if ~contains(mode,["V1","V4","IT"])
    chanmsk = ones(size(chanMat),'logical');
else
    chanmsk = zeros(size(chanMat),'logical');
    if contains(mode,"IT")
        chanmsk = chanmsk | ((chan_arr<=32));
    end
    if contains(mode,"V4")
        chanmsk = chanmsk | ((chan_arr>=49));
    end
    if contains(mode,"V1")
        chanmsk = chanmsk | ((chan_arr<49) & (chan_arr>32));
    end
end
end

function [figh, mean_corr_gen, sem_corr_gen, mean_corr_nat, sem_corr_nat] = ...
	plotStatsTraj(Stat_vec,statname,savestr,gen_idx_seq,nat_idx_seq,block_arr,figh)
% Plot trajectory of certain stats through out the evolution exp.
% Stat_vec: A vector, same length as trial num. This func will sort the
%   stat into each block and generated and natural. Then the scores will be
%   displayed as errorbar plot. The function will automatically clip out the last 
%   generation since it's usually incomplete. 
% block_arr: The block numbering, same length as Stat_vec. 
% 
% Rely on global variable `figdir` and `explabel`. 
if nargin<7, figh = figure; else, clf(figh); end
global figdir explabel
mean_corr_gen = cellfun(@(idx)mean(Stat_vec(idx)),gen_idx_seq(1:end-1));
sem_corr_gen = cellfun(@(idx)sem(Stat_vec(idx)),gen_idx_seq(1:end-1));
mean_corr_nat = cellfun(@(idx)mean(Stat_vec(idx)),nat_idx_seq(1:end-1));
sem_corr_nat = cellfun(@(idx)sem(Stat_vec(idx)),nat_idx_seq(1:end-1));
Cord = colororder;
% set(0,'CurrentFigure',figh)
figure(figh);set(figh,'pos',[1000  558  560  420])
row_gen = cell2mat(reshape(gen_idx_seq(1:end-1),[],1));
row_nat = cell2mat(reshape(nat_idx_seq(1:end-1),[],1));
scatter(block_arr(row_gen),Stat_vec(row_gen),25,'MarkerEdgeAlpha',0.6);hold on
scatter(block_arr(row_nat),Stat_vec(row_nat),25,'MarkerEdgeAlpha',0.6);
errorbar([],mean_corr_gen,sem_corr_gen,'color',[Cord(1,:),0.6],'LineWidth',1.5)
errorbar([],mean_corr_nat,sem_corr_nat,'color',[Cord(2,:),0.6],'LineWidth',1.5)
ylabel(statname);
xlabel("Generation");
legend(["Gen","Nat"])
title(compose("%s Evol Trajectory\n%s",statname,explabel),'interpreter','none')
saveas(figh,fullfile(figdir,savestr+"_scoretraj.png"))
saveas(figh,fullfile(figdir,savestr+"_scoretraj.pdf"))
end
% %% Correlation of vector repr
% corr_arr = corr(vec_reprs,vec_targs);
% mean_corr_gen = cellfun(@(idx)mean(corr_arr(idx)),gen_idx_seq);
% sem_corr_gen = cellfun(@(idx)sem(corr_arr(idx)),gen_idx_seq);
% mean_corr_nat = cellfun(@(idx)mean(corr_arr(idx)),nat_idx_seq);
% sem_corr_nat = cellfun(@(idx)sem(corr_arr(idx)),nat_idx_seq);
% Cord = colororder;
% figure;
% scatter(block_arr(row_gen),corr_arr(row_gen),25,'MarkerEdgeAlpha',0.6);hold on
% scatter(block_arr(row_nat),corr_arr(row_nat),25,'MarkerEdgeAlpha',0.6);
% errorbar([],mean_corr_gen,sem_corr_gen,'color',[Cord(1,:),0.6],'LineWidth',1.5)
% errorbar([],mean_corr_nat,sem_corr_nat,'color',[Cord(2,:),0.6],'LineWidth',1.5)
% ylabel("Z Repr Correlation");
% xlabel("Generation");
% legend(["Gen","Nat"])
% title("Correlation of single trial repr and target repr")
%%
function [thread_msks,thread_num,row_gen,row_nat,gen_idx_seq,nat_idx_seq,block_arr,block_num] =...
    sort_trial_into_masks(Trials)
% preprocessing step: sort the image names / trial into masks.
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
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
end
%%
% formattedMatPath = "S:\Data-Ephys-MAT\Alfa-30032021-012_formatted";
% load(formattedMatPath);
% lfps = lfps(1:end-1,:,:);
% rasters = rasters(1:end-1,:,:);
% meta.spikeID = meta.spikeID(1:end-1);
% % meta.stimuli = 
% save(formattedMatPath,'Trials','lfps','meta','rasters')