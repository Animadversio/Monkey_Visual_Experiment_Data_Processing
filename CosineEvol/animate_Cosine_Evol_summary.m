function animate_Cosine_Evol_summary(CStats, movh)
% imgscore_vec, scorestr, reprMat, targVec, chanArr, ...
%    gen_idx_seq, nat_idx_seq, imageName, stimpath, targimg, Animal, ExempN
% Generate overall movie summary for each experiment: 
%   Showing the population vector, target image, evolved image and scoring. 
%   Adapt from `summary_movie` but better interface 
if nargin <=11, ExempN=4; end % best 3 vs worst 3.
if nargin <=1, movh=figure; else, movh = figure(movh); end
for iTr = 1:numel(CStats)
CStat = CStats(iTr);
Animal = CStat.Animal;
figdir = CStat.meta.figdir;
stimpath = CStat.meta.stimuli;
explabel = CStat.meta.explabel;
scorestr = CStat.targ.score_mode{1}; % assume only one thread. 
gen_idx_seq = CStat.stim.gen_idx_seq;
nat_idx_seq = CStat.stim.nat_idx_seq;
block_num = size(gen_idx_seq,2);
% load the target image
targnm = CStat.targ.target_cfg{1}{2};
repr_dir = fileparts(CStat.targ.repr_path); % parent dir of repr path. 
targfullpath = map2fullpath({targnm},repr_dir);
targfullpath = targfullpath{1};
targimg = imread(targfullpath);
% Find the scores
imgscore_vec = CStat.score.offline_vec;

v = VideoWriter(fullfile(figdir,compose('%s_PopEvol_summary_%s',Animal,scorestr)));
v.FrameRate = 2; open(v);
set(movh,'pos',[187  171  1045  660]);clf(movh)
T=tiledlayout(movh,2,3,'Padding','compact','TileSp','compact');
ax0=nexttile(T,1);ax1=nexttile(T,2);ax2=nexttile(T,3);ax3=nexttile(T,4,[1,3]);
title(T,explabel,'interp','none')

% Panel for Target img
set(movh,'CurrentAxes',ax0)
imshow(targimg)
title(ax0,compose("Target"))
Cord = colororder;
clrseq = brewermap(block_num-1,'Spectral'); % from red to blue / purple-ish
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
YLIM=ylim();L = plot([0,0],YLIM,'-.r','LineWidth',1,'HandleVis','off'); % progress bar 
xlabel("Generation");ylabel(scorestr,'interp','none');ylim(YLIM)
legend('Gen','Nat','Location','best','AutoUpdate','off');box('off')
set(movh,'CurrentAxes',ax3)

% Sort out Population activity pattern
mask = CStat.targ.maskVec{1}; % Fmask, get rid of the non-selective units. 
targVec = CStat.targ.targetActVec{1};
reprMat = CStat.resp.evoke_trials - CStat.resp.bslmean;
norm_reprMat = (CStat.resp.evoke_trials - CStat.resp.bslmean - CStat.targ.meanActVec{1}) ./ CStat.targ.stdActVec{1};
norm_targVec = (CStat.targ.targetActVec{1} - CStat.targ.meanActVec{1})./ CStat.targ.stdActVec{1};
% [V1msk, V4msk, ITmsk] = get_areamask(S(iTr).units.spikeID, array_layout)
[targmsk, targ_area] = parse_mode2maskVec(CStat.targ.score_mode, ...
                        CStat.meta.array_layout, CStat.units.spikeID);

targVec = norm_targVec;
reprMat = norm_reprMat;

% [sortchan, sortidx] = sort(chanArr);
% sortchanX = getChanX(sortchan); % X posistion to plot each chanel

chanX = 1:sum(mask); % contiguous X position 
unit_name = CStat.units.unit_name_arr(mask);
unit_num  = CStat.units.unit_num_arr(mask);
chan_num  = CStat.units.spikeID(mask);
unit_area = area_map(chan_num, CStat.meta.array_layout);
targmsk_2plot = targmsk(mask);
[chanX_targ_broke, targVec_targ_broke] = prepare_broken_target(chanX(targmsk_2plot), targVec(mask&targmsk));
reprBlkMat = cell2mat(cellfun(@(idx)mean(reprMat(mask,idx),2),gen_idx_seq,'uni',0));
reprBlkMat_sem = cell2mat(cellfun(@(idx)sem(reprMat(mask,idx),2),gen_idx_seq,'uni',0));
reprBlkMat_nat = cell2mat(cellfun(@(idx)mean(reprMat(mask,idx),2),nat_idx_seq,'uni',0));
TGTL = plot(ax3, chanX, targVec(mask),'color',[0,0,0,0.6],'LineWidth',1.5); hold on
plot(ax3, chanX_targ_broke, targVec_targ_broke,'color',[0,0,0,0.6],'LineWidth',3.0); hold on
xticks(chanX); xticklabels(unit_name); 
YLIM3 = ylim(ax3); xlim([0,max(chanX)+1])
if CStat.meta.array_layout == "Beto_new"
V4V1sep = 0.5 + sum(unit_area=="V4");
V1ITsep = 0.5 + sum(unit_area=="V4") + sum(unit_area=="V1");
else
V4V1sep = 0.5 + sum(unit_area=="IT");
V1ITsep = 0.5 + sum(unit_area=="IT") + sum(unit_area=="V1");
end
plot([V4V1sep, V4V1sep],YLIM3,'-.r','LineWidth',0.5,'HandleVis','off') % IT V1 sep
plot([V1ITsep, V1ITsep],YLIM3,'-.r','LineWidth',0.5,'HandleVis','off') % V1 V4 sep 

xlabel("Channels");ylabel("Activation");
title(compose("PopRepr Evol (Dot line: natural, Solid line: generated, Black: Target Repr)"),'interpreter','none')
% Start iteration through blocks
for blocki = 1:block_num-1
L.XData = [blocki,blocki];
legend(ax2,'Gen','Nat','Location','southeast');
% evol pattern
plot(ax3,chanX,reprBlkMat(:,blocki),...
     'color',[clrseq(blocki,:),0.5],'LineWidth',1.5) % sortidx
plot(ax3,chanX,reprBlkMat_nat(:,blocki),':',...
     'color',[clrseq(blocki,:),0.3],'LineWidth',1.5) % sortidx
% uistack(TGTL,'top'); % move target to the top layer?
% Curate best few images
[maxScore, maxIdx] = sort(imgscore_vec(gen_idx_seq{blocki}),'descend');
score_mean = mean_score_gen(blocki); % mean(imgscore_vec(gen_idx_seq{blocki}));
score_sem = sem_score_gen(blocki); % sem(imgscore_vec(gen_idx_seq{blocki}));
title(ax2,compose('Block %d scores mean %.2f sem %.2f',blocki,score_mean,score_sem),'interp','none')    
Exemplar_score = [maxScore(1:ExempN)];%, maxScore(end-ExempN+1:end)
Exemplar = CStat.imageName(gen_idx_seq{blocki}([maxIdx(1:ExempN)]));
for i = 1:ExempN    
    imshow(imread(fullfile(stimpath,Exemplar(i)+".bmp")),'Parent',ax1)
    title(ax1,compose('%s scores Cur %.2f',scorestr,Exemplar_score(i)),'interp','none')
    pause(0.1)
    Fs = getframe(movh);
    writeVideo(v,Fs);
end %  examplar loop 
end %  blocks loop
close(v);
end
end

function [chanX_broke, vec_broke] = prepare_broken_target(chanX_targ, vec_targ)
% Example 
%   [chanX_broke, vec_broke] = prepare_broken_target([1,2,5,6,8], [10,5,56,4,5]);
% 
%   assert(isequaln(chanX_broke, [1     2   nan     5     6   nan     8]));
%   assert(isequaln(vec_broke, [10     5   nan    56     4   nan     5]));
delta = diff(chanX_targ);
if all(delta == 1)
chanX_broke = chanX_targ;
vec_broke = vec_targ;
else
breakidx = find(delta ~= 1);
chanX_broke = [];
vec_broke = [];
csr = 1;
for i = 1:numel(breakidx)
    idx = breakidx(i);
    chanX_broke = [chanX_broke, reshape(chanX_targ(csr:idx),1,[]), nan];
    vec_broke = [vec_broke, reshape(vec_targ(csr:idx),1,[]), nan];
    csr = breakidx(i) + 1;
end
chanX_broke = [chanX_broke, reshape(chanX_targ(csr:end),1,[])];
vec_broke = [vec_broke, reshape(vec_targ(csr:end),1,[])];
end
end