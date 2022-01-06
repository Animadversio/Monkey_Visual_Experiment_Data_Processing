% Compare BigGAN and FC6 prototpyes
% Updated July 13
%%
mat_dir = "O:\Mat_Statistics"; 
saveroot = "O:\Evol_BigGAN_FC6_cmp"; 
figdir = fullfile(saveroot,"summary");
Animal = "Both"; Set_Path;
if Animal == "Both"
A = load(fullfile(mat_dir, "Alfa" + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats');
B = load(fullfile(mat_dir, "Beto" + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats');
BFEStats = [A.BFEStats;B.BFEStats];
clear A B
else
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
end
%%
score_traces = cell(numel(BFEStats),2);
block_traces = cell(numel(BFEStats),2);
wdw = [51:200]; bslwdw = [1:50];
for iTr = 1:numel(BFEStats)
ui = 1;%BFEStats(iTr).evol.unit_in_pref_chan(1);
activ_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:),[2])),BFEStats(iTr).evol.psth,"uni",0); % cell arr, vector of score per gen. 
score_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:),[2]))-mean(psth(ui,bslwdw,:),[2,3]), BFEStats(iTr).evol.psth,"uni",0);
block_col = cellfun(@(idx) BFEStats(iTr).evol.block_arr(idx), BFEStats(iTr).evol.idx_seq,"uni",0); % cell arr, vector of block id per gen. 
assert(contains(BFEStats(iTr).evol.space_names{1},"fc6"))
assert(contains(BFEStats(iTr).evol.space_names{2},"BigGAN"))
for GANi = 1:2
score_traces{iTr,GANi} = cat(1,score_col{GANi,1:end-1});
% zscore_traces{iTr,GANi} = zscores_tsr(pref_chan_id, row_gen & thread_msks{threadi});
block_traces{iTr,GANi} = cat(1,block_col{GANi,1:end-1});
end
end
%% Newer API verions
[score_m_traj_col, block_traj_col, score_m_traj_extrap_col, block_traj_extrap_col, ...
   extrap_mask_col, sucsmsk_end, sucsmsk_max, tval_end_arr, pval_end_arr, tval_max_arr, pval_max_arr] =...
   optim_traj_process(block_traces, score_traces, ["FC6", "BigGAN"], "max12", 35);
%% 
% sucsmsk_end
% sucsmsk_max
prefchan_arr = arrayfun(@(S)S.evol.pref_chan(1),BFEStats');
area_arr = arrayfun(@area_map,prefchan_arr);
anim_arr = arrayfun(@(S)S.Animal,BFEStats');
V4msk = area_arr=="V4";
ITmsk = area_arr=="IT";
Amsk = anim_arr=="Alfa";
Bmsk = anim_arr=="Beto";
expids = find(sucsmsk_end);
%% Show the paired images interactively
figh=figure('pos',[897   527   693   407]);
for iTr = expids(:)'
    stimpath = BFEStats(iTr).meta.stimuli;
    imageName = BFEStats(iTr).imageName;
    ephysFN = BFEStats(iTr).meta.ephysFN;
    imgpos = BFEStats(iTr).evol.imgpos(1,:);
    imgsize = BFEStats(iTr).evol.imgsize(1);
    pref_chan = BFEStats(iTr).evol.pref_chan(1);
    expstr = compose("%s PrefCh %02d  [%.1f %.1f] %.1f deg",ephysFN,pref_chan,imgpos(1),imgpos(2),imgsize(1));
    ui = BFEStats(iTr).evol.unit_in_pref_chan(1);
    activ_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:),[2])),BFEStats(iTr).evol.psth,"uni",0); % cell arr, vector of score per gen. 
    score_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:),[2]))-mean(psth(ui,bslwdw,:),[2,3]), BFEStats(iTr).evol.psth,"uni",0);
    block_col = cellfun(@(idx) BFEStats(iTr).evol.block_arr(idx), BFEStats(iTr).evol.idx_seq,"uni",0); % cell arr, vector of block id per gen. 
    [mxScore_FC, mxidx_FC] = max(score_col{1,end-1});
    [mxScore_BG, mxidx_BG] = max(score_col{2,end-1});
    mScore_FC = mean(score_col{1,end-1}); sScore_FC = sem(score_col{1,end-1});
    mScore_BG = mean(score_col{2,end-1}); sScore_BG = sem(score_col{2,end-1});
    bestimgidx_FC = BFEStats(iTr).evol.idx_seq{1,end-1}(mxidx_FC);
    bestimgidx_BG = BFEStats(iTr).evol.idx_seq{2,end-1}(mxidx_BG);
    bestimg_FC = imread(fullfile(stimpath,imageName{bestimgidx_FC}+".bmp"));
    bestimg_BG = imread(fullfile(stimpath,imageName{bestimgidx_BG}+".bmp"));
    imshow(imtile({bestimg_FC,bestimg_BG}));
    title(compose("%s   Gen%d FC6 %.1f+-%.1f BigGAN %.1f+-%.1f ",expstr,size(score_col,2)-1,mScore_FC,sScore_FC,mScore_BG,sScore_BG))
    pause
end
%%

saveallform([outdir],"ProtoCmp_Exp_"+ephysFN,figh)