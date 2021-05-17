%% Evol_Cosine_summary
saveroot = "O:\Evol_Cosine";
mkdir(saveroot)
refcoldir = "N:\Stimuli\2020-CosineEvol\RefCollection";
targmap = get_refimg_map(refcoldir);
%% Annotation for the cosineImage
annot_tab = readtable('CosineImage_annotate.csv','Delimiter',',');
target_annot_map = containers.Map(annot_tab.imgName, annot_tab.Annotation);

%% Summary the evolution trajectory
subfdrs = string(ls(saveroot+"\*Alfa*"));
% figure;montage(fullfile(saveroot,deblank(subfdrs(6:end)),"online_scoretraj.png"),...
%         'Size', [5,8], 'BorderSize', 4,'ThumbnailSize',[656   875])
mtg = imtile(fullfile(saveroot,deblank(subfdrs(6:end)),"online_scoretraj.png"),'GridSize', [4,9],'BorderSize', 4);
imwrite(mtg,fullfile(saveroot,"CosineSummary","CosineTraj_summary.png"))
%% Collect statistics of these fitting and tile into table.
%% Assert one thread.
MAXUNUM = 4;
for iTr = 50%:numel(meta_new)
	meta = meta_new{iTr};
	rasters = rasters_new{iTr};
	Trials = Trials_new{iTr};
	imageName = string(Trials.imageName);
	stimpath = meta.stimuli;
    targ_matnm = Trials.TrialRecord.User.target_cfg{1}{1};
    targ_imgnm = Trials.TrialRecord.User.target_cfg{1}{2};
    score_mode = Trials.TrialRecord.User.score_mode{1};
    GANstr = Trials.TrialRecord.User.space_cfg{1}{1};
    imgpos = Trials.TrialRecord.User.evoConfiguration{1,2};
    imgsize = Trials.TrialRecord.User.evoConfiguration{1,3};
    optimstr = Trials.TrialRecord.User.evoConfiguration{1,5};
    [thread_msks,thread_num,row_gen,row_nat,gen_idx_seq,nat_idx_seq,block_arr,block_num] =...
    		sort_trial_into_masks(Trials); 
    assert(thread_num==1)
    targD = load(fullfile(stimpath,targ_matnm)); 
    %% Different scores
    scores_rec = Trials.TrialRecord.User.scores_record; % recorded score from online Exp 
    natscores_rec = Trials.TrialRecord.User.natscores_record; % recorded score for natural images
    % Characterize Successfulness of evol 
    [expS, S] = evol_score_stats({scores_rec, natscores_rec}, "ol", gen_idx_seq,nat_idx_seq);
    %% Get activations 
    chan_arr = [1:64]';
	unit_arr = 0:MAXUNUM;
	chanMat = repmat(chan_arr,1,MAXUNUM+1); % const array
	unitMat = repmat(unit_arr,64,1); % const array
    
    act_mat = squeeze(mean(rasters(:,51:200,:),2));
	bsl_mat = squeeze(mean(rasters(:,1:50,:),[2,3])); % baseline matrix for each channel
	act_tsr = nan(64,MAXUNUM+1,size(rasters,3));
	bslMat = nan(64,MAXUNUM+1);
	for iCh = 1:size(rasters,1)
	chani = meta.spikeID(iCh);
	uniti = meta.unitID(iCh);
	act_tsr(chani,uniti+1,:) = act_mat(iCh,:);
	bslMat(chani,uniti+1) = bsl_mat(iCh);
    end
	for iThr = 1:thread_num % fetch target information.
	meanActMat = Trials.TrialRecord.User.meanActMat{iThr}; % normalize info
	stdActMat = Trials.TrialRecord.User.stdActMat{iThr}; % normalize info
	targetActMat = Trials.TrialRecord.User.targetActMat{iThr}; % target info
	FMat = Trials.TrialRecord.User.maskMat{iThr}; % Mask in the ReprTsr file, per F significance
	end
	objMask = parse_mode2mask(mode); % Parse out the chanXunit mask used to evaluate objective. (area specific)
	maskMat = objMask & FMat & ~isnan(targetActMat); % Channel Selection & F significant

    %% Other offline scores, noise ceiling of the scores. 
    [expS, S] = evol_score_stats(L1_vec, "L1", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(MSE_vec, "MSE", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(corr_vec, "corr", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(L1All_vec, "L1_All", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(MSE_All_vec, "MSE_All", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(corr_All_vec, "corr_All", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(MSE_IT_vec, "MSE_IT", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(corr_IT_vec, "corr_IT", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(MSE_V1V4_vec, "MSE_V1V4", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(corr_V1V4_vec, "corr_V1V4", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(MSE_V4_vec, "MSE_V4", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(corr_V4_vec, "corr_V4", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(MSE_V1_vec, "MSE_V1", gen_idx_seq,nat_idx_seq,expS,S);
    [expS, S] = evol_score_stats(corr_V1_vec, "corr_V1", gen_idx_seq,nat_idx_seq,expS,S);
    %%
    expS.meta = meta;
    expS.imageName = imageName;
    % S = struct();
    % pfx = "ol";
    % if ~strcmp(pfx,"ol"), S.(pfx+"_vec") = score_vec; end
    % S.(pfx+"_gen_avg") = cellfun(@mean,scores_rec);
    % S.(pfx+"_gen_sem") = cellfun(@sem, scores_rec);
    % S.(pfx+"_nat_avg") = cellfun(@mean,natscores_rec);
    % S.(pfx+"_nat_sem") = cellfun(@sem, natscores_rec);
    % Fgen = anova_cells(scores_rec);
    % S.(pfx+"_gen_F") = Fgen.F;
    % S.(pfx+"_gen_F_P") = Fgen.F_P;
    % S.(pfx+"_gen_F_df") = Fgen.STATS.df;
    % Fnat = anova_cells(natscores_rec);
    % S.(pfx+"_nat_F") = Fnat.F;
    % S.(pfx+"_nat_F_P") = Fnat.F_P;
    % S.(pfx+"_nat_F_df") = Fnat.STATS.df;
    % % [cval,pval] = corr([1:block_num]',OLscore_avg,'type','spearman');
    % % S.(pfx+"_gen_corr") = cval;
    % % S.(pfx+"_gen_corr_P") = pval;
    % [~,P,~,TST] = ttest2(cell2mat(scores_rec(end-1:end)), cell2mat(scores_rec(1:2)));% maybe I should use 2:3
    % S.(pfx+"_gen_last_T") = TST.tstat;
    % S.(pfx+"_gen_last_T_P") = P;
    % [maxscore,max_idx] = max(movmean(OLscore_avg,3));
    % if max_idx == numel(scores_rec), max_idx = max_idx - 1; end
    % [~,P,~,TST] = ttest2(cell2mat(scores_rec(max_idx:max_idx+1)), cell2mat(scores_rec(1:2)));
    % S.(pfx+"_gen_best_T") = TST.tstat;
    % S.(pfx+"_gen_best_T_P") = P;
    % S.(pfx+"_bestscore") = maxscore;

    %% Other offline scores, noise ceiling of the scores. 
    
    %%
    
    %% Compare image similarity. 
    
    % Image Similarity under some mask. 
    
    % 
end
%%

parse_mode2mask("V1ITV4")
function [D, S] = evol_score_stats(scorevec,pfx,gen_idx_seq,nat_idx_seq,D,S)
    if nargin == 1, prefix=""; end
    if nargin <= 4, S=struct(); D=struct(); end
    if ~strcmp(pfx,"ol"), % if not online score, then sort
        D.(pfx+"_vec") = scorevec; 
        scores_rec = cellfun(@(idx)scorevec(idx),gen_idx_seq,'uni',0);
        natscores_rec = cellfun(@(idx)scorevec(idx),nat_idx_seq,'uni',0);
    else, % if it's online score, then no need to sort/ 
        scores_rec = scorevec{1};
        natscores_rec = scorevec{2};
    end
    score_avg = cellfun(@mean,scores_rec);
    D.(pfx+"_gen_avg") = cellfun(@mean,scores_rec);
    D.(pfx+"_gen_sem") = cellfun(@sem, scores_rec);
    D.(pfx+"_nat_avg") = cellfun(@mean,natscores_rec);
    D.(pfx+"_nat_sem") = cellfun(@sem, natscores_rec);
    Fgen = anova_cells(scores_rec);
    S.(pfx+"_gen_F") = Fgen.F;
    S.(pfx+"_gen_F_P") = Fgen.F_P;
    S.(pfx+"_gen_F_df") = Fgen.STATS.df;
    Fnat = anova_cells(natscores_rec);
    S.(pfx+"_nat_F") = Fnat.F;
    S.(pfx+"_nat_F_P") = Fnat.F_P;
    S.(pfx+"_nat_F_df") = Fnat.STATS.df;
    % [cval,pval] = corr([1:block_num]',OLscore_avg,'type','spearman');
    % S.(pfx+"_gen_corr") = cval;
    % S.(pfx+"_gen_corr_P") = pval;
    [~,P,~,TST] = ttest2(cell2mat(scores_rec(end-1:end)), cell2mat(scores_rec(1:2)));% maybe I should use 2:3
    S.(pfx+"_gen_last_T") = TST.tstat;
    S.(pfx+"_gen_last_T_P") = P;
    [maxscore,max_idx] = max(movmean(score_avg,3));
    if max_idx == numel(scores_rec), max_idx = max_idx - 1; end
    [~,P,~,TST] = ttest2(cell2mat(scores_rec(max_idx:max_idx+1)), cell2mat(scores_rec(1:2)));
    S.(pfx+"_gen_best_T") = TST.tstat;
    S.(pfx+"_gen_best_T_P") = P;
    S.(pfx+"_bestscore") = maxscore;
    D.(pfx+"_stat") = S;
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