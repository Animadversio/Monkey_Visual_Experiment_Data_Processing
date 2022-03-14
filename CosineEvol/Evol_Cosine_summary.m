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
saveroot = "O:\Evol_Cosine";
matdir = "O:\Mat_Statistics";
array_layout = "Alfa";
CosSummary = [];
CosStats = [];
for iTr = 4:numel(meta_new)
    meta = meta_new{iTr};
    rasters = rasters_new{iTr};
    Trials = Trials_new{iTr};
    if isempty(meta), fprintf("Cosine Exp %02d loading not successful\n",iTr);continue; end
    fprintf("Cosine Exp %02d, %s\n",iTr,meta.ephysFN)
    imageName = string(Trials.imageName);
    stimpath = meta.stimuli;
    targ_matnm = Trials.TrialRecord.User.target_cfg{1}{1};
    targ_imgnm = Trials.TrialRecord.User.target_cfg{1}{2};
    targ_annot = target_annot_map(targ_imgnm);
    score_mode = Trials.TrialRecord.User.score_mode{1};
    GANstr = Trials.TrialRecord.User.space_cfg{1}{1};
    imgpos = Trials.TrialRecord.User.evoConfiguration{1,2};
    imgsize = Trials.TrialRecord.User.evoConfiguration{1,3};
    optimstr = Trials.TrialRecord.User.evoConfiguration{1,5};
    [thread_msks,thread_num,row_gen,row_nat,gen_idx_seq,nat_idx_seq,block_arr,block_num] =...
            sort_trial_into_masks(Trials); 
    assert(thread_num==1) % current process is compatible with single thread population evolution. 
    targD = load(fullfile(stimpath,targ_matnm)); 
    %% Create dir for stimuli and title labels
    stimparts = split(meta.stimuli,"\");
    expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
    fdrnm = compose("%s-%s", stimparts{end-1}, score_mode);
    % fdrnm = compose("%s-%s-Chan%02d-1",datestr(expday,'yyyy-mm-dd'), Animal, pref_chan(1));
    figdir = fullfile(saveroot, fdrnm);
    if exist(figdir,'dir'),warning("%s figure directory exist! Beware",figdir);end
    mkdir(figdir)
    explabel = stimparts{end-1};
    explabel = explabel + compose(" %s targ: %s (%s)\n%s (%s)",score_mode,targ_imgnm,targ_annot,GANstr,optimstr);
    explabel = explabel + compose("   Pos [%.1f,%.1f] Size %.1f deg",imgpos,imgsize);
    %% Different scores
    scores_rec = Trials.TrialRecord.User.scores_record; % recorded score from online Exp 
    if isfield(Trials.TrialRecord.User,"natscores_record")
    natscores_rec = Trials.TrialRecord.User.natscores_record; % recorded score for natural images
    else
    natscores_rec
    end
    %% Get activations 
    chan_arr = [1:64]';
    unit_arr = 0:MAXUNUM;
    chanMat = repmat(chan_arr,1,MAXUNUM+1); % const array
    unitMat = repmat(unit_arr,64,1); % const array
    % Arrange act mat into tensor format. 
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
    evk_tsr = act_tsr - bslMat;
    for iThr = 1:thread_num % fetch target information.
    meanActMat = Trials.TrialRecord.User.meanActMat{iThr}; % normalize info
    stdActMat = Trials.TrialRecord.User.stdActMat{iThr}; % normalize info
    targetActMat = Trials.TrialRecord.User.targetActMat{iThr}; % target info
    FMat = Trials.TrialRecord.User.maskMat{iThr}; % Mask in the ReprTsr file, per F significance
    end
    baseMask = FMat & ~isnan(targetActMat); % basic mask for valid entries
    [objMask, targ_area] = parse_mode2mask(score_mode,array_layout); % Parse out the chanXunit mask used to evaluate objective. (area specific)
    objMask = objMask & FMat & ~isnan(targetActMat); % Same Channel Selection as Evol & F significant

    % Characterize Successfulness of evol 
    [expS, S] = evol_score_stats({scores_rec, natscores_rec}, "ol", gen_idx_seq,nat_idx_seq);
    %% Other offline scores, noise ceiling of the scores. 
    % Same areas as online (use objMask), but use different criterion
    L1_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, objMask, "L1",array_layout);
    [expS, S] = evol_score_stats(L1_vec, "L1", gen_idx_seq,nat_idx_seq,expS,S);
    MSE_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, objMask, "MSE",array_layout);
    [expS, S] = evol_score_stats(MSE_vec, "MSE", gen_idx_seq,nat_idx_seq,expS,S);
    corr_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, objMask, "corr",array_layout);
    [expS, S] = evol_score_stats(corr_vec, "corr", gen_idx_seq,nat_idx_seq,expS,S);
    
    % Different areas, different criterions.
    L1_All_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "L1_All",array_layout);
    [expS, S] = evol_score_stats(L1_All_vec, "L1_All", gen_idx_seq,nat_idx_seq,expS,S);
    MSE_All_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "MSE_All",array_layout);
    [expS, S] = evol_score_stats(MSE_All_vec, "MSE_All", gen_idx_seq,nat_idx_seq,expS,S);
    corr_All_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "corr_All",array_layout);
    [expS, S] = evol_score_stats(corr_All_vec, "corr_All", gen_idx_seq,nat_idx_seq,expS,S);
    MSE_IT_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "MSE_IT",array_layout);
    [expS, S] = evol_score_stats(MSE_IT_vec, "MSE_IT", gen_idx_seq,nat_idx_seq,expS,S);
    corr_IT_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "corr_IT",array_layout);
    [expS, S] = evol_score_stats(corr_IT_vec, "corr_IT", gen_idx_seq,nat_idx_seq,expS,S);
    MSE_V1V4_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "MSE_V1V4",array_layout);
    [expS, S] = evol_score_stats(MSE_V1V4_vec, "MSE_V1V4", gen_idx_seq,nat_idx_seq,expS,S);
    corr_V1V4_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "corr_V1V4",array_layout);
    [expS, S] = evol_score_stats(corr_V1V4_vec, "corr_V1V4", gen_idx_seq,nat_idx_seq,expS,S);
    MSE_V4_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "MSE_V4",array_layout);
    [expS, S] = evol_score_stats(MSE_V4_vec, "MSE_V4", gen_idx_seq,nat_idx_seq,expS,S);
    corr_V4_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "corr_V4",array_layout);
    [expS, S] = evol_score_stats(corr_V4_vec, "corr_V4", gen_idx_seq,nat_idx_seq,expS,S);
    MSE_V1_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "MSE_V1",array_layout);
    [expS, S] = evol_score_stats(MSE_V1_vec, "MSE_V1", gen_idx_seq,nat_idx_seq,expS,S);
    corr_V1_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "corr_V1",array_layout);
    [expS, S] = evol_score_stats(corr_V1_vec, "corr_V1", gen_idx_seq,nat_idx_seq,expS,S);
    %%
    expStat = struct();
    expStat.scores = expS;
    expStat.meta = meta;
    expStat.imageName = imageName;
    expStat.meta.figdir = figdir;
    expStat.meta.fdrnm = fdrnm;
    
    expStat.evol.targ_imgnm = targ_imgnm;
    expStat.evol.targ_matnm = targ_matnm;
    expStat.evol.targ_area = targ_area;
    expStat.evol.score_mode = score_mode;
    expStat.evol.GANstr = GANstr;
    expStat.evol.imgpos = imgpos;
    expStat.evol.imgsize = imgsize;
    expStat.evol.optimstr = optimstr;
    expStat.evol.explabel = explabel; % added
    expStat.evol.thread_msks = thread_msks;
    expStat.evol.thread_num = thread_num;
    expStat.evol.row_gen = row_gen;
    expStat.evol.row_nat = row_nat;
    expStat.evol.gen_idx_seq = gen_idx_seq;
    expStat.evol.nat_idx_seq = nat_idx_seq;
    expStat.evol.block_arr = block_arr;
    expStat.evol.block_num = block_num;
    expStat.activ.act_tsr = act_tsr;
    expStat.activ.bslMat = bslMat;

    expStat.targ.meanMat = meanActMat; 
    expStat.targ.stdMat = stdActMat; 
    expStat.targ.targetMat = targetActMat; 
    expStat.targ.FMask = FMat; 
    expStat.targ.baseMask = baseMask; 

    S.Animal = Animal;
    S.ephysFN = meta.ephysFN;
    S.fdrnm = fdrnm;
    S.targ_imgnm = targ_imgnm;
    S.targ_annot = target_annot_map(targ_imgnm);
    S.targ_area = targ_area;
    S.targ_unitN = nansum(objMask,'all');
    S.score_mode = score_mode;
    S.GANstr = GANstr;
    S.imgpos = imgpos;
    S.imgsize = imgsize;
    S.optimstr = optimstr;
    save(fullfile(figdir,"expStat.mat"), 'expStat')
    save(fullfile(figdir,"expSummary.mat"), 'S')
    CosSummary = [CosSummary,S];
    CosStats = [CosStats,expStat];
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
save(fullfile(matdir, Animal+"_CosStats.mat"), "CosSummary", "CosStats")
%%
CosSumTab = struct2table(CosSummary); 
writetable(CosSumTab, fullfile(matdir, Animal+"_CosineSummary.csv"))
writetable(CosSumTab, fullfile(saveroot, "summary", Animal+"_CosineSummary.csv"))
%%
CosSumTab = readtable(fullfile(saveroot, "summary", Animal+"_CosineSummary.csv"));
%% create masks for comparison. 
fc6msk = contains(CosSumTab.GANstr,'fc6');
BGmsk = contains(CosSumTab.GANstr,'BigGAN');
SG2msk = contains(CosSumTab.GANstr,'StyleGAN2');
MSEmsk = contains(CosSumTab.score_mode,"MSE");
dotmsk = contains(CosSumTab.score_mode,"dot");
L1msk = contains(CosSumTab.score_mode,"L1");
CCmsk = contains(CosSumTab.score_mode,"corr");
ITmsk = contains(CosSumTab.targ_area,"IT");
V1msk = contains(CosSumTab.targ_area,"V1");
I4msk = contains(CosSumTab.targ_area,"V4");
%%
figdir = fullfile(saveroot, "summary");
h = stripe_plot(CosSumTab, "corr_bestscore", {fc6msk,BGmsk,SG2msk}, ["FC6","BigGAN","StyleGAN"], ...
                    "all Exp", "GAN_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
h = stripe_plot(CosSumTab, "MSE_bestscore", {fc6msk,BGmsk,SG2msk}, ["FC6","BigGAN","StyleGAN"], ...
                    "all Exp", "GAN_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
%% 
h = stripe_plot(CosSumTab, "corr_bestscore", {MSEmsk,CCmsk,dotmsk,L1msk}, ["MSE","CC","dot","L1"], ...
                    "all Exp", "Obj_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
h = stripe_plot(CosSumTab, "MSE_bestscore", {MSEmsk,CCmsk,dotmsk,L1msk}, ["MSE","CC","dot","L1"], ...
                    "all Exp", "Obj_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);

%%
% cos_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, objMask, "cosine_V1V4",array_layout);
% cos_vec = scorePopulationResponse_vec(evk_tsr, targetActMat, meanActMat, stdActMat, baseMask, "corr_IT",array_layout);
% figure;plot(block_arr(row_gen), cos_vec(row_gen), '.') % Fast visualization of a stats
% parse_mode2mask("V1ITV4","Alfa")
function [D, S] = evol_score_stats(scorevec,pfx,gen_idx_seq,nat_idx_seq,D,S)
    % This will clip the last generation in gen_idx_seq, nat_idx_seq by default. 
    if nargin == 1, prefix=""; end
    if nargin <= 4, S=struct(); D=struct(); end
    if ~strcmp(pfx,"ol"), % if not online score, then sort
        scorevec = reshape(scorevec,[],1);% use column vector by default.
        D.(pfx+"_vec") = scorevec; 
        scores_rec = cellfun(@(idx)scorevec(idx),reshape(gen_idx_seq(1:end-1),[],1),'uni',0);
        natscores_rec = cellfun(@(idx)scorevec(idx),reshape(nat_idx_seq(1:end-1),[],1),'uni',0);
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
% 
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