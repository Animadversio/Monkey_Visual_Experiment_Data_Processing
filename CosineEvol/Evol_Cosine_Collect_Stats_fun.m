%% Evol_PCCosine_Collect_Stats_fun
function CosStats = Evol_Cosine_Collect_Stats_fun(meta_new, rasters_new, Trials_new)
% The objective evolve towards a target vector instead of a pattern. 
% 
% dataroot = "N:\Stimuli\2022-PCCosineEvol";
% saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_PCCosine";
dataroot = "N:\Stimuli\2020-CosineEvol";
saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_Cosine";
saveroot = "E:\OneDrive - Harvard University\Evol_Cosine";
for iTr = 1:numel(meta_new)
fprintf("Process Experiment %d\n",iTr)
% lfps = lfps_new{Triali};
meta = meta_new{iTr};
rasters = rasters_new{iTr};
Trials = Trials_new{iTr};
% Parse out animal name. 
if contains(meta.ephysFN,["Alfa","ALfa"]), Animal = "Alfa";
elseif contains(meta.ephysFN,["Beto"]), Animal = "Beto";
end

imageName = string(Trials.imageName);
target_cfg = Trials.TrialRecord.User.target_cfg;
score_mode = Trials.TrialRecord.User.score_mode;
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);

%% Now summarize key information for analysis
stimparts = split(meta.stimuli,"\");
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
fdrnm = compose("%s-%s", stimparts{end-1}, score_mode{1});
reffdrnm = compose("%s-%s-refRepr",datestr(expday,"yyyy-mm-dd"),Animal);
% fdrnm = compose("%s-%s-Chan%02d-1",datestr(expday,'yyyy-mm-dd'), Animal, pref_chan(1));
figdir = fullfile(saveroot, fdrnm);
if exist(figdir,'dir'),warning("%s figure directory exist! Beware",figdir);end
mkdir(figdir)
% Expi = iTr;% FIXME! Expi should be decoded from the ExpRecord.
CosStats(iTr).Animal = Animal;
CosStats(iTr).imageName = string(Trials.imageName);
CosStats(iTr).meta = meta;
CosStats(iTr).meta.expday = expday; 
CosStats(iTr).meta.fdrnm = fdrnm; 
CosStats(iTr).meta.figdir = figdir;
CosStats(iTr).meta.refdir = fullfile(dataroot, reffdrnm);
assert(exist(CosStats(iTr).meta.refdir,'dir'))
if strcmp(Animal,"Beto") ...
    && (expday > datetime(2021,9,1)) 
% after the Beto re-implantation. 
array_layout = "Beto_new";
else
array_layout = string(Animal);
end
CosStats(iTr).meta.array_layout = array_layout;

%% Collect Information and labels for each unit 
assert(isfield(meta,"unitID"))
unit_num_arr = meta.unitID;
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
activ_msk = unit_num_arr~=0;
CosStats(iTr).units.unit_name_arr = unit_name_arr;
CosStats(iTr).units.unit_num_arr = unit_num_arr;
CosStats(iTr).units.spikeID = meta.spikeID;
CosStats(iTr).units.activ_msk = activ_msk;

%% Settings of evolution for each thread. 
assert(thread_num == 1)
explabel = stimparts{end-1};
CosStats(iTr).evol.thread_num = thread_num;
for iThr = 1:thread_num
spacenm = Trials.TrialRecord.User.space_cfg{iThr}{1};
optimnm = Trials.TrialRecord.User.evoConfiguration{iThr,5};
imgpos = Trials.TrialRecord.User.evoConfiguration{iThr,2};
imgsize = Trials.TrialRecord.User.evoConfiguration{iThr,3};
explabel = explabel + compose(" %s (trg: %s)\n%s (%s)",score_mode{iThr},target_cfg{iThr}{2},spacenm,optimnm);
explabel = explabel + compose("   Pos [%.1f,%.1f] Size %.1f deg",imgpos,imgsize);

CosStats(iTr).evol.spacenm{iThr} = spacenm;
CosStats(iTr).evol.optimnm{iThr} = optimnm;
CosStats(iTr).evol.imgpos{iThr} = imgpos;
CosStats(iTr).evol.imgsize{iThr} = imgsize;
end

CosStats(iTr).meta.explabel = explabel; % Long string description of the experiments.

%% Masks for separating images 
[thread_msks,thread_num,row_gen,row_nat,gen_idx_seq,nat_idx_seq,block_arr,block_num] =...
    sort_trial_into_masks(Trials); 

CosStats(iTr).stim.gen_msk = row_gen;
CosStats(iTr).stim.nat_msk = row_nat;
CosStats(iTr).stim.thread_msks = thread_msks;
CosStats(iTr).stim.gen_idx_seq = gen_idx_seq;
CosStats(iTr).stim.nat_idx_seq = nat_idx_seq;
CosStats(iTr).stim.block_arr = block_arr;

%% Population activity for evolved images throughout evolution. 
evkwdw = [51:200];  bslwdw = [1:40];
bslmean = squeeze(mean(rasters(:, bslwdw, :),[2,3])); %average baseline activity 
bslsem = squeeze(sem(mean(rasters(:, bslwdw, :),[2]),3)); % variability between trials in baseline. 
evoke_trials  = squeeze(mean(rasters(:, evkwdw, :),[2])); % evoked firing rate
evkbsl_trials = evoke_trials - bslmean; % evk - bsl score using overall baseline 

CosStats(iTr).resp.evoke_trials = evoke_trials;
CosStats(iTr).resp.bslmean = bslmean;
CosStats(iTr).resp.bslsem = bslsem;
CosStats(iTr).resp.rspcol = cellfun(@(idx) evoke_trials(:,idx), gen_idx_seq,'uni',0);
CosStats(iTr).resp.rspcol_nat = cellfun(@(idx) evoke_trials(:,idx), nat_idx_seq,'uni',0);

%% Population activity for natural images 

%% Information of the target pattern and scoring method. 

meanActMat = Trials.TrialRecord.User.meanActMat; % cell array of meanActMat
stdActMat = Trials.TrialRecord.User.stdActMat;
targetActMat = Trials.TrialRecord.User.targetActMat;
% targetActStd = Trials.TrialRecord.User.targetActStd;
maskMat = Trials.TrialRecord.User.maskMat;

CosStats(iTr).targ.target_cfg = target_cfg;
CosStats(iTr).targ.score_mode = score_mode;
CosStats(iTr).targ.maskMat = maskMat;
CosStats(iTr).targ.meanActMat = meanActMat;
CosStats(iTr).targ.stdActMat = stdActMat;
CosStats(iTr).targ.targetActMat = targetActMat;
% CosStats(iTr).targ.targetActStd = targetActStd;
% convert the matrix repr of mask to vector for visualziation 
maskVec = {};
meanActVec = {};
stdActVec = {};
targetActVec = {};
for iThr = 1:thread_num
for iCh = 1:size(rasters,1)
chani = meta.spikeID(iCh);
uniti = meta.unitID(iCh);
if isnan(uniti), 
    fprintf("encounter nan unit %d, ignore and regarded as 0\n",uniti);
    uniti = 0;
end
maskVec{iThr}(iCh,1) = maskMat{iThr}(chani,uniti+1);
meanActVec{iThr}(iCh,1) = meanActMat{iThr}(chani,uniti+1);
stdActVec{iThr}(iCh,1) = stdActMat{iThr}(chani,uniti+1);
targetActVec{iThr}(iCh,1) = targetActMat{iThr}(chani,uniti+1);
end
end
if any(isnan(stdActVec{1}(maskVec{1}))) || ...
   any(isnan(meanActVec{1}(maskVec{1}))) || ...
   any(isnan(targetActVec{1}(maskVec{1})))
   fprintf("Un expected nan in mean, std, target encountered\n")
   keyboard;
   maskVec{1} = maskVec{1} & ~isnan(meanActVec{1}) ...
       & ~isnan(stdActVec{1}) & ~isnan(targetActVec{1});
end
assert(~any(isnan(stdActVec{1}(maskVec{1}))))
assert(~any(isnan(meanActVec{1}(maskVec{1}))))
assert(~any(isnan(targetActVec{1}(maskVec{1}))))

CosStats(iTr).targ.maskVec = maskVec;
CosStats(iTr).targ.meanActVec = meanActVec;
CosStats(iTr).targ.stdActVec = stdActVec;
CosStats(iTr).targ.targetActVec = targetActVec;
% Source of the representation matrix, and load it in. 
CosStats(iTr).targ.repr_path = fullfile(CosStats(iTr).meta.refdir, target_cfg{1}{1});
CosStats(iTr).targ.repr_D = load(CosStats(iTr).targ.repr_path);

%% Scores for population images
scores_rec = Trials.TrialRecord.User.scores_record; % recorded score from online Exp 
if isfield("natscores_record",Trials.TrialRecord.User)
    natscores_rec = Trials.TrialRecord.User.natscores_record; % recorded score for natural images
else
    natscores_rec = {}; % Temporary solution is to put empty cell array there. 
end
% genN_rec = arrayfun(@(b)b*ones(size(scores_rec{b})),[1:numel(scores_rec)]','uni',0); % format gen number as the score
% scores_vec = cell2mat(scores_rec);
% genN_vec = cell2mat(genN_rec);

CosStats(iTr).score.online_rec = scores_rec; 
CosStats(iTr).score.online_rec_ref = natscores_rec; 

[V1msk, V4msk, ITmsk] = get_areamask(meta.spikeID, array_layout);
[targmsk, targ_area] = parse_mode2maskVec(score_mode, array_layout, meta.spikeID);
M = struct("allmsk", ones(size(V1msk),'logical'), ...
    "V1msk", V1msk, "V4msk", V4msk, "ITmsk", ITmsk,"targmsk",targmsk);

CosStats(iTr).score.offline_vec = scorePopulationVec(evkbsl_trials, ...
    targetActVec{1}, meanActVec{1}, stdActVec{1}, maskVec{1}&targmsk, score_mode);

for area = ["V1","V4","IT","targ","all"]
for dist_fun = ["cosine","corr","dot"]
CosStats(iTr).score.(dist_fun+"_"+area) = scorePopulationVec(evkbsl_trials, ...
    targetActVec{1}, meanActVec{1}, stdActVec{1}, maskVec{1}&M.(area+"msk"), dist_fun);
end
end
% for dist_fun = ["cosine","corr","dot"]
% CosStats(iTr).score.(dist_fun+"_targ") = scorePopulationVec_direction(evkbsl_trials, ...
%     targetActVec{1}, meanActVec{1}, stdActVec{1}, maskVec{1}&targmsk, dist_fun);
% end

% BFEStats(iTr).evol.idx_seq = gen_idx_seq;
% BFEStats(iTr).evol.psth = gen_psth_col;
% BFEStats(iTr).evol.rspmat = gen_rspmat_col;
% BFEStats(iTr).ref.idx_seq = nat_idx_seq;
% BFEStats(iTr).ref.psth = nat_psth_col;
% BFEStats(iTr).ref.rspmat = nat_rspmat_col;


EvolStat = CosStats(iTr);
save(fullfile(figdir,"EvolStat.mat"),'EvolStat')
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

