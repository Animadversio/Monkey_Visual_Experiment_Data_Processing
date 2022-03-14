%% Evol Cosine PCA Analysis 

Animal = "Beto";Set_Path;
% "220218", "220216"
ftrrows = find(contains(ExpRecord.expControlFN,["220218"]) | ...
               contains(ExpRecord.expControlFN,["220221"]) );
%   contains(ExpRecord.expControlFN,["generate_BigGAN_PCcosine"]));
% & contains(ExpRecord.expControlFN,["210218","210408","210413","210415"]));
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows, Animal, false);
%% Process the RF experiments in a modular fashion. 
RFS_col = RF_Calc_Stats_fun(meta_new([1,7]), rasters_new([1,7]), Trials_new([1,7]));
%% Process the Selectivity experiments in a modular fashion. 
Sel_S = selectivity_Collect_Stats_fun(meta_new([2,8]), rasters_new([2,8]), Trials_new([2,8]));
%% Process the RF experiments in a modular fashion. 


%% 
array_layout = "Beto_new";
saveroot = "O:\Evol_PCCosine";
MAXUNUM = 4;
figh = figure(1);
for iTr = [3:6,9:12]%3:numel(meta_new)%numel(meta_new)-6:
meta = meta_new{iTr};
rasters = rasters_new{iTr};
Trials = Trials_new{iTr};
imageName = string(Trials.imageName);
stdActMat = Trials.TrialRecord.User.stdActMat';
meanActMat = Trials.TrialRecord.User.meanActMat';
targetActMat = Trials.TrialRecord.User.targetActMat';
maskMat = Trials.TrialRecord.User.maskMat';
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

%% Form the activation matrix and tensor
act_mat = squeeze(mean(rasters(:,51:200,:),2));
% bsl_mat = squeeze(mean(rasters(:,1:50,:),[2,3])); % baseline matrix for each channel
bsl_mat = squeeze(mean(rasters(:,1:50,:),[2,3])); % baseline matrix for each channel
act_tsr = nan(MAXUNUM+1,64,size(rasters,3));
bslMat = nan(MAXUNUM+1,64);
for iCh = 1:size(rasters,1)
chani = meta.spikeID(iCh);
uniti = meta.unitID(iCh);
act_tsr(uniti+1,chani,:) = act_mat(iCh,:);
bslMat(uniti+1,chani) = bsl_mat(iCh);
end

%% Masks for calculating distance to target 
chan_arr = [1:64]';
unit_arr = 0:MAXUNUM;
chanMat = repmat(chan_arr,1,MAXUNUM+1)';
unitMat = repmat(unit_arr,64,1)';
for iThr = 1:thread_num
stdActMat = Trials.TrialRecord.User.stdActMat{iThr}';
meanActMat = Trials.TrialRecord.User.meanActMat{iThr}';
targetActMat = Trials.TrialRecord.User.targetActMat{iThr}';
FMat = Trials.TrialRecord.User.maskMat{iThr}'; % Mask in the ReprTsr file, F significance

mode = score_mode{iThr};
targetName = target_cfg{iThr}{2};
% targimg = imread(targmap(targetName)); % get the target image
objMask = parse_mode2mask(mode,array_layout)'; % Parse out the chanXunit mask used to evaluate objective. (area specific)
maskMat = objMask & FMat & ~isnan(targetActMat); % Channel Selection & F significant


%% Note to baseline subtract! or correlation will be bad.
vec_reprs = reshape(act_tsr(repmat(maskMat,1,1,size(act_tsr,3))),[],size(act_tsr,3)) - bslMat(maskMat); 
vec_targs = targetActMat(maskMat);
vec_norm_mean = meanActMat(maskMat);
vec_norm_std = stdActMat(maskMat);
vec_reprs_norm = (vec_reprs-vec_norm_mean)./vec_norm_std; % Z scored repr 
vec_targs_norm = (vec_targs-vec_norm_mean)./vec_norm_std; % Z scored target
vecchan = chanMat(maskMat);
vecunit = unitMat(maskMat);

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
saveallform(figdir,"online_scoretraj",figh0,["png","pdf"])

%% Print the population evolution pattern 
figh = popuEvol(vec_reprs,3*vec_targs,vecchan,vecunit,gen_idx_seq,nat_idx_seq,block_arr,figh);
ylabel("Raw Activation")
saveallform(figdir,"RawAct_Pattern_Evolution",figh,["png","pdf"])
figh = popuEvol(vec_reprs_norm,vec_targs,vecchan,vecunit,gen_idx_seq,nat_idx_seq,block_arr,figh);
ylabel("Z Activation")
saveallform(figdir,"NormAct_Pattern_Evolution",figh,["png","pdf"])

%% Image Evolution 
% Trials.imageName(row_gen)
best_imnam_col = strings();
for geni = 1:numel(gen_idx_seq)
imgnm = Trials.imageName(gen_idx_seq{geni}(1));
best_imnam_col(geni) = fullfile(meta.stimuli,imgnm+".bmp");
end
%%
figh = figure('pos',[511    72   970   900]);
montage(best_imnam_col)
set(figh,'pos',[511    72   970   900])
saveallform(figdir,"Image_Evol_per_gen",figh,["jpg","pdf"])
%% Offline computed scores
figh2 = figure('pos',[511    72   970   900]);
scoreframe_imgs = score_frame_image_arr(best_imnam_col(1:numel(scores_rec))',cellfun(@mean, scores_rec));
montage(scoreframe_imgs)
set(figh2,'pos',[511    72   970   900])
saveallform(figdir,"Image_Evol_per_gen_score_framed",figh2,["jpg","pdf"])
end


%% fetch population responses 


%% fetch target vectors 


end
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

function figh = popuEvol(reprMat,targVec,chanArr,unitId,gen_idx_seq,nat_idx_seq,block_arr,figh)
% This creates figure of population response pattern evolving.
global explabel
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

function imageEvolution()

end