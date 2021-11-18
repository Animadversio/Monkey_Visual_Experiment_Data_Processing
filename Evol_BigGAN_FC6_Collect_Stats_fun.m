function BFEStats = Evol_BigGAN_FC6_Collect_Stats_fun(meta_new, rasters_new, Trials_new)
% Compress the information relevant to the preferred channel evolution into a small structure
% Create a folder containing that structure.
saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_BigGAN_FC6_cmp"; 
for iTr = 1:numel(meta_new)
meta = meta_new{iTr};
if contains(meta.ephysFN,["Alfa","ALfa"]), Animal = "Alfa";
elseif contains(meta.ephysFN,["Beto"]), Animal = "Beto";end
rasters = rasters_new{iTr};
% lfps = lfps_new{Triali};
Trials = Trials_new{iTr};
fprintf("Processing experiment %s %s\n",meta.ephysFN, meta.expControlFN);
%% Record Basic Infos
Expi = iTr;% FIXME! Expi should be decoded from the ExpRecord.
BFEStats(iTr).Animal = Animal;
% BFEStats(Expi).Expi = Expi; 
BFEStats(iTr).imageName = Trials.imageName;
BFEStats(iTr).meta = meta;
if ~isfield(Trials,"TrialRecord") || isempty(Trials.TrialRecord), fprintf("Missing Trial Record\n");continue; end
pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
imgpos = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,2));
imgsize = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,3));
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
if isfield(meta,"unitID")
unit_num_arr = meta.unitID;
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
activ_msk = unit_num_arr~=0;
else
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
end
% note if the evolved channel is marked as null 0 then this can fail! 
try % important to continue going when facing issues.
pref_chan_id = zeros(1, thread_num);
for i = 1:thread_num % note here we use the unit_num_arr which exclude the null channels.
    pref_chan_id(i) = find(meta.spikeID==pref_chan(i) & ... % the id in the raster and lfps matrix 
                    unit_num_arr==unit_in_pref_chan(i)); % match for unit number
end
assert(pref_chan_id(1)==pref_chan_id(2))
catch err % err is an MException struct, save the err in a log file and continue. 
    fprintf('Error message:\n%s\n',err.message);
    fprintf('Error trace:\n%s\n',err.getReport);
    disp(meta)
    continue
end
BFEStats(iTr).units.pref_chan = pref_chan;
BFEStats(iTr).units.unit_name_arr = unit_name_arr;
BFEStats(iTr).units.unit_num_arr = unit_num_arr;
BFEStats(iTr).units.activ_msk = activ_msk;
BFEStats(iTr).units.spikeID = meta.spikeID;
BFEStats(iTr).units.pref_chan_id = pref_chan_id;
BFEStats(iTr).units.unit_in_pref_chan = unit_in_pref_chan;

if size(Trials.TrialRecord.User.evoConfiguration,2) == 5
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, strrep(string(Trials.TrialRecord.User.evoConfiguration{i,end}),"_"," ")]; % use this substitution for better printing effect
end
elseif size(Trials.TrialRecord.User.evoConfiguration,2) == 4 % old fashioned config
    Optim_names = repmat("CMAES", 1,thread_num);
end
if isfield(Trials.TrialRecord.User, "space")
    BFEStats(iTr).evol.space_names = Trials.TrialRecord.User.space; % space names and settings
    BFEStats(iTr).evol.space_cfg = Trials.TrialRecord.User.space_cfg; 
else % handle some early experiments which have a different structure. 
    if strcmp(meta.ephysFN,'Alfa-06072020-002')
    BFEStats(iTr).evol.space_names = {'fc6';'BigGAN'};% space names and settings
    BFEStats(iTr).evol.space_cfg = {{'fc6',[]};{'BigGAN',[]}};
    elseif strcmp(meta.ephysFN,'Alfa-06072020-003')
    BFEStats(iTr).evol.space_names = {'fc6';'BigGAN_class'};% space names and settings
    BFEStats(iTr).evol.space_cfg = {{'fc6',[]};{'BigGAN_class',[]}};
    else,keyboard;
    end
end
BFEStats(iTr).evol.optim_names = Optim_names;
BFEStats(iTr).evol.thread_num = thread_num;
BFEStats(iTr).evol.imgpos = imgpos;
BFEStats(iTr).evol.imgsize = imgsize;
BFEStats(iTr).evol.unit_in_pref_chan = unit_in_pref_chan; 
BFEStats(iTr).evol.pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,1));
% Sort the image names into threads, gens, nats and blocks
% seperate the thread natural images and generated images 
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
assert(sum(cellfun(@sum, thread_msks)) == length(imgnm)) % images comes from all these threads
BFEStats(iTr).stim.gen_msk = row_gen;
BFEStats(iTr).stim.nat_msk = row_nat;
BFEStats(iTr).stim.thread_msks = thread_msks;
%
block_arr = cell2mat(Trials.block);
block_list = min(block_arr):max(block_arr);
block_num = length(block_list);
color_seq = brewermap(block_num, 'spectral');
BFEStats(iTr).color_seq = color_seq;
BFEStats(iTr).evol.block_arr = block_arr;
BFEStats(iTr).evol.block_n = block_num;
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
% Collect trial by trial PSTH for the preferred channesl, for ploting and visualization. 
gen_psth_col = cellfun(@(idx) rasters(pref_chan_id(1), :, idx), gen_idx_seq, 'Uni', 0);
nat_psth_col = cellfun(@(idx) rasters(pref_chan_id(1), :, idx), nat_idx_seq, 'Uni', 0);
% Collect trial by trial response for the every channels, for ANOVA
gen_rspmat_col = cellfun(@(idx) squeeze(mean(rasters(:, 51:200, idx),[2])), gen_idx_seq, 'Uni', 0);
nat_rspmat_col = cellfun(@(idx) squeeze(mean(rasters(:, 51:200, idx),[2])), nat_idx_seq, 'Uni', 0);
BFEStats(iTr).evol.idx_seq = gen_idx_seq;
BFEStats(iTr).evol.psth = gen_psth_col;
BFEStats(iTr).evol.rspmat = gen_rspmat_col;
BFEStats(iTr).ref.idx_seq = nat_idx_seq;
BFEStats(iTr).ref.psth = nat_psth_col;
BFEStats(iTr).ref.rspmat = nat_rspmat_col;
% %% BigGAN images loading
% BGmsk = row_gen & thread_msks{2};
% BGimgnms = imgnm(BGmsk);
% BGimgs = arrayfun(@(nm)imread(fullfile(meta.stimuli, nm+".bmp")),BGimgnms,"Uni",0);
% %% imall = cat(4,BGimgs{:,:});
% distmat_evBG = D.distmat_B(cat(4,BGimgs{:,:})); % 240 sec for 950 images
% toc
% BFEStats(Expi).evol.distmat_BG = distmat_evBG;
%% Create the folder and file at last
stimparts = split(meta.stimuli,"\");
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
% fdrnm = compose("%s-%s-Chan%02d-1",datestr(expday,'yyyy-mm-dd'), Animal, pref_chan(1));
fdrnm = compose("%s-Chan%02d", stimparts{end-1}, pref_chan(1));
figdir = fullfile(saveroot, fdrnm);
BFEStats(iTr).meta.fdrnm = fdrnm; 
BFEStats(iTr).meta.figdir = figdir;
if exist(figdir),warning("%s figure directory exist! Beware",figdir);end
mkdir(figdir)
EvolStat = BFEStats(iTr);
save(fullfile(figdir,"EvolStat.mat"),'EvolStat')
fprintf("ExpStats saved to %s EvolStat.mat!\n",figdir)
end
end