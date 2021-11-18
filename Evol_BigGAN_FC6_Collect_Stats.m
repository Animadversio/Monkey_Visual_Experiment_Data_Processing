saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_BigGAN_FC6_cmp"; 
setMatlabTitle("BigGAN_FC6 compare");

Animal = "Beto"; Set_Path; 
ftrrows = find(contains(ExpRecord.expControlFN,"generate_") & ...
              (ExpRecord.Exp_collection=="BigGAN_fc6" |...
               ExpRecord.Exp_collection=="BigGAN_FC6" |...
               ExpRecord.Exp_collection=="BigGAN_FC6_CMAHess"|...
               ExpRecord.Exp_collection=="BigGAN_Hessian"));
disp(ExpRecord(ftrrows,:))
%%
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows([43,45,46:end]), Animal, false);%43,45
%7,47;7,47,86:

%%
eptymsk = cellfun(@isempty,meta_new);
meta_new(eptymsk)=[];
rasters_new(eptymsk)=[];
Trials_new(eptymsk)=[];
%% New Function
BFEStats = Evol_BigGAN_FC6_Collect_Stats_fun(meta_new(1:end), rasters_new(1:end), Trials_new(1:end));
%%

%% Evol_BigGAN_FC6_Collect_Stats.m
BFEStats = repmat(struct(),1,numel(meta_new));
for Triali = 1:numel(meta_new)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
% lfps = lfps_new{Triali};
Trials = Trials_new{Triali};
fprintf("Processing experiment %s %s\n",meta.ephysFN, meta.expControlFN);
%% Record Basic Infos
Expi = Triali;% FIXME! Expi should be decoded from the ExpRecord.
BFEStats(Expi).Animal = Animal;
BFEStats(Expi).Expi = Expi; 
BFEStats(Expi).imageName = Trials.imageName;
BFEStats(Expi).meta = meta;

pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
imgpos = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,2));
imgsize = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,3));
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% note if the evolved channel is marked as null 0 then this can fail! 
pref_chan_id = zeros(1, thread_num);
for i = 1:thread_num % note here we use the unit_num_arr which exclude the null channels.
    pref_chan_id(i) = find(meta.spikeID==pref_chan(i) & ... % the id in the raster and lfps matrix 
                    unit_num_arr==unit_in_pref_chan(i)); % match for unit number
end
BFEStats(Expi).units.pref_chan = pref_chan;
BFEStats(Expi).units.unit_name_arr = unit_name_arr;
BFEStats(Expi).units.unit_num_arr = unit_num_arr;
BFEStats(Expi).units.activ_msk = activ_msk;
BFEStats(Expi).units.spikeID = meta.spikeID;
BFEStats(Expi).units.pref_chan_id = pref_chan_id;
BFEStats(Expi).units.unit_in_pref_chan = unit_in_pref_chan;

if size(Trials.TrialRecord.User.evoConfiguration,2) == 5
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, strrep(string(Trials.TrialRecord.User.evoConfiguration{i,end}),"_"," ")]; % use this substitution for better printing effect
end
elseif size(Trials.TrialRecord.User.evoConfiguration,2) == 4 % old fashioned config
    Optim_names = repmat("CMAES", 1,thread_num);
end
BFEStats(Expi).evol.space_names = Trials.TrialRecord.User.space; % space names and settings
BFEStats(Expi).evol.space_cfg = Trials.TrialRecord.User.space_cfg; 
BFEStats(Expi).evol.optim_names = Optim_names;
BFEStats(Expi).evol.thread_num = thread_num;
BFEStats(Expi).evol.imgpos = imgpos;
BFEStats(Expi).evol.imgsize = imgsize;
BFEStats(Expi).evol.unit_in_pref_chan = unit_in_pref_chan; 
BFEStats(Expi).evol.pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,1));
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
BFEStats(Expi).stim.gen_msk = row_gen;
BFEStats(Expi).stim.nat_msk = row_nat;
BFEStats(Expi).stim.thread_msks = thread_msks;
%
block_arr = cell2mat(Trials.block);
block_list = min(block_arr):max(block_arr);
block_num = length(block_list);
color_seq = brewermap(block_num, 'spectral');
BFEStats(Expi).color_seq = color_seq;
BFEStats(Expi).evol.block_arr = block_arr;
BFEStats(Expi).evol.block_n = block_num;
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
gen_psth_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), gen_idx_seq, 'Uni', 0);
nat_psth_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), nat_idx_seq, 'Uni', 0);
% Collect trial by trial response for the every channels, for ANOVA
gen_rspmat_col = cellfun(@(idx) squeeze(mean(rasters(:, 51:200, idx),[2])), gen_idx_seq, 'Uni', 0);
nat_rspmat_col = cellfun(@(idx) squeeze(mean(rasters(:, 51:200, idx),[2])), nat_idx_seq, 'Uni', 0);
BFEStats(Expi).evol.idx_seq = gen_idx_seq;
BFEStats(Expi).evol.psth = gen_psth_col;
BFEStats(Expi).evol.rspmat = gen_rspmat_col;
BFEStats(Expi).ref.idx_seq = nat_idx_seq;
BFEStats(Expi).ref.psth = nat_psth_col;
BFEStats(Expi).ref.rspmat = nat_rspmat_col;
% %% BigGAN images loading
% BGmsk = row_gen & thread_msks{2};
% BGimgnms = imgnm(BGmsk);
% BGimgs = arrayfun(@(nm)imread(fullfile(meta.stimuli, nm+".bmp")),BGimgnms,"Uni",0);
% %% imall = cat(4,BGimgs{:,:});
% distmat_evBG = D.distmat_B(cat(4,BGimgs{:,:})); % 240 sec for 950 images
% toc
% BFEStats(Expi).evol.distmat_BG = distmat_evBG;

stimparts = split(meta.stimuli,"\");
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
% fdrnm = compose("%s-%s-Chan%02d-1",datestr(expday,'yyyy-mm-dd'), Animal, pref_chan(1));
fdrnm = compose("%s-Chan%02d", stimparts{end-1}, pref_chan(1));
figdir = fullfile(saveroot, fdrnm);
mkdir(figdir)
end
%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
save(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 

%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
for Animal = ["Alfa","Beto"]
BFEStats = [];
fdrlist = strtrim(string(ls(saveroot)));
for i = 1:numel(fdrlist)
    fdrnm = fdrlist(i);
    if contains(fdrnm,Animal) && fdrnm~="2020-07-20-Beto-01"
    D = load(fullfile(saveroot,fdrnm,"EvolStat.mat"),'EvolStat');
    BFEStats = [BFEStats;D.EvolStat];
    end
end
fprintf("Found %d BigGAN FC6 paired experiment for %s\n",numel(BFEStats),Animal)
save(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
end