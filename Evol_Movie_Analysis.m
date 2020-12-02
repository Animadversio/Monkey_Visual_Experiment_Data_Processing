%% Evol_Movie
Animal = "Alfa"; Set_Path;
ftr = find(contains(ExpRecord.expControlFN,'201111'));
ExpRecord(ftr,:)
[meta_new, rasters_new, ~, Trials_new] = loadExperiments(ftr, Animal);
%%
saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_Movie";
Triali = 5; % 3:numel(meta_new)
rasters = rasters_new{Triali};
meta = meta_new{Triali};
Trials = Trials_new{Triali};
wdw = meta.rasterWindow;
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
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
unit_num_arr = meta.unitID;
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
%%
clrseq = brewermap(numel(psth_col),'Spectral');
for iCh = 1:numel(meta.spikeID)
    psth_col = cellfun(@(idx)mean(rasters(iCh,:,idx),3), gen_idx_seq, 'uni', false);
    psthsem_col = cellfun(@(idx)std(rasters(iCh,:,idx),1,3)/sqrt(numel(idx)), gen_idx_seq, 'uni', false);
    figure(1);clf;hold on
    for i =1:numel(psth_col)-1
    shadedErrorBar(-249:500,psth_col{i},psthsem_col{i},'lineProps',{'color',[clrseq(i,:),0.7],'LineWidth',1.5})
    end
    xlim([-50,400])
    title(unit_name_arr(iCh))
    hold off
    pause
end