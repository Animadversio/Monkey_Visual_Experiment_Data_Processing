% StyleGAN 
Animal="Alfa";Set_Path;
ftr = find(contains(ExpRecord.ephysFN,"21102020") & contains(ExpRecord.Exp_collection, "StyleGAN_Evol") );
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr,Animal);
%%
Triali = 1;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
%
pref_chan = Trials.TrialRecord.User.prefChan;
assert(all(pref_chan==pref_chan(1)))
pref_chan=pref_chan(1);
saveroot = "E:\OneDrive - Washington University in St. Louis\StyleGAN_evol";
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
fdrnm = compose("%s-%s-Chan%02d-1",datestr(expday,'yyyy-mm-dd'), Animal, pref_chan(1));
figdir = fullfile(saveroot, fdrnm);
mkdir(figdir)
%% 
pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
imgpos = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,2));
imgsize = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,3));
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
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
act_mat = squeeze(mean(rasters(:,51:200,:),2));
%%
gen_clr = {'red','magenta'}; nat_clr = {'green','cyan'};
figure;
iCh=10;
for iCh = 1:numel(meta.spikeID)
cla(gca, 'reset')
gen_act_mean = cellfun(@(idx)mean(act_mat(iCh,idx),2), gen_idx_seq(:,1:end-1));
gen_act_sem = cellfun(@(idx)std(act_mat(iCh,idx),1,2)/sqrt(length(idx)), gen_idx_seq(:,1:end-1));
nat_act_mean = cellfun(@(idx)mean(act_mat(iCh,idx),2), nat_idx_seq(:,1:end-1));
nat_act_sem = cellfun(@(idx)std(act_mat(iCh,idx),1,2)/sqrt(length(idx)), nat_idx_seq(:,1:end-1));
% plot(block_list(1:end-1), gen_act_mean')
for thri = 1:thread_num
   shadedErrorBar(block_list(1:end-1), gen_act_mean(thri,:), gen_act_sem(thri,:),'lineProps',{'Color', gen_clr{thri}})
end
for thri = 1:thread_num
   shadedErrorBar(block_list(1:end-1), nat_act_mean(thri,:), nat_act_sem(thri,:),'lineProps',{'Color', nat_clr{thri}})
end
xlabel("Generation");ylabel("Evoked Firing Rate");box off
title(unit_name_arr(iCh))
pause
end