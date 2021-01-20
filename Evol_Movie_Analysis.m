%% Evol_Movie
Animal = "Alfa"; Set_Path;
% ftr = find(contains(ExpRecord.expControlFN,'201111'));
% ftr = find(contains(ExpRecord.ephysFN, "Alfa-09122020-002") | ...
%            contains(ExpRecord.ephysFN, "Alfa-09122020-004") | ...
%            contains(ExpRecord.ephysFN, "Alfa-10122020-002")) ;


ftr = find(contains(ExpRecord.ephysFN, ...
    ["Alfa-11112020-006", "Alfa-11112020-007","Alfa-23112020-003", "Alfa-25112020-002", "Alfa-25112020-003", "Alfa-03122020-002", "Alfa-03122020-003", "Alfa-09122020-002", "Alfa-09122020-004", "Alfa-10122020-002"]));
ExpRecord(ftr,:)
[meta_new, rasters_new, ~, Trials_new] = loadExperiments(ftr, Animal);
%%
saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_Movie";
for Triali = 1:2 % 3:numel(meta_new)
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
fprintf("Going to create folder %s\n",fdrnm)
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
% separate nat, gen
imgnm = Trials.imageName;
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm) & ...% first 2 characters are not digits
        cellfun(@(c) ~contains(c(end-4:end), "_nat"), imgnm) ; % last 4 characters are not `_nat`
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the threads 1 and 2 (and maybe thread 3 4)
thread_msks = cell(1, thread_num);
for threadi = 1:thread_num
    msk = contains(imgnm, compose("thread%03d", threadi - 1));
    thread_msks{threadi} = msk; % store masks in a structure for the ease to iterate
end
assert(sum(cellfun(@sum, thread_msks)) == length(imgnm))
% separate block, nat, gen
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
bslwdw = [-50:45] - wdw(1);
actwdw = [51:400] - wdw(1);
meanscore_syn = cellfun(@(idx) mean(rasters(:,actwdw,idx),[2,3]) - mean(rasters(:,bslwdw,idx),[2,3]), gen_idx_seq', 'uni', 0);
semscore_syn = cellfun(@(idx) squeeze(std(mean(rasters(:,actwdw,idx),2),1,3)) / sqrt(numel(idx)), gen_idx_seq', 'uni', 0);
meanscore_nat = cellfun(@(idx) mean(rasters(:,actwdw,idx),[2,3]) - mean(rasters(:,bslwdw,idx),[2,3]), nat_idx_seq', 'uni', 0);
semscore_nat = cellfun(@(idx) squeeze(std(mean(rasters(:,actwdw,idx),2),1,3)) / sqrt(numel(idx)), nat_idx_seq', 'uni', 0);
meanscore_syn = cell2mat(reshape(meanscore_syn,[],size(gen_idx_seq',1),size(gen_idx_seq',2)));
semscore_syn = cell2mat(reshape(semscore_syn,[],size(gen_idx_seq',1),size(gen_idx_seq',2)));
meanscore_nat = cell2mat(reshape(meanscore_nat,[],size(nat_idx_seq',1),size(nat_idx_seq',2)));
semscore_nat = cell2mat(reshape(semscore_nat,[],size(nat_idx_seq',1),size(nat_idx_seq',2)));
%% Plot PSTH evolving
Exp_lab_str = compose("Movie Evol %s %s %s",meta.ephysFN, Trials.TrialRecord.User.space{1},...
        Trials.TrialRecord.User.evoConfiguration{1,end});
threadi = 1;
EvolStats = [];
clrseq = brewermap(block_num,'Spectral');
for iCh = 1:numel(meta.spikeID)
    psth_col = cellfun(@(idx)rasters(iCh,:,idx), gen_idx_seq, 'uni', false);
    [S,sumstr] = compute_evol_stats(psth_col,actwdw,bslwdw);
    EvolStats = [EvolStats;S];
    psthavg_col = cellfun(@(idx)mean(rasters(iCh,:,idx),3), gen_idx_seq, 'uni', false);
    psthsem_col = cellfun(@(idx)std(rasters(iCh,:,idx),1,3)/sqrt(numel(idx)), gen_idx_seq, 'uni', false);
    figure(1);clf;
    ax1=subplot('position',[0.05,0.09,0.59,0.76]);hold on
    for i =1:block_num-1
    shadedErrorBar(-249:500,psthavg_col{threadi,i},psthsem_col{threadi,i},'lineProps',{'color',[clrseq(i,:),0.7],'LineWidth',1.5})
    end
    xlabel("Time from Onset");xlim([-50,400]);vline(0.0)
    title([unit_name_arr(iCh),sumstr])
    hold off
    ax2=subplot('position',[.66,0.09,0.32,0.76]);hold on
    shadedErrorBar(block_list(1:end-1), meanscore_syn(iCh, 1:end-1, threadi), semscore_syn(iCh, 1:end-1, threadi),...
        'lineprops',{'Color',[0,0,0,0.7]},'transparent',1,'patchSaturation',0.045)
    shadedErrorBar(block_list(1:end-1), meanscore_nat(iCh, 1:end-1, threadi), semscore_nat(iCh, 1:end-1, threadi),...
        'lineprops',{'Color',[0,1,0,0.7]},'transparent',1,'patchSaturation',0.045)
    xlabel("generations"); axis tight
    legend(["Generated img","Natural img"],'Location',"Best")
    title([Exp_lab_str, unit_name_arr(iCh)])%[Exp_label_str, compose('Generation averaged score, channel %s', unit_name_arr{channel_j}), compose("Optimizer %s", Optim_names(threadi))])
    saveas(1,fullfile(figdir,compose("Evol_psth_traj_%s.png",unit_name_arr(iCh))))
    saveas(1,fullfile(figdir,compose("Evol_psth_traj_%s.pdf",unit_name_arr(iCh))))
%     saveas(1,fullfile(figdir,compose("Evol_psth_traj_%s.pdf",unit_name_arr(iCh))))
end
save(fullfile(figdir,'EvolStats.mat'),'EvolStats');
end
%%
function [S,sumstr] = compute_evol_stats(psth_col,actwdw,bslwdw)
% window = [51:200];
assert(size(psth_col,1)==1);
blockn = size(psth_col,2);
block_mean_score = cellfun(@(psth)mean(psth(1,actwdw,:),'all')-mean(psth(1,bslwdw,:),'all'),psth_col(1:blockn-1));
[~,peakBlock] = max(block_mean_score);
if peakBlock == blockn-1, peakBlock = blockn-2; end % avoid the last generation
endspsths = cell2mat(reshape(psth_col(blockn-2:blockn-1),1,1,[]));
peakpsths = cell2mat(reshape(psth_col(peakBlock:peakBlock+1),1,1,[]));
initpsths = cell2mat(reshape(psth_col(1:2),1,1,[]));

initacts = squeeze(mean(initpsths(1,actwdw,:),2));
peakacts = squeeze(mean(peakpsths(1,actwdw,:),2));
endsacts = squeeze(mean(endspsths(1,actwdw,:),2));

[~,P,CI,STATS] = ttest2(peakacts, initacts);
S.t_pk = STATS.tstat;
S.t_P_pk = P;
S.DA_pk = (mean(peakacts) - mean(initacts));
S.DAOA_pk = (mean(peakacts) - mean(initacts)) / mean(initacts);
[~,P,CI,STATS] = ttest2(endsacts, initacts);
S.t_ed = STATS.tstat;
S.t_P_ed = P;
S.DA_ed = (mean(endsacts) - mean(initacts));
S.DAOA_ed = (mean(endsacts) - mean(initacts)) / mean(initacts);
sumstr = compose("peak - init: dAct %.1f dA/A0 %.2f t=%.2f(%.1e)\n",S.DA_pk,S.DAOA_pk,S.t_pk,S.t_P_pk)+...
        compose("end - init: dAct %.1f dA/A0 %.2f t=%.2f(%.1e)",S.DA_ed,S.DAOA_ed,S.t_ed,S.t_P_ed);
S.sumstr = sumstr;
end