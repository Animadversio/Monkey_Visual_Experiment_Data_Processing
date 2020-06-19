function reportStats = calc_tune_stats(psth_cells,wdw_vect)
if nargin == 1
    wdw_vect = [];
    MovF = false;
else
    MovF = true;
end
% PSTH Chacteristics F to baseline, t of baseline 
reportStats = struct();
groupsize = cellfun(@(psth) size(psth,3), psth_cells);
indices = reshape(1:numel(psth_cells),size(psth_cells));
idx_vect = arrayfun(@(L, idx) idx*ones(L,1), groupsize, indices, 'UniformOutput', false);

act_wdw = 51:200; bsl_wdw=1:50;
score_vect = cellfun(@(psth)squeeze(mean(psth(1,act_wdw,:),2)),psth_cells,'UniformOutput',false);
basel_vect = cellfun(@(psth)squeeze(mean(psth(1,bsl_wdw,:),2)),psth_cells,'UniformOutput',false);
idx_vect = cell2mat(reshape(idx_vect,[],1));
score_vect = cell2mat(reshape(score_vect,[],1));
basel_vect = cell2mat(reshape(basel_vect,[],1));
[P,ANOVATAB,STATS] = anova1(score_vect,idx_vect,'off');
reportStats.F = ANOVATAB{2,5};
reportStats.F_P = P;
[P,ANOVATAB,STATS] = anova1(basel_vect,idx_vect,'off');
reportStats.F_bsl = ANOVATAB{2,5};
reportStats.F_P_bsl = P;
[H,P,CI,STATS] = ttest2(score_vect,basel_vect);
reportStats.T = STATS.tstat;
reportStats.t_P = P;
if MovF
F_wdw = [];
F_P_wdw = [];
for fi = 1:size(wdw_vect,1)
movscore_vect = cellfun(@(psth)squeeze(mean(psth(1,wdw_vect(fi,1):wdw_vect(fi,2),:),2)),psth_cells,'UniformOutput',false);
movscore_vect = cell2mat(reshape(movscore_vect,[],1));
[P,ANOVATAB,STATS] = anova1(movscore_vect,idx_vect,'off');
F_wdw(end+1) = ANOVATAB{2,5};
F_P_wdw(end+1) = P;
end
reportStats.F_wdw = F_wdw;
reportStats.F_P_wdw = F_P_wdw;
end
end