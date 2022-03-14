function visualize_Cosine_score_traj(S,figh)
for iTr = 1:numel(S)
    if nargin==1, figh=figure('pos',[472   554   490  440]);
    else, figh = figure(figh);set(figh,'pos',[472   554   490  440]);
    end 
    plotStatsTraj_S(S(iTr).score.offline_vec,...
        'Target Score',"score_offline",S(iTr),figh);
    % other score scheme could be added. 
end
end

function [figh, mean_corr_gen, sem_corr_gen, mean_corr_nat, sem_corr_nat] = ...
    plotStatsTraj_S(stat_vec,statname,savestr,S,figh)
% Plot trajectory of certain stats through out the evolution exp. (with
% simpler interface using the structure S. 
% Stat_vec: A vector, same length as trial num. This func will sort the
%   stat into each block and generated and natural. Then the scores will be
%   displayed as errorbar plot. The function will automatically clip out the last 
%   generation since it's usually incomplete. 
% S: Sturcture of an PCCosine Evolution  
% 
if nargin<5, figh = figure; else, clf(figh); end
gen_idx_seq = S.stim.gen_idx_seq;
nat_idx_seq = S.stim.nat_idx_seq;
block_arr = S.stim.block_arr;
mean_corr_gen = cellfun(@(idx)mean(stat_vec(idx)),gen_idx_seq(1:end-1));
sem_corr_gen = cellfun(@(idx)sem(stat_vec(idx)),gen_idx_seq(1:end-1));
mean_corr_nat = cellfun(@(idx)mean(stat_vec(idx)),nat_idx_seq(1:end-1));
sem_corr_nat = cellfun(@(idx)sem(stat_vec(idx)),nat_idx_seq(1:end-1));
row_gen = cell2mat(reshape(gen_idx_seq(1:end-1),[],1));
row_nat = cell2mat(reshape(nat_idx_seq(1:end-1),[],1));
Cord = colororder;

figure(figh);set(figh,'pos',[472   554   490  440]);hold on
shadedErrorBar([],mean_corr_gen,sem_corr_gen,'patchSaturation',0.4,'lineprops',{'color',[Cord(1,:),0.6],'LineWidth',2.5})
shadedErrorBar([],mean_corr_nat,sem_corr_nat,'patchSaturation',0.4,'lineprops',{'color',[Cord(2,:),0.6],'LineWidth',2.5})
ylabel(statname,'FontSize',12);
xlabel("Generation",'FontSize',12);
legend(["Gen","Nat"],'FontSize',10,'Location','Best')
title(compose("%s Evol Trajectory\n%s",statname,S.meta.explabel),'interpreter','none')
saveallform(S.meta.figdir,savestr+"_scoretraj",figh)
clf(figh);hold on
scatter(block_arr(row_gen),stat_vec(row_gen),25,'MarkerEdgeAlpha',0.5);
scatter(block_arr(row_nat),stat_vec(row_nat),25,'MarkerEdgeAlpha',0.5);
shadedErrorBar([],mean_corr_gen,sem_corr_gen,'patchSaturation',0.4,'lineprops',{'color',[Cord(1,:),0.6],'LineWidth',2.5})
shadedErrorBar([],mean_corr_nat,sem_corr_nat,'patchSaturation',0.4,'lineprops',{'color',[Cord(2,:),0.6],'LineWidth',2.5})
ylabel(statname,'FontSize',12);
xlabel("Generation",'FontSize',12);
legend(["Gen","Nat"],'FontSize',10,'Location','Best')
title(compose("%s Evol Trajectory\n%s",statname,S.meta.explabel),'interpreter','none')
saveallform(S.meta.figdir,savestr+"_scoretraj_w_pnts",figh)
end