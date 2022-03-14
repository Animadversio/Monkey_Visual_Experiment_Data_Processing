Animal = "Beto";Set_Path;
% "220218", "220216", "220221"
% ftrrows = find(contains(ExpRecord.expControlFN,["220218"]) | ...
%                contains(ExpRecord.expControlFN,["220221"]) );
%   contains(ExpRecord.expControlFN,["generate_BigGAN_PCcosine"]));
saveroot = "O:\Evol_PCCosine";
ftrrows = find(contains(ExpRecord.expControlFN,...
    ["220225","220119"]));
%     ["220121","220209","220211","220216","220218","220221"]));
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows, Animal, false);
%%
bhvfns = ExpRecord.expControlFN(ftrrows);%{cell2mat(meta_new).expControlFN}';
%% Process the RF experiments in a modular fashion.
rfidx = contains(bhvfns,"rfMapper");
RFS_col = RF_Calc_Stats_fun(meta_new(rfidx), rasters_new(rfidx), Trials_new(rfidx));
%%  Process RF data and get the masks saved to disk. 
for iRF = 1:numel(RFS_col)
    RFStat = RFS_col(iRF);
    maskS = RFStats_indiv_chan_gen_mask(RFStat);
    expdir = fullfile(saveroot, compose("%s-%s-RF",datestr(RFStat.meta.datetime,"yyyy-mm-dd"),RFStat.Animal));
    mkdir(expdir)
    save(fullfile(expdir,'RFStat.mat'),'RFStat')
    save(fullfile(expdir,'maskStat.mat'),'maskS')
end
%%

%% Process the Selectivity experiments in a modular fashion. 
selidx = contains(bhvfns,"selectivity_basic");
SelS_col = selectivity_Collect_Stats_fun(meta_new(selidx), rasters_new(selidx), Trials_new(selidx));
%% Visualize response distribution of all channels 
visusalize_resp_distri_allchan(SelS_col);
%% Visualize correlation structure and PCA spectrum 
...
%%
for i = 1:numel(SelS_col)
    seldir = fullfile(saveroot,SelS_col(i).meta.fdrnm);
    SelS_col(i).meta.figdir = seldir; ReprStat = SelS_col(i);
    mkdir(seldir);save(fullfile(seldir,"ReprStat.mat"),'ReprStat')
    fprintf("Selectivity Exp stats saved to %s\n",fullfile(seldir,"ReprStat.mat"))
end
%% Process the RF experiments in a modular fashion. 
evoidx = contains(bhvfns,"generate_BigGAN_PCcosine");
CosStats = Evol_PCCosine_Collect_Stats_fun(meta_new(evoidx), rasters_new(evoidx), Trials_new(evoidx));
%% Visualize these
visualize_PCCosine_PopEvol(CosStats(16:19),5);
visualize_PCCosine_score_traj(CosStats(16:19),10);
visualize_PCCosine_imageEvol(CosStats(16:19),11,12);

function visusalize_resp_distri_allchan(SelS_col)
for iTr = 1:numel(SelS_col)
for iCh = 1:numel(SelS_col(iTr).units.spikeID)
chan = SelS_col(iTr).units.spikeID(iCh);
unit = SelS_col(iTr).units.unit_num_arr(iCh);
label = SelS_col(iTr).units.unit_name_arr(iCh);
figh1 = vis_selectivity_static(SelS_col(iTr),true,chan,unit,'top',25, 1);
figh2 = vis_selectivity_static(SelS_col(iTr),true,chan,unit,'bottom',25, 2);
figh3 = vis_selectivity_static(SelS_col(iTr),false,chan,unit,'top',25, 3);
saveallform(SelS_col(iTr).meta.figdir,compose("resp_dist_%s-top",label),figh1,["jpg"])%,"pdf"
saveallform(SelS_col(iTr).meta.figdir,compose("resp_dist_%s-bottom",label),figh2,["jpg"])%,"pdf"
saveallform(SelS_col(iTr).meta.figdir,compose("resp_dist_%s-top_unsort",label),figh3,["jpg"])%,"pdf"
end
end
end

function visualize_PCCosine_imageEvol(S,figh,figh2)
for iTr = 1:numel(S)
    if nargin==1,  
        figh = figure('pos',[511    72   970   900]);
        figh2 = figure('pos',[511    72   1030   900]);
    else, 
        figh = figure(figh);set(figh,'pos',[511    72   970   900]);
        figh2 = figure(figh2);set(figh2,'pos',[511    72   1030   900]);
    end 
    gen_idx_seq = S(iTr).stim.gen_idx_seq;
    best_imnam_col = strings();
    for geni = 1:numel(gen_idx_seq)
        imgnm = S(iTr).imageName(gen_idx_seq{geni}(1));
        best_imnam_col(geni) = fullfile(S(iTr).meta.stimuli,imgnm+".bmp");
    end
    %%
    figh = figure(figh);
    montage(best_imnam_col(1:end-1))
    set(figh,'pos',[511    72   970   900])
    title(compose("%s",S(iTr).meta.explabel),'interpreter','none','FontSize',12)
    saveallform(S(iTr).meta.figdir,"Image_Evol_per_gen",figh,["jpg","pdf"])
    
    %% Offline computed scores
    scores_rec = cellfun(@(idx)S(iTr).score.offline_vec(idx),gen_idx_seq,'uni',0);
    meanscores = cellfun(@mean, scores_rec);
    figh2 = figure(figh2);
    [scoreframe_imgs, Clim] = score_frame_image_arr(best_imnam_col(1:end-1),...
                                            meanscores(1:end-1));
    montage(scoreframe_imgs)
    caxis([Clim]); cb = colorbar(); cb.Label.Interpreter = 'none';
    cb.Label.String = S(iTr).targ.score_mode{1}; cb.Label.FontSize = 12;
    set(figh2,'pos',[511    72   1030   900])
    title(compose("%s",S(iTr).meta.explabel),'interpreter','none','FontSize',12)
    saveallform(S(iTr).meta.figdir,"Image_Evol_per_gen_score_framed",figh2,["jpg","pdf"])

end
end

function visualize_PCCosine_score_traj(S,figh)
for iTr = 1:numel(S)
    if nargin==1, figh=figure('pos',[472   554   490  440]);
    else, figh = figure(figh);set(figh,'pos',[472   554   490  440]);
    end 
    plotStatsTraj_S(S(iTr).score.offline_vec,...
        'Target Score',"score_offline",S(iTr),figh);
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

function visualize_PCCosine_PopEvol(S,figh)
for iTr = 1:numel(S)
    if nargin==1, figh=figure('pos',[472   499   895   445]);
    else, figh = figure(figh);end 
    reprMat = S(iTr).resp.evoke_trials - S(iTr).resp.bslmean;
    norm_reprMat = (S(iTr).resp.evoke_trials - S(iTr).resp.bslmean - S(iTr).targ.meanActVec{1}) ./ S(iTr).targ.stdActVec{1};
    mask = S(iTr).targ.maskVec{1};
    targVec = S(iTr).targ.targetActVec{1};
    % [V1msk, V4msk, ITmsk] = get_areamask(S(iTr).units.spikeID, array_layout)
    [targmsk, targ_area] = parse_mode2maskVec(S(iTr).targ.score_mode, ...
                            S(iTr).meta.array_layout, S(iTr).units.spikeID);
    
    popul_evol_plot(reprMat, targVec, mask & targmsk, S(iTr), figh)
    title(compose("PopRepr Evol in %s\n%s",targ_area,S(iTr).meta.explabel),'interpreter','none')
    ylabel("Evoked - bsl rate")
    saveallform(S(iTr).meta.figdir,"RawAct_Pattern_TargArea_Evolution",figh)

    popul_evol_plot(norm_reprMat, targVec, mask, S(iTr), figh)
    title(compose("PopRepr Evol\n%s",S(iTr).meta.explabel),'interpreter','none')
    ylabel("zscore (Evoked - bsl rate)");
    saveallform(S(iTr).meta.figdir,"NormAct_Pattern_All_Evolution",figh)
    
    popul_evol_plot(norm_reprMat, targVec, mask & targmsk, S(iTr), figh)
    title(compose("PopRepr Evol in %s\n%s",targ_area,S(iTr).meta.explabel),'interpreter','none')
    ylabel("zscore (Evoked - bsl rate)");
    saveallform(S(iTr).meta.figdir,"NormAct_Pattern_TargArea_Evolution",figh)
    end

end

function popul_evol_plot(norm_reprMat, targVec, maskidx, S, figh)
    reprBlkMat = cell2mat(cellfun(@(idx)mean(norm_reprMat(:,idx),2),S.stim.gen_idx_seq,'uni',0));
    reprBlkMat_sem = cell2mat(cellfun(@(idx)sem(norm_reprMat(:,idx),2),S.stim.gen_idx_seq,'uni',0));
    reprBlkMat_nat = cell2mat(cellfun(@(idx)mean(norm_reprMat(:,idx),2),S.stim.nat_idx_seq,'uni',0));
    blockN = max(S.stim.block_arr);
    clrseq = brewermap(blockN-1,'Spectral');%from red to blue/purple-ish
    % [sortedchan, sortidx] = sort(chanArr);%sortidx to make sure the channel follows from 1:64
    figure(figh);clf; hold on;set(figh,'pos',[472   499   895   445])
    for i =1:blockN-1
        plot(reprBlkMat(maskidx,i),'color',[clrseq(i,:),0.5],'LineWidth',1.5) % sortidx
        plot(reprBlkMat_nat(maskidx,i),':','color',[clrseq(i,:),0.3],'LineWidth',1.5) % sortidx
        %shadedErrorBar(-249:500,psthavg_col{threadi,i},psthsem_col{threadi,i},'lineProps',{'color',[clrseq(i,:),0.7],'LineWidth',1.5})
    end
    plot(targVec(maskidx),'color',[0,0,0,0.75],'LineWidth',2.5);
    channum_msk = S.units.spikeID(maskidx);
    ITV1sep = sum(channum_msk<=32)+0.5;
    V1V4sep = sum(channum_msk<=16)+0.5;
    vline([ITV1sep,V1V4sep],'-.r')
    xticks(1:sum(maskidx))
    xticklabels(S.units.unit_name_arr(maskidx))
    xlim([0,sum(maskidx)+1])
    xlabel("Channels");ylabel("Activation");
end
