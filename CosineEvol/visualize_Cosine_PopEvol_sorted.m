function visualize_Cosine_PopEvol_sorted(S,figh)
for iTr = 1:numel(S)
    if nargin==1, figh=figure('pos',[472   499   895   445]);
    else, figh = figure(figh);end 
    reprMat = S(iTr).resp.evoke_trials - S(iTr).resp.bslmean;
    norm_reprMat = (S(iTr).resp.evoke_trials - S(iTr).resp.bslmean - S(iTr).targ.meanActVec{1}) ./ S(iTr).targ.stdActVec{1};
    mask = S(iTr).targ.maskVec{1};
    targVec = S(iTr).targ.targetActVec{1};
    normtargVec = (S(iTr).targ.targetActVec{1} - S(iTr).targ.meanActVec{1})./ S(iTr).targ.stdActVec{1};
    [V1msk, V4msk, ITmsk] = get_areamask(S(iTr).units.spikeID, S(iTr).meta.array_layout);
    [targmsk, targ_area] = parse_mode2maskVec(S(iTr).targ.score_mode, ...
                            S(iTr).meta.array_layout, S(iTr).units.spikeID);
    
    popul_evol_plot_sorted_perarea(reprMat, targVec, {mask & targmsk & V1msk, ...
                              mask & targmsk & V4msk, mask & targmsk & ITmsk}, ...
                                    {"V1","V4","IT"}, S(iTr), figh)
    title(compose("PopRepr Evol in %s\n%s",targ_area,S(iTr).meta.explabel),'interpreter','none')
    ylabel("Evoked - bsl rate")
    saveallform(S(iTr).meta.figdir,"RawAct_Pattern_TargArea_Evolution_sorted_perarea",figh)

    popul_evol_plot_sorted_perarea(norm_reprMat, normtargVec, {mask & V1msk, mask & V4msk, mask & ITmsk}, ...
                                    {"V1","V4","IT"}, S(iTr), figh)
    title(compose("PopRepr Evol\n%s",S(iTr).meta.explabel),'interpreter','none')
    ylabel("zscore (Evoked - bsl rate)");
    saveallform(S(iTr).meta.figdir,"NormAct_Pattern_All_Evolution_sorted_perarea",figh)
    
    popul_evol_plot_sorted_perarea(norm_reprMat, normtargVec, {mask & targmsk & V1msk, ...
                                       mask & targmsk & V4msk, mask & targmsk & ITmsk}, ...
                                    {"V1","V4","IT"}, S(iTr), figh)
    title(compose("PopRepr Evol in %s\n%s",targ_area,S(iTr).meta.explabel),'interpreter','none')
    ylabel("zscore (Evoked - bsl rate)");
    saveallform(S(iTr).meta.figdir,"NormAct_Pattern_TargArea_Evolution_sorted_perarea",figh)

    popul_evol_plot_sorted(reprMat, targVec, mask & targmsk, S(iTr), figh)
    title(compose("PopRepr Evol in %s\n%s",targ_area,S(iTr).meta.explabel),'interpreter','none')
    ylabel("Evoked - bsl rate")
    saveallform(S(iTr).meta.figdir,"RawAct_Pattern_TargArea_Evolution_sorted",figh)

    popul_evol_plot_sorted(norm_reprMat, normtargVec, mask, S(iTr), figh)
    title(compose("PopRepr Evol\n%s",S(iTr).meta.explabel),'interpreter','none')
    ylabel("zscore (Evoked - bsl rate)");
    saveallform(S(iTr).meta.figdir,"NormAct_Pattern_All_Evolution_sorted",figh)
    
    popul_evol_plot_sorted(norm_reprMat, normtargVec, mask & targmsk, S(iTr), figh)
    title(compose("PopRepr Evol in %s\n%s",targ_area,S(iTr).meta.explabel),'interpreter','none')
    ylabel("zscore (Evoked - bsl rate)");
    saveallform(S(iTr).meta.figdir,"NormAct_Pattern_TargArea_Evolution_sorted",figh)

    end

end

function popul_evol_plot_sorted(norm_reprMat, targVec, maskidx, S, figh)
    reprBlkMat = cell2mat(cellfun(@(idx)mean(norm_reprMat(:,idx),2),S.stim.gen_idx_seq,'uni',0));
%     reprBlkMat_sem = cell2mat(cellfun(@(idx)sem(norm_reprMat(:,idx),2),S.stim.gen_idx_seq,'uni',0));
    reprBlkMat_nat = cell2mat(cellfun(@(idx)mean(norm_reprMat(:,idx),2),S.stim.nat_idx_seq,'uni',0));
    blockN = max(S.stim.block_arr);
    [~, sortidx] = sort(targVec(maskidx),'descend');
    clrseq = brewermap(blockN-1,'Spectral');%from red to blue/purple-ish
    % [sortedchan, sortidx] = sort(chanArr);%sortidx to make sure the channel follows from 1:64
    figure(figh);clf; hold on;set(figh,'pos',[472   499   895   445])
    for i =1:blockN-1
        reprVecMsk = reprBlkMat(maskidx,i);
        reprVecMsk_nat = reprBlkMat_nat(maskidx,i);
        plot(reprVecMsk(sortidx),'color',[clrseq(i,:),0.5],'LineWidth',1.5) % sortidx
        plot(reprVecMsk_nat(sortidx),':','color',[clrseq(i,:),0.3],'LineWidth',1.5) % sortidx
        %shadedErrorBar(-249:500,psthavg_col{threadi,i},psthsem_col{threadi,i},'lineProps',{'color',[clrseq(i,:),0.7],'LineWidth',1.5})
    end
    targVecMsk = targVec(maskidx);
    plot(targVecMsk(sortidx),'color',[0,0,0,0.75],'LineWidth',2.5);
    xticks(1:sum(maskidx))
    unit_arr_msk = S.units.unit_name_arr(maskidx);
    xticklabels(unit_arr_msk(sortidx));
    xlim([0,sum(maskidx)+1])
    xlabel("Channels");ylabel("Activation");
end

function popul_evol_plot_sorted_perarea(norm_reprMat, targVec, maskidxs, masklabels, S, figh)
    reprBlkMat = cell2mat(cellfun(@(idx)mean(norm_reprMat(:,idx),2),S.stim.gen_idx_seq,'uni',0));
%     reprBlkMat_sem = cell2mat(cellfun(@(idx)sem(norm_reprMat(:,idx),2),S.stim.gen_idx_seq,'uni',0));
    reprBlkMat_nat = cell2mat(cellfun(@(idx)mean(norm_reprMat(:,idx),2),S.stim.nat_idx_seq,'uni',0));
    blockN = max(S.stim.block_arr);
    chan_id_col = {};
    chan_id_cat = [];
    xtick_col = {};
    for mi = 1:numel(maskidxs)
        maskidx = maskidxs{mi};
        chan_ids = find(maskidx);
        [~, sortidx] = sort(targVec(maskidx),'descend');
        chan_id_col{mi} = chan_ids(sortidx);
        xtick_col{mi} = numel(chan_id_cat)+1:numel(chan_id_cat)+numel(chan_ids);
        chan_id_cat = [chan_id_cat; chan_ids(sortidx)];
    end
    clrseq = brewermap(blockN-1,'Spectral');%from red to blue/purple-ish
    % [sortedchan, sortidx] = sort(chanArr);%sortidx to make sure the channel follows from 1:64
    figure(figh);clf; hold on;set(figh,'pos',[472   499   895   445])
    for i =1:blockN-1
%         plot(reprBlkMat(chan_id_cat,i),'color',[clrseq(i,:),0.5],'LineWidth',1.5) % sortidx
%         plot(reprBlkMat_nat(chan_id_cat,i),':','color',[clrseq(i,:),0.3],'LineWidth',1.5) % sortidx
        for mi = 1:numel(maskidxs)
            plot(xtick_col{mi},reprBlkMat(chan_id_col{mi},i),'color',[clrseq(i,:),0.5],'LineWidth',1.5) % sortidx
            plot(xtick_col{mi},reprBlkMat_nat(chan_id_col{mi},i),':','color',[clrseq(i,:),0.3],'LineWidth',1.5) % sortidx
        end
    end
%     plot(targVec(chan_id_cat),'color',[0,0,0,0.75],'LineWidth',2.5);
    for mi = 1:numel(maskidxs)
        plot(xtick_col{mi},targVec(chan_id_col{mi}),'color',[0,0,0,0.75],'LineWidth',2.5);
    end
    xticks(1:length(chan_id_cat))
    xticklabels(S.units.unit_name_arr(chan_id_cat));
    xlim([0,length(chan_id_cat)+1])
    xlabel("Channels");ylabel("Activation");
end