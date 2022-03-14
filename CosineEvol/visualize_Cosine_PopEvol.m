function visualize_Cosine_PopEvol(S,figh)
for iTr = 1:numel(S)
    if nargin==1, figh=figure('pos',[472   499   895   445]);
    else, figh = figure(figh);end 
    reprMat = S(iTr).resp.evoke_trials - S(iTr).resp.bslmean;
    norm_reprMat = (S(iTr).resp.evoke_trials - S(iTr).resp.bslmean - S(iTr).targ.meanActVec{1}) ./ S(iTr).targ.stdActVec{1};
    mask = S(iTr).targ.maskVec{1};
    targVec = S(iTr).targ.targetActVec{1};
    normtargVec = (S(iTr).targ.targetActVec{1} - S(iTr).targ.meanActVec{1})./ S(iTr).targ.stdActVec{1};
    % [V1msk, V4msk, ITmsk] = get_areamask(S(iTr).units.spikeID, array_layout)
    [targmsk, targ_area] = parse_mode2maskVec(S(iTr).targ.score_mode, ...
                            S(iTr).meta.array_layout, S(iTr).units.spikeID);
    
    popul_evol_plot(reprMat, targVec, mask & targmsk, S(iTr), figh)
    title(compose("PopRepr Evol in %s\n%s",targ_area,S(iTr).meta.explabel),'interpreter','none')
    ylabel("Evoked - bsl rate")
    saveallform(S(iTr).meta.figdir,"RawAct_Pattern_TargArea_Evolution",figh)

    popul_evol_plot(norm_reprMat, normtargVec, mask, S(iTr), figh)
    title(compose("PopRepr Evol\n%s",S(iTr).meta.explabel),'interpreter','none')
    ylabel("zscore (Evoked - bsl rate)");
    saveallform(S(iTr).meta.figdir,"NormAct_Pattern_All_Evolution",figh)
    
    popul_evol_plot(norm_reprMat, normtargVec, mask & targmsk, S(iTr), figh)
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
    % Plot the channel and unit number underneath. 
    channum_msk = S.units.spikeID(maskidx);
    ITV1sep = sum(channum_msk<=32)+0.5;
    V1V4sep = sum(channum_msk<=16)+0.5;
    vline([ITV1sep,V1V4sep],'-.r')
    xticks(1:sum(maskidx))
    xticklabels(S.units.unit_name_arr(maskidx))
    xlim([0,sum(maskidx)+1])
    xlabel("Channels");ylabel("Activation");
end