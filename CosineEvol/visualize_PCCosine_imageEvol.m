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
        [maxscore,maxidx] = max(S(iTr).score.offline_vec(gen_idx_seq{geni}));
        imgnm = S(iTr).imageName(gen_idx_seq{geni}(maxidx)); % updated on Mar.8th
        % imgnm = S(iTr).imageName(gen_idx_seq{geni}(1));
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