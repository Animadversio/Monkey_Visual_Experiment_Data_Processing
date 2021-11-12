Animal = "Alfa";Set_Path;
ftrrows = find(contains(ExpRecord.expControlFN,["generate_BigGAN_cosine"]));%&contains(ExpRecord.expControlFN,["210218","210408","210413","210415"]));%...
%%
refcoldir = "N:\Stimuli\2020-CosineEvol\RefCollection";
targmap = get_refimg_map(refcoldir);
annot_tab = readtable('CosineImage_annotate.csv','Delimiter',',');
target_annot_map = containers.Map(annot_tab.imgName, annot_tab.Annotation);
%%
global figdir
Animal = "Alfa";Set_Path;
saveroot = "O:\Evol_Cosine";
matdir = "O:\Mat_Statistics";
figdir = fullfile(saveroot, "summary");
%%
load(fullfile(matdir, Animal+"_CosStats.mat"), "CosSummary", "CosStats")
load(fullfile(matdir, Animal+"_CosImMatchStat.mat"), "CosImgMatch", "CosImgMatchSummary")
%%
CosSumTab = readtable(fullfile(matdir, Animal+"_CosineSummary.csv"));
CosImMatchTab = readtable(fullfile(matdir, Animal+"_CosineImMatchSummary.csv"));
% CosSumTab = readtable(fullfile(saveroot, "summary", Animal+"_CosineSummary.csv"));
%%
fc6msk = contains(CosSumTab.GANstr,'fc6');
BGmsk = contains(CosSumTab.GANstr,'BigGAN');
SG2msk = contains(CosSumTab.GANstr,'StyleGAN2');
MSEmsk = contains(CosSumTab.score_mode,"MSE");
dotmsk = contains(CosSumTab.score_mode,"dot");
L1msk = contains(CosSumTab.score_mode,"L1");
CCmsk = contains(CosSumTab.score_mode,"corr");
ITmsk = contains(CosSumTab.targ_area,"IT");
V1msk = contains(CosSumTab.targ_area,"V1");
I4msk = contains(CosSumTab.targ_area,"V4");
%%
h = stripe_plot(CosImMatchTab, "corr_IT_gen_cc", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_plot(CosImMatchTab, "corr_gen_cc", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_plot(CosImMatchTab, "MSE_IT_gen_cc", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_plot(CosImMatchTab, "corr_IT_gen_cc", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);

%%
sucsmsk = CosSumTab.ol_gen_last_T_P<0.01;
[~,P,~,TST] = ttest(CosImMatchTab.corr_IT_gen_cc(sucsmsk));
%%

h = stripe_plot(CosImMatchTab, "corr_IT_gen_cc", {sucsmsk}, ["online success"], ...
        "all sessions (online score incr p<0.01)", "SucsExp", {},'MarkerEdgeAlpha',0.9);


%% Image Distance Evolution curve
Corrsucsmsk = CosSumTab.corr_gen_best_T_P<0.001;
MSEsucsmsk = CosSumTab.MSE_gen_best_T_P<0.001;
%% Single image distance 
blockvec_col = arrayfun(@(Expi)CosStats(Expi).evol.block_arr(CosStats(Expi).evol.row_gen),1:numel(CosStats),'uni',0);
imdistvec_col = arrayfun(@(Expi)CosImgMatch(Expi).squ.fulldistvec(CosStats(Expi).evol.row_gen),1:numel(CosStats),'uni',0);
[imdist_m,imdist_s,blockvec] = sort_scoreblock(cell2mat(blockvec_col(:)), cell2mat(imdistvec_col(:)));
figure;
shadedErrorBar(blockvec, imdist_m, imdist_s)
xlim([0,40])
%% block average image distance for each session. 
[imdist_m_col,imdist_s_col,blockarr_col] = arrayfun(@(Expi)sort_scoreblock(CosStats(Expi).evol.block_arr(CosStats(Expi).evol.row_gen),...
											CosImgMatch(Expi).squ.fulldistvec(CosStats(Expi).evol.row_gen)),[1:numel(CosStats)]','uni',0);
%%
[imdist_avg_m,imdist_avg_s,blockvec_avg] = sort_scoreblock(cell2mat(blockarr_col(Corrsucsmsk)'), cell2mat(imdist_m_col(Corrsucsmsk)'));
% [imdist_avg_m,imdist_avg_s,blockvec_avg] = sort_scoreblock(cell2mat(blockvec_col(MSEsucsmsk)'), cell2mat(imdistvec_col(MSEsucsmsk)'));
figure;
shadedErrorBar(blockvec_avg, imdist_avg_m, imdist_avg_s)
xlim([0,40]);xlabel("Generation");ylabel("ImDist squ")
saveallform(figdir, "CorrSuccessEvol_ImDistCurv")
%% 
corr_gen_avg_col = arrayfun(@(Expi)CosStats(Expi).scores.corr_gen_avg,[1:numel(CosStats)]','uni',0);
imdist_m_nlast_col = cellfun(@(imdst)imdst(1:end-1)',imdist_m_col,'uni',0);
%%
%%
corr_gen_col = arrayfun(@(Expi)CosStats(Expi).scores.corr_IT_vec(CosStats(Expi).evol.row_gen),[1:numel(CosStats)]','uni',0);
imdist_col = arrayfun(@(Expi)CosImgMatch(Expi).squ.fulldistvec(CosStats(Expi).evol.row_gen),[1:numel(CosStats)]','uni',0);
%
figure;
scatter(cell2mat(corr_gen_col(Corrsucsmsk)),cell2mat(imdist_col(Corrsucsmsk)),9,'MarkerEdgeAlpha',0.2)
ylabel("ImDist to Target");xlabel("Correlation with Target vector")
[cval,pval]=corr(cell2mat(corr_gen_col(Corrsucsmsk)),cell2mat(imdist_col(Corrsucsmsk)));
title(compose("Correlation of Image Matching vs Neural Matching\n for Success Exps Corr %.3f (p %.1e)",cval,pval))
saveallform(figdir, "ImDistNeuroMatch_sucs_corr")
%%
figure;
scatter(cell2mat(corr_gen_avg_col(Corrsucsmsk)),cell2mat(imdist_m_nlast_col(Corrsucsmsk)))
ylabel("ImDist to Target");xlabel("Correlation with Target vector")




%% Redo the explabels 
for Expi = 1:numel(CosStats)
targ_annot = target_annot_map(CosStats(Expi).evol.targ_imgnm);
CosStats(Expi).evol.targ_annot = targ_annot;
stimparts = split(CosStats(Expi).meta.stimuli,"\");
explabel = stimparts{end-1};
explabel = explabel + compose(" %s targ: %s (%s)\n%s (%s)",CosStats(Expi).evol.score_mode,CosStats(Expi).evol.targ_imgnm,CosStats(Expi).evol.targ_annot,CosStats(Expi).evol.GANstr,CosStats(Expi).evol.optimstr);
explabel = explabel + compose("   Pos [%.1f,%.1f] Size %.1f deg",CosStats(Expi).evol.imgpos,CosStats(Expi).evol.imgsize);
CosStats(Expi).evol.explabel = explabel;
end
%%

















%%
save(fullfile(matdir, Animal+"_CosStats.mat"), "CosSummary", "CosStats")

function [score_m,score_s,blockvec] = sort_scoreblock(blockarr,scorearr)
% sort an array of scores according to the block array labels. compute the
% mean and std for each block. 
% really useful function to summarize multiple evolution trajectories into
% a mean one. 
blockvec = min(blockarr):max(blockarr);
score_m = [];score_s = [];
for blocki = min(blockarr):max(blockarr)
    score_m(blocki) = mean(scorearr(blockarr==blocki));
    score_s(blocki) = sem(scorearr(blockarr==blocki));
end
end