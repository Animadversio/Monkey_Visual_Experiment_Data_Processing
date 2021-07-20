%% Evol_Cosine_ImMatch.m
%  Measure image distance from target to each evolved image. 
D = torchImDist();
%% 
refcoldir = "N:\Stimuli\2020-CosineEvol\RefCollection";
targmap = get_refimg_map(refcoldir);
%%
batchsize = 50;blur_sigma=1;
% CosImgMatch = repmat(struct(),1,numel(CosStats));
for Expi = 14:numel(CosStats)
	targimg = imread(targmap(CosStats(Expi).evol.targ_imgnm));
	if size(targimg,3)==4, targimg = targimg(:,:,1:3); % RGBA case
    elseif ndims(targimg)==2, targimg = repmat(targimg,1,1,3); end % grayscale case
	targimg_rsz = imresize(targimg, [256,256]);
	tic
    fprintf("Cosine Exp %d\n",Expi)
    stimpath = CosStats(Expi).meta.stimuli; %ExpRecord.stimuli{241}; % Manif Expi 11
    evol_imgnms = CosStats(Expi).imageName(CosStats(Expi).evol.row_gen);
    imds = imageDatastore(stimpath);
    imds.ReadSize = batchsize;imds.ReadFcn = @readResizeImage; % StyleGAN image needs resize
    distvec = []; fnm_vec = [];
    while imds.hasdata()
	    [imgbatch, imginfo] = imds.read();
	    imgbatch = cell2mat(reshape(imgbatch,1,1,1,[]));
	    if strcmp(CosImgMatch(Expi).evol.GANstr,'fc6')
	    imgbatch_blur = imgaussfilt(imgbatch,blur_sigma);% blure out the artifacts
	    else,
	    imgbatch_blur=imgbatch;
    end
    distvec_bat = D.distance(imgbatch_blur, targimg_rsz)';
    distvec = [distvec; distvec_bat]; % column vector
    fnm_vec = [fnm_vec;imginfo.Filename];
    end
    [~,imgnm_vec,sfxs] = arrayfun(@(fn)fileparts(fn),string(fnm_vec));
    CosImgMatch(Expi).fullnm_vec = fnm_vec;
    CosImgMatch(Expi).imgnm_vec = imgnm_vec;
    CosImgMatch(Expi).squ.distvec = distvec;
    fprintf("Finish evol dist %.1fs\t",toc)
    %% Distance to reference images. 
    refimgMap = get_refimg_map(stimpath,true);
    refimds = imageDatastore(cellstr(refimgMap.values));
    refimds.ReadSize = batchsize;refimds.ReadFcn = @readResizeImage;
    while refimds.hasdata()
	    [imgbatch, imginfo] = refimds.read();
	    imgbatch = cell2mat(reshape(imgbatch,1,1,1,[]));
	    imgbatch_blur = imgaussfilt(imgbatch,blur_sigma);% blure out the artifacts
	    refdistvec = D.distance(imgbatch_blur, targimg_rsz)';
    end
    refnm_vec = [imginfo.Filename];
    [~,refimgnm_vec,sfxs] = arrayfun(@(fn)fileparts(fn),string(refnm_vec));
    CosImgMatch(Expi).reffullnm_vec = refnm_vec;
    CosImgMatch(Expi).refimgnm_vec = refimgnm_vec;
    CosImgMatch(Expi).squ.refdistvec = refdistvec;
    fprintf("Finish ref dist %.1fs\t",toc)
    %% Build map from imagnm to the distance to target. 
    imgnm2dist = containers.Map(imgnm_vec,distvec);
    for iThr = 1:CosStats(Expi).evol.thread_num
        for i = 1:numel(refimgnm_vec) % add the mapping for ref image to the collection. 
        imgnm2dist(refimgnm_vec(i)+compose("_thread%03d_nat",iThr-1)) = refdistvec(i);
        end
    end
    % Sort out the distance vector using this map
    fulldistvect = cellfun(@(nm)imgnm2dist(nm),CosStats(Expi).imageName);
    CosImgMatch(Expi).squ.fulldistvec = fulldistvect;
    %%
    metricnm = "squ";
    figure(10);clf;hold on;set(10,'pos',[1000    478   620    500]);
    scatter(CosStats(Expi).evol.block_arr(CosStats(Expi).evol.row_gen),fulldistvect(CosStats(Expi).evol.row_gen),16,'o')
    scatter(CosStats(Expi).evol.block_arr(CosStats(Expi).evol.row_nat),fulldistvect(CosStats(Expi).evol.row_nat),36,'.')
    legend(["gen","nat"]);xlabel("generations");ylabel(compose("ImDist to Target (%s)",metricnm))
    title(CosStats(Expi).evol.explabel,'interp','none')
    saveallform(CosStats(Expi).meta.figdir,compose("ImMatch_%s_OptimCurve.png",metricnm),10)
    %%
    ExpImMatch = CosImgMatch(Expi);
    save(fullfile(CosStats(Expi).meta.figdir, "ExpImMatchStat.mat"), 'ExpImMatch')
    fprintf("Finish vis save %.1fs\n",toc)
    %%
    % [codes_all, img_ids, code_geni] = load_codes_all(stimpath, 1); % each row is an evolved code
	% fprintf("Finish loading codes %.1fs\t",toc)
	% fprintf("(%d gens)\t",size(avg_codes,1))
	% fprintf("Finish visualize %.1fs\t",toc)
	% D = D.select_metric("squeeze");
	% distmat_s = D.distmat_B(avg_imgs);
    % D = D.select_metric("alex");
    % distmat_a = D.distmat_B(avg_imgs);
end
save(fullfile(matdir, Animal+"_CosImMatchStat.mat"), "CosImgMatch")
%%
%% Re-Visualize the evolution of image similarities
for Expi = 1:numel(CosStats)
    metricnm = "squ";
    fulldistvect = CosImgMatch(Expi).squ.fulldistvec;
    figure(10);clf;hold on;set(10,'pos',[1000    478   620    500]);
    scatter(CosStats(Expi).evol.block_arr(CosStats(Expi).evol.row_gen),fulldistvect(CosStats(Expi).evol.row_gen),16,'o')
    scatter(CosStats(Expi).evol.block_arr(CosStats(Expi).evol.row_nat),fulldistvect(CosStats(Expi).evol.row_nat),36,'.')
    legend(["gen","nat"]);xlabel("generations");ylabel(compose("ImDist to Target (%s)",metricnm))
    title(CosStats(Expi).evol.explabel,'interp','none')
    saveallform(CosStats(Expi).meta.figdir,compose("ImMatch_%s_OptimCurve",metricnm),10)
end
%% Summarize statistics, how good is the correlation between the MSE in neural vs image. 
CosImgMatchSummary = []; % repmat(struct(),1,numel(CosStats));
for Expi = 1:numel(CosStats)
row_gen = CosStats(Expi).evol.row_gen;
row_nat = CosStats(Expi).evol.row_nat;
gen_idx_seq = CosStats(Expi).evol.gen_idx_seq;
nat_idx_seq = CosStats(Expi).evol.nat_idx_seq;

[D,S] = ImMatching_stats(CosImgMatch(Expi).squ.fulldistvec,"imdist_squ",gen_idx_seq,nat_idx_seq);
for critname = ["L1", "MSE", "corr", "L1_All", "MSE_All", "corr_All", "MSE_IT", "corr_IT", "MSE_V1V4", "corr_V1V4", "MSE_V4", "corr_V4", "MSE_V1", "corr_V1"]
[cval,pval] = corr(CosStats(Expi).scores.(critname+"_vec"), CosImgMatch(Expi).squ.fulldistvec,'type','spearman');
[cval_gen,pval_gen] = corr(CosStats(Expi).scores.(critname+"_vec")(row_gen), CosImgMatch(Expi).squ.fulldistvec(row_gen),'type','spearman');
[cval_nat,pval_nat] = corr(CosStats(Expi).scores.(critname+"_vec")(row_nat), CosImgMatch(Expi).squ.fulldistvec(row_nat),'type','spearman');
% ,'type','pearson')
S.(critname+"_cc") = cval;
S.(critname+"_cc_P") = pval;
S.(critname+"_gen_cc") = cval_gen;
S.(critname+"_gen_cc_P") = pval_gen;
S.(critname+"_nat_cc") = cval_nat;
S.(critname+"_nat_cc_P") = pval_nat;
end
CosImgMatchSummary = [CosImgMatchSummary, S];
end
save(fullfile(matdir, Animal+"_CosImMatchStat.mat"), "CosImgMatch", "CosImgMatchSummary")
%%
CosImMatchTab = struct2table(CosImgMatchSummary); 
writetable(CosImMatchTab, fullfile(matdir, Animal+"_CosineImMatchSummary.csv"))
writetable(CosImMatchTab, fullfile(saveroot, "summary", Animal+"_CosineImMatchSummary.csv"))
%%
CosSumTab = readtable(fullfile(saveroot, "summary", Animal+"_CosineSummary.csv"));
% create masks for comparison. 
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
% h = stripe_minor_plot(CosImMatchTab, "normAUS_bsl", {Fmsk&drvmsk,Fmsk&~drvmsk}, ["Driver","NonDriver"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
%            "all chan ANOVA P<1E-3", "drv_cmp_anim_sep", {[1,2]}, 'color', 'MarkerEdgeAlpha',0.3);
h = stripe_plot(CosImMatchTab, "imdist_squ_bestscore", {fc6msk,BGmsk,SG2msk}, ["FC6","BigGAN","StyleGAN2"], ...
                   "all sessions", "GANsep_cmp", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',1);
%%
h = stripe_plot(CosImMatchTab, "imdist_squ_gen_23best_T", {fc6msk,BGmsk,SG2msk}, ["FC6","BigGAN","StyleGAN2"], ...
                   "all sessions", "GANsep_cmp", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',1);
%%
% h = stripe_plot(CosImMatchTab, "imdist_squ_gen_23best_T", {fc6msk,BGmsk,SG2msk}, ["FC6","BigGAN","StyleGAN2"], ...
%                    "all sessions", "GANsep_cmp", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',1);
h = stripe_plot(CosImMatchTab, "imdist_squ_gen_23best_T", {MSEmsk,CCmsk,dotmsk,L1msk}, ["MSE","CC","dot","L1"], ...
                    "all sessions", "Objsep_cmp", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_plot(CosImMatchTab, "imdist_squ_gen_23best_T", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_plot(CosImMatchTab, "imdist_squ_gen_23last_T", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_plot(CosImMatchTab, "imdist_squ_gen_last_T", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_plot(CosImMatchTab, "imdist_squ_gen_best_T", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);

%%
h = stripe_plot(CosImMatchTab, "corr_gen_cc", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_plot(CosImMatchTab, "MSE_gen_cc", {MSEmsk,CCmsk}, ["MSE","CC"], ...
                    "all sessions", "Objsep_msecc_cmp", {[1,2]},'MarkerEdgeAlpha',0.9);
%%
h = stripe_minor_plot(CosImMatchTab, "MSE_gen_cc", {MSEmsk,CCmsk}, ["MSE","CC"],...
                   {fc6msk,BGmsk,SG2msk}, ["FC6","BigGAN","StyleGAN2"],...
                   "all sessions", "Obj_GANsep_msecc_cmp", {[1,2]}, 'marker', 'MarkerEdgeAlpha',1);
%%


function data = readResizeImage(filename)
% Turn off warning backtrace before calling imread
onState = warning('off', 'backtrace');
c = onCleanup(@() warning(onState));
data = imread(filename);
if size(data,3)==4, data = data(:,:,1:3); % RGBA case
elseif ndims(data)==2, data = repmat(data,1,1,3); % grayscale case
end 
if ~all(size(data,[1,2])==[256,256])
    data = imresize(data,[256,256]);
end
end

function [D,S] = ImMatching_stats(scorevec,pfx,gen_idx_seq,nat_idx_seq)
    % This will clip the last generation in gen_idx_seq, nat_idx_seq by default. 
    if nargin == 1, prefix="imdist"; end
    if nargin <= 4, S=struct(); D=struct(); end

    scorevec = reshape(scorevec,[],1);% use column vector by default.
    D.(pfx+"_vec") = scorevec; 
    scores_rec = cellfun(@(idx)scorevec(idx),reshape(gen_idx_seq(1:end-1),[],1),'uni',0);
    natscores_rec = cellfun(@(idx)scorevec(idx),reshape(nat_idx_seq(1:end-1),[],1),'uni',0);

    score_avg = cellfun(@mean,scores_rec);
    D.(pfx+"_gen_avg") = cellfun(@mean,scores_rec);
    D.(pfx+"_gen_sem") = cellfun(@sem, scores_rec);
    D.(pfx+"_nat_avg") = cellfun(@mean,natscores_rec);
    D.(pfx+"_nat_sem") = cellfun(@sem, natscores_rec);
    Fgen = anova_cells(scores_rec);
    S.(pfx+"_gen_F") = Fgen.F;
    S.(pfx+"_gen_F_P") = Fgen.F_P;
    S.(pfx+"_gen_F_df") = Fgen.STATS.df;
    % [cval,pval] = corr([1:block_num]',OLscore_avg,'type','spearman');
    % S.(pfx+"_gen_corr") = cval;
    % S.(pfx+"_gen_corr_P") = pval;
    [~,P,~,TST] = ttest2(cell2mat(scores_rec(end-1:end)), cell2mat(scores_rec(1:2)));% maybe I should use 2:3
    S.(pfx+"_gen_last_T") = TST.tstat;
    S.(pfx+"_gen_last_T_P") = P;
    [~,P,~,TST] = ttest2(cell2mat(scores_rec(end-1:end)), cell2mat(scores_rec(2:3)));% maybe I should use 2:3
    S.(pfx+"_gen_23last_T") = TST.tstat;
    S.(pfx+"_gen_23last_T_P") = P;
    [minscore,min_idx] = min(movmean(score_avg,3));
    if min_idx == numel(scores_rec), min_idx = min_idx - 1; end
    [~,P,~,TST] = ttest2(cell2mat(scores_rec(min_idx:min_idx+1)), cell2mat(scores_rec(1:2)));
    S.(pfx+"_gen_best_T") = TST.tstat;
    S.(pfx+"_gen_best_T_P") = P;
    [~,P,~,TST] = ttest2(cell2mat(scores_rec(min_idx:min_idx+1)), cell2mat(scores_rec(2:3)));
    S.(pfx+"_gen_23best_T") = TST.tstat;
    S.(pfx+"_gen_23best_T_P") = P;

    S.(pfx+"_bestscore") = minscore;
    S.(pfx+"_lastscore") = score_avg(end);
    D.(pfx+"_stat") = S;
end