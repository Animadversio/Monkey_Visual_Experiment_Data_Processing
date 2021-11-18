%% Summary and Additional Plotting script for Manif_Map_SpatialRelation
figdir = "E:\OneDrive - Washington University in St. Louis\Manif_PopStats";
figdir = "E:\OneDrive - Washington University in St. Louis\CortiDistCorr\summary";
Set_Path;
%% Combine struct of two monkeys
CortiDisCmb = [];
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir, Animal+'_CortiDisCorr.mat'), 'CortiDisCorr')
CortiDisCmb = [CortiDisCmb, CortiDisCorr];
end
fprintf("Loaded tuning map similarity matrices from %d experiments\n",numel(CortiDisCmb))
%% Collect into table
CortiDistSumm = [];
CD_vect_Sum = [];
for i=1:numel(CortiDisCmb)
CortiDistSumm(i).Expi = CortiDisCmb(i).Expi;
CortiDistSumm(i).Animal = CortiDisCmb(i).Animal;
CortiDistSumm(i).prefchan = CortiDisCmb(i).units.pref_chan;
for nm = string(fieldnames(CortiDisCmb(i).cc))' % different masks separations
CortiDistSumm(i).(nm) = CortiDisCmb(i).cc.(nm);
end
for nm = string(fieldnames(CortiDisCmb(i).avg))'
CortiDistSumm(i).(nm) = CortiDisCmb(i).avg.(nm);
end
for nm = string(fieldnames(CortiDisCmb(i).avgsph))'
CortiDistSumm(i).("sph_"+nm) = CortiDisCmb(i).avgsph.(nm);
end
for nm = string(fieldnames(CortiDisCmb(i).res))'
CortiDistSumm(i).("res_"+nm) = CortiDisCmb(i).res.(nm);
end
for nm = string(fieldnames(CortiDisCmb(i).sgtr))'
CortiDistSumm(i).("sgtr_"+nm) = CortiDisCmb(i).sgtr.(nm);
end
CD_vects = add_corr_dist_vect(CortiDisCmb(i));
CD_vect_Sum = [CD_vect_Sum;CD_vects];
end
CortTab = struct2table(CortiDistSumm);


%%
dvec_sph_IT_F_all = cell2mat(arrayfun(@(C)C.dvec_sph_IT_F,CD_vect_Sum,'uni',0));
ccvec_sph_IT_F_all = cell2mat(arrayfun(@(C)C.ccvec_sph_IT_F,CD_vect_Sum,'uni',0));
[cval,pval]=corr(ccvec_sph_IT_F_all,dvec_sph_IT_F_all,'type','spearman');
fprintf("Pair Corr ~ Cortical Distance, Correlation %.3d(%.1e) n=%d\n",cval,pval,numel(ccvec_sph_IT_F_all))
DfChmsk = dvec_sph_IT_F_all>0;
[cval,pval]=corr(ccvec_sph_IT_F_all(DfChmsk),dvec_sph_IT_F_all(DfChmsk),'type','spearman');
fprintf("Pair Corr ~ Cortical Distance, Different Channel, Correlation %.3d(%.1e) n=%d\n",cval,pval,sum(DfChmsk))
%%
dvec_sph_V4_F_all = cell2mat(arrayfun(@(C)C.dvec_sph_V4_F,CD_vect_Sum,'uni',0));
dvec_sph_V1_F_all = cell2mat(arrayfun(@(C)C.dvec_sph_V1_F,CD_vect_Sum,'uni',0));
ccvec_sph_V4_F_all = cell2mat(arrayfun(@(C)C.ccvec_sph_V4_F,CD_vect_Sum,'uni',0));
ccvec_sph_V1_F_all = cell2mat(arrayfun(@(C)C.ccvec_sph_V1_F,CD_vect_Sum,'uni',0));
fprintf("V1 array pair corr %.3f (std %.3f n=%d)\n",mean(ccvec_sph_V1_F_all),std(ccvec_sph_V1_F_all),numel(ccvec_sph_V1_F_all))
fprintf("V4 array pair corr %.3f (std %.3f n=%d)\n",mean(ccvec_sph_V4_F_all),std(ccvec_sph_V4_F_all),numel(ccvec_sph_V4_F_all))
fprintf("IT array pair corr %.3f (std %.3f n=%d)\n",mean(ccvec_sph_IT_F_all),std(ccvec_sph_IT_F_all),numel(ccvec_sph_IT_F_all))
[cval,pval]=corr(ccvec_sph_IT_F_all,dvec_sph_IT_F_all,'type','spearman');
fprintf("IT array Pair Corr ~ Cortical Distance, Correlation %.3f(%.1e) n=%d\n",cval,pval,numel(ccvec_sph_IT_F_all))
DfChmsk = dvec_sph_IT_F_all>0;
[cval,pval]=corr(ccvec_sph_IT_F_all(DfChmsk),dvec_sph_IT_F_all(DfChmsk),'type','spearman');
fprintf("IT array Pair Corr ~ Cortical Distance, Different Channel, Correlation %.3f(%.1e) n=%d\n",cval,pval,sum(DfChmsk))
[cval,pval]=corr(ccvec_sph_V4_F_all,dvec_sph_V4_F_all,'type','spearman');
fprintf("V4 array Pair Corr ~ Cortical Distance, Correlation %.3f(%.1e) n=%d\n",cval,pval,numel(ccvec_sph_V4_F_all))
DfChmsk = dvec_sph_V4_F_all>0;
[cval,pval]=corr(ccvec_sph_V4_F_all(DfChmsk),dvec_sph_V4_F_all(DfChmsk),'type','spearman');
fprintf("V4 array Pair Corr ~ Cortical Distance, Different Channel, Correlation %.3f(%.1e) n=%d\n",cval,pval,sum(DfChmsk))
[cval,pval]=corr(ccvec_sph_V1_F_all,dvec_sph_V1_F_all,'type','spearman');
fprintf("V1 array Pair Corr ~ Cortical Distance, Correlation %.3f(%.1e) n=%d\n",cval,pval,numel(ccvec_sph_V1_F_all))
DfChmsk = dvec_sph_V1_F_all>0;
[cval,pval]=corr(ccvec_sph_V1_F_all(DfChmsk),dvec_sph_V1_F_all(DfChmsk),'type','spearman');
fprintf("V1 array Pair Corr ~ Cortical Distance, Different Channel, Correlation %.3f(%.1e) n=%d\n",cval,pval,sum(DfChmsk))

%%
corrIT_P = arrayfun(@(C)C.avg.IT_P, CortiDisCmb);
corrIT = arrayfun(@(C)C.avg.IT, CortiDisCmb);
corrV4_P = arrayfun(@(C)C.avg.V4_P, CortiDisCmb);
corrV4 = arrayfun(@(C)C.avg.V4, CortiDisCmb);
prefchan_arr = arrayfun(@(C)C.units.pref_chan, CortiDisCmb);
%%
corrsphIT_P = arrayfun(@(C)C.avgsph.IT_P, CortiDisCmb);
corrsphIT = arrayfun(@(C)C.avgsph.IT, CortiDisCmb);
corrsphV4_P = arrayfun(@(C)C.avgsph.V4_P, CortiDisCmb);
corrsphV4 = arrayfun(@(C)C.avgsph.V4, CortiDisCmb);

%%
Alfamsk = CortTab.Animal=="Alfa";
Betomsk = CortTab.Animal=="Beto";
ITdrvmsk = CortTab.prefchan < 33;
V4drvmsk = CortTab.prefchan > 48;
V1drvmsk = CortTab.prefchan <= 48 & CortTab.prefchan >=33;
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V1_F",{[],Alfamsk,Betomsk},["All","Alfa","Beto"])
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V4_F",{[],Alfamsk,Betomsk},["All","Alfa","Beto"])
sum_vect_summary(CD_vect_Sum,"ccvec_sph_IT_F",{[],Alfamsk,Betomsk},["All","Alfa","Beto"])
%%
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V1_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["All","ITdriver","V4driver","V1driver"])
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V4_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["All","ITdriver","V4driver","V1driver"])
sum_vect_summary(CD_vect_Sum,"ccvec_sph_IT_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["All","ITdriver","V4driver","V1driver"])
%%
fprintf("Pairwise correlation in successful evolution exps combined\n")
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V1_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["All","ITdriver","V4driver","V1driver"], sucs_msk)
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V4_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["All","ITdriver","V4driver","V1driver"], sucs_msk)
sum_vect_summary(CD_vect_Sum,"ccvec_sph_IT_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["All","ITdriver","V4driver","V1driver"], sucs_msk)

fprintf("Pairwise correlation in successful evolution exps in 2 monkeys\n")
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V1_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["A All","A ITdriver","A V4driver","A V1driver"], sucs_msk&Alfamsk)
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V1_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["B All","B ITdriver","B V4driver","B V1driver"], sucs_msk&Betomsk)
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V4_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["A All","A ITdriver","A V4driver","A V1driver"], sucs_msk&Alfamsk)
sum_vect_summary(CD_vect_Sum,"ccvec_sph_V4_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["B All","B ITdriver","B V4driver","B V1driver"], sucs_msk&Betomsk)
sum_vect_summary(CD_vect_Sum,"ccvec_sph_IT_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["A All","A ITdriver","A V4driver","A V1driver"], sucs_msk&Alfamsk)
sum_vect_summary(CD_vect_Sum,"ccvec_sph_IT_F",{[],ITdrvmsk,V4drvmsk,V1drvmsk},["B All","B ITdriver","B V4driver","B V1driver"], sucs_msk&Betomsk)

%%
figure; scatter(dvec_sph_IT_F_all,ccvec_sph_IT_F_all,'markeredgealpha',0.3)
%%
h=figure; clf;
histogram(ccvec_sph_V1_F_all);hold on%randn(numel(ccvec_sph_V1_F_all),1)*0.1,
histogram(ccvec_sph_V4_F_all)
histogram(ccvec_sph_IT_F_all)
xlabel("Tuning Correlation");ylabel("Pair Count")
saveallform(figdir,"Pair_cc_Comp",h)
% savenm = compose("Pair_cc_Comp");
% savefig(h,fullfile(figdir,savenm+".fig"))
% saveas(h,fullfile(figdir,savenm+".png"))
% saveas(h,fullfile(figdir,savenm+".pdf"))
%%
figure;
histogram(CortTab.V4DfCh(CortTab.prefchan>48),10)
%%
figure;
histogram(CortTab.V1DfCh(CortTab.prefchan>32),10)
%%
figure;
histogram(CortTab.ITDfCh(CortTab.prefchan<33&sucs_msk),10)
%%
figure;
histogram(CortTab.V4DfCh(CortTab.prefchan>48&sucs_msk),10)
%%
tabA = readtable(fullfile(mat_dir, "Alfa"+"_EvolTrajStats.csv"));
tabB = readtable(fullfile(mat_dir, "Beto"+"_EvolTrajStats.csv"));
ExpTab_cmb = cat(1, tabA, tabB); % load hte table for Evol Traject successfulness
sucs_msk = (ExpTab_cmb.t_p_initend<1E-3)&(ExpTab_cmb.t_p_initmax<1E-3);
%%
figure;
histogram(CortTab.IT_F(CortTab.prefchan<33&CortTab.IT_F_df>100))%sucs_msk&
%%
figure;
histogram(CortTab.ITDfCh_F(CortTab.prefchan<33&CortTab.ITDfCh_F_df>100))
%%
figure;
histogram(CortTab.V4_F(CortTab.prefchan>48&CortTab.V4_F_df>100))
%%
figure;
histogram(CortTab.V1_F(CortTab.prefchan<49&CortTab.prefchan>32&CortTab.V1_F_df>100))

%%
% S=corr_vec_U_msk(CortiDisCmb(1))
ccvec_col = arrayfun(@corr_vec_U_msk,CortiDisCmb);
diary(fullfile(figdir, "corr_pool_summary.txt"))
fprintf("\nFor All\n")
fprintf("Tuning map similarity based on spherical correlation of trial avg tuning maps\n"+...
        "for well modulated unit pairs ANOVA p<0.01\n")
V1_F_cc_avgsph_pool = cell2mat(arrayfun(@(C)C.V1_F_cc_avgsph, ccvec_col,'uni',0)');
V4_F_cc_avgsph_pool = cell2mat(arrayfun(@(C)C.V4_F_cc_avgsph, ccvec_col,'uni',0)');
IT_F_cc_avgsph_pool = cell2mat(arrayfun(@(C)C.IT_F_cc_avgsph, ccvec_col,'uni',0)');
summarize_corr_vec(V1_F_cc_avgsph_pool,"V1_F")
summarize_corr_vec(V4_F_cc_avgsph_pool,"V4_F")
summarize_corr_vec(IT_F_cc_avgsph_pool,"IT_F")
%
fprintf("Tuning map similarity based on Single Trial correlation\n"+...
        "for well modulated unit pairs ANOVA p<0.01\n")
V1_F_cc_res_pool = cell2mat(arrayfun(@(C)C.V1_F_cc_sgtr, ccvec_col,'uni',0)');
V4_F_cc_res_pool = cell2mat(arrayfun(@(C)C.V4_F_cc_sgtr, ccvec_col,'uni',0)');
IT_F_cc_res_pool = cell2mat(arrayfun(@(C)C.IT_F_cc_sgtr, ccvec_col,'uni',0)');
summarize_corr_vec(V1_F_cc_res_pool,"V1_F")
summarize_corr_vec(V4_F_cc_res_pool,"V4_F")
summarize_corr_vec(IT_F_cc_res_pool,"IT_F")

fprintf("Noise correlation: Single Trial correlation of Residue response of individual trials\n"+...
        "for well modulated unit pairs ANOVA p<0.01\n")
V1_F_cc_res_pool = cell2mat(arrayfun(@(C)C.V1_F_cc_res, ccvec_col,'uni',0)');
V4_F_cc_res_pool = cell2mat(arrayfun(@(C)C.V4_F_cc_res, ccvec_col,'uni',0)');
IT_F_cc_res_pool = cell2mat(arrayfun(@(C)C.IT_F_cc_res, ccvec_col,'uni',0)');
summarize_corr_vec(V1_F_cc_res_pool,"V1_F")
summarize_corr_vec(V4_F_cc_res_pool,"V4_F")
summarize_corr_vec(IT_F_cc_res_pool,"IT_F")
diary off
winopen(fullfile(figdir, "corr_pool_summary.txt"))
function summarize_corr_vec(corrvec,label,doSEM)
if nargin<=2,doSEM=false;end
if doSEM
fprintf("%s corr: %.3f+-%.3f, N=%d\n",label,mean(corrvec),sem(corrvec),numel(corrvec))
else
fprintf("%s corr: %.3f+-%.3f, N=%d\n",label,mean(corrvec),std(corrvec),numel(corrvec))
end
end
function S=corr_vec_U_msk(CortiDisStat)
areavec = reshape(arrayfun(@area_map,CortiDisStat.units.spikeID),[],1);
V1msk = areavec=="V1";%(CortiDisStat.units.spikeID <=48) & (CortiDisStat.units.spikeID>=33);
V4msk = areavec=="V4";%CortiDisStat.units.spikeID > 48;
ITmsk = areavec=="IT";%CortiDisStat.units.spikeID < 33;
Fmsk = struct2table(CortiDisStat.FStats).F_P < 0.01;
Fmsk = reshape(Fmsk,[],1);
ltrigmat = tril(ones(numel(CortiDisStat.units.spikeID),'logical'),-1);
vecmsks = {Fmsk&V1msk,Fmsk&V4msk,Fmsk&ITmsk};
msklabels = ["V1_F","V4_F","IT_F"];
for mi = 1:numel(vecmsks)
    vecmsk = vecmsks{mi};
    msklab = msklabels(mi);
    matmsks{mi} = (vecmsk&vecmsk')&ltrigmat;
end
S = struct();
for prefx = ["sgtr","res","avg","avgsph"]
    for mi = 1:numel(matmsks)
        S.(msklabels(mi)+"_cc_"+prefx) = CortiDisStat.(prefx+"_corrmat")(matmsks{mi});
    end
end
% (V4msk&V4msk')&ltrigmat
% (ITmsk&ITmsk')&ltrigmat
end

function S = add_corr_dist_vect(CortiDisStat)
% CortiDisStat: one element from the saved `CortiDisCorr` object. 
V1msk = (CortiDisStat.units.spikeID <=48) & (CortiDisStat.units.spikeID>=33);
V4msk = CortiDisStat.units.spikeID > 48;
ITmsk = CortiDisStat.units.spikeID < 33;
Fmsk = struct2table(CortiDisStat.FStats).F_P < 0.01;
S = [];
S = extract_corr_dist(CortiDisStat.avgsph_corrmat, CortiDisStat.cortexDismat, Fmsk, ...
	{[], V1msk, V4msk, ITmsk}, ["sph_all_F", "sph_V1_F", "sph_V4_F", "sph_IT_F"], S);
S = extract_corr_dist(CortiDisStat.avg_corrmat, CortiDisStat.cortexDismat, Fmsk, ...
	{[], V1msk, V4msk, ITmsk}, ["all_F", "V1_F", "V4_F", "IT_F"], S);
end

function S = extract_corr_dist(corrmat, distmat, commonmsk, masks, entries, S)
% Compute correlation of tuning similarity and distance, under certain
% masks.
% entries: suffices to add the name of the entries.
if nargin < 6, S = struct(); end
% if nargin < 5, entries = labels; end
if size(commonmsk,2)==1 || size(commonmsk,1)==1 % common mask is a boolean vector.
    commonmsk = reshape(commonmsk,[],1) & reshape(commonmsk,1,[]) & tril(ones(size(corrmat),'logical'),-1);
elseif isempty(commonmsk)
    commonmsk = tril(ones(size(corrmat),'logical'),-1);%ones(size(corrmat),'logical');
else, % common mask is a boolean matrix; we will intersect it with the lower triangular matrix.
    commonmsk = commonmsk & tril(ones(size(corrmat),'logical'),-1);
end
for mi = 1:numel(masks)
msk = masks{mi};
if size(msk,2)==1 || size(msk,1)==1
finalmsk = reshape(msk,[],1) & commonmsk & reshape(msk,1,[]);
corrvec = reshape(corrmat(finalmsk),[],1);
distvec = reshape(distmat(finalmsk),[],1);
elseif isempty(msk)
corrvec = reshape(corrmat(commonmsk),[],1);
distvec = reshape(distmat(commonmsk),[],1);
else
corrvec = reshape(corrmat(msk&commonmsk),[],1);
distvec = reshape(distmat(msk&commonmsk),[],1);
end
nanmsk = ~isnan(corrvec)&~isnan(distvec);
S.("ccvec_"+entries(mi)) = corrvec(nanmsk);
S.("dvec_"+entries(mi)) = distvec(nanmsk);
end
end

function sum_vect_summary(CD_vect_Sum,vecnm,msks,labels,commonmsk)
% msks: cell array of boolean masks to select experiments to compute stats
%       from 
% commonmsk: boolean masks common to all msks. (exp inclusion)
if nargin < 5, commonmsk = ones(size(CD_vect_Sum),'logical'); end
for mi = 1:numel(msks)
msk = msks{mi};
if isempty(msk), msk=ones(size(CD_vect_Sum),'logical'); end
concat_vec = cell2mat(arrayfun(@(C)C.(vecnm),CD_vect_Sum(msk&commonmsk),'uni',0));
fprintf("%s, %s %.3f (std %.3f n=%d)\n",labels(mi),vecnm,mean(concat_vec),std(concat_vec),numel(concat_vec))
end
end