%% Manif_Map_Stat_Pop_Synopsis
%  Final script to generate the statistics for the figures in paper. 
%  Useful to replicate stats, check criterion and so. 
global figdir
Animal="Both"; Set_Path;
%% Population Kent fitting data 
poptabdir = "O:\Manif_Fitting\popstats";
figdir = "O:\Manif_Fitting\summary";
alfafittab = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
betofittab = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"),'Format','auto');
FitTab = [alfafittab;betofittab];
FitTab.area = arrayfun(@area_map, FitTab.chan);
%% Create masks for the population. 
valmsk = (FitTab.unitnum>0) & ~((FitTab.Animal=="Alfa") & (FitTab.Expi==10));
Fmsk = FitTab.F_P<0.001;
V1msk = FitTab.area == "V1";
V4msk = FitTab.area == "V4";
ITmsk = FitTab.area == "IT";
Alfamsk = FitTab.Animal == "Alfa";
Betomsk = FitTab.Animal == "Beto";
PC12msk = FitTab.space==1;
PC49msk = FitTab.space==2;
PCRNmsk = FitTab.space==3;
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
drivermsk = zeros(size(FitTab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(FitTab.Animal{i})(FitTab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (FitTab.unitnum(i) == driver_unit) & (FitTab.chan(i) == FitTab.prefchan(i));
end
%% 
violin_plot_masks(FitTab, "R2", {V1msk&Fmsk&valmsk,V4msk&Fmsk&valmsk,ITmsk&Fmsk&valmsk}, ["V1","V4","IT"],...
     {drivermsk&Fmsk, ~drivermsk&Fmsk}, ["Driver","Non-Driver"],'showData',true)
ylim([0, 1]);title("Fitting R2 Across Area (ANOVA P<0.001)")
%%
violin_plot_masks(FitTab, "R2", {Fmsk&valmsk}, ["Fsig"],...
     {drivermsk,~drivermsk}, ["Driver","Non-Driver"],'showData',true)
ylim([0, 1]);title("Fitting R2 Across Area (ANOVA P<0.001)")
saveallform(figdir,"KentFit_R2_drv_sep_cmp")
%%
xscatter_minor_plot(FitTab, "kappa", "beta", {V1msk&Fmsk,V4msk&Fmsk,ITmsk&Fmsk}, ["V1","V4","IT"],...
     {drivermsk&Fmsk, ~drivermsk&Fmsk}, ["Driver","Non-Driver"], "Kent Fitting Shape Param Distr", false, figdir, "MarkerEdgeAlpha", 0.2)
xlim([0,25]);ylim([0,20]);axis equal
saveallform(figdir, "kappa-beta_area_drv_sep_cmp")
%%

h = stripe_plot(FitTab, "kappa", {PC12msk&drivermsk,PC49msk&drivermsk,PCRNmsk&drivermsk}, ["PC12","PC4950","RND12"],...
     "Param Width comparison", "PCspace_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.65);
%% 
varnm = "beta";
PC12pairmsk = PC12msk&drivermsk& (FitTab.Expi<=10 & FitTab.Animal=="Beto");
PC49pairmsk = PC49msk&drivermsk;%&(FitTab.Expi<=10 & FitTab.Animal=="Animal");
PCRNpairmsk = PCRNmsk&drivermsk;%&(FitTab.Expi<=10 & FitTab.Animal=="Animal");
paired_stripe_plot({FitTab.(varnm)(PC12pairmsk),FitTab.(varnm)(PC49pairmsk),FitTab.(varnm)(PCRNpairmsk)}, ...
     ["PC12","PC4950","RND12"], {ITmsk(PC12pairmsk),V4msk(PC12pairmsk)}, ["IT","V4"])%,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non-Parametric Statistic of the Population
nonpardir = "O:\Manif_NonParam\summary";
NonParTab = readtable(fullfile(nonpardir,"Both"+"_Popul_NonParamWidth.csv"),'Format','auto');
figdir = nonpardir;
%% Create masks for the population. (Basically same as above)
Fmsk = NonParTab.F_P<0.001;
V1msk = NonParTab.area == "V1";
V4msk = NonParTab.area == "V4";
ITmsk = NonParTab.area == "IT";
Alfamsk = NonParTab.Animal == "Alfa";
Betomsk = NonParTab.Animal == "Beto";
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
drivermsk = zeros(size(NonParTab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(NonParTab.Animal{i})(NonParTab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (NonParTab.unitnum(i) == driver_unit) & (NonParTab.chan(i) == NonParTab.prefchan(i));
end
%% Stripe Comparison plot Non-param tuning width
% here, the output if to figdir
h = stripe_plot(NonParTab, "normAUS_bsl", {V1msk&Fmsk,V4msk&Fmsk,ITmsk&Fmsk}, ["V1","V4","IT"],...
     "NonParam Width comparison", "area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.65);
%%
% figure;
% violinplot(NonParTab.normAUS_bsl(Fmsk),NonParTab.area(Fmsk),'showData',true,'GroupOrder',{'V1','V4','IT'})
%%
violin_plot_masks(NonParTab, "normAUS_bsl", {V1msk&Fmsk,V4msk&Fmsk,ITmsk&Fmsk}, ["V1","V4","IT"],...
     {drivermsk,~drivermsk}, ["Driver","Non-Driver"],'showData',true)
ylim([0, 2*pi]);title("Non-Parametric Tuning Width Across Area (ANOVA P<0.001)")
saveallform(nonpardir,"normVUS_bsl_cmp_area_drv_sep")
%% 
violin_plot_masks(NonParTab, "normAUS_bsl", {V1msk&Fmsk,V4msk&Fmsk,ITmsk&Fmsk}, ["V1","V4","IT"],...
     {Alfamsk,Betomsk}, ["Alfa","Beto"],'showData',true)
ylim([0, 2*pi]);title("Non-Parametric Tuning Width Across Area,Animal (ANOVA P<0.001)")
saveallform(nonpardir,"normVUS_bsl_cmp_area_anim_sep")

%% Popluation Smoothness Comparison 
tabdir = "O:\Manif_MapSmooth\popstats";
figdir = "O:\Manif_MapSmooth\summary";
% Prepare the tables and the structure.
alfa_SmthTab = readtable(fullfile(tabdir,"Alfa"+"_Exp_all_SmoothStat.csv"),'Format','auto');
beto_SmthTab = readtable(fullfile(tabdir,"Beto"+"_Exp_all_SmoothStat.csv"),'Format','auto');
SmthTab = [alfa_SmthTab; beto_SmthTab];
SmthTab.area = arrayfun(@area_map, SmthTab.chan);
%% Driver has "smoother" tuning maps than non driver. (or shuffling has larger effect on disrupting the landscape of Driver landscape) 
violin_plot_masks(SmthTab, "D1E_Dpr", {Fmsk&valmsk}, ["Fsig"],...
     {drivermsk,~drivermsk}, ["Driver","Non-Driver"],'showData',true)
title("Tuning Smoothness Driver-None-Driver (ANOVA P<0.001)")
saveallform(figdir,"D1E_Dpr_cmp_drv_sep")
%
violin_plot_masks(SmthTab, "TVE_Dpr", {Fmsk&valmsk}, ["Fsig"],...
     {drivermsk,~drivermsk}, ["Driver","Non-Driver"],'showData',true)
title("Tuning Smoothness Driver-None-Driver (ANOVA P<0.001)")
saveallform(figdir,"TVE_Dpr_cmp_drv_sep")
%% 
violin_plot_masks(SmthTab, "TVE_Dpr", {V1msk&Fmsk&valmsk,V4msk&Fmsk&valmsk,ITmsk&Fmsk&valmsk}, ["V1","V4","IT"],...
     {drivermsk,~drivermsk}, ["Driver","Non-Driver"],'showData',true)
title("Tuning Smoothness Across Area,Driver-None (ANOVA P<0.001)")
saveallform(figdir,"TVE_Dpr_cmp_area_drv_sep")
%% 
violin_plot_masks(SmthTab, "D1E_Dpr", {V1msk&Fmsk&valmsk,V4msk&Fmsk&valmsk,ITmsk&Fmsk&valmsk}, ["V1","V4","IT"],...
     {drivermsk,~drivermsk}, ["Driver","Non-Driver"],'showData',true)
title("Tuning Smoothness Across Area,Driver-None (ANOVA P<0.001)")
saveallform(figdir,"D1E_Dpr_cmp_area_drv_sep")
%%

