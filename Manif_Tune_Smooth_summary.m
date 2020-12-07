% This script is the final summary figure maker for population analysis of
% Tuning Smoothness. 
%%
global sumdir
Animal="Alfa"; Set_Path;
mat_dir = "O:\Mat_Statistics";
tabdir = "O:\Manif_MapSmooth\popstats";
sumdir = "O:\Manif_MapSmooth\summary";
mkdir(sumdir)
% load(fullfile(Matdir, "Beto_ManifPopDynamics.mat"))
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
%%
Animal="Alfa"; 
alfa_StatsTab_sum = readtable(fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));
Animal="Beto"; 
beto_StatsTab_sum = readtable(fullfile(tabdir,compose("%s_Exp_all_SmoothStat.csv",Animal)));
StatsTab_sum = [alfa_StatsTab_sum;beto_StatsTab_sum];
%% Useful masks for analysis
drivermsk = (StatsTab_sum.chan==StatsTab_sum.prefchan);
tunemsk = (StatsTab_sum.F_P<1E-6);
V1msk = (StatsTab_sum.chan<=48 & StatsTab_sum.chan>=33);
V4msk = (StatsTab_sum.chan>48);
ITmsk = (StatsTab_sum.chan<33);
Alfamsk = (StatsTab_sum.Animal=="Alfa");
Betomsk = (StatsTab_sum.Animal=="Beto");
%%
[~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk),StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk));
%%
figure(3);clf;hold on 
histogram(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk),'norm','pdf')
histogram(StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk),'norm','pdf')
[~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk),StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk));
N1 = sum(drivermsk & tunemsk);
N2 = sum(~drivermsk & tunemsk);
ylabel("Prob Density")
xlabel("d'( Dirichlet Energy vs Shuffling Ctrl )")
title(compose("Comparison of Gap of Dirichlet Energy for Driver\nand Non-Driver tuned channels in all Experiments\n"+...
        "t=%.2f (p=%.1e (%d))",TSTAT.tstat,P,TSTAT.df))
legend(["Driver", "Non-Drivers"])
saveas(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp_pdf.png"))
saveas(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp_pdf.pdf"))
savefig(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp_pdf.fig"))
%%
figure(3);hold on 
histogram(StatsTab_sum.D1E_Dpr(drivermsk & tunemsk),'norm','count')
histogram(StatsTab_sum.D1E_Dpr(~drivermsk & tunemsk),'norm','count')
ylabel("Prob Density")
xlabel("d'( Dirichlet Energy vs Shuffling Ctrl )")
title(compose("Comparison of Gap of Dirichlet Energy for Driver\nand Non-Driver channels in all Experiments"))
legend(["Driver", "Non-Drivers"])
saveas(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp.png"))
saveas(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp.pdf"))
savefig(3,fullfile(sumdir,"D1E_dpr_driver_non_cmp.fig"))
%% 
figure;
scatter(StatsTab_sum.chan(Betomsk & tunemsk), StatsTab_sum.D1E_Dpr(Betomsk & tunemsk))
xlim([0,65])
%%
hist_cmp_plot(StatsTab_sum, {drivermsk & tunemsk, ~drivermsk & tunemsk},...
        ["Driver", "Non-Drivers"],"All Tune","driver_non","pdf");
    %%
hist_cmp_plot(StatsTab_sum, {Alfamsk&drivermsk & tunemsk, Alfamsk&~drivermsk & tunemsk},...
        ["Driver", "Non-Drivers"],"Alfa Tune","driver_non_alfa","pdf");
hist_cmp_plot(StatsTab_sum, {Betomsk&drivermsk & tunemsk, Betomsk&~drivermsk & tunemsk},...
        ["Driver", "Non-Drivers"],"Beto Tune","driver_non_beto","pdf");
function h = hist_cmp_plot(StatsTab_sum, masks, labels, titstr, namestr, normstr)
global sumdir
h = figure;clf;hold on 
for i = 1:numel(masks)
histogram(StatsTab_sum.D1E_Dpr(masks{i}),20,'norm',normstr)
end
[~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(masks{1}),StatsTab_sum.D1E_Dpr(masks{2}));
Ns = cellfun(@sum,masks);
% N2 = sum(masks{2});
ylabel(normstr)
xlabel("d'( Dirichlet Energy vs Shuffling Ctrl )")
title(compose("Comparison of Gap of Dirichlet Energy for\n %s channels in %s Experiments\n"+...
        "t=%.2f (p=%.1e (%d))",join(labels),titstr,TSTAT.tstat,P,TSTAT.df))
legend(compose("%s(%d)",labels',Ns'))%["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("D1E_dpr_%s_cmp_%s.png", namestr, normstr)))
saveas(h,fullfile(sumdir,compose("D1E_dpr_%s_cmp_%s.pdf", namestr, normstr)))
savefig(h,fullfile(sumdir,compose("D1E_dpr_%s_cmp_%s.fig", namestr, normstr)))
end