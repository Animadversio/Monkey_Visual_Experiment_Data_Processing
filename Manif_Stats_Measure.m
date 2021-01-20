%% This file is written to use the stats extracted from the formatted mat file 
%  and further do Kent Function fitting, analysis and plotting for the key units. 
%  and to export tables and draw figures.
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'))
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
%%
addpath D:\Github\Fit_Spherical_Tuning
addpath e:/Github_Projects/Fit_Spherical_Tuning
ft = fittype( @(theta, phi, psi, kappa, beta, A, x, y) KentFunc(theta, phi, psi, kappa, beta, A, x, y), ...
    'independent', {'x', 'y'},'dependent',{'z'});
%% Fit the Kent Distribution. Collect statistics
tic
Kent_stats = repmat(struct(), 1, length(Stats));
for i = 1:length(Stats)
% there can be multiple units, or multiple subspace, so we have this cell
% array
stats_grid = cell(numel(Stats(i).units.pref_chan_id), Stats(i).manif.subsp_n);
scr_stats_grid = cell(numel(Stats(i).units.pref_chan_id), Stats(i).manif.subsp_n);
subsp_str = ['PC23',"PC4950",'RND12'];
subsp_axis = ["PC2","PC3";"PC49","PC50";"RND1","RND2"];
for subsp_i = 1:Stats(i).manif.subsp_n
imgnm_grid = cellfun(@(idx) unique([Stats(i).imageName(idx)]), Stats(i).manif.idx_grid{subsp_i});
psths = Stats(i).manif.psth{subsp_i};
fprintf(subsp_str{subsp_i}+"\n")
for uniti = 1:1:numel(Stats(i).units.pref_chan_id)
channel_j = Stats(i).units.pref_chan_id(uniti);
chan_label_str = sprintf("%s Exp%d Channel %s", Stats(i).Animal, Stats(i).Expi, ...
            Stats(i).units.unit_name_arr{channel_j});

mean_score_map = cellfun(@(psth) mean(psth(uniti, 51:200, :),[2,3])-mean(psth(uniti, 1:40, :),[2,3]), psths);
[Parameter, gof] = fit_Kent(mean_score_map);
scr_stats_grid{uniti, subsp_i} = gof;
% V = coeffvalues(Parameter);
% CI = confint(Parameter);
% param_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
%                 V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
mean_act_map = cellfun(@(psth) mean(psth(uniti, 51:200, :),[2,3]), psths);
[Parameter, gof] = fit_Kent(mean_act_map);
stats_grid{uniti, subsp_i} = gof;
V = coeffvalues(Parameter);
CI = confint(Parameter);
param_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
                V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
fprintf(chan_label_str+"\n "+param_str+"\n")
end
end
Kent_stats(i).act_fit = stats_grid;
Kent_stats(i).scr_fit = scr_stats_grid;
end
toc % took 6.88s to run through! not intensive
%%
% save('D:\Alfa_Manif_Kent_Fit.mat','Kent_stats')
% save(fullfile(savepath, 'Alfa_Manif_Kent_Fit.mat'),'Kent_stats')
savepath = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Manif_SUHash\summary";
% savepath = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Manif_SUHash\summary";
save(compose("E:\\%s_Manif_Kent_Fit.mat", Animal),'Kent_stats')
save(fullfile(savepath, compose("%s_Manif_Kent_Fit.mat", Animal)),'Kent_stats')
%% Formulate it in a table format
subsp_str = ['PC23',"PC4950",'RND12'];
Kstat = repmat(struct("kappa", 0, "CI", [0,0], "Expi", 0, "unit", 0, "subspace", "PC23"), 0,0);
Fit_strs = repmat("",numel(Stats),3);
csr = 0;
for Expi = 1:numel(Stats)
 % index for SU and Ha within each struct
for si = 1:size(Kent_stats(Expi).act_fit, 2) % space idx
for ui = 1:size(Kent_stats(Expi).act_fit, 1) % unit idx
csr = csr + 1;
Kstat(csr).kappa = Kent_stats(Expi).act_fit{ui,si}.coef(4);% consider subspace PC23 now
Kstat(csr).CI = Kent_stats(Expi).act_fit{ui,si}.confint(:,4);
Kstat(csr).theta = Kent_stats(Expi).act_fit{ui,si}.coef(1);% consider subspace PC23 now
Kstat(csr).CI_th = Kent_stats(Expi).act_fit{ui,si}.confint(:,1);
Kstat(csr).phi = Kent_stats(Expi).act_fit{ui,si}.coef(2);% consider subspace PC23 now
Kstat(csr).CI_ph = Kent_stats(Expi).act_fit{ui,si}.confint(:,2);
Kstat(csr).psi = Kent_stats(Expi).act_fit{ui,si}.coef(3);% consider subspace PC23 now
Kstat(csr).CI_ps = Kent_stats(Expi).act_fit{ui,si}.confint(:,3);
Kstat(csr).beta = Kent_stats(Expi).act_fit{ui,si}.coef(5);% consider subspace PC23 now
Kstat(csr).CI_bt = Kent_stats(Expi).act_fit{ui,si}.confint(:,5);
Kstat(csr).A = Kent_stats(Expi).act_fit{ui,si}.coef(6);% consider subspace PC23 now
Kstat(csr).CI_A = Kent_stats(Expi).act_fit{ui,si}.confint(:,6);
Kstat(csr).unit = ui;
Kstat(csr).Expi = Expi;
Kstat(csr).subspace = subsp_str(si);
Kstat(csr).chan = Stats(Expi).units.pref_chan;
Kstat(csr).chan_id = Stats(Expi).units.pref_chan_id(ui);
Kstat(csr).chan_name = Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui));
Kstat(csr).R2 = Kent_stats(Expi).act_fit{ui,si}.rsquare;
end
Fit_strs(Expi, si) = sprintf("Ch%02d %s R2", Stats(Expi).units.pref_chan, subsp_str(si));
for ui = 1:size(Kent_stats(Expi).act_fit, 1)
    U_str = sprintf("\n %d: %.3f", Kent_stats(Expi).act_fit{ui, si}.rsquare);
    Fit_strs(Expi, si) = Fit_strs(i) + U_str;
end
end
end
%%
KstatTab = struct2table(Kstat); % Table is really good for doing statistics (filtering and ordering)
save(fullfile(savepath, compose("%s_Kstats.mat", Animal)), 'Kstat', 'Fit_strs')
writetable(KstatTab, fullfile(savepath, compose("%s_Kstats.csv", Animal)))
%%
%%
% KstatTab = readtable(fullfile(savepath, compose("%s_Kstats.csv", Animal)));
KstatTab = readtable(fullfile(savepath, compose("%s_Kstats_drive.csv", Animal)));
%% Do some basic statistics and plotting
% prefchan_arr = arrayfun(@(S) S.chan, Kstat);
% unit_arr = arrayfun(@(S) S.unit, Kstat);
% R2_arr = arrayfun(@(S) S.R2, Kstat);
% IT_kappa = arrayfun(@(S) S.kappa, Kstat(prefchan_arr<=32 & unit_arr == 1 & R2_arr > 0.5))
% V1_kappa = arrayfun(@(S) S.kappa, Kstat(prefchan_arr<=48 & prefchan_arr>=33 & unit_arr == 1 & R2_arr > 0.5))
% V4_kappa = arrayfun(@(S) S.kappa, Kstat(prefchan_arr>=49 & unit_arr == 1 & R2_arr > 0.5))
IT_kappa = KstatTab.kappa(KstatTab.chan<=32 & KstatTab.subspace=="PC23" & KstatTab.R2 >0.5)'
V1_kappa = KstatTab.kappa(KstatTab.chan<=48 & KstatTab.chan>=33 & KstatTab.subspace=="PC23" & KstatTab.R2 >0.5)'
V4_kappa = KstatTab.kappa(KstatTab.chan>=49 & KstatTab.subspace=="PC23" & KstatTab.R2 >0.5)'
ttest2([IT_kappa, V4_kappa], V1_kappa)
%%
IT_kappa = KstatTab.kappa(KstatTab.chan<=32 & KstatTab.R2 >0.5)'
V1_kappa = KstatTab.kappa(KstatTab.chan<=48 & KstatTab.chan>=33 & KstatTab.R2 >0.5)'
V4_kappa = KstatTab.kappa(KstatTab.chan>=49 & KstatTab.R2 >0.5)'
ttest2([IT_kappa, V4_kappa], V1_kappa)
%%
SUmsk = contains(KstatTab.chan_name,"A");
if Animal == "Alfa"
    SUmsk(KstatTab.Expi == 10 & KstatTab.unit == 2) = true;
    SUmsk(KstatTab.Expi == 10 & KstatTab.unit == 1) = false;
end
ITmsk = KstatTab.chan<=32 & KstatTab.R2 >0.5;
V1msk = KstatTab.chan<=48 & KstatTab.chan>=33 & KstatTab.R2 >0.5;
V4msk = KstatTab.chan>=49 & KstatTab.R2 >0.5;
figure(2);clf;hold on;set(2,'position',[680   354   439   624])
scatter(1*ones(sum(V1msk),1), KstatTab.kappa(V1msk), (5*(1+SUmsk(V1msk))).^2)
scatter(2*ones(sum(V4msk),1), KstatTab.kappa(V4msk), (5*(1+SUmsk(V4msk))).^2)
scatter(3*ones(sum(ITmsk),1), KstatTab.kappa(ITmsk), (5*(1+SUmsk(ITmsk))).^2)
xlim([0.5,3.5]);xticks([1, 2, 3])
ylabel("kappa",'FontSize',14)
title([Animal, "Manifold Tuning Peakedness Comparison", "Large o-SU, Small o-Hash"],'FontSize',14)
set(gca,'xticklabels',["V1", "V4", "IT"],'fontsize',14)
%%
figure(3);clf;hold on;set(3,'position',[680   354   439   624])
errorbar(1+0.1*randn(sum(V1msk),1), KstatTab.kappa(V1msk), ...
         KstatTab.kappa(V1msk) - cellfun(@(c)c(1),KstatTab.CI(V1msk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(V1msk)) - KstatTab.kappa(V1msk), ...%positive error
         'o','Color',[     0    0.4470    0.7410])
errorbar(2+0.1*randn(sum(V4msk),1), KstatTab.kappa(V4msk), ...
         KstatTab.kappa(V4msk) - cellfun(@(c)c(1),KstatTab.CI(V4msk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(V4msk)) - KstatTab.kappa(V4msk), ...%positive error
         'o','Color',[0.8500    0.3250    0.0980])
errorbar(3+0.1*randn(sum(ITmsk),1), KstatTab.kappa(ITmsk), ...
         KstatTab.kappa(ITmsk) - cellfun(@(c)c(1),KstatTab.CI(ITmsk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(ITmsk)) - KstatTab.kappa(ITmsk), ...%positive error
         'o','Color',[0.9290    0.6940    0.1250])
xlim([0.5,3.5]);xticks([1, 2, 3])
ylabel("kappa",'FontSize',14)
title([Animal, "Manifold Tuning Peakedness Comparison", "Large o-SU, Small o-Hash"],'FontSize',14)
set(gca,'xticklabels',["V1", "V4", "IT"],'fontsize',14)
%xticklabels(["V1", "V4", "IT"],'FontSize',14)
%%
saveas(2, fullfile(savepath, compose("%s_Manif_kappa.png",Animal)))
savefig(2, fullfile(savepath, compose("%s_Manif_kappa.fig",Animal)))
saveas(3, fullfile(savepath, compose("%s_Manif_kappa_werrbar.png",Animal)))
savefig(3, fullfile(savepath, compose("%s_Manif_kappa_werrbar.fig",Animal)))
%% Plot Figures for the driver units only. 
KstatTab2 = struct2table(Kstat_drive); % Table is really good for doing statistics (filtering and ordering)
save(fullfile(savepath, compose("%s_Kstats_drive.mat", Animal)), 'Kstat_drive')
writetable(KstatTab2, fullfile(savepath, compose("%s_Kstats_drive.csv", Animal)))
KstatTab = KstatTab2;

SUmsk = contains(KstatTab.chan_name,"A");
if Animal == "Alfa"
    SUmsk(KstatTab.Expi == 10 & KstatTab.unit == 2) = true;
    SUmsk(KstatTab.Expi == 10 & KstatTab.unit == 1) = false;
end
ITmsk = KstatTab.chan<=32 & KstatTab.R2 >0.5;
V1msk = KstatTab.chan<=48 & KstatTab.chan>=33 & KstatTab.R2 >0.5;
V4msk = KstatTab.chan>=49 & KstatTab.R2 >0.5;
figure(2);clf;hold on;set(2,'position',[680   354   439   624])
scatter(1*ones(sum(V1msk),1), KstatTab.kappa(V1msk), (5*(1+SUmsk(V1msk))).^2)
scatter(2*ones(sum(V4msk),1), KstatTab.kappa(V4msk), (5*(1+SUmsk(V4msk))).^2)
scatter(3*ones(sum(ITmsk),1), KstatTab.kappa(ITmsk), (5*(1+SUmsk(ITmsk))).^2)
xlim([0.5,3.5]);xticks([1, 2, 3])
ylabel("kappa",'FontSize',14)
title([Animal, "Manifold Tuning Peakedness Comparison", "Large o-SU, Small o-Hash","Driver Unit Only"],'FontSize',14)
set(gca,'xticklabels',["V1", "V4", "IT"],'fontsize',14)
%%
figure(3);clf;hold on;set(3,'position',[680   354   439   624])
errorbar(1+0.1*randn(sum(V1msk),1), KstatTab.kappa(V1msk), ...
         KstatTab.kappa(V1msk) - cellfun(@(c)c(1),KstatTab.CI(V1msk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(V1msk)) - KstatTab.kappa(V1msk), ...%positive error
         'o','Color',[     0    0.4470    0.7410])
errorbar(2+0.1*randn(sum(V4msk),1), KstatTab.kappa(V4msk), ...
         KstatTab.kappa(V4msk) - cellfun(@(c)c(1),KstatTab.CI(V4msk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(V4msk)) - KstatTab.kappa(V4msk), ...%positive error
         'o','Color',[0.8500    0.3250    0.0980])
errorbar(3+0.1*randn(sum(ITmsk),1), KstatTab.kappa(ITmsk), ...
         KstatTab.kappa(ITmsk) - cellfun(@(c)c(1),KstatTab.CI(ITmsk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(ITmsk)) - KstatTab.kappa(ITmsk), ...%positive error
         'o','Color',[0.9290    0.6940    0.1250])
xlim([0.5,3.5]);xticks([1, 2, 3])
ylabel("kappa",'FontSize',14)
title([Animal, "Manifold Tuning Peakedness Comparison", "Large o-SU, Small o-Hash","Driver Unit Only"],'FontSize',14)
set(gca,'xticklabels',["V1", "V4", "IT"],'fontsize',14)
%xticklabels(["V1", "V4", "IT"],'FontSize',14)
%%
saveas(2, fullfile(savepath, compose("%s_Manif_kappa_driver.png",Animal)))
savefig(2, fullfile(savepath, compose("%s_Manif_kappa_driver.fig",Animal)))
saveas(3, fullfile(savepath, compose("%s_Manif_kappa_werrbar_driver.png",Animal)))
savefig(3, fullfile(savepath, compose("%s_Manif_kappa_werrbar_driver.fig",Animal)))
%% among all the units 
fprintf("Comparison pooling all experiment and units\n")
[H,P,CI,tstat]=ttest2(KstatTab.kappa(ITmsk), KstatTab.kappa(V1msk));
fprintf("IT~V1: t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
[H,P,CI,tstat]=ttest2(KstatTab.kappa(ITmsk), KstatTab.kappa(V4msk));
fprintf("IT~V4: t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
[H,P,CI,tstat]=ttest2(KstatTab.kappa(V1msk), KstatTab.kappa(V4msk));
fprintf("V4~V1: t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
%% just among Hash units 
fprintf("Comparison pooling all experiment, Hash units only\n")
[H,P,CI,tstat]=ttest2(KstatTab.kappa(ITmsk & ~SUmsk), KstatTab.kappa(V1msk & ~SUmsk));
fprintf("IT~V1: t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
[H,P,CI,tstat]=ttest2(KstatTab.kappa(ITmsk & ~SUmsk), KstatTab.kappa(V4msk & ~SUmsk));
fprintf("IT~V4: t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
[H,P,CI,tstat]=ttest2(KstatTab.kappa(V1msk & ~SUmsk), KstatTab.kappa(V4msk & ~SUmsk));
fprintf("V4~V1: t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
%% just among SU units 
fprintf("Comparison pooling all experiment, Single units only\n")
[H,P,CI,tstat]=ttest2(KstatTab.kappa(ITmsk & SUmsk), KstatTab.kappa(V1msk & SUmsk));
fprintf("IT~V1: t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
[H,P,CI,tstat]=ttest2(KstatTab.kappa(ITmsk & SUmsk), KstatTab.kappa(V4msk & SUmsk));
fprintf("IT~V4: t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
[H,P,CI,tstat]=ttest2(KstatTab.kappa(V1msk & SUmsk), KstatTab.kappa(V4msk & SUmsk));
fprintf("V4~V1: t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
% Comparison pooling all experiment, Hash units only
% IT~V1: t=4.449E+00, p=6.990E-05(39)
% IT~V4: t=1.384E+00, p=1.754E-01(34)
% V4~V1: t=-6.179E+00, p=1.839E-06(25)
% Comparison pooling all experiment and units
% IT~V1: t=4.340E+00, p=6.014E-05(56)
% IT~V4: t=2.548E+00, p=1.365E-02(55)
% V4~V1: t=-6.033E+00, p=1.453E-06(29)
%% SU-HU paired comparison plot (vertically comparison)
SUmsk = contains(KstatTab.chan_name,"A");
HUmsk = contains(KstatTab.chan_name,"B");
if Animal == "Alfa"
    SUmsk(KstatTab.Expi == 10 & KstatTab.unit == 2) = true;
    SUmsk(KstatTab.Expi == 10 & KstatTab.unit == 1) = false;
end
SUmsk = find(SUmsk); HUmsk = find(HUmsk);
i_list = 1:length(SUmsk);
figure(6);clf;hold on
errorbar(i_list, KstatTab.kappa(SUmsk), ...
         KstatTab.kappa(SUmsk) - cellfun(@(c)c(1),KstatTab.CI(SUmsk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(SUmsk)) - KstatTab.kappa(SUmsk), ...%positive error
         'o','Color',[     0    0.4470    0.7410])
errorbar(i_list, KstatTab.kappa(HUmsk), ...
         KstatTab.kappa(HUmsk) - cellfun(@(c)c(1),KstatTab.CI(HUmsk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(HUmsk)) - KstatTab.kappa(HUmsk), ...%positive error
         'o','Color',[0.8500    0.3250    0.0980])
R2_str = arrayfun(@(ch,SU,HA)compose("Ch%02d R2\nSU:%.3f\nHA%.3f",ch,SU,HA),...
    KstatTab.chan(SUmsk),KstatTab.R2(SUmsk),KstatTab.R2(HUmsk)); % form the annotation text by array function
text(i_list, 3.5 * ones(1,length(SUmsk)), R2_str, 'FontSize', 11)
ylim([-1.5,4]);xlim([0.5,max(i_list)+0.5])
xlabel("Exp num");ylabel("kappa (Kentfun)");legend({"SU", "Hash"})
xticks(i_list);xtickangle(-30)
xticklabels(string(arrayfun(@(exp,sp)compose('%d-%s',exp,sp), KstatTab.Expi(SUmsk), KstatTab.subspace(SUmsk))))
set(gca,'fontsize',12)
title([Animal + " Manifold Exp Kappa Parameter Comparison","SU vs Hash, 95% CI on errorbar"])
%%
saveas(6,fullfile(savepath,Animal+"_Manif_SUHash_kappa.png"))
saveas(6,fullfile(savepath,Animal+"_Manif_SUHash_kappa.fig"))
%% Horizontal comparison
jitter = sort(randn(1,length(SUmsk)) * 0.05);
figure(7);clf;hold on;set(7,'position',[712   362   493   616])
errorbar(1 * ones(1,length(SUmsk)) + jitter, KstatTab.kappa(SUmsk), ...
         KstatTab.kappa(SUmsk) - cellfun(@(c)c(1),KstatTab.CI(SUmsk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(SUmsk)) - KstatTab.kappa(SUmsk), ...%positive error
         'o','Color',[     0    0.4470    0.7410])
errorbar(2 * ones(1,length(SUmsk)) + jitter, KstatTab.kappa(HUmsk), ...
         KstatTab.kappa(HUmsk) - cellfun(@(c)c(1),KstatTab.CI(HUmsk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(HUmsk)) - KstatTab.kappa(HUmsk), ...%positive error
         'o','Color',[0.8500    0.3250    0.0980])
plot([1, 2]'*ones(1,length(SUmsk)) + jitter, [KstatTab.kappa(SUmsk),KstatTab.kappa(HUmsk)]','color',[0.5,0.5,0.5])
ylim([-1,4]);xlim([0.8,2.2])
ylabel("kappa (Kentfun)");
xticks([1,2]);xtickangle(-0);xticklabels({"SU", "Hash"})
set(gca,'fontsize',14)
[P,H,STATS] = signrank(KstatTab.kappa(SUmsk),KstatTab.kappa(HUmsk));
fprintf("SU~Hash (paired SignRank): z=%.3f, p=%.3E\n", STATS.zval, P)
title([Animal + " Manifold Exp Kappa Parameter Comparison","SU vs Hash, 95% CI on errorbar",compose("SignRank: p=%.2E",P)])
%%
saveas(7,fullfile(savepath,Animal+"_Manif_SUHash_kappa_horiz.png"))
savefig(7,fullfile(savepath,Animal+"_Manif_SUHash_kappa_horiz.fig"))
%%
[H,P,CI,tstat]=ttest(KstatTab.kappa(SUmsk),KstatTab.kappa(HUmsk));
fprintf("SU~Hash (paired t): t=%.3E, p=%.3E(%d)\n", tstat.tstat, P, tstat.df)
[P,H,STATS] = signrank(KstatTab.kappa(SUmsk),KstatTab.kappa(HUmsk));
fprintf("SU~Hash (paired SignRank): z=%.3f, p=%.3E\n", STATS.zval, P)

%% Plot the comparison and stats for subspaces (N.S) different (Only Beto has this. )
PC23msk = find(contains(KstatTab.subspace,"PC23") & KstatTab.Expi <= 10);
PC49msk = find(contains(KstatTab.subspace,"PC4950") & KstatTab.Expi <= 10);
RND1msk = find(contains(KstatTab.subspace,"RND12") & KstatTab.Expi <= 10);
jitter = sort(randn(1,length(PC23msk)) * 0.05);
figure(8);clf;hold on;set(7,'position',[712   362   493   616])
errorbar(1 * ones(1,length(PC23msk)) + jitter, KstatTab.kappa(PC23msk), ...
         KstatTab.kappa(PC23msk) - cellfun(@(c)c(1),KstatTab.CI(PC23msk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(PC23msk)) - KstatTab.kappa(PC23msk), ...%positive error
         'o')
errorbar(2 * ones(1,length(PC49msk)) + jitter, KstatTab.kappa(PC49msk), ...
         KstatTab.kappa(PC49msk) - cellfun(@(c)c(1),KstatTab.CI(PC49msk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(PC49msk)) - KstatTab.kappa(PC49msk), ...%positive error
         'o')
errorbar(3 * ones(1,length(RND1msk)) + jitter, KstatTab.kappa(RND1msk), ...
         KstatTab.kappa(RND1msk) - cellfun(@(c)c(1),KstatTab.CI(RND1msk)), ...%negative error
         cellfun(@(c)c(2),KstatTab.CI(RND1msk)) - KstatTab.kappa(RND1msk), ...%positive error
         'o')
plot([1, 2, 3]'*ones(1,length(PC23msk)) + jitter, [KstatTab.kappa(PC23msk),KstatTab.kappa(PC49msk),KstatTab.kappa(RND1msk)]','color',[0.5,0.5,0.5])
ylim([-1.5,6]);xlim([0.7,3.3])
ylabel("kappa (Kentfun)");
xticks([1,2,3]);xtickangle(-0);xticklabels({"PC12", "PC4950", "RND12"})
[P,ANOVATAB,Fstat] = anova1([KstatTab.kappa(PC23msk),KstatTab.kappa(PC49msk),KstatTab.kappa(RND1msk)],...
                        ["PC23","PC4950","RND12"],'off');
set(gca,'fontsize',14)
title([Animal + " Manifold Exp Kappa Parameter Comparison","Different subspaces, 95% CI on errorbar",compose("p=%.3f, F=%.2f",P,ANOVATAB{2,5})])
%%
saveas(8,fullfile(savepath,Animal+"_Manif_subspace_kappa_horiz.png"))
savefig(8,fullfile(savepath,Animal+"_Manif_subspace_kappa_horiz.fig"))
%%
[P,ANOVATAB,Fstat] = anova1([KstatTab.kappa(PC23msk),KstatTab.kappa(PC49msk),KstatTab.kappa(RND1msk)],...
                        ["PC23","PC4950","RND12"],'off')
%% Plot the summary Scatter and line 
Expi_col = [1,3,4,5,8,9,10,11,12,13,15,16,17,18,19,20,21,22];
KSU_vec = [];
KSU_CI = [];
KHash_vec = [];
KHash_CI = [];
for i = 1:numel(Expi_col)
Expi = Expi_col(i);
SUi = 1; Hai = 2; % index for SU and Ha within each struct
if Expi == 10 && Animal=='Alfa', SUi = 2; Hai =1; end % At Exp 10 SU label swapped....
KSU_vec = [KSU_vec, Kent_stats(Expi).act_fit{SUi}.coef(4)];
KSU_CI = [KSU_CI, Kent_stats(Expi).act_fit{SUi}.confint(:,4)];

KHash_vec = [KHash_vec, Kent_stats(Expi).act_fit{Hai}.coef(4)];
KHash_CI = [KHash_CI, Kent_stats(Expi).act_fit{Hai}.confint(:,4)];
end
%%
i_list = 1:numel(Expi_col);
figure(6);clf;hold on
errorbar(i_list, KSU_vec, KSU_vec-KSU_CI(1,:), KSU_CI(2,:)-KSU_vec,'o')
errorbar(i_list, KHash_vec, KHash_vec-KHash_CI(1,:), KHash_CI(2,:)-KHash_vec,'o')
axis equal tight
xlabel("Exp num");ylabel("kappa (Kentfun)")
xticks(i_list)
xticklabels(Expi_col)

figure(7);clf;hold on
errorbar(KSU_vec, KHash_vec, ...
        KSU_vec-KSU_CI(1,:), KSU_CI(2,:)-KSU_vec,...
        KHash_vec-KHash_CI(1,:), KHash_CI(2,:)-KHash_vec,'o')
xlabel("kappa (SU)");ylabel("kappa (Hash)")
line([0,2],[0,2])
axis equal
%%
figure(4);clf;hold on
Expi_arr = [29, 28, 31, 30, 32:45];
imgsize = arrayfun(@(c)c.evol.imgsize, EStats(Expi_arr));%28:45
Exp_pairs = reshape(Expi_arr, 2, [])';
idx_pairs = [];labels = {};
for pairi = 1:size(Exp_pairs,1)%28:45
    Expi = Exp_pairs(pairi,:);
    imgsize = arrayfun(@(c)c.evol.imgsize, EStats(Expi));
    [idx,~] = find(KstatTab.Expi==Expi);
    idx_pairs = [idx_pairs; idx'];
    labels = [labels, arrayfun(@(expi,sz)compose("Exp%02d %d deg",expi,sz), Expi, imgsize)];
    errorbar(pairi + [0, 0.5], KstatTab.kappa(idx),...
        KstatTab.kappa(idx) - KstatTab.CI_1(idx),...
        KstatTab.CI_2(idx) - KstatTab.kappa(idx),'o-')
end
R2_str = arrayfun(@(ch,deg1,deg3)compose("Ch%02d R2\ndeg1:%.3f\ndeg3:%.3f",ch,deg1,deg3),...
    KstatTab.chan(idx_pairs(:,1)),KstatTab.R2(idx_pairs(:,1)),KstatTab.R2(idx_pairs(:,2))); % form the annotation text by array function
text(1:size(idx_pairs,1), 0.9 * ones(1,size(idx_pairs,1)), R2_str, 'FontSize', 11)
xlim([0.5,10]);ylim([-0.2,1.1])
xticks(1:0.5:9.5);xticklabels(labels)
[H,P,CI,tstat] = ttest(KstatTab.kappa(idx_pairs(:,1)),KstatTab.kappa(idx_pairs(:,2)))
ylabel("kappa");
title(["Compring peakedness of Manifold Exp",compose("with 1 deg and 3 deg image size t=%.2f(%.3f)",tstat.tstat, P)])
%%
errorbar(kap_SU, kapHash)
%% Plot it with spherical plot 
ang_step = 18;
theta_arr = -90:ang_step:90;
phi_arr   = -90:ang_step:90;
[phi_grid, theta_grid] = meshgrid(theta_arr, phi_arr);
figure(3)
ax1 = subplot(121);
sphere_plot(ax1, theta_grid, phi_grid, mean_act_map);
ylabel("RND1(theta)");zlabel("RND2(phi)");cl = caxis;
title([chan_label_str, "Tuning map on RND 1 2 subspace"])%, stat_str])
ax2 = subplot(122);
Kent_sphere_plot(ax2, coeffvalues(Parameter));
ylabel("RND1(theta)");zlabel("RND2(phi)");caxis(cl);
title([chan_label_str, "Kent Fit Tuning map on RND 1 2 subspace", param_str])
%%

%%
mean_act = nanmean(score_col{uniti}, 3 );
[max_score, max_idx ] = max(mean_act,[],'all','linear');
[r_idx,c_idx] = ind2sub(size(mean_act), max_idx);
ang_PC2 = ang_step * (-5 + r_idx -1);
ang_PC3 = ang_step * (-5 + c_idx -1);
[ang_dist, cos_dist] = dist_grid_on_sphere(ang_PC2, ang_PC3);
subplot(1,numel(pref_chan_id),uniti)
scatter(ang_dist(:), mean_act(:))
[b,bint,~,~,stats] = regress(mean_act(:), [ang_dist(:), ones(numel(cos_dist), 1)]);

