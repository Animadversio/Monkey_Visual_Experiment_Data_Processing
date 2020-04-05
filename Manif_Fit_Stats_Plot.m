%% Manif_Fit_Stats_Plot
%% Plot the summary Scatter and line For SU Hash Comparison
Expi_col = [1,3,4,5,8,9,10,11,12,13,15,16,17,18,19,20,21,22];
KSU_vec = [];
KSU_CI = [];
KHash_vec = [];
KHash_CI = [];
Fit_strs = repmat("",1,numel(Expi_col));
for i = 1:numel(Expi_col)
Expi = Expi_col(i);
SUi = 1; Hai = 2; % index for SU and Ha within each struct
if Expi == 10 && Animal=='Alfa', SUi = 2; Hai =1; end % At Exp 10 SU label swapped....
KSU_vec = [KSU_vec, Kent_stats(Expi).act_fit{SUi}.coef(4)];
KSU_CI = [KSU_CI, Kent_stats(Expi).act_fit{SUi}.confint(:,4)];

KHash_vec = [KHash_vec, Kent_stats(Expi).act_fit{Hai}.coef(4)];
KHash_CI = [KHash_CI, Kent_stats(Expi).act_fit{Hai}.confint(:,4)];
Fit_strs(i) = sprintf("Ch%02d R2\nSU %.3f\nHa %.3f", ...
        Stats(Expi).units.pref_chan,...
        Kent_stats(Expi).act_fit{SUi}.rsquare,...
        Kent_stats(Expi).act_fit{Hai}.rsquare);
end
%%
savepath = "E:\Manif_SUHash\summary";
i_list = 1:numel(Expi_col);set(6,'position', [17         227        1678         511])
figure(6);clf;hold on
errorbar(i_list, KSU_vec, KSU_vec-KSU_CI(1,:), KSU_CI(2,:)-KSU_vec,'o','LineWidth',2,'CapSize',10)
errorbar(i_list, KHash_vec, KHash_vec-KHash_CI(1,:), KHash_CI(2,:)-KHash_vec,'o','LineWidth',2,'CapSize',10)
xlabel("Exp num");ylabel("kappa")
xticks(i_list)
xticklabels(Expi_col)
text(i_list, 3.5 * ones(1,numel(Expi_col)), Fit_strs)
legend({"SU", "Hash"})
axis equal tight
title(["Alfa Manifold Exp Kappa parameter comparison for SU vs Hash", "(0.95 CI on errorbar)"])

figure(7);clf;hold on;set(7,'position',[583   172   621   500])
errorbar(KSU_vec, KHash_vec, ...
        KSU_vec-KSU_CI(1,:), KSU_CI(2,:)-KSU_vec,...
        KHash_vec-KHash_CI(1,:), KHash_CI(2,:)-KHash_vec,'o','CapSize',10)
xlabel("kappa (SU)");ylabel("kappa (Hash)")
line([0,2],[0,2])
axis equal
title(["Alfa Manifold Exp Kappa parameter comparison for SU vs Hash", "(0.95 CI on errorbar)"])
%%
saveas(6, fullfile(savepath, "Alfa_Manif_SUHash_kappa.jpg"))
savefig(6, fullfile(savepath, "Alfa_Manif_SUHash_kappa"),'compact')
saveas(7, fullfile(savepath, "Alfa_Manif_SUHash_kappa_scat.jpg"))
savefig(7, fullfile(savepath, "Alfa_Manif_SUHash_kappa_scat"),'compact')
%% Plot everything together
Kstat = repmat(struct("kappa", 0, "CI", [0,0], "Expi", 0, "unit", 0), 0,0);
Fit_strs = repmat("",1,numel(Stats));
csr = 0;
for Expi = 1:numel(Stats)
 % index for SU and Ha within each struct
for ui = 1:size(Kent_stats(Expi).act_fit, 1)
csr = csr + 1;
Kstat(csr).kappa = Kent_stats(Expi).act_fit{ui,1}.coef(4);% consider subspace PC23 now
Kstat(csr).CI = Kent_stats(Expi).act_fit{ui,1}.confint(:,4);
Kstat(csr).theta = Kent_stats(Expi).act_fit{ui,1}.coef(1);% consider subspace PC23 now
Kstat(csr).CI_th = Kent_stats(Expi).act_fit{ui,1}.confint(:,1);
Kstat(csr).phi = Kent_stats(Expi).act_fit{ui,1}.coef(2);% consider subspace PC23 now
Kstat(csr).CI_ph = Kent_stats(Expi).act_fit{ui,1}.confint(:,2);
Kstat(csr).psi = Kent_stats(Expi).act_fit{ui,1}.coef(3);% consider subspace PC23 now
Kstat(csr).CI_ps = Kent_stats(Expi).act_fit{ui,1}.confint(:,3);
Kstat(csr).beta = Kent_stats(Expi).act_fit{ui,1}.coef(5);% consider subspace PC23 now
Kstat(csr).CI_bt = Kent_stats(Expi).act_fit{ui,1}.confint(:,5);
Kstat(csr).A = Kent_stats(Expi).act_fit{ui,1}.coef(6);% consider subspace PC23 now
Kstat(csr).CI_A = Kent_stats(Expi).act_fit{ui,1}.confint(:,6);
Kstat(csr).unit = ui;
Kstat(csr).Expi = Expi;
Kstat(csr).chan = Stats(Expi).units.pref_chan;
Kstat(csr).chan_id = Stats(Expi).units.pref_chan_id(ui);
Kstat(csr).R2 = Kent_stats(Expi).act_fit{ui,1}.rsquare;
end
Fit_strs(Expi) = sprintf("Ch%02d R2", Stats(Expi).units.pref_chan);
for ui = 1:size(Kent_stats(Expi).act_fit, 1)
    U_str = sprintf("\n %d: %.3f", Kent_stats(Expi).act_fit{ui}.rsquare);
    Fit_strs(Expi) = Fit_strs(i) + U_str;
end
end
%%
KstatTab = struct2table(Kstat);
save(fullfile(savepath, "Alfa_Kstats.mat"), 'Kstat')
writetable(KstatTab, fullfile(savepath, "Alfa_Kstats.csv"))
%%
prefchan_arr = arrayfun(@(S) S.chan, Kstat);
unit_arr = arrayfun(@(S) S.unit, Kstat);
R2_arr = arrayfun(@(S) S.R2, Kstat);
IT_kappa = arrayfun(@(S) S.kappa, Kstat(prefchan_arr<=32 & unit_arr == 1 & R2_arr > 0.5))
V1_kappa = arrayfun(@(S) S.kappa, Kstat(prefchan_arr<=48 & prefchan_arr>=33 & unit_arr == 1 & R2_arr > 0.5))
V4_kappa = arrayfun(@(S) S.kappa, Kstat(prefchan_arr>=49 & unit_arr == 1 & R2_arr > 0.5))
ttest2([IT_kappa, V4_kappa], V1_kappa)