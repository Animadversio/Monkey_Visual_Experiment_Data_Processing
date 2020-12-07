% Make final figures out of this
Animal="Both";Set_Path;
tabdir = "O:\Manif_Fitting\Kent_summary";
mat_dir = "O:\Mat_Statistics";
sumdir = tabdir;
alfatab = readtable(fullfile(tabdir,"Alfa_Kstats.csv"));
betotab = readtable(fullfile(tabdir,"Beto_Kstats.csv"));
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
%%
alltab = [];
Animal_tab = array2table(repmat("Alfa",size(alfatab,1),1),'VariableNames',{'Animal'});
alltab = [alltab; Animal_tab, alfatab];
Animal_tab = array2table(repmat("Beto",size(betotab,1),1),'VariableNames',{'Animal'});
alltab = [alltab; Animal_tab, betotab];
%%
drivermsk = zeros(size(validmsk));
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(alltab.Animal(i))(alltab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = driver_unit == alltab.unit(i);
end
%%
validmsk = ~((alltab.Animal=="Alfa")&(alltab.Expi==10));
Alfamsk = (alltab.Animal=="Alfa");
Betomsk = (alltab.Animal=="Beto");
V1msk = (alltab.chan<=48 & alltab.chan>=33);
V4msk = (alltab.chan>48);
ITmsk = (alltab.chan<33);
%%
figure;hold on
histogram(alltab.R2(validmsk&Alfamsk),20,'FaceAlpha',0.6)
histogram(alltab.R2(validmsk&Betomsk),20,'FaceAlpha',0.6)

%% Show only the real driver units in the channel 
%  Compare the R2 histogram
hist_plot(alltab, "R2", {validmsk&Alfamsk, validmsk&Betomsk},["Alfa","Beto"],...
    "driver valid","anim_sep","count")
hist_plot(alltab, "R2", {validmsk&V1msk, validmsk&V4msk, validmsk&ITmsk},["V1","V4","IT"],...
    "driver valid","area_sep","count")
%%
hist_plot(alltab, "R2", {validmsk& drivermsk& Alfamsk, validmsk& drivermsk& Betomsk},["Alfa","Beto"],...
    "pure driver valid","drv_anim_sep","count")
hist_plot(alltab, "R2", {validmsk& drivermsk& V1msk, validmsk& drivermsk& V4msk, validmsk& drivermsk& ITmsk},["V1","V4","IT"],...
    "pure driver valid","drv_area_sep","count")
%%
hist_plot(alltab, "R2", {validmsk& drivermsk},["All Driver"],...
    "pure driver valid","drv_all","count")
function h = hist_plot(StatsTab_sum, statname, masks, labels, titstr, savestr, normstr)
global sumdir
h = figure;clf;hold on 
for i = 1:numel(masks)
histogram(StatsTab_sum.(statname)(masks{i}),20,'norm',normstr)
medval = median(StatsTab_sum.(statname)(masks{i}));
vline(medval,"-.",labels(i)+" "+num2str(medval))
end
% [~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(masks{1}),StatsTab_sum.D1E_Dpr(masks{2}));
Ns = cellfun(@sum,masks);
ylabel(normstr)
xlabel(statname) % "d'( Dirichlet Energy vs Shuffling Ctrl )"
title(compose("Comparison of %s for\n %s channels in %s Experiments",statname,join(labels),titstr))
        %"\nt=%.2f (p=%.1e (%d))",statname,join(labels),titstr,TSTAT.tstat,P,TSTAT.df))
legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.png", statname, savestr, normstr)))
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.pdf", statname, savestr, normstr)))
savefig(h,fullfile(sumdir,compose("%s_%s_cmp_%s.fig", statname, savestr, normstr)))
end

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
