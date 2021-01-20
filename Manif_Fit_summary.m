% Make final figures out of the Kent fittings
Animal="Both";Set_Path;
tabdir = "O:\Manif_Fitting\Kent_summary";
poptabdir = "O:\Manif_Fitting\popstats";
mat_dir = "O:\Mat_Statistics";
global sumdir
sumdir = "O:\Manif_Fitting\summary";
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
alfatab = readtable(fullfile(tabdir,"Alfa_Kstats.csv"));
betotab = readtable(fullfile(tabdir,"Beto_Kstats.csv"));
preftab = [];
Animal_tab = array2table(repmat("Alfa",size(alfatab,1),1),'VariableNames',{'Animal'});
preftab = [preftab; Animal_tab, alfatab];
Animal_tab = array2table(repmat("Beto",size(betotab,1),1),'VariableNames',{'Animal'});
preftab = [preftab; Animal_tab, betotab];
%% 
drivermsk = zeros(size(validmsk)); % Masks of real driver units
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(alltab.Animal(i))(alltab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = driver_unit == alltab.unit(i);
end
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
hist_plot(alltab, "R2", {validmsk& drivermsk},["All Driver"],...
    "pure driver valid","drv_all","count")

%% Same thing for all the channels (Popu)
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat.csv"));
poptab = [alfatab_pop;betotab_pop];
%%
validmsk = (poptab.unitnum>0) & ~((poptab.Animal=="Alfa") & (poptab.Expi==10));
Alfamsk = (poptab.Animal=="Alfa");
Betomsk = (poptab.Animal=="Beto");
V1msk = (poptab.chan<=48 & poptab.chan>=33);
V4msk = (poptab.chan>48);
ITmsk = (poptab.chan<33);
drivermsk = zeros(size(poptab,1),1); % Masks of real driver units
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end
prefchmsk = poptab.chan==poptab.prefchan;
%% Hist comparison of Kappa value
msk = validmsk & drivermsk & poptab.R2 > 0.4;
h = hist_plot(poptab, "kappa", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "all", "count");
h = hist_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "area_sep", "count");
h = hist_plot(poptab, "kappa", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "anim_sep", "count");
        %%
msk = validmsk & drivermsk & poptab.R2 > 0.4;
h = stripe_plot(poptab, "kappa", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "all");
h = stripe_plot(poptab, "kappa", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "anim_sep");
h = stripe_plot(poptab, "kappa", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "area_sep",{[3,1],[2,1],[3,2]});
        %%
msk = validmsk & drivermsk & poptab.R2 > 0.4;
h = stripe_plot(poptab, "beta", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "all");
h = stripe_plot(poptab, "beta", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "anim_sep");
h = stripe_plot(poptab, "beta", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "area_sep",{[3,1],[2,1],[3,2]});
%%
nondrvmsk = validmsk & ~prefchmsk & poptab.R2 > 0.4 & poptab.F_P < 1E-5;
h = stripe_plot(poptab, "kappa", {nondrvmsk}, ["all driver"], ...
            "All Exp (non pref chan, R2>0.4, ANOVA P<1E-5)", "nonpref");
h = stripe_plot(poptab, "kappa", {nondrvmsk&Alfamsk, nondrvmsk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (non pref chan, R2>0.4, ANOVA P<1E-5)", "nonpref_anim_sep");
h = stripe_plot(poptab, "kappa", {nondrvmsk&V1msk, nondrvmsk&V4msk, nondrvmsk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (non pref chan, R2>0.4, ANOVA P<1E-5)", "nonpref_area_sep",{[3,1],[2,1],[3,2]});
%%
nondrvmsk = validmsk & ~prefchmsk & poptab.F_P < 1E-5;
h = hist_plot(poptab, "R2", {nondrvmsk}, ["all driver"], ...
            "All Exp (non pref chan, ANOVA P<1E-5)", "nonpref", "count");
h = hist_plot(poptab, "R2", {nondrvmsk&Alfamsk, nondrvmsk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (non pref chan, ANOVA P<1E-5)", "nonpref_anim_sep", "count");
h = hist_plot(poptab, "R2", {nondrvmsk&V1msk, nondrvmsk&V4msk, nondrvmsk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (non pref chan, ANOVA P<1E-5)", "nonpref_area_sep", "count",{[3,1],[2,1],[3,2]});
%% Hist comparison of Beta value
h = hist_plot(poptab, "beta", {msk}, ["all driver"], ...
            "All Exp (driver, R2>0.4)", "all", "count");
h = hist_plot(poptab, "beta", {msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"], ...
            "All Exp (driver, R2>0.4)", "area_sep", "count");
h = hist_plot(poptab, "beta", {msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"], ...
            "All Exp (driver, R2>0.4)", "anim_sep", "count");
%% Scatter of the center location of Kent.
msk = validmsk & drivermsk & poptab.R2 > 0.4;
xscatter_plot(poptab,"phi","theta",{msk},["Driver (R2>0.4)"],"All Exp (driver, R2>0.4)", "all")
xscatter_plot(poptab,"phi","theta",{msk&V1msk, msk&V4msk, msk&ITmsk}, ["V1","V4","IT"],...
    "All Exp (driver, R2>0.4)", "area_sep")
xscatter_plot(poptab,"phi","theta",{msk&Alfamsk, msk&Betomsk}, ["Alfa","Beto"],...
    "All Exp (driver, R2>0.4)", "anim_sep")

function h = xscatter_plot(StatsTab_sum, statname1, statname2, masks, labels, titstr, savestr)
global sumdir
h = figure;clf;hold on;set(h,'pos',[1106         327         560         526]);
% statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
scatter(StatsTab_sum.(statname1)(masks{i}),StatsTab_sum.(statname2)(masks{i}))
% medval = median(StatsTab_sum.(statname)(masks{i}));
% vline(medval,"-.",labels(i)+" "+num2str(medval))
end
if any(strcmp(statname1,["phi","theta"])), xlim([-pi/2,pi/2]); end
if any(strcmp(statname2,["phi","theta"])), ylim([-pi/2,pi/2]); 
    pbaspect([1 1 1]);daspect([1 1 1]);end
% FStat = anova_cells(statscol);
% [~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(masks{1}),StatsTab_sum.D1E_Dpr(masks{2}));
title_str = compose("Scatter of %s - %s for\n %s channels %s",...
    statname1,statname2,join(labels),titstr);
%"\nt=%.2f (p=%.1e (%d))",statname,join(labels),titstr,TSTAT.tstat,P,TSTAT.df))
Ns = cellfun(@sum,masks);
xlabel(statname1);ylabel(statname2) % "d'( Dirichlet Energy vs Shuffling Ctrl )"
title(title_str)
legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_xscat_%s.png", statname1, statname2, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_xscat_%s.pdf", statname1, statname2, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_xscat_%s.fig", statname1, statname2, savestr)))
end

function h = stripe_plot(StatsTab_sum, statname, masks, labels, titstr, savestr, Tpairs)
if nargin<7, Tpairs = {}; end
global sumdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
xjitter = 0.15*randn(numel(statscol{i}),1);
scatter(i+xjitter, statscol{i})
% medval = median(StatsTab_sum.(statname)(masks{i}));
% vline(medval,"-.",labels(i)+" "+num2str(medval))
end
xticks(1:numel(masks));xticklabels(labels)
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
Ns = cellfun(@sum,masks);
ylabel(statname)
% xlabel(statname) % "d'( Dirichlet Energy vs Shuffling Ctrl )"
title(title_str)
legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmp.png", statname, savestr)))
saveas(h,fullfile(sumdir,compose("%s_%s_stripcmp.pdf", statname, savestr)))
savefig(h,fullfile(sumdir,compose("%s_%s_stripcmp.fig", statname, savestr)))
end

function h = hist_plot(StatsTab_sum, statname, masks, labels, titstr, savestr, normstr, Tpairs)
if nargin<8, Tpairs = {}; end
global sumdir
h = figure;clf;hold on 
statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
histogram(StatsTab_sum.(statname)(masks{i}),20,'norm',normstr)
medval = median(StatsTab_sum.(statname)(masks{i}));
vline(medval,"-.",labels(i)+" "+num2str(medval))
end
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
Ns = cellfun(@sum,masks);
ylabel(normstr)
xlabel(statname) % "d'( Dirichlet Energy vs Shuffling Ctrl )"
title(title_str)
legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.png", statname, savestr, normstr)))
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.pdf", statname, savestr, normstr)))
savefig(h,fullfile(sumdir,compose("%s_%s_cmp_%s.fig", statname, savestr, normstr)))
end

function h = hist_cmp_plot(StatsTab_sum, statname, masks, labels, titstr, namestr, normstr)
global sumdir
h = figure;clf;hold on 
for i = 1:numel(masks)
histogram(StatsTab_sum.D1E_Dpr(masks{i}),20,'norm',normstr)
end
[~,P,~,TSTAT]=ttest2(StatsTab_sum.D1E_Dpr(masks{1}),StatsTab_sum.D1E_Dpr(masks{2}));
Ns = cellfun(@sum,masks);
ylabel(normstr)
xlabel(statname)
title(compose("Comparison of %s for\n %s channels in %s Experiments\n"+...
        "t=%.2f (p=%.1e (%d))",statname,join(labels),titstr,TSTAT.tstat,P,TSTAT.df))
legend(compose("%s(%d)",labels',Ns'))%["Driver", "Non-Drivers"]
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.png", statname, namestr, normstr)))
saveas(h,fullfile(sumdir,compose("%s_%s_cmp_%s.pdf", statname, namestr, normstr)))
savefig(h,fullfile(sumdir,compose("%s_%s_cmp_%s.fig", statname, namestr, normstr)))
end
