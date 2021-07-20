
%% Manif_NonParametric_Tests
mat_dir = "O:\Mat_Statistics";
nonpardir = "O:\Manif_NonParam\summary";
% measure the non paramtric tuning width or center of mass!
% parameter grid
[phi_grid, theta_grid] = meshgrid(-90:18:90, -90:18:90); 
XX = cosd(theta_grid).* cosd(phi_grid);
YY = sind(theta_grid) .* cosd(phi_grid);
ZZ = sind(phi_grid);
%  Integration Weight Matrix: 
%  Integrate the area of the patch that each node governs
phi1_grid = max(phi_grid - 9, -90) /180 *pi;
phi2_grid = min(phi_grid + 9,  90) /180 *pi;
theta1_grid = max(theta_grid - 9, -90) /180 *pi;
theta2_grid = min(theta_grid + 9,  90) /180 *pi;
Wgrid = abs(sin(phi2_grid) - sin(phi1_grid)).*(theta2_grid - theta1_grid);
%%
poptabdir = "O:\Manif_Fitting\popstats";
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"));
poptab = [alfatab_pop;betotab_pop];
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
drivermsk = zeros(size(poptab,1),1); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end

%% Big Loop Across 2 monkeys
tic
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')

idxlist = find(drivermsk & poptab.Animal==Animal);
NPStats = [];
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
Expi = Tab.Expi; spi = Tab.space; ui = Tab.unitnum; ci = Tab.chan;
for nm = ["Animal", "Expi", "unitstr", "unitnum", "chan", "prefchan", "space","F","F_P"] % copy stats.
S.(nm) = Tab.(nm);
end
S.area = area_map(Tab.chan);
% if Tab.chan < 33
% S.area = "IT";
% elseif Tab.chan < 49
% S.area = "V1";  
% else
% S.area = "V4"; 
% end

iCh = find((MapVarStats(Expi).units.spikeID==ci & MapVarStats(Expi).units.unit_num_arr==ui));
S.iCh = iCh;
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2), MapVarStats(Expi).manif.act_col{spi}, 'uni',1);
[maxAct, maxIdx] = max(actmap_mean,[],'all','linear');
S.maxAct = maxAct;
bslmat = cell2mat(reshape(MapVarStats(Expi).manif.bsl_col{spi},1,[]));
bslmean = mean(bslmat(iCh,:));
bslstd = std(bslmat(iCh,:));
thresh = 0.9 * maxAct; 
[CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap_mean, thresh, XX, YY, ZZ, Wgrid);
S.CoMtheta_90 = CoMtheta;
S.CoMphi_90 = CoMphi;
S.CoMrho_90 = CoMrho;
S.Mtheta_90 = Mtheta;
S.Mphi_90 = Mphi;
S.Mrho_90 = Mrho;
% Core of the loop, use the Spherical Integration to get the activation maps. 
integ = SphIntegration(actmap_mean,0,theta_grid,phi_grid,Wgrid);
S.AUS = integ; % Area under tuning map
S.normAUS = integ / maxAct; % Normalized area under tuning map.
integ = SphIntegration(actmap_mean,bslmean,theta_grid,phi_grid,Wgrid);
S.AUS_bsl = integ;
S.normAUS_bsl = integ / (maxAct - bslmean);
NPStats = [NPStats, S];
end
%%
toc
NPStatTab = struct2table(NPStats);
writetable(NPStatTab,fullfile(nonpardir,Animal+"_Driver_NonParamWidth.csv"))
end
% Re-Read the tables and collect them in a huge table. 
NPtab = [];
for Animal = ["Alfa","Beto"] % merge the tabs for 2 monks
NPStatTab = readtable(fullfile(nonpardir,Animal+"_Driver_NonParamWidth.csv"));
NPtab = [NPtab; NPStatTab];
end
% NPtab.theta_max = NPtab.theta_max / 180 * pi; % change unit 
% NPtab.phi_max = NPtab.phi_max / 180 * pi;
writetable(NPtab,fullfile(nonpardir,"Both"+"_Driver_NonParamWidth.csv"))


%% Reload the tabs to plot the final statistics! 
nonpardir = "O:\Manif_NonParam\summary";
NPtab = readtable(fullfile(nonpardir,"Both"+"_Driver_NonParamWidth.csv"));
global figdir
figdir = "O:\Manif_NonParam\summary";
msk = NPtab.F_P < 1E-3;
%% Creat masks for stats
validmsk = ~ (NPtab.Animal=="Alfa" & NPtab.Expi==10);
Fmsk = NPtab.F_P < 1E-3;
Alfamsk = (NPtab.Animal=="Alfa");
Betomsk = (NPtab.Animal=="Beto");
V1msk = (NPtab.chan<=48 & NPtab.chan>=33);
V4msk = (NPtab.chan>48);
ITmsk = (NPtab.chan<33);
%%
msk = Fmsk & validmsk;
h = stripe_plot(NPtab, "normAUS_bsl", {V1msk & msk, V4msk & msk, ITmsk & msk}, ["V1", "V4", "IT"], ...
	"All Exp (driver, ANOVA P<1E-3)", "area_sep", {[3,1],[2,1],[3,2]});
h = stripe_plot(NPtab, "normAUS", {V1msk & msk, V4msk & msk, ITmsk & msk}, ["V1", "V4", "IT"], ...
	"All Exp (driver, ANOVA P<1E-3)", "area_sep", {[3,1],[2,1],[3,2]});
h = stripe_minor_plot(NPtab, "normAUS", {V1msk & msk, V4msk & msk, ITmsk & msk}, ["V1", "V4", "IT"], ...
    {Alfamsk & msk, Betomsk & msk}, ["Alfa", "Beto"], ...
	"All Exp (driver, ANOVA P<1E-3)", "area_anim_sep", {[3,1],[2,1],[3,2]});
h = stripe_minor_plot(NPtab, "normAUS_bsl", {V1msk & msk, V4msk & msk, ITmsk & msk}, ["V1", "V4", "IT"], ...
    {Alfamsk & msk, Betomsk & msk}, ["Alfa", "Beto"], ...
	"All Exp (driver, ANOVA P<1E-3)", "area_anim_sep", {[3,1],[2,1],[3,2]});
%%
diary(fullfile(figdir,'WidthProg_Stats.log'))
fprintf("Test areal progression of statistics normAUS_bsl :\n")
test_progression(NPtab, "normAUS_bsl", {V1msk & msk, V4msk & msk, ITmsk & msk});
diary off
%% Functions for mass production figures 
msk = poptab.F_P < 1E-3&drivermsk;
Alfamsk = (poptab.Animal=="Alfa");
Betomsk = (poptab.Animal=="Beto");
V1msk = (poptab.chan<=48 & poptab.chan>=33);
V4msk = (poptab.chan>48);
ITmsk = (poptab.chan<33);
% plotTuneCenterMultMsk(poptab,{msk&Alfamsk, msk&Betomsk},["Alfa","Beto"],"","","Kent fit",'Kent_Fsig_anim_sep')
% plotTuneCenterMultMsk(poptab,{msk&V1msk, msk&V4msk, msk&ITmsk},["V1","V4","IT"],"","","Kent fit",'Kent_Fsig_area_sep')
plotTuneCenterMultMsk(poptab,{msk},["driver"],"","","Kent fit",'Kent_Fsig')

function lm = test_progression(tab,varnm,msks,labels)
var_vec = []; idx_vec = [];
for i = 1:numel(msks)
msk = msks{i};
var_arr = tab.(varnm)(msk);
idx_arr = i*ones(numel(var_arr),1);
var_vec = [var_vec; var_arr];
idx_vec = [idx_vec; idx_arr];
end
[cval,pval] = corr(idx_vec, var_vec, 'Type', 'Spearman');
fprintf("Spearman Correlation %.3f(%.1e) df=%d\n",cval,pval,numel(var_vec)-1)
lm = fitlm(idx_vec, var_vec);
fprintf("Linear Regression t=%.3f(p=%.1e) R2=%.3f\n",...
    lm.Coefficients.tStat(2),lm.Coefficients.pValue(2),lm.Rsquared.Adjusted)
disp(lm)
end

function plotTuneCenterMultMsk(NPtab,msks,labels,prefix,suffix,titlenote,savestr)
global figdir
if nargin==4, titlenote=""; end
titlestr = compose("Non-Parametric Estimate (%s) of Tuning Peak\n",titlenote);
h=figure;set(h,'pos',[1000         253         617         725])
for mi = 1:numel(msks)
msk = msks{mi};
scatter(NPtab.(prefix+"theta"+suffix)(msk), NPtab.(prefix+"phi"+suffix)(msk)); hold on
titlestr = titlestr + compose("%s: ",labels{mi});
[~,P,CI,TST] = ttest2(NPtab.(prefix+"theta"+suffix)(msk),NPtab.(prefix+"phi"+suffix)(msk));
titlestr = titlestr + compose("Paired Theta - Phi: CI=[%.2f-%.2f] t=%.2f(p=%.1e) df=%d\n",CI(1),CI(2),TST.tstat,P,TST.df);
[~,P,CI,TST] = ttest(NPtab.(prefix+"theta"+suffix)(msk));
titlestr = titlestr + compose("Theta distrib: CI=[%.2f-%.2f] t=%.2f(p=%.1e) df=%d\n",CI(1),CI(2),TST.tstat,P,TST.df);
[~,P,CI,TST] = ttest(NPtab.(prefix+"phi"+suffix)(msk));
titlestr = titlestr + compose("Phi distrib: CI=[%.2f-%.2f] t=%.2f(p=%.1e) df=%d\n",CI(1),CI(2),TST.tstat,P,TST.df);
end
xline(0.0,'-.');yline(0.0,'-.');
title(titlestr);legend(labels);
xlabel("PC2(Theta)");ylabel("PC3(Phi)")
axis equal;xlim([-pi/2,pi/2]);ylim([-pi/2,pi/2])
savenm = compose("PeakLoc_%s%s_%s",prefix,suffix,savestr);
saveas(h,fullfile(figdir,savenm+".png"));
saveas(h,fullfile(figdir,savenm+".pdf"));
savefig(h,fullfile(figdir,savenm+".fig"));
end

function plotTuneCenter(NPtab,msk,prefix,suffix,titlenote,savestr)
global figdir
if nargin==4, titlenote=""; savestr=""; end
% suffix = "_P90"; prefix = "CoM";
h=figure;set(h,'pos',[1000         253         617         725])
% scatter(mod(NPtab.CoMtheta_P90 + pi/2,pi)-pi/2, NPtab.CoMphi_P90)
scatter(NPtab.(prefix+"theta"+suffix)(msk), NPtab.(prefix+"phi"+suffix)(msk));
[~,P,CI,TST] = ttest2(NPtab.(prefix+"theta"+suffix)(msk),NPtab.(prefix+"phi"+suffix)(msk));
titlestr = compose("Non-Parametric Estimate (%s) of Tuning Peak\n",titlenote);
titlestr = titlestr + compose("Paired Theta - Phi: CI=[%.2f-%.2f] t=%.2f(p=%.1e) df=%d\n",CI(1),CI(2),TST.tstat,P,TST.df);
[~,P,CI,TST] = ttest(NPtab.(prefix+"theta"+suffix)(msk));
titlestr = titlestr + compose("Theta distrib: CI=[%.2f-%.2f] t=%.2f(p=%.1e) df=%d\n",CI(1),CI(2),TST.tstat,P,TST.df);
[~,P,CI,TST] = ttest(NPtab.(prefix+"phi"+suffix)(msk));
titlestr = titlestr + compose("Phi distrib: CI=[%.2f-%.2f] t=%.2f(p=%.1e) df=%d\n",CI(1),CI(2),TST.tstat,P,TST.df);
title(titlestr)
xline(0.0,'-.');yline(0.0,'-.');
xlabel("PC2(Theta)");ylabel("PC3(Phi)")
axis equal;xlim([-pi/2,pi/2]);ylim([-pi/2,pi/2])
savenm = compose("PeakLoc_%s%s_%s",prefix,suffix,savestr);
saveas(h,fullfile(figdir,savenm+".png"));
saveas(h,fullfile(figdir,savenm+".pdf"));
savefig(h,fullfile(figdir,savenm+".fig"));
end

function [CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap, thresh, XX, YY, ZZ, Weight)
% actmap: is a 2d array, 
% XX, YY, ZZ: are 2d array same shape as actmap, corresponding to it entry by entry. 
% thresh: is a scaler that you want to compute Center of Mass or Mean that activation pass.
% Weight: A 2d array of weights which is useful for 2d integration on
%    hemisphere
msk = actmap > thresh;
Xvec = XX(msk); Yvec = YY(msk); Zvec = ZZ(msk); Wvec = Weight(msk);
Actvec = actmap(msk);
CoMvec = [sum(Xvec .* Actvec .* Wvec), sum(Yvec .* Actvec .* Wvec), sum(Zvec .* Actvec .* Wvec)] / sum(Actvec .* Wvec); 
Mvec = [sum(Xvec.* Wvec), sum(Yvec.* Wvec), sum(Zvec.* Wvec)]/sum(Wvec);
% [TH,PHI,R] = cart2sph(CoMvec(),CoMvec(),CoMvec())
[CoMtheta, CoMphi, CoMrho] = cart2sph(CoMvec(1),CoMvec(2),CoMvec(3));
[Mtheta, Mphi, Mrho] = cart2sph(Mvec(1),Mvec(2),Mvec(3));
end

function [Integ,meanInteg] = SphIntegration(actmap, baseline, theta_grid, phi_grid, Weight)
% actmap: is a 2d array, 
% XX, YY, ZZ: are 2d array same shape as actmap, corresponding to it entry by entry. 
% thresh: is a scaler that you want to compute Center of Mass or Mean that activation pass.
% Weight: A 2d array of weights which is useful for 2d integration on
%    hemisphere
Integ = mean((actmap - baseline) .* Weight,'all');
meanInteg = sum((actmap - baseline) .* Weight,'all') / sum(Weight,'all');
% Xvec = XX(msk); Yvec = YY(msk); Zvec = ZZ(msk); Wvec = Weight(msk);
% Actvec = actmap(msk);
% CoMvec = [sum(Xvec .* Actvec .* Wvec), sum(Yvec .* Actvec .* Wvec), sum(Zvec .* Actvec .* Wvec)] / sum(Actvec .* Wvec); 
% Mvec = [sum(Xvec.* Wvec), sum(Yvec.* Wvec), sum(Zvec.* Wvec)]/sum(Wvec);
end

function h = stripe_plot(tab, statname, masks, labels, titstr, savestr, Tpairs)
if nargin<7, Tpairs = {}; end
global figdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
statscol = cellfun(@(M)tab.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
xjitter = 0.15*randn(numel(statscol{i}),1);
scatter(i+xjitter, statscol{i})
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
saveas(h,fullfile(figdir,compose("%s_%s_stripcmp.png", statname, savestr)))
saveas(h,fullfile(figdir,compose("%s_%s_stripcmp.pdf", statname, savestr)))
savefig(h,fullfile(figdir,compose("%s_%s_stripcmp.fig", statname, savestr)))
end

function h = stripe_minor_plot(tab, statname, masks, labels, minormsks, minorlabels, titstr, savestr, Tpairs)
if nargin<7, Tpairs = {}; end
marker = 'o*x^v';
global figdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
labelcol = []; Ns = [];
for i = 1:numel(masks)
    for j = 1:numel(minormsks)
    M = masks{i}&minormsks{j};
    xjitter = 0.15*randn(sum(M),1);
    scatter(i+xjitter, tab.(statname)(M),'Marker',marker(j))
    labelcol = [labelcol, labels(i)+minorlabels(j)];
    Ns = [Ns, sum(M)];
    end
end
xticks(1:numel(masks));xticklabels(labels)
ylabel(statname)
statscol = cellfun(@(M)tab.(statname)(M), masks,'uni',0);
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e df=%d)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P,FStat.STATS.df);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
title(title_str)
% Ns = cellfun(@sum,masks);
% legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
legend(compose("%s(%d)",labelcol',Ns'),'Location','best') % ["Driver", "Non-Drivers"]
saveas(h,fullfile(figdir,compose("%s_%s_stripcmpmarker.png", statname, savestr)))
saveas(h,fullfile(figdir,compose("%s_%s_stripcmpmarker.pdf", statname, savestr)))
savefig(h,fullfile(figdir,compose("%s_%s_stripcmpmarker.fig", statname, savestr)))
end

function addField()
end
% %% Visualize
% msk = NPtab.F_P<1E-3;
% suffix = "_P90"; prefix = "CoM";
% h=figure;set(h,'pos',[1000         253         617         725])
% % scatter(mod(NPtab.CoMtheta_P90 + pi/2,pi)-pi/2, NPtab.CoMphi_P90)
% scatter(NPtab.(prefix+"theta"+suffix)(msk), NPtab.(prefix+"phi"+suffix)(msk));
% xline(0.0,'-.');yline(0.0,'-.');
% xlabel("PC2(Theta)");ylabel("PC3(Phi)")
% axis equal;xlim([-pi/2,pi/2]);ylim([-pi/2,pi/2])
% [~,P,CI,TST] = ttest2(NPtab.(prefix+"theta"+suffix)(msk),NPtab.(prefix+"phi"+suffix)(msk));
% titlestr = compose("Non-Parametric Estimate of Tuning Peak\n");
% titlestr = titlestr + compose("Paired Theta - Phi: CI=[%.2f-%.2f] t=%.2f(p=%.1e) df=%d\n",CI(1),CI(2),TST.tstat,P,TST.df);
% [~,P,CI,TST] = ttest(NPtab.(prefix+"theta"+suffix)(msk));
% titlestr = titlestr + compose("Theta distrib: CI=[%.2f-%.2f] t=%.2f(p=%.1e) df=%d\n",CI(1),CI(2),TST.tstat,P,TST.df);
% [~,P,CI,TST] = ttest(NPtab.(prefix+"phi"+suffix)(msk));
% titlestr = titlestr + compose("Phi distrib: CI=[%.2f-%.2f] t=%.2f(p=%.1e) df=%d\n",CI(1),CI(2),TST.tstat,P,TST.df);
% disp(titlestr)
% title(titlestr)