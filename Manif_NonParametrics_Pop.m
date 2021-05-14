%% Manif_NonParametric_Tests
%  Majorly get the tuning center by some ways like center of mass!
%  And estimate the tuning width by sth like "volume" under the surface.
%  Population version 
nonpardir = "O:\Manif_NonParam\summary";
[phi_grid, theta_grid] = meshgrid(-90:18:90, -90:18:90); 
XX = cosd(theta_grid) .* cosd(phi_grid);
YY = sind(theta_grid) .* cosd(phi_grid);
ZZ = sind(phi_grid);
%% Integration Weight Matrix
%  Integrate the area of the patch that each node governs
phi1_grid = max(phi_grid - 9, -90) /180 *pi;
phi2_grid = min(phi_grid + 9,  90) /180 *pi;
theta1_grid = max(theta_grid - 9,  -90) /180 *pi;
theta2_grid = min(theta_grid + 9,   90) /180 *pi;
Wgrid = abs(sin(phi2_grid) - sin(phi1_grid)).*(theta2_grid - theta1_grid);

%% Get population data on alfa and beto
% Set_Path; 
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
poptabdir = "O:\Manif_Fitting\popstats";
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"));
poptab = [alfatab_pop;betotab_pop];
load(fullfile(mat_dir,"Alfa"+"_Evol_stats.mat"),'EStats')
EStats_all.Alfa = EStats;
load(fullfile(mat_dir,"Beto"+"_Evol_stats.mat"),'EStats')
EStats_all.Beto = EStats;
% Masks of real driver units instead of using the first. 
drivermsk = zeros(size(poptab,1),1); 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end

%% Collect statistics across all exps
tic
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
idxlist = find(poptab.Animal==Animal);
NPStats = [];
for i = 1:numel(idxlist) % Loop through all entries (Exp-iCh) for this animal
Tab = poptab(idxlist(i),:);
Expi = Tab.Expi; spi = Tab.space; ui = Tab.unitnum; ci = Tab.chan;
for nm = ["Animal", "Expi", "unitstr", "unitnum", "chan", "prefchan", "space","F","F_P"] % copy stats.
S.(nm) = Tab.(nm);
end
iCh = find((MapVarStats(Expi).units.spikeID==ci & MapVarStats(Expi).units.unit_num_arr==ui));
S.iCh = iCh;
S.driver = drivermsk(idxlist(i));
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2), MapVarStats(Expi).manif.act_col{spi}, 'uni',1);
[maxAct, maxIdx] = max(actmap_mean,[],'all','linear');
[ri,ci] = ind2sub(size(actmap_mean),maxIdx);
S.phi_max = phi_grid(ri,ci);
S.theta_max = theta_grid(ri,ci);
S.maxAct = maxAct;
bslmat = cell2mat(reshape(MapVarStats(Expi).manif.bsl_col{spi},1,[]));
bslmean = mean(bslmat(iCh,:));
bslstd = std(bslmat(iCh,:));
thresh = 0.9 * maxAct; 
S.bslmean = bslmean;
S.bslstd = bslstd;
S.CoMthresh = thresh;
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
integ = SphIntegration(actmap_mean,bslmean,theta_grid,phi_grid,Wgrid,bslmean); % activity should be larger than bslmean to be taken into account
S.AUS_bsl = integ; % Area under tuning map - bsl
S.normAUS_bsl = integ / (maxAct - bslmean); % Normalized area under tuning map - bsl.
NPStats = [NPStats, S];
end
toc % Around 21 sec for an animal.
NPStatTab = struct2table(NPStats);
writetable(NPStatTab,fullfile(nonpardir,Animal+"_Pop_NonParamStat.csv"))
end
%
NPtab = [];
for Animal = ["Alfa","Beto"] % merge the tabs for 2 monks
NPStatTab = readtable(fullfile(nonpardir,Animal+"_Pop_NonParamStat.csv"));
NPtab = [NPtab; NPStatTab];
end
NPtab.theta_max = NPtab.theta_max / 180 * pi; % change unit 
NPtab.phi_max = NPtab.phi_max / 180 * pi;
writetable(NPtab,fullfile(nonpardir,"Both"+"_Pop_NonParamStat.csv"))
%%
nonpardir = "O:\Manif_NonParam\summary";
NPtab = readtable(fullfile(nonpardir,"Both"+"_Pop_NonParamStat.csv"));
%%
global figdir
figdir = "O:\Manif_NonParam\summary";
Fmsk = NPtab.F_P < 1E-3;
Alfamsk = (NPtab.Animal=="Alfa");
Betomsk = (NPtab.Animal=="Beto");
V1msk = (NPtab.chan<=48 & NPtab.chan>=33);
V4msk = (NPtab.chan>48);
ITmsk = (NPtab.chan<33);
drvmsk = NPtab.driver;
%%
h = stripe_minor_plot(NPtab, "normAUS_bsl", {drvmsk&Fmsk&V1msk,drvmsk&Fmsk&V4msk,drvmsk&Fmsk&ITmsk}, ["V1","V4","IT"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
					"all chan ANOVA P<1E-3", "drv_area_anim_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.7);
h = stripe_minor_plot(NPtab, "AUS_bsl", {drvmsk&Fmsk&V1msk,drvmsk&Fmsk&V4msk,drvmsk&Fmsk&ITmsk}, ["V1","V4","IT"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
					"all chan ANOVA P<1E-3", "drv_area_anim_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.7);     
%%
title_str = testProgression(NPtab, "normAUS_bsl", {drvmsk&Fmsk&V1msk,drvmsk&Fmsk&V4msk,drvmsk&Fmsk&ITmsk},...
    ["V1","V4","IT"], "area", "");
%%
h = stripe_plot(NPtab, "normAUS_bsl", {Fmsk&drvmsk,Fmsk&~drvmsk}, ["Driver","NonDriver"], ...
                    "all chan ANOVA P<1E-3", "drv_cmp", {[1,2]},'MarkerEdgeAlpha',0.3);
h = stripe_minor_plot(NPtab, "normAUS_bsl", {Fmsk&drvmsk,Fmsk&~drvmsk}, ["Driver","NonDriver"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
					"all chan ANOVA P<1E-3", "drv_cmp_anim_sep", {[1,2]},'MarkerEdgeAlpha',0.3);
%%
h = stripe_plot(NPtab, "AUS_bsl", {Fmsk&drvmsk,Fmsk&~drvmsk}, ["Driver","NonDriver"], ...
                    "all chan ANOVA P<1E-3", "drv_cmp", {[1,2]},'MarkerEdgeAlpha',0.3);
h = stripe_minor_plot(NPtab, "AUS_bsl", {Fmsk&drvmsk,Fmsk&~drvmsk}, ["Driver","NonDriver"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
					"all chan ANOVA P<1E-3", "drv_cmp_anim_sep", {[1,2]},'MarkerEdgeAlpha',0.3);
%%
% h = stripe_plot(NPtab, "normAUS_bsl", {Fmsk&drvmsk,Fmsk&~drvmsk}, ["Driver","NonDriver"], "all chan ANOVA P<1E-3", "drv_cmp", {[1,2]});
h = stripe_minor_plot(NPtab, "normAUS_bsl", {Fmsk&V1msk,Fmsk&V4msk,Fmsk&ITmsk}, ["V1","V4","IT"], {drvmsk, ~drvmsk}, ["Driver", "NonDriver"],...
					"all chan ANOVA P<1E-3", "drv_cmp_area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.2);
h = stripe_minor_plot(NPtab, "AUS_bsl", {Fmsk&V1msk,Fmsk&V4msk,Fmsk&ITmsk}, ["V1","V4","IT"], {drvmsk, ~drvmsk}, ["Driver", "NonDriver"],...
					"all chan ANOVA P<1E-3", "drv_cmp_area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.2);     
%%
h = stripe_minor_plot(NPtab, "normAUS_bsl", {Fmsk&V1msk,Fmsk&V4msk,Fmsk&ITmsk}, ["V1","V4","IT"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
					"all chan ANOVA P<1E-3", "anim_area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.2);
h = stripe_minor_plot(NPtab, "AUS_bsl", {Fmsk&V1msk,Fmsk&V4msk,Fmsk&ITmsk}, ["V1","V4","IT"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
					"all chan ANOVA P<1E-3", "anim_area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.2);     
%%
msk = ~drvmsk & Fmsk;
h = stripe_minor_plot(NPtab, "normAUS_bsl", {msk&V1msk,msk&V4msk,msk&ITmsk}, ["V1","V4","IT"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
					"non driver chan ANOVA P<1E-3", "nondrv_anim_area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.2);
h = stripe_minor_plot(NPtab, "AUS_bsl", {msk&V1msk,msk&V4msk,msk&ITmsk}, ["V1","V4","IT"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
					"non driver chan ANOVA P<1E-3", "nondrv_anim_area_sep", {[1,2],[2,3],[1,3]},'MarkerEdgeAlpha',0.2);     

                
%%
plotTuneCenter(NPtab,Fmsk,"CoM","_P90","CoM 90ptile",'Fsig')
% plotTuneCenter(NPtab,msk,"M","_P90","Mean 90ptile",'Fsig')
plotTuneCenter(NPtab,Fmsk,"CoM","_90","CoM 90% max",'Fsig')
plotTuneCenter(NPtab,Fmsk,"","_max","max",'Fsig')
%%
plotTuneCenterMultMsk(NPtab,{Fmsk&V1msk, Fmsk&V4msk, Fmsk&ITmsk},["V1","V4","IT"],"CoM","_P90","CoM 90ptile",'Fsig_area_sep')
plotTuneCenterMultMsk(NPtab,{Fmsk&V1msk, Fmsk&V4msk, Fmsk&ITmsk},["V1","V4","IT"],"CoM","_90","CoM 90% max",'Fsig_area_sep')
plotTuneCenterMultMsk(NPtab,{Fmsk&Alfamsk, Fmsk&Betomsk},["Alfa","Beto"],"CoM","_P90","CoM 90ptile",'Fsig_anim_sep')
plotTuneCenterMultMsk(NPtab,{Fmsk&Alfamsk, Fmsk&Betomsk},["Alfa","Beto"],"CoM","_90","CoM 90% max",'Fsig_anim_sep')

%% Functions for mass production figures 
% plotTuneCenterMultMsk(poptab,{msk&Alfamsk, msk&Betomsk},["Alfa","Beto"],"","","Kent fit",'Kent_Fsig_anim_sep')
% plotTuneCenterMultMsk(poptab,{msk&V1msk, msk&V4msk, msk&ITmsk},["V1","V4","IT"],"","","Kent fit",'Kent_Fsig_area_sep')
plotTuneCenterMultMsk(poptab,{Fmsk},["driver"],"","","Kent fit",'Kent_Fsig')
%% Computational Utilities
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

function [Integ,meanInteg] = SphIntegration(actmap, baseline, theta_grid, phi_grid, Weight, thresh)
if nargin == 5, thresh = baseline; end
% actmap: is a 2d array, 
% XX, YY, ZZ: are 2d array same shape as actmap, corresponding to it entry by entry. 
% thresh: is a scaler that you want to compute Center of Mass or Mean that activation pass.
% Weight: A 2d array of weights which is useful for 2d integration on
%    hemisphere
mask = actmap > thresh;
Integ = sum((actmap - baseline) .* Weight .* mask,'all');
meanInteg = sum((actmap - baseline) .* Weight .* mask,'all') / sum(Weight.* mask,'all');
% Xvec = XX(msk); Yvec = YY(msk); Zvec = ZZ(msk); Wvec = Weight(msk);
% Actvec = actmap(msk);
% CoMvec = [sum(Xvec .* Actvec .* Wvec), sum(Yvec .* Actvec .* Wvec), sum(Zvec .* Actvec .* Wvec)] / sum(Actvec .* Wvec); 
% Mvec = [sum(Xvec.* Wvec), sum(Yvec.* Wvec), sum(Zvec.* Wvec)]/sum(Wvec);
end
function addField()
end

function title_str=testProgression(tab, statname, masks, labels, sepvarnm, titstr)
title_str = compose("Test progression of %s ~ %s, for %s\n",statname, sepvarnm, titstr);
valuevec = [];
predvec = [];
for i = 1:numel(masks)
    value = tab.(statname)(masks{i});
    valuevec = [valuevec; value];
    predvec = [predvec; ones(numel(value),1)*(i-1)];
    title_str = title_str+compose("%s %.3f+-%.3f (n=%d)\t",labels(i),mean(value),sem(value),numel(value));
end
[cval,pval] = corr(predvec, valuevec, 'Type', 'Spearman');
title_str = title_str+compose("\nSpearman correlation %.3f(%.1e) n=%d",cval,pval,numel(predvec));
lm = fitlm(predvec, valuevec);
[F_p,F_tbl] = anova1(valuevec,predvec,'off');
Fval = F_tbl{2,5}; F_df = F_tbl{4,3};
anovastr = compose("\nANOVA F=%.3f p=%.1e(df=%d,%d)\n",Fval,F_p,F_tbl{2,3},F_tbl{3,3});
lmstr = compose("Linear Regres %s = %.3f + %s * %.3f \n Intercept %.3f+-%.3f, Slope %.3f+-%.3f\n Slope!=0: T=%.1f P=%.1e\n Fit Rsquare=%.3f",...
                statname, lm.Coefficients.Estimate(1), sepvarnm, lm.Coefficients.Estimate(2),...
                lm.Coefficients.Estimate(1), lm.Coefficients.SE(1), ...
                lm.Coefficients.Estimate(2), lm.Coefficients.SE(2), ...
                lm.Coefficients.tStat(2), lm.Coefficients.pValue(2), ...
                lm.Rsquared.Ordinary);
title_str = title_str + anovastr+lmstr;
% disp(lm)
disp(title_str);
end

function h = stripe_plot(tab, statname, masks, labels, titstr, savestr, Tpairs, varargin)
if nargin<7, Tpairs = {}; end
global figdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
statscol = cellfun(@(M)tab.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
mean_i = mean(statscol{i});
sem_i = sem(statscol{i});
N_i = numel(statscol{i});
legstr = compose("%s %.3f(%.3f) (%d)",labels(i),mean_i,sem_i,N_i);
xjitter = 0.15*randn(N_i,1);
scatter(i+xjitter, statscol{i},'DisplayName',legstr,varargin{:});
end
xticks(1:numel(masks)); xticklabels(labels)
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
Ns = cellfun(@sum,masks);
ylabel(statname,'interpreter','none'); title(title_str,'interpreter','none'); % xlabel(statname) 
% legend(compose("%s(%d)",labels',Ns')) % ["Driver", "Non-Drivers"]
legend('Location','best')
saveas(h,fullfile(figdir,compose("%s_%s_stripcmp.png", statname, savestr)))
saveas(h,fullfile(figdir,compose("%s_%s_stripcmp.pdf", statname, savestr)))
savefig(h,fullfile(figdir,compose("%s_%s_stripcmp.fig", statname, savestr)))
end

function h = stripe_minor_plot(tab, statname, masks, labels, minormsks, minorlabels, titstr, savestr, Tpairs, varargin)
if nargin<7, Tpairs = {}; end
marker = 'o*x^v';
global figdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
labelcol = []; Ns = [];
for i = 1:numel(masks)
    for j = 1:numel(minormsks)
    M = masks{i}&minormsks{j};
    statcol = tab.(statname)(M);
    mean_i = mean(statcol);
	sem_i = sem(statcol);
    N_i = sum(M);
    xjitter = 0.15*randn(N_i,1);
	legstr = compose("%s %.3f(%.3f) (%d)",labels(i)+minorlabels(j),mean_i,sem_i,N_i);
    scatter(i+xjitter,statcol,'Marker',marker(j),'DisplayName',legstr,varargin{:})
    labelcol = [labelcol, labels(i)+minorlabels(j)];
    Ns = [Ns, sum(M)];
    end
end
xticks(1:numel(masks));xticklabels(labels)
ylabel(statname,'interpreter','none')

statscol = cellfun(@(M)tab.(statname)(M), masks,'uni',0);
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e df=%d)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P,FStat.STATS.df);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
title(title_str,'interpreter','none')
% Ns = cellfun(@sum,masks);
% legend(compose("%s(%d)",labelcol',Ns'),'Location','best') % ["Driver", "Non-Drivers"]
legend('Location','best')%
saveas(h,fullfile(figdir,compose("%s_%s_stripcmpmarker.png", statname, savestr)))
saveas(h,fullfile(figdir,compose("%s_%s_stripcmpmarker.pdf", statname, savestr)))
savefig(h,fullfile(figdir,compose("%s_%s_stripcmpmarker.fig", statname, savestr)))
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