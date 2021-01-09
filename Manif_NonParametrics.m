%% Manif_NonParametric_Tests
nonpardir = "O:\Manif_NonParam\summary";

[phi_grid, theta_grid] = meshgrid(-90:18:90, -90:18:90); 
XX = cosd(theta_grid).* cosd(phi_grid);
YY = sind(theta_grid) .* cosd(phi_grid);
ZZ = sind(phi_grid);
%% Integration Weight Matrix
phi1_grid = max(phi_grid - 9, -90) /180 *pi;
phi2_grid = min(phi_grid + 9,  90) /180 *pi;
theta1_grid = max(theta_grid - 9,  -90) /180 *pi;
theta2_grid = min(theta_grid + 9,   90) /180 *pi;
Wgrid = abs(sin(phi2_grid) - sin(phi1_grid)).*(theta2_grid - theta1_grid);
%%
poptabdir = "O:\Manif_Fitting\popstats";
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"));
poptab = [alfatab_pop;betotab_pop];
%%
tic
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')

idxlist = find(drivermsk & poptab.Animal==Animal);
NPStats = [];
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
Expi = Tab.Expi; spi = Tab.space; ui = Tab.unitnum; ci = Tab.chan;
for nm = ["Animal", "Expi", "unitstr", "unitnum", "chan", "prefchan", "space","F","F_P"]
S.(nm) = Tab.(nm);
end
iCh = find((MapVarStats(Expi).units.spikeID==ci & MapVarStats(Expi).units.unit_num_arr==ui));
S.iCh = iCh;
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2), MapVarStats(Expi).manif.act_col{spi}, 'uni',1);
[maxAct, maxIdx] = max(actmap_mean,[],'all','linear');
[ri,ci] = ind2sub(size(actmap_mean),maxIdx);
S.phi_max = phi_grid(ri,ci);
S.theta_max = theta_grid(ri,ci);
S.maxAct = maxAct;
thresh = 0.9 * maxAct; 
[CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap_mean, thresh, XX, YY, ZZ, Wgrid);
S.CoMtheta_90 = CoMtheta;
S.CoMphi_90 = CoMphi;
S.CoMrho_90 = CoMrho;
S.Mtheta_90 = Mtheta;
S.Mphi_90 = Mphi;
S.Mrho_90 = Mrho;
thresh = 0.75 * maxAct; 
[CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap_mean, thresh, XX, YY, ZZ, Wgrid);
S.CoMtheta_75 = CoMtheta;
S.CoMphi_75 = CoMphi;
S.CoMrho_75 = CoMrho;
S.Mtheta_75 = Mtheta;
S.Mphi_75 = Mphi;
S.Mrho_75 = Mrho;
thresh = prctile(actmap_mean,90,'all') ;
[CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap_mean, thresh, XX, YY, ZZ, Wgrid);
S.CoMtheta_P90 = CoMtheta;
S.CoMphi_P90 = CoMphi;
S.CoMrho_P90 = CoMrho;
S.Mtheta_P90 = Mtheta;
S.Mphi_P90 = Mphi;
S.Mrho_P90 = Mrho;
thresh = prctile(actmap_mean,75,'all') ;
[CoMvec, CoMtheta, CoMphi, CoMrho, Mvec, Mtheta, Mphi, Mrho] = CentofMass(actmap_mean, thresh, XX, YY, ZZ, Wgrid);
S.CoMtheta_P75 = CoMtheta;
S.CoMphi_P75 = CoMphi;
S.CoMrho_P75 = CoMrho;
S.Mtheta_P75 = Mtheta;
S.Mphi_P75 = Mphi;
S.Mrho_P75 = Mrho;
NPStats = [NPStats, S];
end
%%
toc
NPStatTab = struct2table(NPStats);
writetable(NPStatTab,fullfile(nonpardir,Animal+"_Driver_NonParamStat.csv"))
end
%
NPtab = [];
for Animal = ["Alfa","Beto"]
NPStatTab = readtable(fullfile(nonpardir,Animal+"_Driver_NonParamStat.csv"));
NPtab = [NPtab; NPStatTab];
end
NPtab.theta_max = NPtab.theta_max / 180 * pi;
NPtab.phi_max = NPtab.phi_max / 180 * pi;
writetable(NPtab,fullfile(nonpardir,"Both"+"_Driver_NonParamStat.csv"))

%%
global figdir
figdir = "O:\Manif_NonParam\summary";
msk = NPtab.F_P < 1E-3;
plotTuneCenter(NPtab,msk,"CoM","_P90","CoM 90ptile",'Fsig')
% plotTuneCenter(NPtab,msk,"M","_P90","Mean 90ptile",'Fsig')
plotTuneCenter(NPtab,msk,"CoM","_90","CoM 90% max",'Fsig')
plotTuneCenter(NPtab,msk,"","_max","max",'Fsig')
%%
msk = NPtab.F_P < 1E-3;
Alfamsk = (NPtab.Animal=="Alfa");
Betomsk = (NPtab.Animal=="Beto");
V1msk = (NPtab.chan<=48 & NPtab.chan>=33);
V4msk = (NPtab.chan>48);
ITmsk = (NPtab.chan<33);
plotTuneCenterMultMsk(NPtab,{msk&V1msk, msk&V4msk, msk&ITmsk},["V1","V4","IT"],"CoM","_P90","CoM 90ptile",'Fsig_area_sep')
plotTuneCenterMultMsk(NPtab,{msk&V1msk, msk&V4msk, msk&ITmsk},["V1","V4","IT"],"CoM","_90","CoM 90% max",'Fsig_area_sep')
plotTuneCenterMultMsk(NPtab,{msk&Alfamsk, msk&Betomsk},["Alfa","Beto"],"CoM","_P90","CoM 90ptile",'Fsig_anim_sep')
plotTuneCenterMultMsk(NPtab,{msk&Alfamsk, msk&Betomsk},["Alfa","Beto"],"CoM","_90","CoM 90% max",'Fsig_anim_sep')

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