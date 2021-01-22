%% Visualize all the manifold tuning maps and the fittings like a montage!
%  This script plot the maps on sphere. 
Animal="Both";Set_Path;
tabdir = "O:\Manif_Fitting\Kent_summary";
poptabdir = "O:\Manif_Fitting\popstats";
addpath D:\Github\Fit_Spherical_Tuning
addpath e:\Github_Projects\Fit_Spherical_Tuning
alfatab = readtable(fullfile(tabdir,"Alfa_Kstats.csv"));
betotab = readtable(fullfile(tabdir,"Beto_Kstats.csv"));
preftab = [];
Animal_tab = array2table(repmat("Alfa",size(alfatab,1),1),'VariableNames',{'Animal'});
preftab = [preftab; Animal_tab, alfatab];
Animal_tab = array2table(repmat("Beto",size(betotab,1),1),'VariableNames',{'Animal'});
preftab = [preftab; Animal_tab, betotab];

%%
ang_step = 18;
[theta_grid2, phi_grid2] = meshgrid(-180:18:180, -90:18:90);
theta_grid2d = theta_grid2/180*pi;
phi_grid2d = phi_grid2/180*pi;

figure(23);
T=tiledlayout(5,10,'TileSpacing','compact','padding','compact');
idxlist = find(drivermsk & Alfamsk);
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
fitval = KentFunc(Tab.theta, Tab.phi, Tab.psi, Tab.kappa, Tab.beta, Tab.A, reshape(theta_grid2d,[],1), reshape(phi_grid2d,[],1));
fitval = reshape(fitval,size(theta_grid2d));
ax = nexttile(T,i);
sphere_plot(ax,theta_grid2d,phi_grid2d,fitval);axis off
title(compose("%s Exp %d Ch%s R2 %.3f\n[the,phi]=[%.1f,%.1f] psi=%.1f\nkappa=%.2f beta=%.2f A=%.1f",...
	 Tab.Animal{1}, Tab.Expi, Tab.unitstr{1}, Tab.R2, Tab.theta/pi*180, Tab.phi/pi*180, Tab.psi/pi*180, Tab.kappa, Tab.beta, Tab.A),'FontSize',10);
end
figure(24);
T=tiledlayout(6,11,'TileSpacing','compact','padding','compact');
idxlist = find(drivermsk & Betomsk);
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
fitval = KentFunc(Tab.theta, Tab.phi, Tab.psi, Tab.kappa, Tab.beta, Tab.A, reshape(theta_grid2d,[],1), reshape(phi_grid2d,[],1));
fitval = reshape(fitval,size(theta_grid2d));
ax = nexttile(T,i);
sphere_plot(ax,theta_grid2d,phi_grid2d,fitval);axis off
title(compose("%s Exp %d Ch%s R2 %.3f\n[the,phi]=[%.1f,%.1f] psi=%.1f\nkappa=%.2f beta=%.2f A=%.1f",...
	 Tab.Animal{1}, Tab.Expi, Tab.unitstr{1}, Tab.R2, Tab.theta/pi*180, Tab.phi/pi*180, Tab.psi/pi*180, Tab.kappa, Tab.beta, Tab.A),'FontSize',10);
end
%%
mtgdir = "O:\Manif_Fitting\montage";
saveas(23,fullfile(mtgdir,"Alfa_fit_montage.png"));
saveas(23,fullfile(mtgdir,"Alfa_fit_montage.pdf"));
savefig(23,fullfile(mtgdir,"Alfa_fit_montage.fig"));
saveas(24,fullfile(mtgdir,"Beto_fit_montage.png"));
saveas(24,fullfile(mtgdir,"Beto_fit_montage.pdf"));
savefig(24,fullfile(mtgdir,"Beto_fit_montage.fig"));

%% Fitting Result with Baseline 
% poptab = readtable(fullfile(tabdir,compose("%s_Exp_all_KentStat_bsl_pole.csv",Animal)));
alfatab_pop = readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"));
betotab_pop = readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"));
poptab = [alfatab_pop;betotab_pop];

ang_step = 18;
[theta_grid2, phi_grid2] = meshgrid(-180:18:180, -90:18:90);
theta_grid2d = theta_grid2/180*pi;
phi_grid2d = phi_grid2/180*pi;

figure(21);
T=tiledlayout(5,10,'TileSpacing','compact','padding','compact');
idxlist = find(drivermsk & Alfamsk);
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
fitval = Tab.bsl + KentFunc(Tab.theta, Tab.phi, Tab.psi, Tab.kappa, Tab.beta, Tab.A, reshape(theta_grid2d,[],1), reshape(phi_grid2d,[],1));
fitval = reshape(fitval,size(theta_grid2d));
ax = nexttile(T,i);
sphere_plot(ax,theta_grid2d,phi_grid2d,fitval);axis off
title(compose("%s Exp %d Ch%s R2 %.3f\n[\\theta,\\phi]=[%.1f,%.1f] \\psi=%.1f\n\\kappa=%.2f \\beta=%.2f A=%.1f bsl=%.1f",...
	 Tab.Animal{1}, Tab.Expi, Tab.unitstr{1}, Tab.R2, Tab.theta/pi*180, Tab.phi/pi*180, Tab.psi/pi*180, Tab.kappa, Tab.beta, Tab.A, Tab.bsl),'FontSize',10);
end
figure(22);
T=tiledlayout(6,11,'TileSpacing','compact','padding','compact');
idxlist = find(drivermsk & Betomsk);
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
fitval = Tab.bsl + KentFunc(Tab.theta, Tab.phi, Tab.psi, Tab.kappa, Tab.beta, Tab.A, reshape(theta_grid2d,[],1), reshape(phi_grid2d,[],1));
fitval = reshape(fitval,size(theta_grid2d));
ax = nexttile(T,i);
sphere_plot(ax,theta_grid2d,phi_grid2d,fitval);axis off
title(compose("%s Exp %d Ch%s R2 %.3f\n[\\theta,\\phi]=[%.1f,%.1f] \\psi=%.1f\n\\kappa=%.2f \\beta=%.2f A=%.1f bsl=%.1f",...
	 Tab.Animal{1}, Tab.Expi, Tab.unitstr{1}, Tab.R2, Tab.theta/pi*180, Tab.phi/pi*180, Tab.psi/pi*180, Tab.kappa, Tab.beta, Tab.A, Tab.bsl),'FontSize',10);
end
%%
mtgdir = "O:\Manif_Fitting\montage";
saveas(21,fullfile(mtgdir,"Alfa_bslfit_montage.png"));
saveas(21,fullfile(mtgdir,"Alfa_bslfit_montage.pdf"));
savefig(21,fullfile(mtgdir,"Alfa_bslfit_montage.fig"));
saveas(22,fullfile(mtgdir,"Beto_bslfit_montage.png"));
saveas(22,fullfile(mtgdir,"Beto_bslfit_montage.pdf"));
savefig(22,fullfile(mtgdir,"Beto_bslfit_montage.fig"));



%%
Animal = "Alfa";
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
%%
ang_step = 18;
[phi_grid, theta_grid] = meshgrid(-90:18:90, -90:18:90); 
% Note this is different from above theta, phi here should match the
% underlying parameter map of actmap

figure(25);
T=tiledlayout(5,10,'TileSpacing','compact','padding','compact');
idxlist = find(drivermsk & Alfamsk);
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
Expi = Tab.Expi; spi = Tab.space; ui = Tab.unitnum; ci = Tab.chan;
iCh = find((MapVarStats(Expi).units.spikeID==ci & MapVarStats(Expi).units.unit_num_arr==ui));
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2),MapVarStats(Expi).manif.act_col{spi},'uni',1);
ax = nexttile(T,i);
sphere_plot(ax,theta_grid,phi_grid,actmap_mean);axis off
title(compose("%s Exp %d Ch%s R2 %.3f\n[the,phi]=[%.1f,%.1f] psi=%.1f\nkappa=%.2f beta=%.2f A=%.1f",...
	 Tab.Animal{1}, Tab.Expi, Tab.unitstr{1}, Tab.R2, Tab.theta/pi*180, Tab.phi/pi*180, Tab.psi/pi*180, Tab.kappa, Tab.beta, Tab.A),'FontSize',10);
end

Animal = "Beto";
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
figure(26);
T=tiledlayout(6,11,'TileSpacing','compact','padding','compact');
idxlist = find(drivermsk & Betomsk);
for i = 1:numel(idxlist)
Tab = poptab(idxlist(i),:);
Expi = Tab.Expi; spi = Tab.space; ui = Tab.unitnum; ci = Tab.chan;
iCh = find((MapVarStats(Expi).units.spikeID==ci & MapVarStats(Expi).units.unit_num_arr==ui));
actmap_mean = cellfun(@(A)...
    mean(A(iCh,:),2),MapVarStats(Expi).manif.act_col{spi},'uni',1);
ax = nexttile(T,i);
sphere_plot(ax,theta_grid,phi_grid,actmap_mean);axis off
title(compose("%s Exp %d Ch%s R2 %.3f\n[the,phi]=[%.1f,%.1f] psi=%.1f\nkappa=%.2f beta=%.2f A=%.1f",...
	 Tab.Animal{1}, Tab.Expi, Tab.unitstr{1}, Tab.R2, Tab.theta/pi*180, Tab.phi/pi*180, Tab.psi/pi*180, Tab.kappa, Tab.beta, Tab.A),'FontSize',10);
end
%%
mtgdir = "O:\Manif_Fitting\montage";
saveas(25,fullfile(mtgdir,"Alfa_real_bsl_montage.png"));
saveas(25,fullfile(mtgdir,"Alfa_real_bsl_montage.pdf"));
savefig(25,fullfile(mtgdir,"Alfa_real_bsl_montage.fig"));
saveas(26,fullfile(mtgdir,"Beto_real_bsl_montage.png"));
saveas(26,fullfile(mtgdir,"Beto_real_bsl_montage.pdf"));
savefig(26,fullfile(mtgdir,"Beto_real_bsl_montage.fig"));

%% Single Sphere Demo. What it looks like?
% [phi_grid2, theta_grid2] = meshgrid(-90:18:90, -180:ang_step:180);%phi_arr
[theta_grid2, phi_grid2] = meshgrid(-180:18:180, -90:18:90);
theta_grid2d = theta_grid2/180*pi;
phi_grid2d = phi_grid2/180*pi;
figure;ax=subplot(1,1,1);
fitval = KentFunc(Tab.theta, Tab.phi, Tab.psi, Tab.kappa, Tab.beta, Tab.A, reshape(theta_grid2d,[],1), reshape(phi_grid2d,[],1));
fitval = reshape(fitval,size(theta_grid2d));
sphere_plot(ax,theta_grid2d,phi_grid2d,fitval);axis off


