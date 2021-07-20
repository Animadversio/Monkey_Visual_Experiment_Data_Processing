%% vis Kent
%  This script is to show how the shape of a Kent function depends on the
%  parameters. Showing a grid of Kents.
addpath D:\Github\Fit_Spherical_Tuning
addpath e:\Github_Projects\Fit_Spherical_Tuning
demodir = "O:\Manuscript_Manifold\Figure2\Kent";
%%
[theta_grid2, phi_grid2] = meshgrid(-180:9:180, -90:9:90);
theta_grid2d = theta_grid2/180*pi;
phi_grid2d = phi_grid2/180*pi;

figure(3);clf;set(3,'pos',[703   47  1050  950])
[kappa_grid, beta_grid] = meshgrid([0, 0.25, 0.5, 1, 2, 4],[0, 0.25, 0.5, 1, 2, 4]);
T = tiledlayout(6,6,'TileSpacing','compact','padding','compact');
for i=1:numel(kappa_grid)
kappa = kappa_grid(i);beta = beta_grid(i);
fitval = KentFunc(0.0, 0.0, pi/4, kappa, beta, 1, reshape(theta_grid2d,[],1), reshape(phi_grid2d,[],1));
fitval = reshape(fitval,size(theta_grid2d));
ax = nexttile(T,i);
sphere_plot(ax,theta_grid2d,phi_grid2d,fitval);axis off
title(compose("\\kappa=%.2f \\beta=%.2f",kappa,beta))
end
title(T,"A=1, b=0, (\theta,\phi)=(0,0), \psi=\pi/4")
demodir = "O:\Manuscript_Manifold\Figure2\Kent";
saveallform(demodir, "KentFuncDemo2",3)

%%
[theta_grid2, phi_grid2] = meshgrid(-90:18:90, -90:18:90); % sparser grid smaller file...
theta_grid2d = theta_grid2/180*pi;
phi_grid2d = phi_grid2/180*pi;
figure(1);clf;set(3,'pos',[703   47  1050  950])
[kappa_grid, beta_grid] = meshgrid([0, 0.25, 0.5, 1, 2, 4],[0, 0.25, 0.5, 1, 2, 4]);
T = tiledlayout(6,6,'TileSpacing','compact','padding','compact');
for i=1:numel(kappa_grid)
kappa = kappa_grid(i);beta = beta_grid(i);
fitval = KentFunc(0.0, 0.0, pi/4, kappa, beta, 1, reshape(theta_grid2d,[],1), reshape(phi_grid2d,[],1));
fitval = reshape(fitval,size(theta_grid2d));
ax = nexttile(T,i);
sphere_plot(ax,theta_grid2d,phi_grid2d,fitval);axis off
title(compose("\\kappa=%.2f \\beta=%.2f",kappa,beta))
end
title(T,"A=1, b=0, (\theta,\phi)=(0,0), \psi=\pi/4")
saveallform(demodir, "KentFuncDemo4",1) % this function can save the vector pdf
