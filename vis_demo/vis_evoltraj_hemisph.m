Animal = "Beto";
Set_Path;mat_dir = "O:\Mat_Statistics";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
for Expi = 11
EStats(Expi).meta.stimuli;
[codes_all, img_ids, generations] = load_codes_all(EStats(Expi).meta.stimuli, 1); %, "all"

end
%%
[basisvec, projcoef, latent, tsquared, explained] = pca(codes_all,'NumComponents',3);
%%
figure;
scatter3(projcoef(:,1),projcoef(:,2),projcoef(:,3),36,generations)
axis equal image
%%
sphnorm = mean(norm_axis(codes_all(generations==max(generations),:),2));
figure;
axis equal image
%%
[PHI,THE] = meshgrid(-90:18:90,-90:18:90);
[XX,YY,ZZ] = sph2cart(THE/180*pi,PHI/180*pi,sphnorm);
figure;surf(XX,YY,ZZ);axis equal;shading flat
xlabel("PC1");ylabel("PC2");zlabel("PC3")
%%
figure(7);hold on
[PHI,THE] = meshgrid(-90:18:90,-90:18:90);
[XX,YY,ZZ] = sph2cart(THE/180*pi,PHI/180*pi,sphnorm);
for theta = -90:18:90
[X,Y,Z] = sph2cart(theta/180*pi,[-90:90]/180*pi,sphnorm);
plot3(X,Y,Z,'Color',[0,0,0,0.5])
[X,Y,Z] = sph2cart([-90:90]/180*pi,theta/180*pi,sphnorm);
if numel(Z) < numel(X), Z = repmat(Z,size(X,1),size(X,2)); end
plot3(X,Y,Z,'Color',[0,0,0,0.5])
end
% colormap("parula")
colorimg = squeeze(ind2rgb(generations',parula(80)));
scatter3(codes_all*basisvec(:,1),codes_all*basisvec(:,2),codes_all*basisvec(:,3),36,colorimg)
surf(XX,YY,ZZ,'FaceAlpha',0.2,'FaceColor',[0.8500 0.3250 0.0980]);%shading flat
axis equal image;
xlabel("PC1");ylabel("PC2");zlabel("PC3")
%%
figdir = "O:\Manuscript_Manifold\Figure1";
savenm = "Beto_Exp11_Evol_Traj_Hemisphere";
savefig(7,fullfile(figdir,savenm+".fig"))
saveas(7,fullfile(figdir,savenm+".png"))
saveas(7,fullfile(figdir,savenm+".pdf"))
% saveas(7,fullfile(figdir,savenm+".png"))
%%
exportgraphics 