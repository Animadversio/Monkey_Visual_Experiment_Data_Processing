%% vis L1 ball 
radius = 2;
ticks = -3:0.1:3;
cent = [0,0,1];
[XX, YY, ZZ] = meshgrid(ticks,ticks,ticks);
V = abs(XX-cent(1)) + abs(YY-cent(2)) + abs(ZZ-cent(3));
figure(2);clf
scatter3(cent(1),cent(2),cent(3),9);
hold on
%%
fv = isosurface(XX,YY,ZZ,V,radius);
p = patch(fv);
isonormals(XX,YY,ZZ,V,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
p.FaceAlpha = 0.35;
daspect([1 1 1]);axis equal
view(3); 
camlight 
lighting flat
% Section with x, y plane
sectRad = radius - sum(abs(cent));
xsect = [ 0 sectRad 0 -sectRad]+cent(1);
ysect = [-sectRad 0 sectRad  0]+cent(2);
psect = patch(xsect, ysect, 'blue','FaceAlpha',0.9);
% x, y plane background
XLIM = xlim();YLIM = ylim();
xplane = [XLIM, XLIM(end:-1:1)];
yplane = reshape([YLIM; YLIM],1,[]);
zplane = -1E-3*ones(1,4);
pplane = patch(xplane,yplane,zplane,'black','FaceAlpha',0.1);
%% 
savefig(2,"L1sphere_plane_section_SYY.fig")