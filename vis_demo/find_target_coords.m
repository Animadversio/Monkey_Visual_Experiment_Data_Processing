cent = [9.5, -5];
X0=cent(1); Y0=cent(2);
Radi = norm(cent);
%%
ang = rad2deg(atan2(cent(2), cent(1)));
%%
Delta_ang = 30;
[X1,Y1] = pol2cart(deg2rad(ang + Delta_ang), Radi);
[X2,Y2] = pol2cart(deg2rad(ang - Delta_ang), Radi);
%
L = 7;
figure;hold on 
rectangle('Position',[X0-L/2, Y0-L/2, L, L])
rectangle('Position',[X1-L/2, Y1-L/2, L, L])
rectangle('Position',[X2-L/2, Y2-L/2, L, L])
axis equal
fprintf("Delta Angle(deg): %f\n",Delta_ang)
fprintf("Origin  %.3f,%.3f\n",X0,Y0)
fprintf("Choic1  %.3f,%.3f\n",X1,Y1)
fprintf("Choic2  %.3f,%.3f\n",X2,Y2)