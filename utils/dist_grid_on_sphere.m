function [ang_dist, cos_dist] = dist_grid_on_sphere(q_ang1, q_ang2)
% compute a grid of spherical distance (angle or 1-cos) on the hemisphere 
% to a vector (q_ang1, q_ang2)
ang_step = 18;
theta = [-5:5] * ang_step;
phi = [-5:5] * ang_step;
[PHI, THETA] = meshgrid(phi, theta);
vectmat = cat(3, cos(THETA/180*pi) .* cos(PHI/180*pi),...
                 sin(THETA/180*pi) .* cos(PHI/180*pi),...
                 sin(PHI/180*pi));
%q_ang1 = 0;q_ang2 = 0;
% create query vector from the angle
q_vect = [cos(q_ang1/180*pi) .* cos(q_ang2/180*pi),...
          sin(q_ang1/180*pi) .* cos(q_ang2/180*pi),...
          sin(q_ang2/180*pi)]; % create a 3d grid on hemisphere
q_vect = reshape(q_vect, 1, 1, []);
cos_arr = sum(vectmat .* q_vect,3);
ang_dist = acos(cos_arr);
cos_dist = 1 - cos_arr;
return