%% Code for bilinear interpolation in nd array
interpi = 6 + data.Theta / pi * 10;
interpj = 6 + data.Phi / pi * 10;
i_grid = [floor(interpi), ceil(interpi)];
j_grid = [floor(interpj), ceil(interpj)];
if interpi == floor(interpi)
    i_W = [1, 0];
else
    i_W = [i_grid(2) - interpi, interpi - i_grid(1)];
end
if interpj == floor(interpj)
    j_W = [1, 0];
else
    j_W = [j_grid(2) - interpj, interpj - j_grid(1)];
end
sum(psth_avg_tsr(:,:,i_grid,j_grid) .* reshape(i_W'*j_W,[1,1,2,2]),[1,3,4])