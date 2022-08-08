% Plot a mountain landscape as a histogram / contour / surface plot.
[XX,YY] = meshgrid(-5:0.1:5,-5:0.1:5);
%%
figdir = "E:\OneDrive - Harvard University\Manuscript_Manifold\Response\Schematics";
%%
height = ampl + 0.1*perlin_smth;
figure();
subplot(121)
surfc(XX,YY,height,'edgecolor','none');
lighting flat
subplot(122)
% [f,xi] = ksdensity(ampl); 
% plot(xi,f);
histogram(height(:))
%%
[XX,YY] = meshgrid(-5:0.1:5,-5:0.1:5);
X0 = -5; Y0 = 5; c = 0.1; sigma = 4;
ampl = exp(-((XX - X0) .^ 2 + (YY - Y0) .^ 2 + c * (XX - X0) .* (YY - Y0))/sigma^2);
%%
n = 101;
m = 101;
im = zeros(n, m);
perlin_im = perlin_noise(im);
perlin_im = perlin_im / std(perlin_im,1,'all');
perlin_smth = imgaussfilt(perlin_im, 3);
% figure; imagesc(perlin_im); colormap gray;
%%
height = 0.2+ ampl + 0.1*perlin_smth;
figure('pos',[400,300,1400,400]);
subplot(131)
surfc(XX,YY,height,'edgecolor','none');hold on 
contour3(XX,YY,height,'LevelList',[0:0.1:1.2],'Color','black')
% lighting flat
lighting gouraud
view([az,el]);
subplot(132)
histogram(height(:))
subplot(133)
contourf(XX,YY,height,'LevelList',[0:0.1:1.2])
saveallform(figdir, "landscape_contour_2")
save(fullfile(figdir,'plotdata2.mat'),'height','XX','YY','ampl','perlin_smth')
%%
function im = perlin_noise(im)

    [n, m] = size(im);
    i = 0;
    w = sqrt(n*m);

    while w > 3
        i = i + 1;
        d = interp2(randn(n, m), i-1, 'spline');
        im = im + i * d(1:n, 1:m);
        w = w - ceil(w/2 - 1);
    end
end