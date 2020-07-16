%% GAN invert batch example
addpath D:\Github\CMAES_optimizer_matlab\geomtry_util
%%
data = load("D:\Generator_DB_Windows\init_population\init_code.mat");
init_code = data.codes; % 30 by 4096
init_imgs = G.visualize(init_code);
%%
MAXSTEPS = 150;
[code_fit_med, img_fit_med, loss_med] = GAN_invert_fun(G, init_imgs, MAXSTEPS);
%%
MAXSTEPS = 300;
[code_fit, img_fit, loss] = GAN_invert_fun(G, init_imgs, MAXSTEPS);
%%
figure; 
subtightplot(1,3,1)
montage(init_imgs)
subtightplot(1,3,2)
montage(img_fit_med)
subtightplot(1,3,3)
montage(img_fit)
%%
corrmat = pdist([init_code;code_fit_med;code_fit], 'correlation');
corrmat = squareform(corrmat);
figure(13);
imagesc(1-corrmat);
axis image; colorbar()
xticks([15,45,75]);xticklabels(["original code", "fit code (150 step)", "fit code (300 step)"])
yticks([15,45,75]);yticklabels(["original code", "fit code (150 step)", "fit code (300 step)"]);ytickangle(30)
title("Code Correlation Between original init code and fit code")
%%
x = cmdscale(corrmat);
labels = [ones(size(init_code,1),1);  ones(size(code_fit_med,1),1).*2; ones(size(code_fit,1),1).*3 ; ones(size(pasu_code,1),1).*4 ];
figure;
gscatter(x(:,1),x(:,2),labels);
%%
corrmat = pdist([init_code;code_fit_med;code_fit;pasu_code], 'correlation');
corrmat = squareform(corrmat);
figure;
imagesc(1-corrmat);axis image
%%

