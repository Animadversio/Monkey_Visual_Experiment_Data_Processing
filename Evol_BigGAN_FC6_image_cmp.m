% Analyze Evolved image in BigGAN and FC6 

setMatlabTitle("BigGAN FC6 image comparison"); 
%% Image Analysis for BigGAN and FC6
D = torchImDist("squeeze", true);
%%
res = 119;
distmap = D.distance(imresize(im1,[res,res]), imresize(im2,[res,res]));
figure(1);clf; T=tiledlayout(1,3,'Padding','compact');
nexttile(1);imshow(im1); 
nexttile(2);imshow(im2); 
nexttile(3);imagesc(distmap);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map")
%%
Animal="Both"; Set_Path;
ftrrows = find(contains(ExpRecord.expControlFN,"generate_") & (ExpRecord.Exp_collection=="BigGAN_fc6" |...
               ExpRecord.Exp_collection=="BigGAN_FC6"));
figdir = "E:\OneDrive - Washington University in St. Louis\Evol_BigGAN_FC6"; 
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows, Animal, false);
%%
stimfdr = "N:\Stimuli\2020-BigGAN\2020-07-22-Beto-01\2020-07-22-10-14-22";
im1 = imread(fullfile(stimfdr, "block048_thread000_gen_gen047_003089.bmp"));
im2 = imread(fullfile(stimfdr, "block048_thread001_gen_gen047_003114.bmp"));
%%
im1 = imread(fullfile("N:\Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55","block043_thread000_gen_gen042_001709.bmp"));
im2 = imread(fullfile("N:\Stimuli\2020-BigGAN\2020-07-20-Beto-01\2020-07-20-12-40-51","block042_thread000_gen_gen041_001046.bmp"));
%%
imfdr = "N:\Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55";
imnames = string(ls(imfdr+"\block043_thread000*"));
im_gen1 = cellfun(@(impath)imread(fullfile(imfdr,impath)),imnames,'Uni',0);
imfdr = "N:\Stimuli\2020-BigGAN\2020-07-20-Beto-01\2020-07-20-12-40-51";
imnames = string(ls(imfdr+"\block042_thread000*"));
im_gen2 = cellfun(@(impath)imread(fullfile(imfdr,impath)),imnames,'Uni',0);
%% Demo of population based distance map averaging procedure.
res = 120;
distmap_col = {};
for i = 1:numel(im_gen1)
    for j = 1:numel(im_gen2)
    im1 = im_gen1{i}; im2 = im_gen2{j};
    distmap_col{i,j} = D.distance(imresize(im1,[res,res]), imresize(im2,[res,res]));
    end
end
distmap_stack = cell2mat(reshape(distmap_col,1,1,[]));
distmap_mean = mean(distmap_stack,3); 
figure;imagesc(distmap_mean);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map")
%% 
for i = 1:numel(im_gen1)
    for j = 1:numel(im_gen2)
    im1 = im_gen1{i}; im2 = im_gen2{j};
    distmap = D.distance(imresize(im1,[res,res]), imresize(im2,[res,res]));
    figure(1);clf; T=tiledlayout(1,3,'Padding','compact');
    nexttile(1);imshow(im1)
    nexttile(2);imshow(im2)
    nexttile(3);imagesc(distmap);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map")
    title(T,compose("FC6 %d BigGAN %d", i, j))
    pause
    end
end
%% Load the images for 2 threads in paired generation, compute the average distance map. 
[imFC_fina, imnmFC_fina] = loadEvolImages("N:\Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55", 0, -2);
[imBG_fina, imnmBG_fina] = loadEvolImages("N:\Stimuli\2020-BigGAN\2020-07-20-Beto-01\2020-07-20-12-40-51", 0, -2);
[distmap_fina, distmap_col_fina] = calc_distmap_avg(D,imFC_fina, imBG_fina);

[imFC_mid, imnmFC_mid] = loadEvolImages("N:\Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55", 0, 21);
[imBG_mid, imnmBG_mid] = loadEvolImages("N:\Stimuli\2020-BigGAN\2020-07-20-Beto-01\2020-07-20-12-40-51", 0, 21);
[distmap_mid, distmap_col_mid] = calc_distmap_avg_subsamp(D, imFC_mid, imBG_mid);
%
[imFC_init, imnmFC_init] = loadEvolImages("N:\Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55", 0, 2);
[imBG_init, imnmBG_init] = loadEvolImages("N:\Stimuli\2020-BigGAN\2020-07-20-Beto-01\2020-07-20-12-40-51", 0, 2);
[distmap_init, distmap_col_init] = calc_distmap_avg(D,imFC_init, imBG_init);
%%
figdir = "E:\OneDrive - Washington University in St. Louis\BigGAN_FC6_Evol_cmp\2020-07-20-Beto-01";
figure; T = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile(1);imagesc(distmap_init);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map");title("First Gen")
nexttile(2);imagesc(distmap_mid);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map");title("Mid Gen")
nexttile(3);imagesc(distmap_fina);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map");title("Final Gen")
%%
figure(13); T = tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
nexttile(3);imagesc(distmap_init);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map");title("First Gen")
nexttile(1);imshow(imFC_init{randsample(numel(imFC_init),1)});title("FC6 GAN")
nexttile(2);imshow(imBG_init{randsample(numel(imBG_init),1)});title("BigGAN")
nexttile(6);imagesc(distmap_mid);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map");title("Mid Gen")
nexttile(4);imshow(imFC_mid{randsample(numel(imFC_mid),1)});
nexttile(5);imshow(imBG_mid{randsample(numel(imBG_mid),1)});
nexttile(9);imagesc(distmap_fina);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map");title("Final Gen")
nexttile(7);imshow(imFC_fina{randsample(numel(imFC_fina),1)});
nexttile(8);imshow(imBG_fina{randsample(numel(imBG_fina),1)});
title(T, compose("2020-07-20 Beto 01 Evol BigGAN FC6 Chan 20"))
saveas(13,fullfile(figdir, "Evol_Img_Similarity.png"))
savefig(13,fullfile(figdir, "Evol_Img_Similarity.fig"))
%%
figure(5); clf;hold on;
histogram(distmap_init(:),0.2:0.01:1,'FaceAlpha',0.5);
histogram(distmap_mid(:),0.2:0.01:1,'FaceAlpha',0.5); 
histogram(distmap_fina(:),0.2:0.01:1,'FaceAlpha',0.5); 
xlabel("Value in Distmap");legend(["Initial Generation", "Middle Generation", "Final Generation"])
title(["Patch Distance Distribution between BigGAN and FC6.","(First, Mid, Last Gen)"])
saveas(5,fullfile(figdir, "Img_Distmap_distribution.png"))
%%
figure;imagesc(distmap_fina);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map")
%% Create that figure Gen by Gen
FCdir = "N:\Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55";
BGdir = "N:\Stimuli\2020-BigGAN\2020-07-20-Beto-01\2020-07-20-12-40-51";
distmap_traj = {};
distmap_col_traj = {};
tic
for geni = 1:42
  [imFC_fina, imnmFC_fina] = loadEvolImages(FCdir, 0, geni);
  [imBG_fina, imnmBG_fina] = loadEvolImages(BGdir, 0, geni);
  [distmap_fina, distmap_col_fina] = calc_distmap_avg_subsamp(D,imFC_fina, imBG_fina, 10);
  distmap_traj{geni} = distmap_fina;
  distmap_col_traj{geni} = distmap_col_fina;
  toc
end
%%
for geni = 1:42
    dist_prctile(geni) = prctile(distmap_traj{geni}(:),5);
    distmap_allpair = cell2mat(reshape(distmap_col_traj{geni},1,1,[]));
    dist_pair_prctile(geni) = prctile(distmap_allpair(:),5);
end
figure;hold on;
plot(dist_prctile);plot(dist_pair_prctile)
legend(["5%-ile in avg distmap", "5%-ile in distmap of all pairs"])
title(["BigGAN FC6 Dual Evol Image Cmp","2020-07-20-Beto-01"])
xlabel("Generation");ylabel("Patch Distance (Squeeze)")
%%
savefig(2,fullfile("E:\OneDrive - Washington University in St. Louis\BigGAN_FC6_ImgCmp\2020-07-20-Beto-01", "Evol_imgsim_traj_2020-07-20-Beto-01.fig"))
saveas(2,fullfile("E:\OneDrive - Washington University in St. Louis\BigGAN_FC6_ImgCmp\2020-07-20-Beto-01", "Evol_imgsim_traj_2020-07-20-Beto-01.png"))
%%
function [im_gen, imnames] = loadEvolImages(stimfolder, thread, block, suffix)
% Util function to load in evolved image of certain block x thread from a
% folder
% block: negative `block` variable means counting from the last block. e.g. -1 means
%        the last block,-2 means the penultimate block
if nargin<=3, suffix=".bmp";end
if block <= 0 % negative block number means counting from the last one
allimgs = string(ls(stimfolder+"\block*thread*gen*"));
lastgeni = regexp(allimgs(end), "block(\d*)_thread\d*_gen", 'tokens');
lastgeni = str2num(lastgeni{1});
blockid = lastgeni + block + 1;
else
blockid = block;
end
imnames = string(ls(stimfolder+compose("\\block%03d_thread%03d_gen*",blockid,thread)));
im_gen = cellfun(@(impath)imread(fullfile(stimfolder,impath)),imnames,'Uni',0);
end

function [distmap_mean, distmap_col] = calc_distmap_avg(D, im_gen1, im_gen2, res)
% Compute the image distance map between each pair of image from 2 cell
% arrays of images. 
if nargin <= 3,  res= 120;  end
distmap_col = {};
for i = 1:numel(im_gen1)
    for j = 1:numel(im_gen2)
    im1 = im_gen1{i}; im2 = im_gen2{j};
    distmap_col{i,j} = D.distance(imresize(im1,[res,res]), imresize(im2,[res,res]));
    end
end
distmap_stack = cell2mat(reshape(distmap_col,1,1,[]));
distmap_mean = mean(distmap_stack,3); 
end

function [distmap_mean, distmap_col] = calc_distmap_avg_subsamp(D, im_gen1, im_gen2, sampleN, res)
% Compute the image distance map between each pair of image from 2 cell
% arrays of images. 
if nargin <= 4,  res= 120;  end
if nargin == 3 %, sampleN = 0;
imgidx_i = 1:numel(im_gen1);
imgidx_j = 1:numel(im_gen2);
else
imgidx_i = randsample(numel(im_gen1),sampleN)';
imgidx_j = randsample(numel(im_gen2),sampleN)';
end
distmap_col = {};
for i = imgidx_i
    for j = imgidx_j
    im1 = im_gen1{i}; im2 = im_gen2{j};
    distmap_col{i,j} = D.distance(imresize(im1,[res,res]), imresize(im2,[res,res]));
    end
end
distmap_stack = cell2mat(reshape(distmap_col,1,1,[]));
distmap_mean = mean(distmap_stack,3); 
end
