% Analyze 
%% Image Analysis for BigGAN and FC6
D = torchImDist("squeeze", true);
% loadEvolImages(thread=1)
% loadEvolImages(thread=2)
% D.spatial=true;
%%
res = 119;
distmap = D.distance(imresize(im1,[res,res]), imresize(im2,[res,res]));
figure(1);clf; T=tiledlayout(1,3,'Padding','compact');
nexttile(1);imshow(im1)
nexttile(2);imshow(im2)
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
%%
res = 120;
distmap_col = {};
for i = 1:numel(im_gen1)
    for j = 1:numel(im_gen2)
    im1 = im_gen1{i}; im2 = im_gen2{j};
    distmap_col{i,j} = D.distance(imresize(im1,[res,res]), imresize(im2,[res,res]));
    end
end
%%
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
%%
[imFC_fina, imnmFC_fina] = loadEvolImages("N:\Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55", 0, -1);
[imBG_fina, imnmBG_fina] = loadEvolImages("N:\Stimuli\2020-BigGAN\2020-07-20-Beto-01\2020-07-20-12-40-51", 0, -1);
[distmap_infa, distmap_col_fina] = calc_distmap_avg(D,imFC_fina, imBG_fina);
%
[imFC_init, imnmFC_init] = loadEvolImages("N:\Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55", 0, 1);
[imBG_init, imnmBG_init] = loadEvolImages("N:\Stimuli\2020-BigGAN\2020-07-20-Beto-01\2020-07-20-12-40-51", 0, 1);
[distmap_init, distmap_col_init] = calc_distmap_avg(D,imFC_init, imBG_init);
%%
figure;imagesc(distmap_infa);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map")
figure;imagesc(distmap_init);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map")
%%
function [im_gen, imnames] = loadEvolImages(stimfolder, thread, block, suffix)
if nargin<=3, suffix=".bmp";end
if block <= 0
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
function [distmap_mean, distmap_col] = calc_distmap_avg(D,im_gen1, im_gen2)
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
end