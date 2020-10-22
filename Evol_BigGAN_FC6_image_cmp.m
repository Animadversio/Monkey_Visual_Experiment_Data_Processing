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