%%
D = torchLPIPS("squeeze", 1); % spatial output
%% Compare Evolved image in BigGAN and FC6 in batch
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
%%
figdir = "E:\OneDrive - Washington University in St. Louis\BigGAN_FC6_ImgCmp"; 
for Expi = 1:numel(BFEStats)
tic
stimfolder = BFEStats(Expi).meta.stimuli;
[imFC_init, imnmFC_init] = loadEvolImages(stimfolder, 0, 1);
[imBG_init, imnmBG_init] = loadEvolImages(stimfolder, 1, 1);
[distmap_init, distmap_col_init] = calc_distmap_avg(D,imFC_init, imBG_init);
toc
midgenN = round(BFEStats(Expi).evol.block_n / 2); 
[imFC_midd, imnmFC_midd] = loadEvolImages(stimfolder, 0, midgenN);
[imBG_midd, imnmBG_midd] = loadEvolImages(stimfolder, 1, midgenN);
[distmap_midd, distmap_col_midd] = calc_distmap_avg(D,imFC_midd, imBG_midd);
toc
[imFC_fina, imnmFC_fina] = loadEvolImages(stimfolder, 0, -2);
[imBG_fina, imnmBG_fina] = loadEvolImages(stimfolder, 1, -2);
[distmap_fina, distmap_col_fina] = calc_distmap_avg(D,imFC_fina, imBG_fina);
toc
save(fullfile(figdir, compose("%s_Exp%02d.mat",Animal,Expi)), 'distmap_init', 'distmap_col_init', ...
	'distmap_midd', 'distmap_col_midd', 'distmap_fina', 'distmap_col_fina', 'midgenN')
stimparts = split(BFEStats(Expi).meta.stimuli,"\");
prefchan_str = BFEStats(Expi).units.unit_name_arr(BFEStats(Expi).units.pref_chan_id(1));
% fdrnm = compose("%s-Chan%02d", stimparts{end-1}, pref_chan(1));
% figdir = fullfile(saveroot, fdrnm);
%%
figure(13); T = tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
nexttile(3);imagesc(distmap_init);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map");title("First Gen")
nexttile(1);imshow(imFC_init{randsample(numel(imFC_init),1)});title("FC6 GAN")
nexttile(2);imshow(imBG_init{randsample(numel(imBG_init),1)});title("BigGAN")
nexttile(6);imagesc(distmap_midd);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map");title("Mid Gen")
nexttile(4);imshow(imFC_midd{randsample(numel(imFC_midd),1)});
nexttile(5);imshow(imBG_midd{randsample(numel(imBG_midd),1)});
nexttile(9);imagesc(distmap_fina);axis image;colormap(flipud(gray));colorbar();xlabel("patch dist map");title("Final Gen")
nexttile(7);imshow(imFC_fina{randsample(numel(imFC_fina),1)});
nexttile(8);imshow(imBG_fina{randsample(numel(imBG_fina),1)});
title(T, compose("%s Evol BigGAN FC6 Chan %s", stimparts{end-1}, prefchan_str))
saveas(13,fullfile(figdir, compose("%s_Exp%02d_chan%s_Evol_Img_Similarity.png",Animal,Expi,prefchan_str)))
savefig(13,fullfile(figdir, compose("%s_Exp%02d_chan%s_Evol_Img_Similarity.fig",Animal,Expi,prefchan_str)))
toc
%%
figure(8); clf;hold on;
histogram(distmap_init(:),0.2:0.01:.8,'FaceAlpha',0.4);
histogram(distmap_midd(:),0.2:0.01:.8,'FaceAlpha',0.4); 
histogram(distmap_fina(:),0.2:0.01:.8,'FaceAlpha',0.4); 
xlabel("Value in Distmap");legend(["Initial Generation", "Middle Generation", "Final Generation"])
title(["Patch Distance Distribution between BigGAN and FC6.","(First, Mid, Last Gen)"])
saveas(8,fullfile(figdir, compose("%s_Exp%02d_chan%s_Distmap_distribution.png",Animal,Expi,prefchan_str)))
savefig(8,fullfile(figdir, compose("%s_Exp%02d_chan%s_Distmap_distribution.fig",Animal,Expi,prefchan_str)))
toc
end
end
%% Redo the plot... did sth wrong.
figdir = "E:\OneDrive - Washington University in St. Louis\BigGAN_FC6_ImgCmp"; 
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
for Expi = 1:numel(BFEStats)
prefchan_str = BFEStats(Expi).units.unit_name_arr(BFEStats(Expi).units.pref_chan_id(1));
load(fullfile(figdir, compose("%s_Exp%02d.mat",Animal,Expi)), 'distmap_init', 'distmap_col_init', ...
	'distmap_midd', 'distmap_col_midd', 'distmap_fina', 'distmap_col_fina', 'midgenN')
figure(8); clf;hold on;
histogram(distmap_init(:),0.2:0.01:.8,'FaceAlpha',0.4);
histogram(distmap_midd(:),0.2:0.01:.8,'FaceAlpha',0.4); 
histogram(distmap_fina(:),0.2:0.01:.8,'FaceAlpha',0.4); 
xlabel("Value in Distmap");legend(["Initial Generation", "Middle Generation", "Final Generation"])
title(compose("Patch Distance Distribution between BigGAN and FC6.%s Exp %d\n(First, Mid, Last Gen)",Animal,Expi))
saveas(8,fullfile(figdir, compose("%s_Exp%02d_chan%s_Distmap_distribution.png",Animal,Expi,prefchan_str)))
savefig(8,fullfile(figdir, compose("%s_Exp%02d_chan%s_Distmap_distribution.fig",Animal,Expi,prefchan_str)))
end
end
%%

mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
figdir = "E:\OneDrive - Washington University in St. Louis\BigGAN_FC6_ImgCmp"; 
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
for Expi = 1:numel(BFEStats)
    tic
    fprintf("Processing %s Exp %d\n",Animal,Expi)
    FCdir = BFEStats(Expi).meta.stimuli;
    BGdir = BFEStats(Expi).meta.stimuli;
    genN = generation_num_from_fdr(BGdir);
    distmap_traj = {};
    distmap_col_traj = {};
    dist_prctile = [];
    dist_pair_prctile = [];
    toc
    for geni = 1:genN-1
      [imFC_fina, imnmFC_fina] = loadEvolImages(FCdir, 0, geni);
      [imBG_fina, imnmBG_fina] = loadEvolImages(BGdir, 0, geni);
      [distmap_fina, distmap_col_fina] = calc_distmap_avg_subsamp(D, imFC_fina, imBG_fina, 10);
      
      distmap_traj{geni} = distmap_fina;
      distmap_col_traj{geni} = cell2mat(reshape(distmap_col_fina,1,1,[])); % distmap_col_fina;
      
      dist_prctile(geni) = prctile(distmap_traj{geni}(:),5);
      dist_pair_prctile(geni) = prctile(distmap_col_traj{geni}(:),5);
    end
    toc
    save(fullfile(figdir, compose("%s_Exp%02d_traj.mat",Animal,Expi)), 'distmap_traj', 'distmap_col_traj', ...
            'dist_prctile', 'dist_pair_prctile', 'genN')
    
    stimparts = split(BFEStats(Expi).meta.stimuli,"\");
    prefchan_str = BFEStats(Expi).units.unit_name_arr(BFEStats(Expi).units.pref_chan_id(1));
    
    figure(10);clf;hold on;
    plot(dist_prctile);plot(dist_pair_prctile)
    legend(["5%-ile in avg distmap", "5%-ile in distmap of all pairs"])
    title(["BigGAN FC6 Dual Evol Image Cmp", stimparts{end-1}])
    xlabel("Generation");ylabel("Patch Distance (Squeeze)")
    saveas(10,fullfile(figdir, compose("%s_Exp%02d_chan%s_Dist_prctile_traj.png",Animal,Expi,prefchan_str)))
    savefig(10,fullfile(figdir, compose("%s_Exp%02d_chan%s_Dist_prctile_traj.fig",Animal,Expi,prefchan_str)))
end
end
%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
figdir = "E:\OneDrive - Washington University in St. Louis\BigGAN_FC6_ImgCmp"; 
for Animal = ["Alfa", "Beto"]
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
for Expi = 1:numel(BFEStats)
    load(fullfile(figdir, compose("%s_Exp%02d_traj.mat",Animal,Expi)), ...
            'dist_prctile', 'dist_pair_prctile')
    stimparts = split(BFEStats(Expi).meta.stimuli,"\");
    prefchan_str = BFEStats(Expi).units.unit_name_arr(BFEStats(Expi).units.pref_chan_id(1));
    
    figure(10);clf;hold on;
    plot(dist_prctile);plot(dist_pair_prctile)
    legend(["5%-ile in avg distmap", "5%-ile in distmap of all pairs"])
    title(["BigGAN FC6 Dual Evol Image Cmp", stimparts{end-1}])
    xlabel("Generation");ylabel("Patch Distance (Squeeze)")
    saveas(10,fullfile(figdir, compose("%s_Exp%02d_chan%s_Dist_prctile_traj.png",Animal,Expi,prefchan_str)))
    savefig(10,fullfile(figdir, compose("%s_Exp%02d_chan%s_Dist_prctile_traj.fig",Animal,Expi,prefchan_str)))
end
end
%%
generation_num_from_fdr("N:\Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55")
%%
function genN = generation_num_from_fdr(stimfolder)
allimgs = string(ls(stimfolder+"\block*thread*gen*"));
lastgeni = regexp(allimgs(end), "block(\d*)_thread\d*_gen", 'tokens');
lastgeni = str2num(lastgeni{1});
genN = lastgeni;
end

function [im_gen, imnames] = loadEvolImages(stimfolder, thread, block, suffix)
% Util function to load in evolved image of certain block x thread from a
% folder
% block: negative `block` variable means counting from the last block. e.g. -1 means
%        the last block,-2 means the penultimate block
% thread: Python numbering convention here, 0 for first thread, 1 for 2nd
%        thread
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

function [distmap_mean, distmap_col] = calc_distmap_avg(D, im_gen1, im_gen2, sampleN, res)
% Compute the image distance map between each pair of image from 2 cell
% arrays of images. 
if nargin <= 4,  res= 120;  end
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
