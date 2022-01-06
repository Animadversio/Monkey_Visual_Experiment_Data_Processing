%% vis Evol ColorFrame Montage for paper. Figure 2B 
%%
Set_Path;
mat_dir = "O:\Mat_Statistics"; 
%%
Animal = "Beto"; 
load(fullfile(mat_dir, Animal+"_Evol_stats.mat"))
load(fullfile(mat_dir, Animal+"_Manif_stats.mat"))
load(fullfile(mat_dir, Animal+"_Manif_stats.mat"))
%% Evolution stimuli
stimdir = EStats(11).meta.stimuli;
%% Calc mean scores.
actcol = cellfun(@(P)squeeze(mean(P(1,51:200,:),[1,2])),EStats(11).evol.psth,'Unif',0);
bslcol = cellfun(@(P)squeeze(mean(P(1,1:51,:),[1,2])),EStats(11).evol.psth,'Unif',0);
actmean = cellfun(@mean, actcol);
actsem = cellfun(@sem, actcol);
bslmean = mean(cat(1,bslcol{:}));
bslsem = sem(cat(1,bslcol{:}));
%% Load up all codes 
[codes_all, img_ids, generations] = load_codes_all(stimdir,1);
%% Get mean codes 
code_mean = arrayfun(@(geni)mean(codes_all(generations==geni,:),1),[min(generations):max(generations)]','uni',0);
code_mean = cat(1,code_mean{:});
%% Visualize the mean images 
G = FC6Generator();
code_mean_imgs = G.visualize(code_mean(:,:));
%% Export the Evolution Image Trajectory
outdir = "E:\OneDrive - Harvard University\Manuscript_Manifold\Figure2";
img_frame = score_frame_image_arr(code_mean_imgs(:,:,:,[1:4:80,80]),actmean([1:4:80,80]),[-0.1150   13.1639]);
mtg = imtile(img_frame, 'GridSize',[3,7],'Thumb',[296,296]);
imwrite(mtg,fullfile(outdir,"Beto_Exp11_EvolImgTraj.png"))
%% 
figure;
montage(img_frame,'Size',[3,7],'Thumb',[296,296])
%% 
mapactcol = cellfun(@(P)squeeze(mean(P(1,51:200,:),[1,2])),Stats(11).manif.psth{1},'Unif',0);
mapbslcol = cellfun(@(P)squeeze(mean(P(1,1:40,:),[1,2])),Stats(11).manif.psth{1},'Unif',0);
mapact_mean = cellfun(@mean,mapactcol);
mapbslmean = cellfun(@mean,mapbslcol);
mapbslmean_all = mean(cat(1,mapbslcol{:}));
mapact_mean - mapbslmean_all;
mapevkmean = mapact_mean - mapbslmean_all;
%% formatting map image name into array
mapstimdir = strrep(Stats(11).meta.stimuli,"\\storage1.ris.wustl.edu\crponce\Active","N:\");
mapimgnms = string(cellfun(@(idx)Stats(11).imageName{idx(1)},Stats(11).manif.idx_grid{1},'Unif',0));
mapimgfns = fullfile(mapstimdir,mapimgnms+".jpg");
%% Load and frame the images, and tile them accordingly
mapimg_frame = score_frame_image_arr(mapimgfns, mapevkmean);
mtg = imtile(mapimg_frame', 'GridSize',[11,11],'Thumb',[296,296]);
imwrite(mtg,fullfile(outdir,"Beto_Exp11_TuningMapMtg.png"))
%%
figh = figure('pos',[2565         248         480         420]);
imagesc(-90:18:90,-90:18:90,mapevkmean);
axis image
xticks(-90:18:90);
yticks(-90:18:90)
colorbar()
xlabel("direction 2");
ylabel("direction 1")
title("Beto Exp11 Tuning Map PrefChan 26")
saveallform(outdir,"Beto_Exp11_TuningMapAct",figh)