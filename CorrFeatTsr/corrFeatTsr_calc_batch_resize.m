% Using the `corrFeatTsr_func` to compute the correlation tensor and store basic statistics into a stats mat!
%  Resize image to get rid of artifaces
net = vgg16;
%%
savedir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
hier_savedir = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_Hierarchy";
%%
Animal="Alfa"; 
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
% Animal="Alfa"; ExpType="Manif"; 

for Animal=["Alfa"]%"Alfa",
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')
% flags = struct("batch",50,"online_compute",0,"load_all_img",1,"shuffleN",100,"resize",true);
ExpType="Evol"; 
flags = struct("batch",60,"online_compute",1,"load_all_img",1,"shuffleN",50,"resize",true);
wdw_vect = [[1, 20] + 10 * [0:18]'; [1,50]+[0:50:150]'; [51,200]];
for Expi = 2:numel(EStats)
imgsize = EStats(Expi).evol.imgsize;
imgpix = imgsize * 40;
flags.imgpix = imgpix;
T00 = tic;
[imgfn, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, flags);
toc(T00)
%%
layernames = ["fc8", "fc7", "fc6", "conv5_3", "conv4_3", "conv3_3"];
for iLayer = 1:length(layernames)
layername = layernames(iLayer);
[cc_tsr, MFeat, StdFeat, cc_refM, cc_refS] = corrFeatTsr_func(imgfn, score_vect, net, layername, flags);
outfn = fullfile(savedir, compose("%s_%s_Exp%d_%s_resize.mat",Animal,ExpType,Expi,layername)); % LW means long window
save(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
end
end
end


%%
Animal = "Beto"; 
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(mat_dir, compose("%s_ImageRepr.mat", Animal)), 'ReprStats')

function dlimg = loadimges(images)
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(imgnm)),[224,224]), images, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[])); 
end

function [imgfullnm_vect, score_vect] = loadManifData(Stats, EStats, Animal, ExpType, Expi, flags)
ui=1;si=1;
assert(EStats(Expi).Animal == Animal && Stats(Expi).Animal == Animal)
fprintf("Processing %s Exp %d pref chan %d\n",ExpType,Expi,EStats(Expi).units.pref_chan)
if ExpType == "Manif"
imgnm_grid = string(cellfun(@(idx) unique(Stats(Expi).imageName(idx)), Stats(Expi).manif.idx_grid{si}));
imgnm_vect = reshape(imgnm_grid, [], 1);
stimpath = Stats(Expi).meta.stimuli;
% imgN=121; 
elseif ExpType == "Evol"
index_vect = cell2mat(EStats(Expi).evol.idx_seq');
imgnm_vect = EStats(Expi).imageName(index_vect); % reshape(imgnm_grid, [], 1);
stimpath = EStats(Expi).meta.stimuli;
stimpath = strrep(stimpath,"N:\Stimuli\2019-Manifold", "E:\Network_Data_Sync\Stimuli\2019-Manifold");
end
imgN = length(imgnm_vect);
tmpfn = ls(fullfile(stimpath, imgnm_vect(1)+"*"));
tmpparts = split(tmpfn,".");suffix = "."+tmpparts{2};
imgfullnm_vect = cellfun(@(imgnm) fullfile(stimpath, imgnm+suffix),imgnm_vect);
% Compute the score(firing rate at different time slices) 
% the time window info in recorded in `wdw_vect` and saved to disk. 
if ExpType == "Manif"
psth_all = cellfun(@(psth) reshape(mean(psth(ui,:,:),[3]),1,1,[]), Stats(Expi).manif.psth{si}, ...
    'UniformOutput', false); % note there is trial averaging here. 
psth_all = reshape(cell2mat(psth_all),imgN,[]);
elseif ExpType == "Evol"
psth_all = squeeze(cell2mat(reshape(EStats(Expi).evol.psth,1,1,[])))'; % imgN by 200
end
% This part could be abbrieviated. 
score_vect = movmean(psth_all,20,2,'Endpoints','discard'); % short time window 20ms average 
score_vect = score_vect(:, 1:10:end); % subsample to decrease redunancy
wdw_vect = [1, 20] + 10 * [0:18]';
score_vect = [score_vect, mean(psth_all(:,1:50),2),mean(psth_all(:,51:100),2),...
    mean(psth_all(:,101:150),2),mean(psth_all(:,151:200),2),mean(psth_all(:,51:200),2)]; % [Trials, nTimeWindows]
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]]; % [nTimeWindows, 2]
end