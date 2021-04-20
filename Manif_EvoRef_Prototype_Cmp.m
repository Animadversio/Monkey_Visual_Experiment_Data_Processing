%% Tightly connected to Manif_Imgs_Dissim_ExamplarImg.m, image version of that
%  Plot the best and worst images in different image spaces in a montage.
Animal = "Both"; Set_Path;
mat_dir = "O:\Mat_Statistics";
%%
figdir = "O:\Manif_Proto_Cmp";
Animal = "Alfa";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
pasu_idx_vec = reshape(Stats(12).ref.pasu_idx_grid',[],1); % reshape into one row. But transpose to make similar object adjacent
pasu_val_msk = ~cellfun(@isempty,pasu_idx_vec);
%% Find Examplar Images 
flag.doEvoRef = true;
P.topN = 3;
P.botN = 3;
for Expi = 34:numel(Stats)
manifdir = Stats(Expi).meta.stimuli;
evodir = EStats(Expi).meta.stimuli;
evorefdir = fileparts(EStats(Expi).meta.stimuli);
gabordir = "N:\Stimuli\2019-Manifold\gabor"; % N:\Stimuli\2019-Parameterized-Shapes\Gabors
pasudir = "N:\Stimuli\2019-Manifold\pasupathy-wg-f-4-ori"; % "N:\Stimuli\2019-Parameterized-Shapes\Pasupathy"
imgcol = {};
imgnmcol = {};
scorecol = [];
semcol = [];
stdcol = [];
titstr_col = [];

bsl_VEC_ALL = [];
% Manifold images. 
si=1; ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_sem_vec = cellfun(@(psth)sem(squeeze(mean(psth(ui,51:200,:),[1,2]))),reshape(Stats(Expi).manif.psth{si},[],1));
score_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).manif.psth{si},[],1),'uni',0)); % single trial vector
bsl_VEC_ALL = [bsl_VEC_ALL; score_bsl_VEC];
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
manifImgnm = cellfun(@(idx)string(Stats(Expi).imageName(idx(1))),Stats(Expi).manif.idx_grid{1});
fprintf("Top N Manifold image: ")
disp(manifImgnm(sortId(1:P.topN)))
fprintf("Bottom N Manifold image: " );
disp(manifImgnm(sortId(end-P.botN+1:end)))
suffix = findSuffix(manifImgnm(1), manifdir);
for ranki = [1:P.topN, numel(sortId)-P.botN+1:numel(sortId)]
	img = loadImgWsuffix(manifImgnm(sortId(ranki)), manifdir, suffix);
	imgcol{end+1} = img;
	imgnmcol{end+1} = manifImgnm(sortId(ranki));
	scorecol(end+1) = score_vec(sortId(ranki));
	semcol(end+1) = score_sem_vec(sortId(ranki));
	stdcol(end+1) = score_std_vec(sortId(ranki));
% 	titstr_col(end+1) = [];
end

if flag.doEvoRef
evoref_vec = cellfun(@(psth)mean(psth(1,51:200,:),'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_std_vec = cellfun(@(psth)std(mean(psth(1,51:200,:),[1,2]),1,'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_sem_vec = cellfun(@(psth)sem(squeeze(mean(psth(1,51:200,:),[1,2]))),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(1,1:45,:),[2])),reshape(EStats(Expi).ref.psth_arr,[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; evoref_bsl_VEC];
[evoref_sortScore,sortId] = sort(evoref_vec,'Descend');
[evoref_maxScore,evoref_maxId] = max(evoref_vec);
fprintf("Top N natural reference image: ")
disp(EStats(Expi).ref.imgnm(sortId(1:P.topN)))
fprintf("Bottom N natural reference image: " );
disp(EStats(Expi).ref.imgnm(sortId(end-P.botN+1:end)))
suffix = findSuffix(EStats(Expi).ref.imgnm(1), evorefdir);
for ranki = [1:P.topN, numel(sortId)-P.botN+1:numel(sortId)]
	img = loadImgWsuffix(EStats(Expi).ref.imgnm(sortId(ranki)), evorefdir);
	imgcol{end+1} = img;
	imgnmcol{end+1} = EStats(Expi).ref.imgnm(sortId(ranki));
	scorecol(end+1) = evoref_vec(sortId(ranki));
	semcol(end+1) = evoref_sem_vec(sortId(ranki));
	stdcol(end+1) = evoref_std_vec(sortId(ranki));
% 	titstr_col(end+1) = [];
end
end

% Pasupathy patches
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
nanmsk = ~pasu_val_msk|isnan(pasu_vec);
pasu_vec(nanmsk) = []; % isnan(pasu_vec) % get rid of non-existing pasupathy images. 
pasu_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
% pasu_sem_vec = cellfun(@(psth)sem(squeeze(mean(psth(ui,51:200,:),[1,2]))),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_sem_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all')/sqrt(size(psth,3)),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_std_vec(nanmsk) = [];
pasu_sem_vec(nanmsk) = [];
pasu_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).ref.pasu_psths',[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; pasu_bsl_VEC];
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
pasuimgnm = cellfun(@(idx)string(unique(Stats(Expi).imageName(idx))),reshape(Stats(Expi).ref.pasu_idx_grid',[],1),'uni',0);
pasuimgnm(nanmsk) = [];
pasuimgnm = string(pasuimgnm);

fprintf("Top N Pasupathy image: ")
disp(pasuimgnm(sortId(1:P.topN)))
fprintf("Bottom N Pasupathy image: " );
disp(pasuimgnm(sortId(end-P.botN+1:end)))
suffix = findSuffix(pasuimgnm(1), pasudir);
for ranki = [1:P.topN, numel(sortId)-P.botN+1:numel(sortId)]
	img = loadImgWsuffix(pasuimgnm(sortId(ranki)), pasudir, suffix);
	imgcol{end+1} = img;
	imgnmcol{end+1} = pasuimgnm(sortId(ranki));
	scorecol(end+1) = pasu_vec(sortId(ranki));
	semcol(end+1) = pasu_sem_vec(sortId(ranki));
	stdcol(end+1) = pasu_std_vec(sortId(ranki));
% 	titstr_col(end+1) = [];
end
end
% Gabor patches
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
nanmsk = isnan(gab_vec);
gab_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
% gab_sem_vec = cellfun(@(psth)sem(squeeze(mean(psth(ui,51:200,:),[1,2]))),reshape(Stats(Expi).ref.gab_psths,[],1));
gab_sem_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all')/sqrt(size(psth,3)),reshape(Stats(Expi).ref.gab_psths,[],1));
gab_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).ref.gab_psths,[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; gab_bsl_VEC];
gaborimgnm = cellfun(@(idx)unique(Stats(Expi).imageName(idx)),reshape(Stats(Expi).ref.gab_idx_grid,[],1),'Uni',0);
gaborimgnm(nanmsk) = [];
gaborimgnm = string(gaborimgnm); 
gab_vec(nanmsk) = [];
gab_std_vec(nanmsk) = [];
gab_sem_vec(nanmsk) = [];
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);

fprintf("Top N Gabor image: ")
disp(gaborimgnm(sortId(1:P.topN)))
fprintf("Bottom N Gabor image: " );
disp(gaborimgnm(sortId(end-P.botN+1:end)))
suffix = findSuffix(gaborimgnm(1), gabordir);
for ranki = [1:P.topN, numel(sortId)-P.botN+1:numel(sortId)]
	img = loadImgWsuffix(gaborimgnm(sortId(ranki)), gabordir, suffix);
	imgcol{end+1} = img;
	imgnmcol{end+1} = gaborimgnm(sortId(ranki));
	scorecol(end+1) = gab_vec(sortId(ranki));
	semcol(end+1) = gab_sem_vec(sortId(ranki));
	stdcol(end+1) = gab_std_vec(sortId(ranki));
% 	titstr_col(end+1) = [];
end
end
bsl_rate = nanmean(bsl_VEC_ALL);
bsl_rate_std = nanstd(bsl_VEC_ALL,1);
bsl_rate_sem = sem(bsl_VEC_ALL,1);
figure(1);
montage(imgcol, 'Size', [4, P.botN+P.topN], 'BorderSize', 5)
figure(2);
T = montageWscore(imgcol, [4, P.botN+P.topN], imgnmcol, scorecol, stdcol, semcol);
savenm = compose("%s_Exp%02d_pref%02d_prototypes",Animal,Expi,Stats(Expi).units.pref_chan);
saveas(1,fullfile(figdir,savenm+".png"))
saveas(1,fullfile(figdir,savenm+".pdf"))
savefig(1,fullfile(figdir,savenm+".fig"))
saveas(2,fullfile(figdir,savenm+"_labels.png"))
saveas(2,fullfile(figdir,savenm+"_labels.pdf"))
savefig(2,fullfile(figdir,savenm+"_labels.fig"))

end
%%

function suffix = findSuffix(name, folder)
% suffix = findSuffix("gab_ori_0.0_0.5_thread000_nat", fileparts(EStats(Expi).meta.stimuli));
nm = ls(fullfile(folder,string(name+"*")));
parts = split(nm,'.');
assert(numel(nm) > 0)
suffix = string(parts{end});
end

function img = loadImgWsuffix(name, folder, suffix)
% img = loadImgWsuffix("gab_ori_0.0_0.5_thread000_nat", fileparts(EStats(Expi).meta.stimuli), "bmp");
if nargin == 2
suffix = findSuffix(name, folder);
end
imgnm = fullfile(folder, name+"."+suffix);
img = imread(imgnm);
end

function T = montageWscore(imgcol, tilesize, imgnmcol, scorecol, stdcol, semcol)
H = tilesize(1); W = tilesize(2);
T = tiledlayout(H, W, 'TileSpacing', 'compact', 'Padding', 'compact');
for imgi = 1:numel(imgcol)
ax = nexttile(T, imgi);
imshow(imgcol{imgi})
if contains(imgnmcol{imgi}, "thread00")
    idx = strfind(imgnmcol{imgi}, "thread00");
    imgnmcol{imgi} = extractBetween(imgnmcol{imgi},1,idx-2);
end
titstr = compose("%s\n%.1f(%.1f)",imgnmcol{imgi},scorecol(imgi),semcol(imgi));
title(titstr, 'Interpreter', 'none')
end
end