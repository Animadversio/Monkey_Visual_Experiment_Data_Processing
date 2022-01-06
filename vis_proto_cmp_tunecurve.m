% Make tuning curve figure and exemplar prototype simplified compare figures
%% Prep data collection
Animal="Both";Set_Path;
py.importlib.import_module("pickle");
for Animal=["Alfa","Beto"]
load(fullfile(mat_dir, Animal+"_ImageRepr.mat"), 'ReprStats');
load(fullfile(mat_dir, Animal+"_Manif_Stats.mat"), 'Stats');
load(fullfile(mat_dir, Animal+"_Evol_Stats.mat"), 'EStats');
load(fullfile(mat_dir, Animal+"_ManifMapVarStats.mat"),'MapVarStats')
MapVarStatsCmb.(Animal) = MapVarStats;
ReprStatsCmb.(Animal) = ReprStats;
StatsCmb.(Animal) = Stats;
EStatsCmb.(Animal) = EStats;
end
%% Fetch prototype image 
modeldir = "O:\corrFeatTsr_FactorVis\models\resnet50_linf8-layer3_NF3_bdr1_Tthresh_3__nobdr_res-robust_CV";
Explist =  {{"Alfa", 36},
			{"Alfa", 34},
			{"Alfa", 1},
			{"Beto", 34},
			{"Beto", 22},
			{"Beto", 27}};
protocol = [];
for k = 1:numel(Explist)
expname = Explist{k};
Animal = expname{1};Expi = expname{2};
protocol(k).manifproto = ReprStatsCmb.(Animal)(Expi).Manif.BestImg;
protocol(k).evolproto = ReprStatsCmb.(Animal)(Expi).Evol.BestImg;
factordata = py.pickle.load(py.open(fullfile(modeldir,compose("%s_Exp%02d_factors.pkl",Animal,Expi)),'rb'));
Hmaps = factordata.Hmaps.double;
ccfactor = factordata.ccfactor.double;
Hmaps_renorm = Hmaps .* reshape(vecnorm(ccfactor,2,1),1,1,[]);
Hmaps_maxnorm = Hmaps_renorm / prctile(Hmaps_renorm,99,'all');
protocol(k).Hmaps_maxnorm = Hmaps_maxnorm;
protocol(k).facttsr_img = factordata.get("tsr_proto").single; 
end
% factordata = py.pickle.load(py.open(fullfile(modeldir,compose("%s_Exp%02d_factors.pkl",Animal,Expi)),'rb'));
save_proto_tiles({protocol.evolproto},"evolproto")
save_proto_tiles({protocol.manifproto},"manifproto")
save_proto_tiles({protocol.Hmaps_maxnorm},"factmodel_Hmaps_")
save_proto_tiles({protocol.facttsr_img},"factmodel_proto")

%%

%% Machinery to plot 1-d Tuning curves
figdir = "O:\Manuscript_Manifold\Figure3\TuningCurve";
Explist =  {{"Alfa", 20}};%,
%             {"Alfa", 4},
% 			{"Alfa", 6},
% 			{"Alfa", 33},
% 			{"Alfa", 26},
% 			{"Beto", 6},
% 			{"Beto", 7},
% 			{"Beto", 20},
% 			{"Beto", 29}};
for k = 1:numel(Explist)
expname = Explist{k};
Animal = expname{1};Expi = expname{2};
% Animal = "Beto";Expi=11;
[manif_imgfps, mapper] = map2fullpath({}, StatsCmb.(Animal)(Expi).meta.stimuli);
PC12idx_grid = StatsCmb.(Animal)(Expi).manif.idx_grid{1};
imageName = StatsCmb.(Animal)(Expi).imageName;
manifimgnms = cellfun(@(I)imageName(I(1)),PC12idx_grid);
manifimgfps = cellfun(@(C)string(mapper(C)),manifimgnms);
manifimgs = cellfun(@(fp)imread(fp),manifimgfps,'Uni',0);
imgpos = EStatsCmb.(Animal)(Expi).evol.imgpos;
imsize_deg = EStatsCmb.(Animal)(Expi).evol.imgsize;
unit_name_arr = StatsCmb.(Animal)(Expi).units.unit_name_arr;
iCh = MapVarStatsCmb.(Animal)(Expi).units.pref_chan_id(1);
explabel = compose("PrefChan %s img %.1f deg pos [%.1f %.1f]",unit_name_arr(iCh),...
    imsize_deg,imgpos(1),imgpos(2));
act_col = MapVarStatsCmb.(Animal)(Expi).manif.act_col{1};
actmap_mean = cellfun(@(A)mean(A(iCh,:),'all'),act_col,'uni',1);
[maxScore,maxIdx] = max(actmap_mean,[],'all','linear');
[row,col] = ind2sub(size(actmap_mean),maxIdx);%row=6;col=6;
assert(actmap_mean(row,col)==maxScore)
Clim = prctile(actmap_mean,[0,100],'all')';
[figh,T] = TuneCurvePlot(manifimgs,actmap_mean,[row,col],Clim);
title(T,compose("%s Exp%02d %s Best image name: %s",Animal,Expi,explabel,manifimgnms{row,col}),'interp','none')
saveallform(figdir,compose("%s_Exp%02d_BestTuneCurv",Animal,Expi),figh,["pdf","png"])
[figh,T] = TuneCurvePlot(manifimgs,actmap_mean,[6,6],Clim);
title(T,compose("%s Exp%02d %s Cent image name: %s",Animal,Expi,explabel,manifimgnms{6,6}),'interp','none')
saveallform(figdir,compose("%s_Exp%02d_CentTuneCurv",Animal,Expi),figh,["pdf","png"])
% figh = figure('pos',[207   512    2200   500]);
% T=tiledlayout(2,1,'tilesp','compact','pad','none');
% frame_img_list = score_frame_image_arr(manifimgs(row,:),actmap_mean(row,:),Clim); % fix PC2, move PC3
% nexttile(1);imshow(imtile(frame_img_list,'GridSize',[1,11]))
% frame_img_list = score_frame_image_arr(manifimgs(:,col),actmap_mean(:,col),Clim); % fix PC3, move PC2
% nexttile(2);imshow(imtile(frame_img_list,'GridSize',[1,11]))
% caxis(Clim);colorbar;
end
%%
%% Copy representative data 
% imgcol = {};
% for expname = Explist'
% Animal = expname{1}{1};Expi = expname{1}{2};
% copyfile(fullfile(modeldir,compose("%s_Exp%02d_summary.png",Animal,Expi)),...
%      fullfile(figdir,compose("%s_Exp%02d_summary.png",Animal,Expi)))
% end
%%
function [figh,T] = TuneCurvePlot(manifimgs,actmap_mean,cent,Clim)
row = cent(1);col = cent(2);
figh = figure('pos',[207   512    2200   500]);
T=tiledlayout(2,1,'tilesp','compact','pad','none');
frame_img_list = score_frame_image_arr(manifimgs(row,:),actmap_mean(row,:),Clim); % fix PC2, move PC3
nexttile(1);imshow(imtile(frame_img_list,'GridSize',[1,11]))
frame_img_list = score_frame_image_arr(manifimgs(:,col),actmap_mean(:,col),Clim); % fix PC3, move PC2
nexttile(2);imshow(imtile(frame_img_list,'GridSize',[1,11]))
caxis(Clim);colorbar;
end
function save_proto_tiles(imgcol,savenm)
figdir = "O:\Manuscript_Manifold\Figure3\protoProg";
h=figure('pos', [1000         616         560         362]);
T=tiledlayout(2,3,'tilesp','compact','pad','none');
for i=1:6
   nexttile(T,i); imshow(imgcol{i}) 
end
saveallform(figdir,compose("%s_collection",savenm),h,["pdf","png"])
end