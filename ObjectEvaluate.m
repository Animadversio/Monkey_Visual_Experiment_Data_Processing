ObjPath = "D:\DL_Projects\Vision\object-proposals-master\objectness-release-v2.2";
cd(ObjPath)
%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
savedir = "E:\OneDrive - Washington University in St. Louis\Repr_with_Mask";

%%
load(fullfile(savedir, 'summary',"Both_MaskCol.mat"),'MasksCol')
ObjMapStats = struct();
%%
tic
thr = 0.7;
for Animal=["Alfa","Beto"]
load(fullfile(mat_dir, Animal+"_ImageRepr.mat"), 'ReprStats');
load(fullfile(mat_dir, Animal+"_Evol_Stats.mat"), 'EStats');
load(fullfile(mat_dir, Animal+"_Evol_ccFtMask.mat"),'ccMskStats')
for Expi = 1:46
    imgsize = EStats(Expi).evol.imgsize;
    imsize_pix = imgsize * 40;
    img = ReprStats(Expi).Evol.BestBlockAvgImg;
    img_rsz = imresize(img, [imsize_pix, imsize_pix]);
    boxes = runObjectness(img,100);
    [objHeatMap, rawmap] = computeObjectnessHeatMap(img,boxes);
    boxes_rsz = runObjectness(img_rsz,100);
    [objHeatMap_rsz, rawmap_rsz] = computeObjectnessHeatMap(img_rsz,boxes_rsz);
    ObjMapStats.(Animal)(Expi).evolblck.rawmap = rawmap;
    ObjMapStats.(Animal)(Expi).evolblck.rawmap = rawmap/100;
    ObjMapStats.(Animal)(Expi).evolblck.heatmap = objHeatMap;
    ObjMapStats.(Animal)(Expi).evolblck.boxes = boxes;
    ObjMapStats.(Animal)(Expi).evolblck.rawmap_rsz = rawmap_rsz;
    ObjMapStats.(Animal)(Expi).evolblck.rawmap_rsz = rawmap_rsz/100;
    ObjMapStats.(Animal)(Expi).evolblck.heatmap_rsz = objHeatMap_rsz;
    ObjMapStats.(Animal)(Expi).evolblck.boxes_rsz = boxes_rsz;
    
    [objHeatMap, rawmap] = computeObjectnessHeatMap(img,boxes(boxes(:,5)>thr,:));
    [objHeatMap_rsz, rawmap_rsz] = computeObjectnessHeatMap(img_rsz,boxes_rsz(boxes_rsz(:,5)>thr,:));
    ObjMapStats.(Animal)(Expi).evolblck.rawmap_thr = rawmap;
    ObjMapStats.(Animal)(Expi).evolblck.rawmap_thr_norm = rawmap/(sum(boxes(:,5)>thr));
    ObjMapStats.(Animal)(Expi).evolblck.heatmap_thr = objHeatMap;
    ObjMapStats.(Animal)(Expi).evolblck.rawmap_rsz_thr = rawmap_rsz;
    ObjMapStats.(Animal)(Expi).evolblck.rawmap_rsz_thr_norm = rawmap_rsz/(sum(boxes_rsz(:,5)>thr));
    ObjMapStats.(Animal)(Expi).evolblck.heatmap_rsz_thr = objHeatMap_rsz;
    toc
end
end
%%
save(fullfile(mat_dir, "Both"+"_ObjMaskStats.mat"),'ObjMapStats')

%%
thr = 0.5;
tic
for Expi = 1:46
    imgsize = EStats(Expi).evol.imgsize;
    imsize_pix = imgsize * 40;
    img = ReprStats(Expi).Evol.BestBlockAvgImg;
    img_rsz = imresize(img, [imsize_pix, imsize_pix]);
    boxes = ObjMapStats.(Animal)(Expi).evolblck.boxes;
    boxes_rsz = ObjMapStats.(Animal)(Expi).evolblck.boxes_rsz;
    [objHeatMap, rawmap] = computeObjectnessHeatMap(img,boxes(boxes(:,5)>thr,:));
    [objHeatMap_rsz, rawmap_rsz] = computeObjectnessHeatMap(img_rsz,boxes_rsz(boxes_rsz(:,5)>thr,:));
    ObjMapStats.(Animal)(Expi).evolblck.rawmap_thr = rawmap;
    ObjMapStats.(Animal)(Expi).evolblck.heatmap_thr = objHeatMap;
    ObjMapStats.(Animal)(Expi).evolblck.rawmap_rsz_thr = rawmap_rsz;
    ObjMapStats.(Animal)(Expi).evolblck.heatmap_rsz_thr = objHeatMap_rsz;
    toc
end
%%
Expi = 3;
img = ReprStats(Expi).Evol.BestBlockAvgImg;
% img_rsz = imresize(img, [imsize_pix, imsize_pix]);
boxes = runObjectness(img,200);
[objHeatMap, rawmap] = computeObjectnessHeatMap(img,boxes);
%%
tab = [];
P=[];P.plotMtg = false;
layername="conv5_3";CMINPrct = 20; CMAXPrct = 75; ALIM = [0, 1];
for Animal=["Alfa","Beto"]
load(fullfile(mat_dir, Animal+"_ImageRepr.mat"), 'ReprStats');
load(fullfile(mat_dir, Animal+"_Evol_Stats.mat"), 'EStats');
load(fullfile(mat_dir, Animal+"_Evol_ccFtMask.mat"),'ccMskStats')
for Expi = 1:numel(EStats)
img = ReprStats(Expi).Evol.BestBlockAvgImg;
imgsize = EStats(Expi).evol.imgsize;
imsize_pix = imgsize * 40;
imgpos = EStats(Expi).evol.imgpos;
prefchan = EStats(Expi).evol.pref_chan;
interpMsk = MasksCol.(Animal)(Expi).(layername).interp;
CLIM = prctile(interpMsk,[CMINPrct,CMAXPrct],'all');
alphaMsk = clip((interpMsk - CLIM(1)) ./ (CLIM(2) - CLIM(1)), 0, 1) * (ALIM(2) - ALIM(1))+ALIM(1);
alphaMsk_rsz = imresize(alphaMsk, [imsize_pix, imsize_pix]);
objHeatMap = ObjMapStats.(Animal)(Expi).evolblck.heatmap_thr;
objHeatMap_rsz = ObjMapStats.(Animal)(Expi).evolblck.heatmap_rsz_thr;
boxes = ObjMapStats.(Animal)(Expi).evolblck.boxes;
boxes_rsz = ObjMapStats.(Animal)(Expi).evolblck.boxes_rsz;
rawmap = ObjMapStats.(Animal)(Expi).evolblck.rawmap_thr;
rawmap_rsz = ObjMapStats.(Animal)(Expi).evolblck.rawmap_rsz_thr;
rawmap_nrm = ObjMapStats.(Animal)(Expi).evolblck.rawmap_thr_norm;
rawmap_rsz_nrm = ObjMapStats.(Animal)(Expi).evolblck.rawmap_rsz_thr_norm;

alphaMsk_comp = 1 - alphaMsk;
% alphaMsk_comp(isnan(alphaMsk_comp)) = 1;
alphaMsk_rsz_comp = 1 - alphaMsk_rsz;
% alphaMsk_rsz_comp(isnan(alphaMsk_rsz_comp)) = 1;
obj_inmsk = nansum(rawmap .* alphaMsk,'all')/nansum(alphaMsk,'all');
obj_outmsk = nansum(rawmap .* alphaMsk_comp,'all')/nansum(alphaMsk_comp,'all');
obj_inmsk_rsz = nansum(rawmap_rsz .* alphaMsk_rsz,'all')/nansum(alphaMsk_rsz,'all');
obj_outmsk_rsz = nansum(rawmap_rsz .* alphaMsk_rsz_comp,'all')/nansum(alphaMsk_rsz_comp,'all');
obj_inmsk_nrm = nansum(rawmap_nrm .* alphaMsk,'all')/nansum(alphaMsk,'all');
obj_outmsk_nrm = nansum(rawmap_nrm .* alphaMsk_comp,'all')/nansum(alphaMsk_comp,'all');
obj_inmsk_rsz_nrm = nansum(rawmap_rsz_nrm .* alphaMsk_rsz,'all')/nansum(alphaMsk_rsz,'all');
obj_outmsk_rsz_nrm = nansum(rawmap_rsz_nrm .* alphaMsk_rsz_comp,'all')/nansum(alphaMsk_rsz_comp,'all');
if P.plotMtg
figure(1);
montage({img,alphaMsk,objHeatMap,objHeatMap_rsz})
title(compose("Exp %d prefchan %d %.1f deg [%.1f %.1f]\nObjectiveness: "+...
	"In msk %.3f Out msk %.3f; Resize In msk %.3f Out msk %.3f \n"+...
	"Normalized Objectiveness: In msk %.3f Out msk %.3f; Resize In msk %.3f Out msk %.3f ",...
	Expi,prefchan,imgsize,imgpos(1),imgpos(2),...
	obj_inmsk, obj_outmsk, obj_inmsk_rsz, obj_outmsk_rsz,...
	obj_inmsk_nrm, obj_outmsk_nrm, obj_inmsk_rsz_nrm, obj_outmsk_rsz_nrm))
pause%(0.2)
else
fprintf(compose("Exp %d prefchan %d %.1f deg [%.1f %.1f]\nObjectiveness: "+...
	"In msk %.3f Out msk %.3f; Resize In msk %.3f Out msk %.3f \n"+...
	"Normalized Objectiveness: In msk %.3f Out msk %.3f; Resize In msk %.3f Out msk %.3f \n",...
	Expi,prefchan,imgsize,imgpos(1),imgpos(2),...
	obj_inmsk, obj_outmsk, obj_inmsk_rsz, obj_outmsk_rsz,...
	obj_inmsk_nrm, obj_outmsk_nrm, obj_inmsk_rsz_nrm, obj_outmsk_rsz_nrm))  
T = struct();
for varnm = ["Animal", "Expi", "prefchan", "imgsize", "imgpos", ...
	"obj_inmsk", "obj_outmsk", "obj_inmsk_rsz", "obj_outmsk_rsz", ...
	"obj_inmsk_nrm", "obj_outmsk_nrm", "obj_inmsk_rsz_nrm", "obj_outmsk_rsz_nrm"]
	T.(varnm) = eval(varnm);
end
T.obj_box_max = max(boxes(:,5));
T.obj_box_rsz_max = max(boxes_rsz(:,5));

tab = [tab;T];
end
end
end
tab = struct2table(tab);
%%
V1msk = tab.prefchan<49&tab.prefchan>32;
V4msk = tab.prefchan>=49&tab.prefchan>32;
ITmsk = tab.prefchan<49&tab.prefchan<=32;
szmsk = tab.imgsize>1.1;
%%
Tsummary(tab,{V1msk,V4msk,ITmsk},["V1","V4","IT"],"obj_inmsk_rsz_nrm","obj_outmsk_rsz_nrm")
T2summary(tab,{V1msk,V4msk,ITmsk},["V1","V4","IT"],{[3,2],[3,1],[2,1]},"obj_inmsk_rsz_nrm")
%%
Tsummary(tab,{V1msk&szmsk,V4msk&szmsk,ITmsk&szmsk},["V1","V4","IT"],"obj_inmsk_rsz_nrm","obj_outmsk_rsz_nrm")
T2summary(tab,{V1msk&szmsk,V4msk&szmsk,ITmsk&szmsk},["V1","V4","IT"],{[3,2],[3,1],[2,1]},"obj_inmsk_rsz_nrm")
%%

Tsummary(tab,{V1msk,V4msk,ITmsk},["V1","V4","IT"],"obj_box_max","obj_outmsk_rsz_nrm")
T2summary(tab,{V1msk,V4msk,ITmsk},["V1","V4","IT"],{[3,2],[3,1],[2,1]},"obj_box_max")
%%
alphaMsk_comp = 1 - alphaMsk;
alphaMsk_comp(isnan(alphaMsk_comp)) = 1;
% alphaMsk_rsz_comp = 1 - alphaMsk_rsz;
% alphaMsk_rsz_comp(isnan(alphaMsk_rsz_comp)) = 1;
obj_inmsk = nansum(rawmap .* alphaMsk,'all')/nansum(alphaMsk,'all');
obj_outmsk = nansum(rawmap .* alphaMsk_comp,'all')/nansum(alphaMsk_comp,'all');

function Tsummary(tab,msks,labels,varnm1,varnm2)
for m = 1:numel(msks)
[~,P,CI,TST] = ttest(tab.(varnm1)(msks{m}),tab.(varnm2)(msks{m}));
fprintf("%s: %s - %s: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",labels(m),varnm1,varnm2,P,TST.tstat,TST.df,CI(1),CI(2));
end
end
function T2summary(tab,msks,labels,pairs,varnm)
for p = 1:numel(pairs)
i = pairs{p}(1); j = pairs{p}(2);
[~,P,CI,TST] = ttest2(tab.(varnm)(msks{i}),tab.(varnm)(msks{j}));
fprintf("%s - %s: P=%.1e t=%.3f(df=%d),CI=[%.1f,%.1f]\n",labels(i),labels(j),P,TST.tstat,TST.df,CI(1),CI(2));
end
end