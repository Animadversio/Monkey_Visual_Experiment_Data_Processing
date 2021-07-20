%% Code to create mask and merge on images from the masks stored in the correlation. 
%  Dup from corrFeatTsr2mask and majorly used to creat montage for all experiment. 

ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto"; Expi = 29;  %thresh = "neg"; sum_method="max"; %"both"
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')
ccMskStats = repmat(struct(),1,length(EStats));
%% Review the ccFeatTsr with respect to the evolution of the psth and images. 
%  Actually this is the newer version of the visualization function
ExpType = "Manif";
doPlot = true; doSave = true;
thresh = "both"; 
for Expi = 40:45
ccMskStats(Expi).(ExpType).units = EStats.units;
% winopen(fullfile(moviedir, compose("%s_Evol_Exp%02d_Best_PSTH.mov.avi",Animal,Expi)))
% winopen(fullfile(moviedir, compose("%s_Manif_Exp%02d_Avg_PSTH.mov.avi",Animal,Expi)))
savedir = fullfile(ccmat_dir,compose("%s_%s_Exp%d",Animal,ExpType,Expi));

for layername = ["conv3_3", "conv4_3", "conv5_3"] %["conv5_3"]%"conv5_3";
% if layername=="conv3_3" && Expi ==40,continue;end
outfn = fullfile(ccmat_dir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername));
load(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
%%
for sum_method=["max","L1","L1signif"]
if thresh == "both"
maskTsr = abs(t_signif_tsr)>=5; % Bi sided thresholding 
elseif thresh == "pos"
maskTsr = t_signif_tsr>=5; % positively correlated thresholding
elseif thresh == "neg"
maskTsr = t_signif_tsr<=-5; % negatively correlated thresholding
end
plotTsr = cc_tsr;
plotTsr(~maskTsr) = 0; 
if sum_method=="L1"
L1plotTsr = squeeze(mean(abs(plotTsr),3));
elseif sum_method=="L1signif"
L1plotTsr = squeeze(sum(abs(plotTsr),3)./(sum(maskTsr,3)));
elseif sum_method=="max"
L1plotTsr = squeeze(max(abs(plotTsr),[],3));
end
%% Calculate the CLIM for each time bin (keep it the same for bins of same length for cmp)
if doPlot
CLIM_arr = zeros(24,2);
CLIM_arr(1:19,:) = CLIM_arr(1:19,:) + prctile(L1plotTsr(:,:,1:19),[2,98],'all')' + [0, 1E-4];
CLIM_arr(20:23,:) = CLIM_arr(20:23,:) + prctile(L1plotTsr(:,:,20:23),[2,98],'all')'+ [0, 1E-4];
CLIM_arr(24,:) = CLIM_arr(24,:) + prctile(L1plotTsr(:,:,24),[2,98],'all')'+ [0, 1E-4];
CLIM_arr(isnan(CLIM_arr)) = 0; % put 0 in the nan place (if there is no activated voxel. then prctile will be all nan)

if doSave
v = VideoWriter(fullfile(savedir,compose("%s_%s_Exp%d_%s_cc_%s-%s.mp4",Animal,ExpType,Expi,layername,thresh,sum_method)));
v.FrameRate = 2;open(v);
end
figure(13);%(12);
% imagesc(mean(abs(cc_tsr(:,:,:,end)),3));colorbar()
for fi = 1:24
imagesc(L1plotTsr(:,:,fi));axis image;caxis(CLIM_arr(fi,:));colorbar();
title(compose("%s %s Exp%d Corr Coef with VGG16 %s\n rate in [%d,%d]ms\n Channel compressed with %s thresholded %s",...
    Animal,ExpType,Expi,layername,wdw_vect(fi,1),wdw_vect(fi,2),thresh,sum_method))
pause(0.2)
if doSave,writeVideo(v,getframe(13));end
drawnow;
end
if doSave,close(v);end
end
ccMskStats(Expi).(ExpType).(layername).(sum_method) = L1plotTsr(:,:,end);
end
end
end
%%
save(fullfile(mat_dir,compose("%s_Evol_ccFtMask.mat",Animal)),'ccMskStats')
%% Reload the computed mask and create mask from it. 
net = vgg16; % load the net to compute the receptive field structure. 
CNNRF = CNN_receptive_field(net);
%% 
Animal = "Beto";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')
load(fullfile(mat_dir, compose("%s_Evol_ccFtMask.mat",Animal)),'ccMskStats')
load(fullfile(mat_dir, compose("%s_ImageRepr.mat",Animal)),'ReprStats')
%% Interpolate and Visualize masks Together with images.
global savedir
savedir = "E:\OneDrive - Washington University in St. Louis\Repr_with_Mask";
mkdir(fullfile(savedir,"summary"))
%% Create mask collection of both monk using the parameters
MasksCol = struct();
P.sum_method = "L1"; % Method to sum up the feature dim
P.sup_boundary = true; % options and params for suppressing boundary activation
P.cut_bdr = struct("conv3_3",3,"conv4_3",2,"conv5_3",1); % how many border elements to cut!
for Animal = ["Alfa","Beto"]
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')
load(fullfile(mat_dir, compose("%s_Evol_ccFtMask.mat",Animal)),'ccMskStats')
for Expi = 1:numel(EStats)
for layername = ["conv3_3","conv4_3","conv5_3"]
MskMap = ccMskStats(Expi).(ExpType).(layername).(P.sum_method);
idx = find(contains(CNNRF.Name,layername));
startpix = CNNRF.start(idx,:);
jumppix = CNNRF.jump(idx,:);
H = size(MskMap,1);W = size(MskMap,2);
[centGridX, centGridY] = meshgrid(startpix(1)+jumppix(1)*(0:H-1),startpix(2)+jumppix(2)*(0:W-1));
if P.sup_boundary 
bdr = P.cut_bdr.(layername);
padMap = nan(size(MskMap));
padMap(bdr+1:end-bdr,bdr+1:end-bdr) = MskMap(bdr+1:end-bdr,bdr+1:end-bdr);
else
padMap = MskMap;
end
% figure(6);
% imagesc(padMap);
% axis image;colorbar();caxis(CLIM)
[interpX, interpY] = meshgrid(linspace(1,224,256), linspace(1,224,256));
interpMsk = griddata(centGridX, centGridY, double(padMap), interpX, interpY);
MasksCol.(Animal)(Expi).(layername).interp = interpMsk;
% MasksCol.(Animal)(Expi).(layername).alpha = alphaMsk;
end
end
end
%%
save(fullfile(savedir,'summary',"Both_MaskCol.mat"),'MasksCol')
%% Separate by animal 
savedir = "E:\OneDrive - Washington University in St. Louis\Repr_with_Mask";
load(fullfile(savedir,'summary',"Both_MaskCol.mat"),'MasksCol')
for Animal =[ "Alfa", "Beto"]
D = load(fullfile(mat_dir, Animal+"_Evol_stats.mat"), 'EStats');
EStats_CMB.(Animal) = D.EStats;
D = load(fullfile(mat_dir, Animal+"_ImageRepr.mat"),'ReprStats');
ReprStats_CMB.(Animal) = D.ReprStats;
end
%% Final parameters to generate mask and tune the alpha value
layername="conv4_3";CMINPrct = 20; CMAXPrct = 75; ALIM = [0.1, 1];
%%
figure(13);clf;set(13,'pos',[335          42        1875         955])
T = tiledlayout(5,10,'Padd','compact','TileSp','none');
figure(14);clf;set(14,'pos',[335          42        1875         955])
T2 = tiledlayout(5,10,'Padd','compact','TileSp','compact');
for Expi = 1:numel(EStats_CMB.(Animal))
prefchan_lab = EStats_CMB.(Animal)(Expi).units.unit_name_arr(EStats_CMB.(Animal)(Expi).units.pref_chan_id);
ReprImg = ReprStats_CMB.(Animal)(Expi).Evol.BestBlockAvgImg; % uint8 format
nexttile(T, Expi)
imshow(ReprImg)
title(compose("Exp%02d, PrefChan %s",Expi,prefchan_lab))
nexttile(T2, Expi)
interpMsk = MasksCol.(Animal)(Expi).(layername).interp;
CLIM = prctile(interpMsk,[CMINPrct,CMAXPrct],'all');
alphaMsk = clip((interpMsk - CLIM(1)) ./ (CLIM(2) - CLIM(1)), 0, 1) * (ALIM(2) - ALIM(1))+ALIM(1);
imshow(double(ReprImg)/255.*alphaMsk)
title(compose("Exp%02d, PrefChan %s",Expi,prefchan_lab))
end
title(T,compose("%s BestBlockAvg",Animal))
title(T2,compose("%s BestBlockAvg Saliency Masked",Animal))
%%
Animals = ["Alfa", "Beto"];
for Ai = 1:2
h(Ai) = figure(Ai);set(h(Ai),'pos',[335          42        1875         955])
T = tiledlayout(h(Ai),5,10,'Padd','compact','TileSp','none');
Animal = Animals(Ai);
for Expi = 1:numel(EStats_CMB.(Animal))
prefchan_lab = EStats_CMB.(Animal)(Expi).units.unit_name_arr(EStats_CMB.(Animal)(Expi).units.pref_chan_id);
ReprImg = ReprStats_CMB.(Animal)(Expi).Evol.BestBlockAvgImg; % uint8 format
nexttile(T, Expi)
imshow(ReprImg)
title(compose("Exp%02d, PrefChan %s",Expi,prefchan_lab))
% nexttile(T2, Expi)
% interpMsk = MasksCol.(Animal)(Expi).(layername).interp;
% CLIM = prctile(interpMsk,[CMINPrct,MAXPrct],'all');
% alphaMsk = clip((interpMsk - CLIM(1)) ./ (CLIM(2) - CLIM(1)), 0, 1) * (ALIM(2) - ALIM(1))+ALIM(1);
% imshow(double(ReprImg)/255.*alphaMsk)
% title(compose("Exp%02d, PrefChan %s",Expi,prefchan_lab))
end
end
%%
tabA = readtable(fullfile(mat_dir, "Alfa"+"_EvolTrajStats.csv"));
tabB = readtable(fullfile(mat_dir, "Beto"+"_EvolTrajStats.csv"));
ExpTab_cmb = cat(1, tabA, tabB); % load hte table for Evol Traject successfulness
%%
sucs_msk = (ExpTab_cmb.t_p_initend<1E-3)&(ExpTab_cmb.t_p_initmax<1E-3);
msk = (ExpTab_cmb.t_p_initmax<1E-2);
V1msk = ExpTab_cmb.pref_chan <=48 & ExpTab_cmb.pref_chan >= 33;
V4msk = ExpTab_cmb.pref_chan <=64 & ExpTab_cmb.pref_chan >= 49;
ITmsk = ExpTab_cmb.pref_chan <=32 & ExpTab_cmb.pref_chan >= 1;
Alfamsk = ExpTab_cmb.Animal == "Alfa";
Betomsk = ExpTab_cmb.Animal == "Beto";
%%
layername="conv4_3";CLIMPrct = [20, 75]; ALIM = [0.1, 1];
[hIT,TIT] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(ITmsk&sucs_msk,:),[6,6],"IT cortex Both (success p<1E-3)",...
    layername,CLIMPrct,ALIM,"Both_IT_Success_"+layername);
[hV4,TV4] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(V4msk&sucs_msk,:),[6,6],"V4 cortex Both (success p<1E-3)",...
    layername,CLIMPrct,ALIM,"Both_V4_Success_"+layername);
[hV1,TV1] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(V1msk&sucs_msk,:),[6,6],"V1 cortex Both (success p<1E-3)",...
    layername,CLIMPrct,ALIM,"Both_V1_Success_"+layername);
%
layername="conv5_3";CLIMPrct = [20, 75]; ALIM = [0.1, 1];
[hIT,TIT] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(ITmsk&sucs_msk,:),[6,6],"IT cortex Both (success p<1E-3)",...
    layername,CLIMPrct,ALIM,"Both_IT_Success_"+layername);
[hV4,TV4] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(V4msk&sucs_msk,:),[6,6],"V4 cortex Both (success p<1E-3)",...
    layername,CLIMPrct,ALIM,"Both_V4_Success_"+layername);
[hV1,TV1] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(V1msk&sucs_msk,:),[6,6],"V1 cortex Both (success p<1E-3)",...
    layername,CLIMPrct,ALIM,"Both_V1_Success_"+layername);
%
layername="conv3_3";CLIMPrct = [20, 75]; ALIM = [0.1, 1];
[hIT,TIT] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(ITmsk&sucs_msk,:),[6,6],"IT cortex Both (success p<1E-3)",...
    layername,CLIMPrct,ALIM,"Both_IT_Success_"+layername);
[hV4,TV4] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(V4msk&sucs_msk,:),[6,6],"V4 cortex Both (success p<1E-3)",...
    layername,CLIMPrct,ALIM,"Both_V4_Success_"+layername);
[hV1,TV1] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(V1msk&sucs_msk,:),[6,6],"V1 cortex Both (success p<1E-3)",...
    layername,CLIMPrct,ALIM,"Both_V1_Success_"+layername);
%% 

[hIT,TIT] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(ITmsk&sucs_msk,:),[6,6],"IT cortex Both (success p<1E-3)",...
    false,"Both_IT_Success");
[hV4,TV4] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(V4msk&sucs_msk,:),[6,6],"V4 cortex Both (success p<1E-3)",...
    false,"Both_V4_Success");
[hV1,TV1] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
    ExpTab_cmb(V1msk&sucs_msk,:),[6,6],"V1 cortex Both (success p<1E-3)",...
    false,"Both_V1_Success");
function [h, T] = montage_proto_msk(EStats_CMB,ReprStats_CMB,MasksCol,...
                ExpTab,mtgsize,titstr,layername,CLIMPrct,ALIM,savename)
global savedir
if islogical(layername) && ~layername && (nargin==7 ||  nargin==8), raw_img=true;
if nargin==8, savename = CLIMPrct;
elseif nargin==7, savename = "Both_Success";end
else, raw_img = false;end
Explist=table2struct(ExpTab);
ExpN = numel(Explist);
mtgsize(1) = int32(ceil(ExpN/mtgsize(2)));
figW = 168 * mtgsize(2);
figH = 175 * mtgsize(1) + 30;
h = figure;set(h,'pos',[600, -50, figW, figH]);
T = tiledlayout(h,mtgsize(1),mtgsize(2),'Padd','none','TileSp','none');
for Ei = 1:numel(Explist)
Expi = Explist(Ei).Expi;
Animal = Explist(Ei).Animal;
prefchan_lab = EStats_CMB.(Animal)(Expi).units.unit_name_arr(...
    EStats_CMB.(Animal)(Expi).units.pref_chan_id);
ReprImg = ReprStats_CMB.(Animal)(Expi).Evol.BestBlockAvgImg; % uint8 format
nexttile(T, Ei)
if raw_img
imshow(ReprImg)
else
interpMsk = MasksCol.(Animal)(Expi).(layername).interp;
CLIM = prctile(interpMsk,CLIMPrct,'all');
alphaMsk = clip((interpMsk - CLIM(1)) ./ (CLIM(2) - CLIM(1)), 0, 1) * (ALIM(2) - ALIM(1))+ALIM(1);
alphaMsk(isnan(interpMsk)) = ALIM(1);
imshow(double(ReprImg)/255.*alphaMsk)
end
title(compose("%s Exp%d, PrfCh%s",Animal(1),Expi,prefchan_lab))
end
if ~raw_img, titstr = titstr+compose(" %s Prctile Thresh %d-%d",layername,CLIMPrct(1),CLIMPrct(2)); 
else, titstr = titstr+" raw"; end
title(T,titstr)
savefig(h,fullfile(savedir,"summary",savename+".fig"))
saveas(h,fullfile(savedir,"summary",savename+".pdf"))
saveas(h,fullfile(savedir,"summary",savename+".png"))
end
% savefig(hIT,fullfile(savedir,"summary","Both_IT_Success.fig"))
% saveas(hIT,fullfile(savedir,"summary","Both_IT_Success.pdf"))
% saveas(hIT,fullfile(savedir,"summary","Both_IT_Success.png"))
% savefig(hV4,fullfile(savedir,"summary","Both_V4_Success.fig"))
% saveas(hV4,fullfile(savedir,"summary","Both_V4_Success.pdf"))
% saveas(hV4,fullfile(savedir,"summary","Both_V4_Success.png"))
% savefig(hV1,fullfile(savedir,"summary","Both_V1_Success.fig"))
% saveas(hV1,fullfile(savedir,"summary","Both_V1_Success.pdf"))
% saveas(hV1,fullfile(savedir,"summary","Both_V1_Success.png"))
%% Refactor the structure of ccMskStats object. 
% ccMskStats_cmb = repmat(struct(),1,length(EStats));
% ExpType = "Evol";
% for Expi = 1:length(EStats)
%     ccMskStats_cmb(Expi).(ExpType).units = ccMskStats(Expi).units;
%     for layername = ["conv3_3", "conv4_3", "conv5_3"]
%         ccMskStats_cmb(Expi).(ExpType).(layername) = ccMskStats(Expi).(layername);
%     end
% end
% ccMskStats = ccMskStats_cmb;

% figure(12);set(12,'position',[566         372        1697         481])
% %%
% % Expi = 29; ExpType = "Evol"; layername = "conv5_3"; 
% sum_method = "L1"; % Method to sum up the feature dim
% sup_boundary = true;bdr = 1; % options and params for suppressing boundary activation
% cut_bdr = struct("conv3_3",3,"conv4_3",2,"conv5_3",1); % how many border elements to cut!
% for Expi = 1:numel(EStats)
% prefchan_lab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
% for layername = ["conv3_3","conv4_3","conv5_3"]
% bdr = cut_bdr.(layername);
% MskMap = ccMskStats(Expi).(ExpType).(layername).(sum_method);
% idx = find(contains(CNNRF.Name,layername));
% startpix = CNNRF.start(idx,:);
% jumppix = CNNRF.jump(idx,:);
% H = size(MskMap,1);W = size(MskMap,2);
% [centGridX, centGridY] = meshgrid(startpix(1)+jumppix(1)*(0:H-1),startpix(2)+jumppix(2)*(0:W-1));
% % ccMskStats(Expi).Evol.conv3_3.L1;
% if sup_boundary 
% padMap = nan(size(MskMap));
% padMap(bdr+1:end-bdr,bdr+1:end-bdr) = MskMap(bdr+1:end-bdr,bdr+1:end-bdr);
% else
% padMap = MskMap;
% end
% CLIM = prctile(padMap,[0,98],'all')';
% % figure(6);
% % imagesc(padMap);
% % axis image;colorbar();caxis(CLIM)
% [interpX, interpY] = meshgrid(linspace(1,224,256), linspace(1,224,256));
% interpMsk = griddata(centGridX, centGridY, double(padMap), interpX, interpY);
% alphaMsk = clip((interpMsk - CLIM(1)) ./ (CLIM(2) - CLIM(1)) *0.9,0,1);
% alphaMsk(isnan(alphaMsk)) = 0;
% %
% ReprImg = ReprStats(Expi).Evol.BestImg; % uint8 format
% figure(12);clf
% subtightplot(1,4,1,0.02,[0.02,0.12],0.05) % this layout has be tuned to be great.
% imshow(ReprImg)
% subtightplot(1,4,2,0.02,[0.02,0.12],0.05)
% imshow(double(ReprImg) / 255 .* alphaMsk)
% subtightplot(1,4,3,0.02,[0.02,0.12],0.05)
% imshow(alphaMsk)
% subtightplot(1,4,4,0.02,[0.02,0.12],0.05)
% imagesc(padMap);
% axis image;colorbar();caxis(CLIM)
% suptitle(compose("\n%s %s Exp%d Corr Coef Mask Channel compressed with %s threshold %s\nPref Chan %s rate in [%d,%d]ms with VGG16 %s, img cent [%.1f, %.1f] size %d deg",...
%     Animal,ExpType,Expi,"bi-sided",sum_method,prefchan_lab,50,200,strrep(layername,"_","-"),EStats(Expi).evol.imgpos(1),EStats(Expi).evol.imgpos(2),EStats(Expi).evol.imgsize))
% pause(0.2)
% saveas(12,fullfile(savedir, compose("%s_%s_Exp%d_%s-%s.jpg",Animal,ExpType,Expi,layername,sum_method)))
% saveas(12,fullfile(savedir, compose("%s_%s_Exp%d_%s-%s.pdf",Animal,ExpType,Expi,layername,sum_method)))
% end
% end