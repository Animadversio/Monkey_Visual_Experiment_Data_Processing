%% Code to create mask and merge on images from the masks stored in the correlation. 
%  Also code to visualize and analyze the mask in below!
ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa"; %Expi = 29;  %thresh = "neg"; sum_method="max"; %"both"
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')
%%
ccMskStats = repmat(struct(),1,length(EStats));
%% Review the ccFeatTsr with respect to the evolution of the psth and images. 
%  Actually this is the newer version of the visualization function
ExpType = "Manif";
doPlot = true; doSave = false;
thresh = "pos"; 
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


%%
% tabA = readtable(fullfile(mat_dir, "Alfa"+"_EvolTrajStats.csv"));
% tabB = readtable(fullfile(mat_dir, "Beto"+"_EvolTrajStats.csv"));
% ExpTab_cmb = cat(1, tabA, tabB); % load hte table for Evol Traject successfulness
% sucs_msk = (ExpTab_cmb.t_p_initend<1E-3)&(ExpTab_cmb.t_p_initmax<1E-3);
%%
Animal = "Beto"; 
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
% load(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')
load(fullfile(mat_dir, compose("%s_ImageRepr.mat", Animal)), 'ReprStats')
sucstab = readtable(fullfile(mat_dir, Animal+"_EvolTrajStats.csv"));
%% Purely view the masks with different compression methods
figdir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr\Compress_Cmp";
for ExpType = ["Evol"]%"Manif",
doPlot = false; doSave = true;
thresh = "pos"; ncol=5;
figure(1);T=tiledlayout(3,ncol,'padd','compact','TileSp','compact');
maplayers = ["conv3_3", "conv4_3", "conv5_3"]; % layers to plot the tensor
methodlist = ["max","L1","L1signif"]; % compression methods to compare
for Expi = 1:31%32:numel(EStats)
ccMskStats(Expi).(ExpType).units = EStats.units;
imgpos = EStats(Expi).evol.imgpos;
imgsize = EStats(Expi).evol.imgsize;
t_end = sucstab.t_initend(Expi);
t_max = sucstab.t_initmax(Expi);
DAOA_end = sucstab.DAOA_initend(Expi);
DAOA_max = sucstab.DAOA_initmax(Expi);
prefchanlab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
% winopen(fullfile(moviedir, compose("%s_Evol_Exp%02d_Best_PSTH.mov.avi",Animal,Expi)))
% winopen(fullfile(moviedir, compose("%s_Manif_Exp%02d_Avg_PSTH.mov.avi",Animal,Expi)))
savedir = fullfile(ccmat_dir,compose("%s_%s_Exp%d",Animal,ExpType,Expi));
for li = 1:numel(maplayers) %["conv5_3"]%"conv5_3";
layername = maplayers(li);
% if layername=="conv3_3" && Expi == 40, continue; end
outfn = fullfile(ccmat_dir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername));
load(outfn,'cc_tsr','MFeat','StdFeat','wdw_vect','cc_refM','cc_refS');
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
%
if thresh == "both"
maskTsr = abs(t_signif_tsr)>=5; % Bi sided thresholding 
elseif thresh == "pos"
maskTsr = t_signif_tsr>=5; % positively correlated thresholding
elseif thresh == "neg"
maskTsr = t_signif_tsr<=-5; % negatively correlated thresholding
end
maskN = squeeze(sum(maskTsr,[1,2,3])); % N of correlated units per
ctmp = cc_tsr(:,:,:,end);
ccmean = mean(ctmp(maskTsr(:,:,:,end)));
ccthresh = min(ctmp(maskTsr(:,:,:,end)),[],'all');
ccmax = max(ctmp,[],'all');
ccprctl = prctile(ctmp,[99,99.9,99.99,99.999],'all');
mskprct = maskN(end) / numel(ctmp);

plotTsr = cc_tsr;
plotTsr(~maskTsr) = 0; 
for mi = 1:numel(methodlist)
sum_method = methodlist(mi);
if sum_method=="L1"
L1plotTsr = squeeze(mean(abs(plotTsr),3));
elseif sum_method=="L1signif"
L1plotTsr = squeeze(sum(abs(plotTsr),3)./(sum(maskTsr,3)));
elseif sum_method=="max"
L1plotTsr = squeeze(max(abs(plotTsr),[],3));
end
nexttile(T,(li-1)*ncol+mi)
imagesc(L1plotTsr(:,:,end));axis image;colorbar();%caxis(CLIM_arr(fi,:));
title(compose("%s, %s compressed\n %d units (%.3f prc), mean cc %.3f",...
    layername,sum_method,maskN(end),mskprct,ccmean))
end
nexttile(T,(li-1)*ncol+4)
histogram(ctmp(:),'EdgeColor','none');box off
xlim([-1,1]);if ~isempty(ccthresh),vline(ccthresh);end%caxis(CLIM_arr(fi,:));
title(compose("%s, cc max %.3f, Prctle 99-99.99:\n [%.3f %.3f %.3f %.3f]",...
    layername,ccmax,ccprctl(1),ccprctl(2),ccprctl(3),ccprctl(4)))
end
title(T,compose("%s %s Exp%d Corr Coef with VGG16, Ch%s %.1f deg [%.1f, %.1f]\nt_{end}=%.2f t_{max}=%.2f DAOA_{end}=%.2f DAOA_{max}=%.2f\nParam: rate in [%d,%d]ms, %s thresholded",...
    Animal,ExpType,Expi,prefchanlab,imgsize,imgpos(1),imgpos(2),t_end,t_max,DAOA_end,DAOA_max,wdw_vect(end,1),wdw_vect(end,2),thresh))
% Show the exemplar images!
nexttile(T,(1-1)*ncol+5)
imshow(ReprStats(Expi).Evol.BestBlockAvgImg);
title("BestBlockAvg image")
nexttile(T,(2-1)*ncol+5)
imshow(ReprStats(Expi).Evol.BestImg);
title("Best Evol Trial image")
nexttile(T,(3-1)*ncol+5)
imshow(ReprStats(Expi).Manif.BestImg);
title("Best Manif image")
if doSave
    savefn = compose("%s_%s_Exp%d_%s_mask_summary",Animal,ExpType,Expi,thresh);
    saveas(1, fullfile(figdir, savefn + ".png"));
    saveas(1, fullfile(figdir, savefn + ".pdf"));
end
end
end
%% Reload the computed mask and create mask from it. 
net = vgg16; % load the net to compute the receptive field structure. 
%%
Animal = "Beto";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')
load(fullfile(mat_dir, compose("%s_Evol_ccFtMask.mat",Animal)),'ccMskStats')
load(fullfile(mat_dir, compose("%s_ImageRepr.mat",Animal)),'ReprStats')
%% Interpolate and Visualize masks Together with images.
savedir = "E:\OneDrive - Washington University in St. Louis\Repr_with_Mask";
figure(12);set(12,'position',[566         372        1697         481])
%%
% Expi = 29; ExpType = "Evol"; layername = "conv5_3"; 
sum_method = "L1"; % Method to sum up the feature dim
sup_boundary = true;bdr = 1; % options and params for suppressing boundary activation
cut_bdr = struct("conv3_3",3,"conv4_3",2,"conv5_3",1); % how many border elements to cut!
for Expi = 1:numel(EStats)
prefchan_lab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
for layername = ["conv3_3","conv4_3","conv5_3"]
bdr = cut_bdr.(layername);
MskMap = ccMskStats(Expi).(ExpType).(layername).(sum_method);
CNNRF = CNN_receptive_field(net);
idx = find(contains(CNNRF.Name,layername));
startpix = CNNRF.start(idx,:);
jumppix = CNNRF.jump(idx,:);
H = size(MskMap,1);W = size(MskMap,2);
[centGridX, centGridY] = meshgrid(startpix(1)+jumppix(1)*(0:H-1),startpix(2)+jumppix(2)*(0:W-1));
% ccMskStats(Expi).Evol.conv3_3.L1;
if sup_boundary 
padMap = nan(size(MskMap));
padMap(bdr+1:end-bdr,bdr+1:end-bdr) = MskMap(bdr+1:end-bdr,bdr+1:end-bdr);
else
padMap = MskMap;
end
CLIM = prctile(padMap,[0,98],'all')';
% figure(6);
% imagesc(padMap);
% axis image;colorbar();caxis(CLIM)
[interpX, interpY] = meshgrid(linspace(1,224,256), linspace(1,224,256));
interpMsk = griddata(centGridX, centGridY, double(padMap), interpX, interpY);
alphaMsk = clip((interpMsk - CLIM(1)) ./ (CLIM(2) - CLIM(1)) *0.9,0,1);
alphaMsk(isnan(alphaMsk)) = 0;
%
figure(12);clf
subtightplot(1,4,1,0.02,[0.02,0.12],0.05)
imshow(ReprStats(Expi).Evol.BestImg)
subtightplot(1,4,2,0.02,[0.02,0.12],0.05)
imshow(double(ReprStats(Expi).Evol.BestImg) / 255 .* alphaMsk)
subtightplot(1,4,3,0.02,[0.02,0.12],0.05)
imshow(alphaMsk)
subtightplot(1,4,4,0.02,[0.02,0.12],0.05)
imagesc(padMap);
axis image;colorbar();caxis(CLIM)
suptitle(compose("\n%s %s Exp%d Corr Coef Mask Channel compressed with %s threshold %s\nPref Chan %s rate in [%d,%d]ms with VGG16 %s, img cent [%.1f, %.1f] size %d deg",...
    Animal,ExpType,Expi,"bi-sided",sum_method,prefchan_lab,50,200,strrep(layername,"_","-"),EStats(Expi).evol.imgpos(1),EStats(Expi).evol.imgpos(2),EStats(Expi).evol.imgsize))
pause(0.2)
saveas(12,fullfile(savedir, compose("%s_%s_Exp%d_%s-%s.jpg",Animal,ExpType,Expi,layername,sum_method)))
saveas(12,fullfile(savedir, compose("%s_%s_Exp%d_%s-%s.pdf",Animal,ExpType,Expi,layername,sum_method)))
end
end
%%

%%
%%

moviedir = "E:\OneDrive - Washington University in St. Louis\Evol_Manif_Movies";
% winopen(fullfile(moviedir, compose("%s_Evol_Exp%02d_Best_PSTH.mov.avi",Animal,Expi)))
winopen(fullfile(moviedir, compose("%s_Manif_Exp%02d_Avg_PSTH.mov.avi",Animal,Expi)))
%%
CNNRF = CNN_receptive_field(net);
idx = find(contains(CNNRF.Name,layername));
startpix = CNNRF.start(idx,:);
jumppix = CNNRF.jump(idx,:);
H = size(cc_tsr,1);W = size(cc_tsr,2);
[centGridX, centGridY] = meshgrid(startpix(1)+jumppix(1)*(0:H-1),startpix(2)+jumppix(2)*(0:W-1));
%%
[interpX, interpY] = meshgrid(linspace(1,224,256), linspace(1,224,256));
interpMsk = griddata(centGridX, centGridY, double(L1plotTsr(:,:,end)), interpX, interpY);
%%
figure;imagesc(interpMsk)

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