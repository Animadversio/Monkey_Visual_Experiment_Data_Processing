ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
Animal = "Beto"; Expi = 29;  %thresh = "neg"; sum_method="max"; %"both"
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
ccMskStats = repmat(struct(),1,length(EStats));
%% Review the ccFeatTsr with respect to the evolution of the psth and images. 
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
%% Calculate the CLIM
if doPlot
CLIM_arr = zeros(24,2);
CLIM_arr(1:19,:) = CLIM_arr(1:19,:) + prctile(L1plotTsr(:,:,1:19),[2,98],'all')' + [0, 1E-4];
CLIM_arr(20:23,:) = CLIM_arr(20:23,:) + prctile(L1plotTsr(:,:,20:23),[2,98],'all')'+ [0, 1E-4];
CLIM_arr(24,:) = CLIM_arr(24,:) + prctile(L1plotTsr(:,:,24),[2,98],'all')'+ [0, 1E-4];
%
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
%
end
end
end
%%
save(fullfile(MatStats_path,compose("%s_Evol_ccFtMask.mat",Animal)),'ccMskStats')
% end
%% Refactor the structure of ccMskStats object. 
ccMskStats_cmb = repmat(struct(),1,length(EStats));
ExpType = "Evol";
for Expi = 1:length(EStats)
    ccMskStats_cmb(Expi).(ExpType).units = ccMskStats(Expi).units;
    for layername = ["conv3_3", "conv4_3", "conv5_3"]
        ccMskStats_cmb(Expi).(ExpType).(layername) = ccMskStats(Expi).(layername);
    end
end
ccMskStats = ccMskStats_cmb;
%%
Animal = "Alfa";
load(fullfile(MatStats_path,compose("%s_Evol_ccFtMask.mat",Animal)),'ccMskStats')
%%
Expi = 29; ExpType = "Evol"; layername = "conv3_3"; sum_method = "L1";
ccMskStats(Expi).(ExpType).(layername).(sum_method)
% ccMskStats(Expi).Evol.conv3_3.L1;


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