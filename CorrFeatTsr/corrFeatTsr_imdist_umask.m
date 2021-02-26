% Test if the reference that are similar to the prototpyes under the masks
% are more activating.
Animal = "Alfa"; Set_Path;%thresh = "neg"; sum_method="max"; %"both"
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(mat_dir, compose("%s_Manif_stats.mat", Animal)), 'Stats')
load(fullfile(mat_dir,compose("%s_Evol_ccFtMask.mat",Animal)),'ccMskStats')
load(fullfile(mat_dir, compose("%s_ImageRepr.mat",Animal)),'ReprStats')
Expi = 3;  
%%
D = torchImDist();
D = D.select_metric("squeeze", true);
% DL = torchLPIPS();
net = vgg16;
CNNRF = CNN_receptive_field(net);
clear net
%%
distMap2Ref = repmat(struct(),numel(EStats),1);
%%
for Expi = 1:numel(EStats)
prefchan_lab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
fprintf("%s Exp %d prefchan %s\n",Animal,Expi,prefchan_lab)
disp_pix_n = EStats(Expi).evol.imgsize * 40;% actual num of pix displayed on screen
protoimg = ReprStats(Expi).Evol.BestBlockAvgImg;
protoimg_disp = imresize(protoimg, [disp_pix_n,disp_pix_n]);
refimgs = {}; refimgs_disp = {};
for i = 1:numel(EStats(Expi).ref.imgnm)
    img = imread(EStats(Expi).ref.impaths(i));
    [H,W,C] = size(img);
	if C==1, img = repmat(img,1,1,3); end
    refimgs{i} = imresize(img, [256,256]);
    refimgs_disp{i} = imresize(img, [disp_pix_n,disp_pix_n]);
end
refimgtsr = cat(4,refimgs{:});
refimgstsr_disp = cat(4,refimgs_disp{:});
distmaps_disp = D.distance(protoimg_disp, refimgstsr_disp);
distmaps_disp = permute(distmaps_disp,[2,3,1]);
distmaps = D.distance(protoimg, refimgtsr);
distmaps = permute(distmaps,[2,3,1]);
ui = EStats(Expi).units.unit_in_pref_chan;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
distMap2Ref(Expi).score_vec = score_vec;
distMap2Ref(Expi).distmaps = distmaps;
distMap2Ref(Expi).distmaps_disp = distmaps_disp;
% figure;montage(distmaps_disp)
%%
for layername = ["conv3_3","conv4_3","conv5_3"]
MskMap = ccMskStats(Expi).("Evol").(layername).("L1signif");
[interpMsk,alphaMsk] = ccArr2mask(MskMap,layername,CNNRF, [50,50]);
msk_dists = mask_w_avgmap(distmaps_disp,alphaMsk);
msk_out_dists = mask_w_avgmap(distmaps_disp,alphaMsk,true);
[cval,pval] = corr(score_vec, msk_dists);
[cval_o,pval_o] = corr(score_vec, msk_out_dists,'row','complete');
fprintf("Im Dist (under Alpha mask of layer %s) %.3f(%.1e)  out dist %.3f(%.1e)\n",...
    layername, cval,pval,cval_o,pval_o)
distMap2Ref(Expi).(layername).alphaMsk = alphaMsk;
distMap2Ref(Expi).(layername).interpMsk = interpMsk;
distMap2Ref(Expi).(layername).cval = cval;
distMap2Ref(Expi).(layername).pval = pval;
end
end
%%
figure;
scatter(ones(numel(distMap2Ref),1),arrayfun(@(R)R.(layername).cval,distMap2Ref))
%%
tabA = readtable(fullfile(mat_dir, "Alfa"+"_EvolTrajStats.csv"));
tabB = readtable(fullfile(mat_dir, "Beto"+"_EvolTrajStats.csv"));
ExpTab_cmb = cat(1, tabA, tabB); % load hte table for Evol Traject successfulness
%%
dist_refact_corr = arrayfun(@(R)R.(layername).cval,distMap2Ref);
%%
figure;
succmsk = ~(tabA.pref_chan>=33&tabA.pref_chan<=48);
scatter(ones(numel(distMap2Ref(succmsk)),1),arrayfun(@(R)R.(layername).cval,distMap2Ref(succmsk)))

function w_avgmap = mask_w_avgmap(mapstack, mask, inverse)
if nargin == 2, inverse = false;end
if inverse
   mask = max(mask,[],'all') - mask; 
end
mapsize = size(mapstack,[1,2]);
mask_rsz = imresize(mask,mapsize);
w_avgmap = squeeze(nansum(mapstack.*mask_rsz,[1,2]) / nansum(mask_rsz,[1,2]));
end
function [interpMsk,alphaMsk] = ccArr2mask(MskMap,layername,CNNRF,CLIMPrct,P)
if nargin < 5, 
% P.sum_method = "L1"; % Method to sum up the feature dim
P.sup_boundary = true; % options and params for suppressing boundary activation
P.cut_bdr = struct("conv3_3",3,"conv4_3",2,"conv5_3",1); % how many border elements to cut!
end
if nargin < 4, CLIMPrct = [20,75]; end
bdr = P.cut_bdr.(layername);
idx = find(contains(CNNRF.Name,layername));
startpix = CNNRF.start(idx,:);
jumppix = CNNRF.jump(idx,:);
H = size(MskMap,1);W = size(MskMap,2);
[centGridX, centGridY] = meshgrid(startpix(1)+jumppix(1)*(0:H-1),startpix(2)+jumppix(2)*(0:W-1));
if P.sup_boundary 
padMap = nan(size(MskMap));
padMap(bdr+1:end-bdr,bdr+1:end-bdr) = MskMap(bdr+1:end-bdr,bdr+1:end-bdr);
else
padMap = MskMap;
end
[interpX, interpY] = meshgrid(linspace(1,224,256), linspace(1,224,256));
interpMsk = griddata(centGridX, centGridY, double(padMap), interpX, interpY);
CLIM = prctile(interpMsk,CLIMPrct,'all');
alphaMsk = clip((interpMsk - CLIM(1)) ./ (CLIM(2) - CLIM(1)),0,1);
alphaMsk(isnan(alphaMsk)) = 0;
end