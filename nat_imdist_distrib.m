% Natural image distance distribution. 
% And the distribution of im distance among GAN images, Gabor imges and Pasu images.
global figroot 
figroot = "E:\OneDrive - Washington University in St. Louis\ImDist_ref";
D = torchImDist("squeeze");
G6 = FC6Generator();
GB = torchBigGAN();
%% Randomly sample natural images
imagenet_root = "N:\Stimuli\imagenet-2012\imagenet12\images";
[imgids,INimgs] = imagenet_sampler(imagenet_root, 500);
ImgNet_imdist = struct();
%% Compute the distance based on different metrics
tic
D = torchImDist("squeeze");
distMat = D.distmat_B(INimgs,50);
ImgNet_imdist.squ = distMat;
D.metric = "L2";
ImgNet_imdist.L2 = D.distmat(INimgs);
D=D.load_net();
D.metric = "FC6_L2";
ImgNet_imdist.FC6 = D.distmat(INimgs);
D.metric = "FC6_corr";
ImgNet_imdist.FC6_corr = D.distmat(INimgs);
%% SSIM can take quite a while 
D = D.select_metric("SSIM"); 
ImgNet_imdist.SSIM = D.distmat(INimgs); 
%%
h = figure();
vis_distmat(ImgNet_imdist.SSIM,"ImageNet_test","SSIM",h)
vis_distmat(ImgNet_imdist.squ,"ImageNet_test","LPIPS_SqueezeNet",h)
vis_distmat(ImgNet_imdist.L2,"ImageNet_test","L2",h)
vis_distmat(ImgNet_imdist.FC6,"ImageNet_test","FC6",h)
vis_distmat(ImgNet_imdist.FC6_corr,"ImageNet_test","FC6_corr",h)
%%
save(fullfile(figroot,'ImageNet_test_dist.mat'),'ImgNet_imdist','imgids')
%%
figure;montage(INimgs(:,:,:,1:49))
saveas(gcf,fullfile(figroot,'ImageNet_test_examples.png'))

%% Compute distance between nearest neighbors
imdist_lpips = ImgNet_imdist.squ;
imdist_lpips = imdist_lpips + diag(nan(1,size(imdist_lpips,1)));
nearneighb_dist = nanmin(imdist_lpips,[],1);
D_m = mean(nearneighb_dist);
D_s = std(nearneighb_dist);
D_prc = prctile(nearneighb_dist,[5,95]);
fprintf("Distance between nearest neighbors\n%s Mean %.3f+-%.3f [5, 95] Percentile [%.3f, %.3f]\n", "squ", D_m, D_s, D_prc(1), D_prc(2))


%% Randomly sample images from the GAN space
%% FC6 GAN space.
imgs = G6.visualize(5*randn(500,4096));
FC6G_imdist = struct();
%%
D = torchImDist("squeeze");
distMat = D.distmat_B(imgs,50); % batch process is faster here
FC6G_imdist.squ = distMat;
%%
D.metric = "L2";
FC6G_imdist.L2 = D.distmat(imgs);
D=D.load_net();
D.metric = "FC6_L2";
FC6G_imdist.FC6 = D.distmat(imgs);
D.metric = "FC6_corr";
FC6G_imdist.FC6_corr = D.distmat(imgs);
%% SSIM can take quite a while 
tic
D = D.select_metric("SSIM"); 
FC6G_imdist.SSIM = D.distmat(imgs); 
toc % 2111 sec for 500 images matrix
%% Visualize distribution and quantify
h = figure();
vis_distmat(FC6G_imdist.squ,"FC6GAN_rand_img","LPIPS_SqueezeNet",h)
vis_distmat(FC6G_imdist.SSIM,"FC6GAN_rand_img","SSIM",h)
vis_distmat(FC6G_imdist.L2,"FC6GAN_rand_img","L2",h)
vis_distmat(FC6G_imdist.FC6,"FC6GAN_rand_img","FC6",h)
vis_distmat(FC6G_imdist.FC6_corr,"FC6GAN_rand_img","FC6_corr",h)
%%
save(fullfile(figroot,'FC6GAN_rand_img_dist.mat'),'FC6G_imdist')
%%
figure;montage(imgs(:,:,:,1:49))
saveas(gcf,fullfile(figroot,'FC6GAN_rand_examples.png'))


%% BigGAN images space. 
BigGAN_imdist = struct();
%%
codes = GB.sample_noise(500);
imgs = GB.visualize_latent(codes);
%%
D = torchImDist("squeeze");
distMat = D.distmat_B(imgs,50); % batch process is faster here
BigGAN_imdist.squ = distMat;
D.metric = "L2";
BigGAN_imdist.L2 = D.distmat(imgs);
D=D.load_net();
D.metric = "FC6_L2";
BigGAN_imdist.FC6 = D.distmat(imgs);
D.metric = "FC6_corr";
BigGAN_imdist.FC6_corr = D.distmat(imgs);
%% SSIM can take quite a while 
tic
D = D.select_metric("SSIM"); 
BigGAN_imdist.SSIM = D.distmat(imgs); 
toc % 2111 sec for 500 images matrix
%%
h = figure();
vis_distmat(BigGAN_imdist.squ,"BigGAN_rand_img","LPIPS_SqueezeNet",h)
vis_distmat(BigGAN_imdist.SSIM,"BigGAN_rand_img","SSIM",h)
vis_distmat(BigGAN_imdist.L2,"BigGAN_rand_img","L2",h)
vis_distmat(BigGAN_imdist.FC6,"BigGAN_rand_img","FC6",h)
vis_distmat(BigGAN_imdist.FC6_corr,"BigGAN_rand_img","FC6_corr",h)
%%
save(fullfile(figroot,'BigGAN_rand_img_dist.mat'),'BigGAN_imdist','codes')
%%
figure;montage(imgs(:,:,:,1:49))
saveas(gcf,fullfile(figroot,'BigGAN_rand_examples.png'))


%% Actual FC6 Images in Manifold Experiments
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
minDvec_all = [];
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir, Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
for Expi = 1:numel(ManifImDistStat)
minDvec = dist2nneighbor(ManifImDistStat(Expi).squ, 1);
% nneighbor_dist_summary(minDvec,"lpips")
minDvec_all = [minDvec_all, minDvec];
end
end
nneighbor_dist_summary(minDvec_all,"lpips")
%%
load(fullfile(mat_dir, "gab_imdist.mat"),'gab_imdist')
load(fullfile(mat_dir, "pasu_imdist.mat"),'pasu_imdist')
%%
fprintf("In the gabor patch space, ")
minDvec = dist2nneighbor(gab_imdist.squ, 0);
nneighbor_dist_summary(minDvec,"lpips")
fprintf("In the Pasupathy shape space, ")
minDvec = dist2nneighbor(pasu_imdist.squ, 0);
nneighbor_dist_summary(minDvec,"lpips")
%% Summary statistics
%% Print the summary string for each image space 
figroot = "E:\OneDrive - Washington University in St. Louis\ImDist_ref";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(figroot,'ImageNet_test_dist.mat'),'ImgNet_imdist')
load(fullfile(figroot,'FC6GAN_rand_img_dist.mat'),'FC6G_imdist')
load(fullfile(figroot,'BigGAN_rand_img_dist.mat'),'BigGAN_imdist')
load(fullfile(mat_dir, "gab_imdist.mat"),'gab_imdist')
load(fullfile(mat_dir, "pasu_imdist.mat"),'pasu_imdist')
%% 
diary(fullfile(figroot,"ref_imdist_stats.log"))
dist_summary(ImgNet_imdist,"ImageNet Random")
dist_summary(FC6G_imdist,"FC6 Random")
dist_summary(BigGAN_imdist,"BigGAN Random")
dist_summary(gab_imdist,"Gabor ")
dist_summary(pasu_imdist,"Pasupathy ")
diary off

function dist_summary(DistStruct,name)
fprintf("%s Image Set\n",name)
for metric = string(fieldnames(DistStruct)')
    Dmat = DistStruct.(metric);
    Dmat = Dmat + diag(nan(1,size(Dmat,1)));
    D_m = nanmean(Dmat,'all');
    D_s = nanstd(Dmat,1,'all');
    D_prc = prctile(Dmat,[5,95],'all');
    D_prc2 = prctile(Dmat,[20,80],'all');
    fprintf("%s Mean %.3f +- %.3f [20, 80] Percentile [%.3f, %.3f] [5, 95] Percentile [%.3f, %.3f]\n", metric, D_m, D_s, D_prc2(1), D_prc2(2), D_prc(1), D_prc(2))
end
end

function nneighbor_dist_summary(mindist, metric)
D_m = nanmean(mindist,'all');
D_s = nanstd(mindist,1,'all');
D_prc = prctile(mindist,[5,95],'all');
D_prc2 = prctile(mindist,[20,80],'all');
fprintf("Distance to nearest neightbor: %s Mean %.3f+-%.3f [5, 95] Percentile [%.3f, %.3f]\n", metric, D_m, D_s, D_prc(1), D_prc(2))
end

function mindist = dist2nneighbor(distmat, manif)
% return vector of distance to nearest neighbor.
if nargin == 1, manif = false; end
if manif
    assert(all(size(distmat)==[121,121]))
    distmat = distmat(11:111,11:111);
end
distmat = distmat + diag(nan(1,size(distmat,1)));
mindist = nanmin(distmat,[],1);
end
% vis_distmat: function description
function vis_distmat(distMat,imspace_str,metric_str,h) % [outputs] = 
    if isempty(h) || nargin < 4, h=figure();end
    global figroot
    distMat_nodiag = distMat + diag(nan(size(distMat,1),1));
    figure(h);
    histogram(distMat_nodiag(:),'norm','prob','edgecolor','none')
    XLIM = xlim;xlim([0,XLIM(2)]);box off
    statstr = compose("mean %.3f med %.3f std %.3f",nanmean(distMat_nodiag,'all'),...
        nanmedian(distMat_nodiag,'all'),nanstd(distMat_nodiag,1,'all'));
    vline(double(nanmean(distMat_nodiag,'all')),'-.',"Mean")
    % vline(nanmedian(distMat_nodiag,'all'),'-.',"Median")
    title([compose("%s %s ImDist",strrep(imspace_str,"_",' '),strrep(metric_str,"_",' ')),statstr])
    xlabel(strrep(metric_str,"_",' '));ylabel("prob")
    saveas(h,fullfile(figroot,compose("%s_%s.png",imspace_str,metric_str)))
    saveas(h,fullfile(figroot,compose("%s_%s.pdf",imspace_str,metric_str)))
end
function [imgids,imgtsr,imgs] = imagenet_sampler(imagenet_root, imgnum, imsize)
if nargin < 3, imsize=[256,256]; end
% imlist = ls(fullfile(imagenet_root,"test\*.JPEG")); % no need for this,
% pattern is fixed.
imgids = randsample(100000,imgnum);
imgpaths = arrayfun(@(ids) fullfile(imagenet_root, "test", compose("ILSVRC2012_test_%08d.JPEG",ids)),imgids);
imgs = arrayfun(@(imgpath) imread(imgpath), imgpaths,'Uni',0);
imgtsr = {};
for i = 1:numel(imgs)
[H,W,C] = size(imgs{i});
if C==1, imgs{i} = repmat(imgs{i},1,1,3); end
if H >= W
pad = (H-W)/2;
crp_rsz_img = imcrop(imgs{i},[0,pad,W,W]);
else
pad = (W-H)/2;
crp_rsz_img = imcrop(imgs{i},[pad,0,H,H]);
end
imgtsr{i} = imresize(crp_rsz_img, imsize);
end
imgtsr = cat(4,imgtsr{:});
end
%
%%
% % [d2first,sortidx] = sort(distMat(:,1));
% % distmat_sort = distMat(sortidx,sortidx);
% % figure; 
% % imagesc(distmat_sort)
% distMat_nodiag = distMat + diag(nan(size(distMat,1),1));
% figure(1);
% histogram(distMat_nodiag(:),'norm','prob','edgecolor','none')
% XLIM = xlim;xlim([0,XLIM(2)]);box off
% statstr = compose("mean %.3f med %.3f std %.3f",nanmean(distMat_nodiag,'all'),nanmedian(distMat_nodiag,'all'),nanstd(distMat_nodiag,1,'all'));
% vline(nanmean(distMat_nodiag,'all'),'-.',"Mean")
% % vline(nanmedian(distMat_nodiag,'all'),'-.',"Median")
% title(["FC6 random images ImDist",statstr])
% xlabel("LPIPS SqueezeNet");ylabel("prob")
% saveas(1,fullfile(figroot,"FC6GAN_squnetdist.png"))