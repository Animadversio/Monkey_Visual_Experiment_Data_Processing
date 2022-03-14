%% Pred w mask

net = alexnet;
py.importlib.import_module("numpy");
py.importlib.import_module("pickle");
%%
modelroot = "E:\OneDrive - Washington University in St. Louis\corrFeatTsr_FactorVis\models";
modeldir = fullfile(modelroot,"resnet50_linf8-layer3_NF3_bdr1_Tthresh_3__nobdr_res-robust_CV");
% py.numpy.load()
%%
figdir = "E:\OneDrive - Harvard University\similar_img_w_mask";
figh1 = 1; figh2=2; figh3=3; figh4=4; figh5=5;
for Animal = ["Alfa", "Beto"]
for Expi = 1:46
Ddict = py.pickle.load(py.open(fullfile(modeldir,compose("%s_Exp%02d_factors.pkl",Animal,Expi)),"rb"));
bestimg = Ddict.get('BestImg').double / 255.0;
Hmaps = Ddict.get('Hmaps').single;
tsr_proto = Ddict.tsr_proto.single;
fact_proto1 = Ddict.fact_protos{1}.single;
fact_proto2 = Ddict.fact_protos{2}.single;
fact_proto3 = Ddict.fact_protos{3}.single;
explabel = compose("%s-Exp%02d",Animal,Expi);
%%
figure(figh1);set(figh1,'pos',[20         401        1825         480]);
subplot(141)
imagesc(bestimg);axis image
title("Best Evolved image");
subplot(142)
imagesc(Hmaps/max(Hmaps,[],'all')); axis image
title("Importance Mask of the three factors");
featmsk = sum(Hmaps,3)/max(sum(Hmaps,3),[],'all');
% imagesc(featmsk); axis image
subplot(143)
alphamsk_rsz = imresize(featmsk,[256,256]);
imagesc(alphamsk_rsz); title("interpolated alpha mask"); axis image
subplot(144)
vmin = 0.2; vmax = 0.8;
alphamsk_clip = min(vmax,max(vmin,alphamsk_rsz));
alphamsk_clip = (alphamsk_clip - vmin) / (vmax - vmin);
maskedimg = (double(bestimg) .* alphamsk_clip);
imshow(maskedimg);title("masked best image")
sgtitle(explabel)

figure(figh2);set(figh2,'pos',[20         401        1825         480]);
subplot(141)
imshow(tsr_proto);title("facttsr model proto")
subplot(142)
imshow(fact_proto1); title("fact1 model proto")
subplot(143)
imshow(fact_proto2); title("fact2 model proto")
subplot(144)
imshow(fact_proto3); title("fact3 model proto")
sgtitle(explabel+" Prototypes of the Factorized Model")
%%
% reprvec_mask = squeeze(activations(net,maskedimg,'fc6'));
% reprvec_unmask = squeeze(activations(net,bestimg,'fc6'));
im2proc = cat(4, bestimg, maskedimg, tsr_proto, fact_proto1, fact_proto2, fact_proto3) * 255;
repr_all = squeeze(activations(net,im2proc,'fc6'));

%%
distmeasure = "correlation";
% cosdist_mask  = pdist2(DS.allFeatures, reprvec_mask',distmeasure);
% cosdist_unmask  = pdist2(DS.allFeatures, reprvec_unmask',distmeasure);
dist_all  = pdist2(DS.allFeatures, repr_all',distmeasure);

unmskbest_dist = dist_all(:, 1);
mskbest_dist = dist_all(:, 2);
tsr_dist = dist_all(:, 3);
fact1_dist = dist_all(:, 4);
fact2_dist = dist_all(:, 5);
fact3_dist = dist_all(:, 6);
%
figure(figh3);clf;set(figh3,'pos', [680   558   560   420])
histogram(unmskbest_dist);hold on
histogram(mskbest_dist)
histogram(tsr_dist)
legend(["with mask","without mask","tsr reevolve"])
xlabel(distmeasure+" distance")
title([explabel,"Distance to natural images"])
%%
% [minval_mask,minidx_mask] = mink(cosdist_mask,25);
% [minval_unmask,minidx_unmask] = mink(cosdist_unmask,25);
[~, minidx_unmskbest] = mink(unmskbest_dist, 25);
[~, minidx_mskbest] = mink(mskbest_dist, 25);
[~, minidx_tsr] = mink(tsr_dist, 25);
[~, minidx_fact1] = mink(fact1_dist, 25);
[~, minidx_fact2] = mink(fact2_dist, 25);
[~, minidx_fact3] = mink(fact3_dist, 25);

%%
figure(figh4);set(figh4,'pos', [0          97        1890         760]);
subplot(121);
montage(DS.allFeatureNames(minidx_mskbest), 'Size', [5, 5],'ThumbnailSize',[256,256])
title("most similar images with Mask")
subplot(122);
montage(DS.allFeatureNames(minidx_unmskbest), 'Size', [5, 5],'ThumbnailSize',[256,256])
title("most similar images without Mask")
sgtitle(explabel)
%
figure(figh5);set(figh5,'pos', [0          97        3300         760]);
subplot(141);
montage(DS.allFeatureNames(minidx_tsr), 'Size', [5, 5],'ThumbnailSize',[256,256])
title("most similar images to prototype with tensor models")
subplot(142);
montage(DS.allFeatureNames(minidx_fact1), 'Size', [5, 5],'ThumbnailSize',[256,256])
title("most similar images to prototype with factor 1")
subplot(143);
montage(DS.allFeatureNames(minidx_fact2), 'Size', [5, 5],'ThumbnailSize',[256,256])
title("most similar images to prototype with factor 2")
subplot(144);
montage(DS.allFeatureNames(minidx_fact3), 'Size', [5, 5],'ThumbnailSize',[256,256])
title("most similar images to prototype with factor 3")
sgtitle(explabel)
saveallform(figdir, explabel+"_proto_mask", figh1, ["png"])
saveallform(figdir, explabel+"_factorProto_mask", figh2, ["png"])
saveallform(figdir, explabel+"_dist_hist_cmp", figh3, ["png"])
saveallform(figdir, explabel+"_closest_img_msk_unmsk", figh4, ["png"])
saveallform(figdir, explabel+"_closest_img_factprotos", figh5, ["png"])
end
end