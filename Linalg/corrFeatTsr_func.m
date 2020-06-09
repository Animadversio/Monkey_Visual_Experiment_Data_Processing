%% Wrapper for the correlation method. 
%  
%  Input: Images (path or ), Response vectors, net, Layer, flag to see if shuffle or not
%  Output: correlation tensor, mean feat value
function [cc_tsr, MFeat, StdFeat, varargout] = corrFeatTsr_func(images, score_vect, net, layername, flags)
Bsz = flags.batch;
if ~isfield(flags,"shuffleN"), flags.shuffleN=100; end
if isstr(images(1)) || isstring(images(1))
    flags.loadIMG = true; 
    imgN = length(images); 
else
    flags.loadIMG = false; 
    imgN = size(images,4); 
    dlimg = images;
end
shuffleN = flags.shuffleN; 
score_shuffle = [score_vect];
for i = 1:shuffleN
    score_shuffle = [score_shuffle, score_vect(randperm(imgN),:)];
end
dummy = activations(net, zeros(224,224,3), layername);
ft_shape = size(dummy); % size(feat_tsr,[1,2,3]);
nfeat = prod(ft_shape); % number of feature predictors in total
ntpnt = size(score_vect,2);

if flags.online_compute % batch computation of correlation coefficient tensor. 
if flags.load_all_img % pre load all the images into dlimg
    imgcol = cellfun(@(imgnm) imresize(imread(fullfile(imgnm)),[224,224]), images, 'UniformOutput', false);
    dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
    clear imgcol
end
%% initialize the "sufficient statistics"
SSqFeat = zeros(nfeat,1,'single');
SFeat = zeros(nfeat,1,'single');
SSqrsp = zeros(ntpnt*(1+shuffleN),1,'single');
Srsp = zeros(ntpnt*(1+shuffleN),1,'single');
InnProd = zeros(nfeat, ntpnt*(1+shuffleN),'single');
%% batch computation of correlation coefficient tensor.
%% Batch by batch
curN = 0;
csr = 1; 
T0_all = tic;
while curN < imgN
    T0 = tic;
    csr = curN + 1; 
    csr_end = min(imgN, csr + Bsz -1);
    rsp = score_shuffle(csr:csr_end,:); % change this into 
    if ~flags.load_all_img && flags.loadIMG
    imgcol = cellfun(@(imgnm) imresize(imread(fullfile(imgnm)),[224,224]), images(csr:csr_end), 'UniformOutput', false);
    dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
    T1 = toc(T0);
    feat_tsr_tmp = activations(net, dlimg, layername);
    T2 = toc(T0);
    else % use the 
    T1 = toc(T0);
    feat_tsr_tmp = activations(net, dlimg(:,:,:,csr:csr_end), layername);
    T2 = toc(T0);
    end
    feat_tsr_tmp = reshape(feat_tsr_tmp, nfeat, []); % [nfeat, batchsize]
    SSqFeat = SSqFeat + sum(feat_tsr_tmp.^2, 2);
    SFeat = SFeat + sum(feat_tsr_tmp, 2);
    SSqrsp = SSqrsp + sum(rsp.^2, 1)'; % summation over samples
    Srsp = Srsp + sum(rsp, 1)';
    InnProd = InnProd + feat_tsr_tmp*rsp;
    curN = csr_end;
    T3 = toc(T0);
    fprintf("Latencies: load img %.1f CNN proc %.1f compute innprod %.1f\n", T1, T2, T3);
end
toc(T0_all)
% Compute mean, and std from the running stats. 
Mrsp = Srsp / curN;
MFeat = SFeat / curN;
Stdrsp = sqrt(SSqrsp / curN - Mrsp.^2);
StdFeat = sqrt(SSqFeat / curN - MFeat.^2);
cc_tsr_ref = (InnProd./curN - MFeat * Mrsp') ./ (StdFeat * Stdrsp');
cc_tsr_ref = single(reshape(cc_tsr_ref, [ft_shape, ntpnt, shuffleN + 1]));
cc_tsr = cc_tsr_ref(:,:,:,:,1);
if shuffleN > 0 
cc_tsr_ref = cc_tsr_ref(:,:,:,:,2:end);
cc_refM = mean(cc_tsr_ref,5);
cc_refS = std(cc_tsr_ref,0,5);
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
varargout{1} = cc_refM;
varargout{2} = cc_refS;
end
MFeat = reshape(MFeat, ft_shape);
StdFeat = reshape(StdFeat, ft_shape);
else % Commpute feature tensor over 3000+ images directly using activations. (faster then doing loop)
T0 = tic;
imgcol = cellfun(@(imgnm) imresize(imread(fullfile(imgnm)),[224,224]), images, 'UniformOutput', false);
dlimg = cell2mat(reshape(imgcol,1,1,1,[]));
T1 = toc(T0);
feat_tsr = activations(net, dlimg, layername, 'MiniBatchSize', Bsz);
T2 = toc(T0);
cc_tsr_ref = corr(reshape(feat_tsr,nfeat,imgN)', score_shuffle);
cc_tsr_ref = single(reshape(cc_tsr_ref, [ft_shape, ntpnt, shuffleN + 1]));
cc_tsr = cc_tsr_ref(:,:,:,:,1);
if shuffleN > 0 
cc_tsr_ref = cc_tsr_ref(:,:,:,:,2:end);
cc_refM = mean(cc_tsr_ref,5);
cc_refS = std(cc_tsr_ref,0,5);
t_signif_tsr = (cc_tsr - cc_refM) ./ cc_refS;
varargout{1} = cc_refM;
varargout{2} = cc_refS;
end
MFeat = mean(feat_tsr,4);
StdFeat = std(feat_tsr,0,4);
T3 = toc(T0);
fprintf("Latencies: load img %.1f CNN proc %.1f compute innprod %.1f\n", T1, T2, T3);
end

end
