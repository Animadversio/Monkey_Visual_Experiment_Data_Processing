%% Regress the codes onto the PSTH mat
n2c_beta = mvregress(BGrspmat(:,31:end)',code4trial(31:end,:));
c2n_beta = mvregress(code4trial(31:end,:),BGrspmat(:,31:end)');
%% Linear Regression
nsamp = size(BGrspmat,2);
n2c_beta = pinv([ones(nsamp,1),BGrspmat']) * code4trial;
c2n_beta = pinv([ones(nsamp,1),code4trial]) * BGrspmat';
%%
[evc_all, eva_all, evc_cls, eva_cls, evc_nos, eva_nos] = loadHessian();
%%
fixnoise = Trials.TrialRecord.User.space_cfg{2}{2};
G = torchBigGAN();
G = G.select_space("class", fixnoise);
%%
imgidx = 20:40;
figure(2);
subplot(121)
montage(G.visualize(code4trial(imgidx,:)))
subplot(122)
montage(G.visualize([ones(numel(imgidx), 1), BGrspmat(:,imgidx)'] * n2c_beta))
%%
ncomp = 50;
[Rsp_L,code_L,Rsp_S,code_S,betaPLS] = plsregress(BGrspmat',code4trial,ncomp);
%%
pred_code = [ones(nsamp,1),BGrspmat'] * betaPLS;
%%
imgidx = 31:60;
figure(3);
subplot(121)
montage(G.visualize(code4trial(imgidx,:)))
subplot(122)
montage(G.visualize(pred_code(imgidx,:)))
%% Decoding using Hessian Data. 
Expi = 11;
% Find the BigGAN Code for Hessian Experiments
data = py.numpy.load(fullfile(HEStats(Expi).meta.stimuli,"Hess_imgs","summary","class_ImDist_root_data.npz"));
codes_all_cls = data.get('vecs_arr').double;
codes_mat_cls = reshape(codes_all_cls, [], 256);
data = py.numpy.load(fullfile(HEStats(Expi).meta.stimuli,"Hess_imgs","summary","noise_ImDist_root_data.npz"));
codes_all_nos = data.get('vecs_arr').double;
codes_mat_nos = reshape(codes_all_nos, [], 256);
ncell = numel(HessBGStats(Expi).units.spikeID);
rspmat_cls = reshape(HessBGStats(Expi).class.resp_mat(:,:,:),ncell,[]);
rspmat_nos = reshape(HessBGStats(Expi).noise.resp_mat(:,:,:),ncell,[]);
%%
nsamp = size(rspmat,2);
rspmat = [rspmat_cls, rspmat_nos];
codes_mat = [codes_mat_cls; codes_mat_nos];
valmsk = find(sum(isnan(rspmat),1)==0); % some image don't have response...repitition.
%% PLS Decoding
ncomp = 80;
[Rsp_L,code_L,Rsp_S,code_S,betaPLS] = plsregress(rspmat(:,valmsk)',codes_mat(valmsk,:),ncomp);
pred_code = [ones(nsamp,1),rspmat'] * betaPLS;
%%
imgidx = 71:100;
figure(3);
subplot(121)
montage(G.visualize_latent(codes_mat(imgidx,:)))
subplot(122)
montage(G.visualize_latent(pred_code(imgidx,:)))
%% OLS Decoding
nvalid = numel(valmsk);
n2c_beta = pinv([ones(nvalid,1),rspmat(:,valmsk)']) * codes_mat(valmsk,:);
c2n_beta = pinv([ones(nvalid,1),codes_mat(valmsk,:)]) * rspmat(:,valmsk)';
pred_code_OLS = [ones(nsamp,1),rspmat'] * n2c_beta;
%%
imgidx = 71:110;
figure(4);
tiledlayout(1,2,'TileSp','compact',"Pad",'compact')
nexttile(1);% subplot(121)
montage(G.visualize_latent(codes_mat(imgidx,:)))
title("Original Image in Hessian Selectivity")
nexttile(2);
montage(G.visualize_latent(pred_code_OLS(imgidx,:)))
title("Generated from Reconstructed code from Neural response")

