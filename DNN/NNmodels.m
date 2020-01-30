%% Demo to use Matlab Deep Learning Framework to fit a "Linear" model to neural response and 
%% Interpret that by using DeepDream to avctivate the fitted model. 
%% This project finally becomes anothor repo https://github.com/Animadversio/Visual_Neuron_Modelling
%% The matlab Deep Network Visualization Framework doesn't work well....Majorly because the DeepDream function doesn't allow using weights to weight the units. 
%% Many internal CNN functions I've tried to modified are copied here. 
%% 2019 Dec. BXW

%% Import Codes

stimPath = "\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191002a\backup_10_02_2019_14_33_18";
codes_fns = string(ls(fullfile(stimPath,"*.mat")));
load(fullfile(stimPath, codes_fns(1)))
%%
codes_all = [];
img_ids = {};
for i = 1:numel(codes_fns)
    load(fullfile(stimPath, codes_fns(i)))
    codes_all = [codes_all ; codes];
    img_ids = [img_ids, ids];
end
img_ids = string(img_ids');
%%
imgnm = string(Trials.imageName);
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
imds = imageDatastore(cellfun(@(nm)fullfile(stimPath,nm+".jpg"),imgnm(row_gen))); % all the evolved images. 
%%
net = alexnet;
%%
Bsize = 20;tic;
activMat = [];
for i=1:Bsize:1600
activTsr = activations(net, imds.subset(i:i+Bsize-1), "conv4");
activMat = [activMat, reshape(activTsr,[],Bsize)];
end
toc
%% Fetch scores
sampsz = size(activMat,2);
psz = size(activMat,1);
%
prefch_idx = find(meta.spikeID==Trials.TrialRecord.User.prefChan);
gen_row_idx = find(row_gen);
scores = rasters(prefch_idx(2),:,gen_row_idx(1:sampsz));
scores = permute(scores, [3,2,1]);
scores_red = mean(scores(:,50:200),2)-mean(scores(:,1:40),2);
%% Directly regress scores onto the activations
tic;
%[Mdl,FitInfo] = fitrlinear(activMat', scores_red ,'Lambda',1E-4,'Regularization','ridge','Learner','svm',...
%    'Solver','dual','Beta',randn(psz,1)/psz*10);
hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus');
[Mdl,FitInfo] = fitrlinear(activMat', scores_red ,...
    'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',hyperopts);
%     'Lambda',1E-4,'Regularization','ridge','Learner','leastsquare',...
%     'Solver','dual','Beta',randn(psz,1)/psz*10);
toc
%% 
predScore = Mdl.predict(activMat') ;
figure;plot(scores_red);hold on ;plot(predScore)
title("Linear fitting ")
xlabel("id of generated Image")
legend(["Neural Score","Fitted Score"])
%%
max(predScore)
sum(Mdl.Beta~=0)
%%
weights = reshape(Mdl.Beta,Tsrsz(1:3));
figure,imagesc(squeeze(mean(abs(weights),3)));axis image;colorbar
title("Conv4 weight heatmap")
%%
figure,plot(squeeze(mean(weights(3:6,4:7,:),[1,2])))
title("Conv4 hotspot channel weight")
%%
FeatVect = squeeze(mean(weights(3:6,4:7,:),[1,2]));
%%
FeatVect = squeeze(mean(weights(5,5,:),[1,2]));
%%
% I = deepDreamImage_extra(net, "conv4", 1:384, 'Weight', FeatVect);
%%
I3 = deepDreamImage_extra(net, "conv4", 2, 'channels',1:length(FeatVect),'Weight', FeatVect/std(FeatVect),... 
        'NumIterations',40);
figure;imshow(I3)