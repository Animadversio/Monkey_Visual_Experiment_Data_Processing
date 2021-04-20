%% Manifold Decoding
Set_Path;
%%
mat_dir = "O:\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
load(fullfile(mat_dir, Animal+"_ManifMapVarStats.mat"), 'MapVarStats')
%% 
%%
[PHI, THETA] = meshgrid(-90:18:90.1,-90:18:90.1);
XX = cosd(PHI).*cosd(THETA);
YY = cosd(PHI).*sind(THETA);
ZZ = sind(PHI);
targ_coord = [reshape(XX,[],1),reshape(YY,[],1),reshape(ZZ,[],1)]; % 121 by 3
%%
Expi = 3;si = 1;
% Fstats
nCh = numel(MapVarStats(Expi).units.spikeID);
actmap_col = arrayfun(@(iCh)cellfun(@(A)...
    squeeze(A(iCh,:)),MapVarStats(Expi).manif.act_col{si},'uni',0),...
    1:nCh,'uni',0);
FStats = cellfun(@anova_cells,actmap_col); 
Fmsk = struct2table(FStats).F_P < 0.001; 
%%
actmean = cellfun(@(A)mean(A,2),MapVarStats(Expi).manif.act_col{si},'uni',0);
actmean = cell2mat(reshape(actmean,1,[]));
actmeanmat = actmean';
actmean_sel = actmeanmat(:,Fmsk);
[zact_sel, chmean, chstd] = zscore(actmean_sel,1);
%% OLS decoding.... very much not successful, generalization really bad
Kfold = 120;
CV = cvpartition(121,'KFold',Kfold);
% CV.NumTestSets
result = repmat(struct(),Kfold,1);
for cvi = 1:Kfold
    trainmsk = CV.training(cvi);
    validmsk = CV.test(cvi);
    nTrain = sum(trainmsk);
    nValid = sum(validmsk);
    Ytrain = targ_coord(trainmsk,:);
    Xtrain = zact_sel(trainmsk,:);
    Yvalid = targ_coord(validmsk,:);
    Xvalid = zact_sel(validmsk,:);
    n2c_beta = pinv([ones(nTrain,1),Xtrain]) * Ytrain;
    coord_fit  = [ones(nTrain,1),Xtrain] * n2c_beta;
    coord_pred = [ones(nValid,1),Xvalid] * n2c_beta;
    
    MSEvalid = mean(sum((Yvalid - coord_pred).^2,2),1); 
    MSEtrain = mean(sum((Ytrain - coord_fit) .^2,2),1); 
    fprintf("CV partition %d MSE valid %.3f, train %.3f\n",cvi,MSEvalid,MSEtrain)
    result(cvi).beta = n2c_beta;
    result(cvi).MSEvalid = MSEvalid;
    result(cvi).MSEtrain = MSEtrain;
    result(cvi).coord_fit = coord_fit;
    result(cvi).coord_pred = coord_pred;
    result(cvi).Ytrain = Ytrain;
    result(cvi).Yvalid = Yvalid;
end
MSEvalid_arr = arrayfun(@(R)R.MSEvalid,result);
MSEtrain_arr = arrayfun(@(R)R.MSEtrain,result);
fprintf("%d fold CV MSE valid %.3f(%.3f), train %.3f(%.3f)\n",Kfold,mean(MSEvalid_arr),sem(MSEvalid_arr),...
    mean(MSEtrain_arr),sem(MSEtrain_arr))
%%
Kfold = 6;
CV = cvpartition(121,'KFold',Kfold);
%%
n_features = sum(Fmsk);
network = MLPNet();
network.AddInputLayer(n_features,false);
network.AddHiddenLayer(10,'leakyrelu',true);
% network.AddHiddenLayer(50,'leakyrelu',true);
network.AddOutputLayer(3,'linear',false);
network.NetParams('rate',0.001,'momentum','adam','lossfun','rmse',...
    'regularization','L2');
network.trainable = true;
network.Summary();
%%
cvi = 2;
n_data = sum(CV.training(cvi));
trainmsk = CV.training(cvi);
validmsk = CV.test(cvi);
Y = targ_coord(trainmsk,:);
X = zact_sel(trainmsk,:);
Y_val = targ_coord(validmsk,:);
X_val = zact_sel(validmsk,:);
%%
decodedir = "O:\Manif_PopDecode";
mkdir(decodedir);
diary(fullfile(decodedir,Animal+"_MLP_fit.log"))
for Expi = 1:numel(MapVarStats)
si = 1;
prefchanlab = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
fprintf("\nProcess Exp %d prefchan %s\n",Expi,prefchanlab)
nCh = numel(MapVarStats(Expi).units.spikeID);
actmap_col = arrayfun(@(iCh)cellfun(@(A)...
    squeeze(A(iCh,:)),MapVarStats(Expi).manif.act_col{si},'uni',0),...
    1:nCh,'uni',0);

FStats = cellfun(@anova_cells,actmap_col); 
Fmsk = struct2table(FStats).F_P < 0.001; 
n_features = sum(Fmsk);
fprintf("Informative Chan num %d\n",n_features)

network = MLPNet();
network.AddInputLayer(n_features,false);
network.AddHiddenLayer(10,'leakyrelu',true);
% network.AddHiddenLayer(50,'leakyrelu',true);
network.AddOutputLayer(3,'linear',false);
network.NetParams('rate',0.001,'momentum','adam','lossfun','rmse',...
    'regularization','L2');
network.trainable = true;
network.Summary();

actmean = cellfun(@(A)mean(A,2),MapVarStats(Expi).manif.act_col{si},'uni',0);
actmean = cell2mat(reshape(actmean,1,[]));
actmeanmat = actmean';
actmean_sel = actmeanmat(:,Fmsk);
Kfold = 10;
CV = cvpartition(121,'KFold',Kfold);
[zact_sel, chmean, chstd] = zscore(actmean_sel,1);
[networks, loss_traces, err_trains, err_valids] = train_regress_net_CV(network,zact_sel,targ_coord,CV);
save(fullfile(decodedir,compose(Animal+"_Manif_Exp%02d_MLPdecode.mat",Expi)), "networks", "loss_traces", "err_trains", "err_valids", "actmean", "actmean_sel", "zact_sel", "FStats", "Fmsk", "chmean", "chstd", "CV")
end
diary off;

function newObj = copyObject(obj)
% https://undocumentedmatlab.com/articles/general-use-object-copy
    try
        % R2010b or newer - directly in memory (faster)
        objByteArray = getByteStreamFromArray(obj);
        newObj = getArrayFromByteStream(objByteArray);
    end
end
function [networks, loss_traces, err_trains, err_valids] = train_regress_net_CV(network,X,Y,CV)
networks = [];
loss_traces = {};
err_trains = {};
err_valids = {};
Kfold = CV.NumTestSets;
for cvi = 1:Kfold
    network.init_weights();
    n_data = sum(CV.training(cvi));
    trainmsk = CV.training(cvi);
    validmsk = CV.test(cvi);
    Y_train = Y(trainmsk,:);
    X_train = X(trainmsk,:);
    Y_val = Y(validmsk,:);
    X_val = X(validmsk,:);
    [network, d_loss, err_train, err_valid] = train_regress_net(network,X_train,Y_train,X_val,Y_val,false);
    networks = [networks,copyObject(network)];
    loss_traces{cvi} = d_loss;
    err_trains{cvi} = err_train;
    err_valids{cvi} = err_valid;
end
end

function [network, d_loss, err_train, err_valid] = train_regress_net(network,X,Y,X_val,Y_val,verbose)
if nargin == 5, verbose = false; end
n_data = size(X,1);
rmserr = 1E4;                          % pre-allocate training accuracy
n_batch = 101;                    % Size of the minibatch
max_epoch = 800;                 % Maximum number of epochs
max_batch_idx = ceil(n_data/n_batch);          % Maximum batch index
max_num_batches = max_batch_idx.*max_epoch;     % Maximum number of batches

% Pre-allocate for epoch and error vectors (for max iteration)
epoch = zeros(1,max_num_batches);
d_loss = epoch;
% ce_test = zeros(max_epoch,1);
% ce_train = zeros(max_epoch,1);
% ce_val = zeros(max_epoch,1);
err_train = nan(max_epoch,1);
err_valid = nan(max_epoch,1);

% Initialize iterator and timer
batch_idx = 1;      % Index to keep track of minibatches
epoch_idx = 1;      % Index to keep track of epochs

target_error = 0.001; % Desired classification accuracy
tic;
while ((epoch(batch_idx)<max_epoch)&&(rmserr>target_error))
    
    % Compute current epoch
    epoch(batch_idx+1) = batch_idx*n_batch/n_data;

    % Randomly sample data to create a minibatch
    rand_ind = randsample(n_data,n_batch);

    % Index into input and output data for minibatch
    X_batch = X(rand_ind,:);    % Sample Input layer
    Y_batch = Y(rand_ind,:);    % Sample Output layer
    
    % Train model
    d_loss(batch_idx+1) = network.training(X_batch,Y_batch)./n_batch;
    
    % Only compute error/classification metrics after each epoch
    if ~(mod(batch_idx,max_batch_idx))
        % Compute error metrics for training, test, and validation set
        [err_train(epoch_idx),ce_train(epoch_idx),~]=network.NetworkError(X,Y,'regression');%classification
        
        [err_valid(epoch_idx),ce_val(epoch_idx),~]=network.NetworkError(X_val,Y_val,'regression');
        eval_time = toc;
%         [~,ce_test(epoch_idx),~]=network.NetworkError(X_test,Y_test,'regression');
%       
        if verbose
        fprintf('\n-----------End of Epoch %i------------\n', epoch_idx);
        fprintf('Loss function: %f \n',d_loss(batch_idx+1));
        fprintf('Validation Set RMSE: %f Training Set RMSE: %f \n',err_valid(epoch_idx),err_train(epoch_idx));
        fprintf('Test Set Evaluation Time: %f s\n\n',eval_time);
        end
        rmserr = err_valid(epoch_idx);
        epoch_idx = epoch_idx+1;    % Update epoch index
    end

    % Update batch index
    batch_idx = batch_idx+1;
end
fprintf('-----------End of Training Epoch %i------------\n', epoch_idx-1);
fprintf('Loss function: %f \n',d_loss(batch_idx));
fprintf('Validation Set RMSE: %f Training Set RMSE: %f \n',err_valid(epoch_idx-1),err_train(epoch_idx-1));
fprintf('Total Training time: %f sec\n',eval_time);
end