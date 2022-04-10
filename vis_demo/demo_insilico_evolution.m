%% Get CNN and GAN    
% Setup GAN
% Download `matlabGANfc6.mat` 
%   https://drive.google.com/file/d/1zJkKKgpk0NJB5qAxJRG0JM573D7kiqLt/view?usp=sharing
%
G = FC6Generator("matlabGANfc6.mat"); % put the path to this matfile here.
net = alexnet(); % Use one unit in AlexNet as in silico model of the neuron
% Choose your unit. 
    % ix, iy are locations on feature map.
    % iChan is the channel. 
    % default one is the "neuron" corresponding to gold fish class from
    % AlexNet.
layer = "fc8"; iChan = 2; ix = 1; iy = 1;                        
%% Set up CMAES optimizer
options = struct("init_sigma",3.0);
optim = CMAES_simple(4096, [], options);
disp(optim.opts)
%% Simple Optimization loop, use origin as init gen
init_z = zeros(1,4096);
codes = init_z;
for iGen = 1:100
    imgs = G.visualize(codes);
    acts = activations(net,imgs,layer);
    scores = squeeze(acts(:,:,iChan,:));
    [codes_new] = optim.doScoring(codes,scores,true);
    codes = codes_new;
end
%%
figure;montage(imgs)


%% Optimization loop, use textures as init gen
% Download `texture_init_code.mat` (not mandatory)
%   https://drive.google.com/file/d/1qXEKP_jtkHqqGoSjOnXzLPqmwinCR3qW/view?usp=sharing
options = struct("init_sigma",3.0);
optim = CMAES_simple(4096, [], options);
disp(optim.opts)
data = load("texture_init_code.mat");
init_z = data.codes;
codes = init_z;
for iGen = 1:100
    imgs = G.visualize(codes);
    acts = activations(net,imgs,layer);
    scores = squeeze(acts(:,:,iChan,:));
    [codes_new] = optim.doScoring(codes,scores,true);
    codes = codes_new;
end




%% Full example, with recording and visualization.
options = struct("init_sigma",3.0);
optim = CMAES_simple(4096, [], options);
disp(optim.opts)
data = load("texture_init_code.mat");
init_z = data.codes;
codes = init_z;
codes_all = [];
scores_all = [];
generations = [];
img_traj = {};
for iGen = 1:100
    imgs = G.visualize(codes);
    acts = activations(net,imgs,layer);
    scores = squeeze(acts(:,:,iChan,:));
    [codes_new] = optim.doScoring(codes,scores,true);
    % record some info for analysis
    codes_all = [codes_all; codes];
    scores_all = [scores_all; scores];
    generations = [generations; ones(numel(scores),1)*iGen];
    img_traj{iGen} = G.visualize(mean(codes,1));
    codes = codes_new;
end
% 
figure;
montage(img_traj)
xlabel("Generation"); ylabel("activation")
figure;
scatter(generations, scores_all)
xlabel("Generation"); ylabel("activation")