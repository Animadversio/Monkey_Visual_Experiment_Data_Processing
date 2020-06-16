function net_un = getUntrainedNet(name)
% Get the untrained network initialized with kaiming initialization for
% ReLU nonlinearity
if nargin == 1
    name = "vgg16";
end
% net = vgg16; % obsolete
net_un = vgg16('Weights', 'ImageNet');
unLayers = net_un.Layers;
idxLayersWithWeights = [2 4 7 9 12 14 16 19 21 23 26 28 30 33 36 39];
% unLayers(1).Mean = net.Layers(1).Mean;
for idx = idxLayersWithWeights
%     unLayers(idx).Weights = net.Layers(idx).Weights;
%     unLayers(idx).Bias = net.Layers(idx).Bias;
    unLayers(idx).Weights = kaiming_Weight(unLayers(idx));
    unLayers(idx).Bias = zeros(size(unLayers(idx).Bias),'single');
end
net_un = SeriesNetwork(unLayers);
end
function Weights = kaiming_Weight(Layer)
if isa(Layer,'nnet.cnn.layer.Convolution2DLayer')
    nin = prod(size(Layer.Weights,[1,2,3]));
    Weights = randn(size(Layer.Weights),'single') * single(sqrt(2.0 / nin));
elseif isa(Layer,'nnet.cnn.layer.FullyConnectedLayer')
    nin = prod(size(Layer.Weights,[2]));
    Weights = randn(size(Layer.Weights),'single') * single(sqrt(2.0 / nin));
end
end
function Layer = kaiming(Layer)
if isa(Layer,'nnet.cnn.layer.Convolution2DLayer')
    nin = prod(size(Layer.Weights,[1,2,3]));
    Layer.Weights = randn(size(Layer.Weights),'single') * single(sqrt(2.0 / nin));
elseif isa(Layer,'nnet.cnn.layer.FullyConnectedLayer')
    nin = prod(size(Layer.Weights,[2]));
    Layer.Weights = randn(size(Layer.Weights),'single') * single(sqrt(2.0 / nin));
end
Layer.Bias = zeros(size(Layer.Bias),'single');
end