function I = deepDreamImage_extra(net, layer, channels, varargin)
%deepDreamImage   Visualize network features using Deep Dream.
%
%  Deep Dream is a deep learning feature visualization technique that
%  synthesizes images that strongly activate network layers. Visualizing
%  these images highlight the image features learned by a network. These
%  images are useful for understanding and diagnosing network behavior.
%
%  I = deepDreamImage(network, layer, channel) returns an image
%  that strongly activates the channels of a layer within the network.
%  network can be a SeriesNetwork object or a DAGNetwork object. If net is
%  a SeriesNetwork object, layer can be a numeric index, a string scalar,
%  or a character vector corresponding to one of the layers.
%  If net is a DAGNetwork, layer can be a string scalar or a character
%  vector. channel must be a scalar or a vector of channel indices.
%  When channel is a vector, the activations of each channel are optimized
%  independently.
%
%  The output image, I, is a sequence of grayscale or truecolor (RGB)
%  images stored in a 4-D array. Images are concatenated along the fourth
%  dimension of I such that the image that maximizes the output of
%  channel(k) is I(:,:,:,k).
%
%  I = deepDreamImage(..., Name, Value) specifies additional name-value
%  pair arguments described below.
%
%  'InitialImage'   The image used to initialize Deep Dream. Use this to
%                   see how an image is modified to maximize network layer
%                   activations. The size of the inital image depends on
%                   the selected layer. For layers towards the end of the
%                   network, the image must be the same size or larger than
%                   the network's input size.
%
%                   Default: If you do not specify an image, the function
%                            uses a random initial image with pixel values
%                            drawn from a normal distribution with mean 0
%                            and standard deviation 1.
%
%  'PyramidLevels'  Number of multiresolution image pyramid levels to use
%                   to generate the output image. Increase the number of
%                   pyramid levels to produce larger output images at the
%                   expense of additional computation. Set the number of
%                   levels to 1 to produce an image the same size as
%                   'InitialImage'.
%
%                   Default: 3
%
%  'PyramidScale'   The scale between each pyramid level. Reduce the
%                   pyramid scale while increasing PyramidLevels to
%                   incorporate more fine-grain details into the
%                   synthesized image. This can help generate more
%                   informative images for layers at the beginning of the
%                   network.
%
%                   Default: 1.4
%
%  'NumIterations'  The number of iterations per pyramid level. Increase
%                   the number of iterations to produce more detailed
%                   images at the expense of additional computation.
%
%                   Default: 10
%
%  'OutputScaling'  The type of scaling to apply to the output image. Valid
%                   values are 'linear' or 'none'. Select 'linear' to scale
%                   output pixel values between [0 1]. The output image
%                   corresponding to each layer channel, I(:,:,:,channel),
%                   is scaled independently. Select 'none' to disable
%                   output scaling.
%
%                   Default: 'linear'
%
%  'Verbose'        Set to true to display progress information.
%
%                   Default: true
%
%  'ExecutionEnvironment'  The execution environment for the network.
%                          Specify what hardware resources to use to run
%                          the optimization. Possible values are:
%
%                          'auto' - Use a GPU if available, otherwise use
%                                   the CPU.
%                          'gpu'  - Use the GPU. To use a GPU, you must
%                                   have Parallel Computing Toolbox(TM),
%                                   and a CUDA-enabled NVIDIA GPU with
%                                   compute capability 3.0 or higher. If a
%                                   suitable GPU is not available, the
%                                   function returns an error.
%                          'cpu'  - Use the CPU.
%
%                           Default: 'auto'
%
% Notes:
% ------
% - This function implements a version of Deep Dream that uses a
%   multi-resolution image pyramid and Laplacian Pyramid Gradient
%   Normalization to generate high-resolution images.
%
% - Selecting ReLU or dropout layers for visualization may not produce
%   useful images because of the effect those layers have on the network
%   gradients.
%
% - To visualize classification layer features, select the last fully
%   connected layer before the classification layer.
%
% Example: Visualize convolutional neural network features
% ----------------------------------------------------------
% % Train a convolutional neural network on digit data.
% [XTrain, TTrain] = digitTrain4DArrayData;
%
% layers = [ ...
%     imageInputLayer([28 28 1])
%     convolution2dLayer(5,20)
%     reluLayer()
%     maxPooling2dLayer(2,'Stride',2)
%     fullyConnectedLayer(10)
%     softmaxLayer()
%     classificationLayer()];
%
% options = trainingOptions('sgdm', 'Plots', 'training-progress');
% net = trainNetwork(XTrain, TTrain, layers, options);
% net.Layers
%
% % Select the last fully-connected layer for visualization. Choosing this
% % layer illustrates what images the network thinks look like digits.
% layer = 'fc'
%
% % Select all ten output channels for visualization.
% channels = 1:10;
%
% % Generate images.
% I = deepDreamImage(net, layer, channels);
%
% % Display the image corresponding to digit 0.
% figure
% imshow(I(:,:,:,1))
%
% See also trainNetwork, SeriesNetwork, alexnet, vgg16, vgg19.

% References
% ----------
% DeepDreaming with TensorFlow :
%    https://github.com/tensorflow/tensorflow/blob/master/tensorflow/examples/tutorials/deepdream/deepdream.ipynb

% Copyright 2016-2019 The MathWorks, Inc.

try
    I = doDeepDreamImage(net, layer, channels, varargin{:});
catch ex
    throw(ex)
end
end

function I = doDeepDreamImage(net, layer, channels, varargin)
[reducedInternalNetwork, params] = iParseInputs(net, layer, channels, varargin{:});

numChannels = reducedInternalNetwork.OutputLayers{1}.ExternalCustomLayer.nWeights;
%numel(...
    %reducedInternalNetwork.OutputLayers{1}.ExternalCustomLayer.Channels);
   % FIXED: Use not the channels number but the groups of weight to define the
   % input!

% Convert learnable parameters to the correct format
GPUShouldBeUsed = nnet.internal.cnn.util.GPUShouldBeUsed( ...
    params.ExecutionEnvironment );

if GPUShouldBeUsed
    executionSettings = struct( ...
        'executionEnvironment', 'gpu', ...
        'useParallel', false );
else
    executionSettings = struct( ...
        'executionEnvironment', 'cpu', ...
        'useParallel', false );
end

reducedInternalNetwork = reducedInternalNetwork.prepareNetworkForTraining( executionSettings );

% Move data to the GPU.
if GPUShouldBeUsed
    X = gpuArray(single(repmat(params.InitialImage, [1,1,1,numChannels])));
else
    X = single(repmat(params.InitialImage, [1,1,1,numChannels]));
end

I = nnet.internal.cnn.visualize.deepDreamImageLaplacianNorm(...
    reducedInternalNetwork, X, ...
    params.NumIterations, ...
    params.PyramidLevels, ...
    params.PyramidScale, ...
    params.TileSize,...
    params.StepSize,...
    params.LaplacianGradNorm, ...
    params.Verbose);

% Scale output for display.
if ~strcmpi(params.OutputScaling, 'none')
    I = iScaleImage(I);
end

% Return data on the host.
I = gather(I);
end

%--------------------------------------------------------------------------
function [reducedInternalNetwork, params] = iParseInputs(net, layer, channels, varargin)

iCheckNetwork(net);

% This will allow us to find the input, and from that the input size
internalNetwork = nnet.internal.cnn.util.externalToInternalNetwork(net);
inputSize = internalNetwork.InputLayers{1}.InputSize;

iValidateDAGIsStringOrCharArray(net, layer);

layerIdx = iValidateNetworkLayerNameOrIndex(net, layer, mfilename);

iValidateLayerHasSingleOutput(net.Layers(layerIdx));

iCheckChannels(net, layer, channels, inputSize);

p = inputParser;
addParameter(p, 'InitialImage',    []);
addParameter(p, 'NumIterations',   10);
addParameter(p, 'PyramidLevels',   3);
addParameter(p, 'PyramidScale',    1.4);
addParameter(p, 'Verbose',         true);
addParameter(p, 'OutputScaling',   'linear');
addParameter(p, 'ExecutionEnvironment', 'auto');
addParameter(p, 'Weights', ones(1,length(channels)));
addParameter(p, 'channels', 1:channels);

parse(p, varargin{:});

userInput = p.Results;

iCheckNumIterations(userInput.NumIterations);

iCheckPyramidLevels(userInput.PyramidLevels);

iCheckPyramidScale(userInput.PyramidScale);

iCheckLogicalParameter(userInput.Verbose, 'Verbose');

outputScaling = iCheckOutputScaling(userInput.OutputScaling);

params.ExecutionEnvironment = iValidateAndReturnExecutionEnvironment( ...
    userInput.ExecutionEnvironment);

% Reduce the network and add an optimizechannelaverage layer at the end
[reducedInternalNetwork, reducedNetwork] = ...
    reduceNetworkForDeepDream_extra(net,...
    layerIdx, userInput.channels, userInput.Weights);  % add the weights to the final layer

params.InitialImage = iValidateAndReturnInitialImage(p, reducedNetwork, layer);

% User visible parameters
params.NumIterations = double(userInput.NumIterations);
params.PyramidLevels = double(userInput.PyramidLevels);
params.PyramidScale  = double(userInput.PyramidScale);
params.Verbose       = logical(userInput.Verbose);
params.OutputScaling = convertStringsToChars(outputScaling);

% Internal parameters
params.TileSize          = max(inputSize(1),inputSize(2));
params.StepSize          = 1;
params.LaplacianGradNorm = true;

end

%--------------------------------------------------------------------------
function iCheckNetwork(net)
validateattributes(net, {'SeriesNetwork','DAGNetwork'}, {'scalar'}, mfilename, 'network', 1);

assertSingleInputSingleOutputNetwork( net );

assertNoSequenceInputLayer( net );

assertNo3DInputLayer( net );
end

%--------------------------------------------------------------------------
function scaling = iCheckOutputScaling(value)
scaling = validatestring(value, {'none', 'linear'}, ...
    mfilename, 'OutputScaling');
end

%--------------------------------------------------------------------------
function iCheckChannels(net, layerName, channels, inputSize)
% net.Layers(layerIndex) must have specified channels

layerOutputSize = iLayerOutputSize(net, layerName, inputSize);

numChannels = layerOutputSize(3);

validateattributes(channels, {'numeric'}, ...
    {'vector', 'nonempty', 'real', 'nonsparse', 'integer', '>=', 1, '<=' numChannels}, ...
    mfilename, 'channel', 3);
end

%--------------------------------------------------------------------------
function iCheckInitialImage(I, sortedNet, layerNumber)

sortedInternalNetwork = nnet.internal.cnn.util.externalToInternalNetwork(sortedNet);
netInputSize = sortedInternalNetwork.InputLayers{1}.InputSize;

validateattributes(I, {'numeric'}, ...
    {'nonempty', 'real', 'nonsparse', 'size', [NaN NaN netInputSize(3)]},...
    mfilename, 'InitialImage');

iCheckInitialImageSize(sortedInternalNetwork, layerNumber, size(I));
end

%--------------------------------------------------------------------------
function iCheckNumIterations(n)
validateattributes(n, {'numeric'}, ...
    {'scalar', 'integer', 'nonempty', 'real', 'positive', 'nonsparse'}, ...
    mfilename, 'NumIterations');
end

%--------------------------------------------------------------------------
function iCheckPyramidLevels(n)
validateattributes(n, {'numeric'}, ...
    {'scalar', 'integer', 'nonempty', 'real', 'positive', 'nonsparse'}, ...
    mfilename, 'PyramidLevels');
end

%--------------------------------------------------------------------------
function iCheckPyramidScale(n)
validateattributes(n, {'numeric'}, ...
    {'scalar', 'nonempty', 'real', '>', 1, 'finite', 'nonsparse'}, ...
    mfilename, 'PyramidScale');
end

%--------------------------------------------------------------------------
function iCheckLogicalParameter(tf,paramName)
validateattributes(tf, {'logical','numeric'},...
    {'nonnan', 'scalar', 'nonempty', 'real','nonsparse'},...
    mfilename,paramName);

end

%--------------------------------------------------------------------------
function iCheckInitialImageSize(internalNetwork,layerNumber,forwardSize)

outputSize = internalNetwork.inferOutputSizesGivenInputSizes({forwardSize}, layerNumber);
outputSize = iExtractFromCell(outputSize, 1:2);

if any(outputSize < 1)
    error(message('nnet_cnn:deepDreamImage:InitialImageNotValidImage'));
end

end

%--------------------------------------------------------------------------
function initialImage = iValidateAndReturnInitialImage(p, reducedNet, layer)
na = nnet.internal.cnn.analyzer.NetworkAnalyzer(reducedNet);
sortedExternalLgraph = na.LayerGraph;
sortedExternalNetwork = assembleNetwork(sortedExternalLgraph);

sortedNames = [na.LayerAnalyzers.Name];

% This check is done because the layer may be an index or a layer name for
% Series networks
if ~isa(layer, 'double')
    layer = find(sortedNames == layer);
end
% If no initial image is provided, then generate one
if ismember('InitialImage', p.UsingDefaults )
    
    receptiveFieldSize = ...
        nnet.internal.cnn.visualize.computeReceptiveFieldSize(...
        reducedNet, layer);
    
    % Generate an initial image of random Gaussian noise
    % with the same mean as the network's ImageInputLayer,
    % so it's representative of images the network expects.
    % If the initial image's mean is too low, the gradients
    % may get zeroed out by the ReLU layers and no optimization
    % will occur.
    if isempty(sortedExternalNetwork.Layers(1).Mean)
        % If the Mean property is empty, add nothing to the normal distribution
        imageMean = 0;
    else
        imageMean = sortedExternalNetwork.Layers(1).Mean;
        
        % Ensure that the mean is taken across the spatial dimensions, as 
        % the receptive field size may not be the same as the input layer 
        % size.
        imageMean = mean(imageMean, [1,2]);
    end
    initialImage = single(randn(receptiveFieldSize) + imageMean);
else
    
    % Else, subtract the mean value of the initial image, to give it
    % mean zero.
    % Load in the initial image, if given
    initialImage = p.Results.InitialImage;
    
    iCheckInitialImage(initialImage, sortedExternalNetwork, layer);
    
    % Cast to precision used by optimization process.
    initialImage = single(initialImage);
    
end

end

%--------------------------------------------------------------------------
function validatedCharArray = iValidateAndReturnExecutionEnvironment(originalCharArray)
validExecutionEnvironments = {'auto', 'gpu', 'cpu'};
validatedCharArray = validatestring( ...
    originalCharArray, validExecutionEnvironments, ...
    'deepDreamImage', 'ExecutionEnvironment');
end

%--------------------------------------------------------------------------
% Scale image such that the max value maps to 1 and min value maps to zero,
% similar to MAT2GRAY. Note that each layer channel is scaled independently.
%--------------------------------------------------------------------------
function B = iScaleImage(I)

numChannels = size(I,4);
B = zeros(size(I), 'like', I);
for i = 1:numChannels
    A = I(:,:,:,i);
    range = [min(A(:)) max(A(:))];
    delta = 1 ./ (range(2) - range(1));
    B(:,:,:,i) = delta * I(:,:,:,i) - range(1) * delta;
end
end

%--------------------------------------------------------------------------
function outputSize = iLayerOutputSize(net, layerName, inputSize)

initialImage = single(randn(inputSize));

acts = net.activations(initialImage, layerName);
outputSize = size(acts);

end

%--------------------------------------------------------------------------
function layerIndex = iValidateNetworkLayerNameOrIndex(...
    net, layerNameOrIndex,fname)

% This will be fed into nnet.internal.cnn.layer.Layer.findLayerByName which
% will output a cellarray of the layer names in the right order
internalLayers = nnet.cnn.layer.Layer.getInternalLayers(net.Layers);

% This will allow us to find the input and output layers of a network
internalNetwork = nnet.internal.cnn.util.externalToInternalNetwork(net);

layerNameOrIndex = convertStringsToChars(layerNameOrIndex); % Convert string to char
if ischar(layerNameOrIndex)
    name = layerNameOrIndex;
    
    [layerIndex, layerNames] = nnet.internal.cnn.layer.Layer.findLayerByName(internalLayers, name);
    
    % first and last layer are not valid.
    firstAndLastLayers = {internalNetwork.OutputLayers{1}.Name,...
        internalNetwork.InputLayers{1}.Name};
    
    % Display a clear error message when user picks the first or last layer
    if any(ismember(firstAndLastLayers, name))
        error(message('nnet_cnn:deepDreamImage:FirstOrLastLayer'));
    end
    
    firstAndLastLayerIdx = ismember(layerNames, firstAndLastLayers);
    layerNames = layerNames(~firstAndLastLayerIdx);
    
    iValidateLayerName( name, layerNames );
    
else
    % first and last layer index is not allowed.
    validateattributes(layerNameOrIndex, {'numeric'},...
        {'positive', 'integer', 'real', 'scalar', '>' 1, '<', numel(internalLayers)}, ...
        fname, 'layer');
    layerIndex = layerNameOrIndex;
end
end

%--------------------------------------------------------------------------
function assertNo3DInputLayer( net )
layer = nnet.internal.cnn.util.findLayers(net, 'nnet.cnn.layer.Image3DInputLayer');
if ~isempty(layer)
    error(message('nnet_cnn:deepDreamImage:NotAvailableFor3D'));
end
end

%--------------------------------------------------------------------------
function assertSingleInputSingleOutputNetwork( net )
% Verify we don't have a multi-input or multi-output network.
if length(net.InputNames) > 1 || length(net.OutputNames) > 1
    error(message("nnet_cnn:deepDreamImage:NotAvailableForMIMO"));
end
end

function assertNoSequenceInputLayer( net )
internalLayers = nnet.internal.cnn.layer.util.ExternalInternalConverter.getInternalLayers( net.Layers );
isRNN = nnet.internal.cnn.util.isRNN( internalLayers );
if isRNN
    error(message('nnet_cnn:deepDreamImage:NotAvailableForRNN'));
end
end

%--------------------------------------------------------------------------
function iValidateDAGIsStringOrCharArray(net, layer)
if isa(net,"DAGNetwork") && ~iIsValidStringOrCharArray(layer)
    error(message('nnet_cnn:deepDreamImage:LayerNameMustBeStringForDAG'));
end
end

%--------------------------------------------------------------------------
function tf = iIsValidStringOrCharArray(x)
tf = nnet.internal.cnn.layer.paramvalidation.isValidStringOrCharArray(x);
end

function iValidateLayerName(layerName, layerNames)
if ~iLayerNameExists(layerName, layerNames)
    error(message('nnet_cnn:deepDreamImage:LayerDoesNotExist',...
        layerName));
end
end

function tf = iLayerNameExists(layerName, layerNames)
tf = any(strcmp(layerNames, layerName));
end

function iValidateLayerHasSingleOutput(layer)
if iLayerHasMultipleOutputs(layer)
    error(message('nnet_cnn:deepDreamImage:MultipleOutputLayer', layer.Name));
end
end

function tf = iLayerHasMultipleOutputs(layer)
if isprop(layer, 'NumOutputs') && layer.NumOutputs > 1
    tf = true;
else
    tf = false; % Custom output layers would fall under this category
end
end

function outOfCell = iExtractFromCell(inCell, spatialDimension)
outOfCell = inCell{1};
if isa(outOfCell, 'cell')
    outOfCell = outOfCell{1};
end
outOfCell = outOfCell(spatialDimension);
end