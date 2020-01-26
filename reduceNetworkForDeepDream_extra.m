function [reducedInternalNet, reducedNet] = reduceNetworkForDeepDream_extra(net, layerIdx, channels, weights)
% Given an external SeriesNetwork or DAGNetwork, the (unsorted) index of a 
% layer, and the channels to optimize, reduces the network so only layers 
% upstream of the specified layer are preserved, then appends an 
% OptimizeChannelAverage layer.

%   Copyright 2019 The MathWorks, Inc.

lgraph = iConvertToLayerGraph(net);

connections = lgraph.Connections;
connections = connections.Variables;

layerName = string(lgraph.Layers(layerIdx).Name);

connections = iRemoveBackslashInOutFromNodeName(connections);

nodeNamesToKeep = iFindLayersToKeep(connections, layerName);
nodeNamesToRemove = connections(~ismember(connections, nodeNamesToKeep));

lgraph = lgraph.removeLayers(nodeNamesToRemove);
%lgraph = lgraph.addLayers(nnet.internal.cnn.visualize.layer.OptimizeChannelAverageLayer(channels));
lgraph = lgraph.addLayers(OptimizeChannelWeightAverageLayer(channels, weights));

lgraph = lgraph.connectLayers(layerName, lgraph.Layers(end).Name);

lgraph = iRemoveAugmentations(lgraph);

if isa(net, "DAGNetwork")
    reducedNet = assembleNetwork(lgraph);
elseif isa(net, "SeriesNetwork")
    reducedNet = assembleNetwork(lgraph.Layers);
end

reducedInternalNet = nnet.internal.cnn.util.externalToInternalNetwork(reducedNet);
end

function connections = iRemoveBackslashInOutFromNodeName(connections)
% This takes in the external connections property of a Network and returns
% the same connections property without the /in's and /out's. This is
% necessary when trying to locate the upstream nodes of a layer in the
% function: iFindLayersToKeep

% Extract whatever is before "/", as layer ports may have variable names

hasSlashInOrOut = contains(connections, "/");

connections(hasSlashInOrOut) = extractBefore(connections(hasSlashInOrOut), "/");
end

function list = iFindLayersToKeep(connections, layerName)
% Finds upstream nodes of a node in a DAGNetwork - starts with the layerName
% which we choose and sees which nodes are upstream of it, and then
% recursively does the same for those nodes it's connected to and so
% on until all of those search paths have finally led to the input node.
% This will end up looking like inputHasBeenReached = [1 1 1] for example,
% so all(inputHasBeenReached) = 1 and we get out of the loop.

list = "";
% The first node is determined by it not having anything upstream of it and
% therefore it will never appear in the second column of the connections
% array property
inputHasBeenReached = ~ismember(connections(:,2),layerName);
while all(inputHasBeenReached) == 0
    
    % Because we are extracting layerNames from a column vector, in some
    % cases a node is connected to multiple upstream connections which
    % returns a vertical array with layerNames, therefore we have to
    % be consistent with concatenating it vertically
    list = [list; layerName];
    inputHasBeenReached = ~ismember(connections(:,2),layerName);
    layerName = connections(~inputHasBeenReached,1);
    
end
list(1) = [];
list = unique(list);

end

function lgraph = iRemoveAugmentations(lgraph)
imageInputIdx = arrayfun(@(l)isa(l,'nnet.cnn.layer.ImageInputLayer'),lgraph.Layers);
for i=find(imageInputIdx')
    oldLayer = lgraph.Layers(i);
    newLayer = imageInputLayer( oldLayer.InputSize, ...
        'Name', oldLayer.Name, ...
        'DataAugmentation', 'none', ...
        'Normalization', oldLayer.Normalization, ...
        'NormalizationDimension', oldLayer.NormalizationDimension );
    for s = {'Mean','StandardDeviation','Min','Max'}
        if ~isempty(oldLayer.(s{1}))
            newLayer.(s{1}) = oldLayer.(s{1});
        end
    end
    lgraph = replaceLayer(lgraph, oldLayer.Name, newLayer);
end
end

function lgraph = iConvertToLayerGraph(net)
if isa(net, 'DAGNetwork')
    lgraph = layerGraph(net);
else
    lgraph = nnet.cnn.LayerGraph(net.Layers);
end
end