function RF = CNN_receptive_field_for_unit(net, layer, unit)
% layer = 'conv4';
% unit = [1, 1];
% RF = CNN_receptive_field_for_unit(net, "pool5", [3,3]);
% Result is formed by [LB; UB] 2 by 2 matrix. 
result = CNN_receptive_field(net);
layerNames = string({net.Layers.Name});
inputSize = net.Layers(1).InputSize;
idx = find(contains(result.Name,layer));
if isempty(idx)
    error("Layer not found.")
end
if length(idx) > 1
    error("Multiple Layers sharing the same name.")
end
if isa(net.Layers(idx), 'nnet.cnn.layer.FullyConnectedLayer') || all(result.feat_tsr_size(idx,1:2)==1)
    LB = [1, 1];
    UB = inputSize(1:2);
else
    feat_map_size = result.feat_tsr_size(idx,1:2);
    isfeasible = all(feat_map_size >= unit) && all(1 <= unit);
    if ~isfeasible
        error("Unit number out of bound, not-feasible");
    end
    LB = (unit - 1) .* result.jump(idx,:) + result.start(idx,:) - result.RF(idx,:) / 2 + 1;
    UB = (unit - 1) .* result.jump(idx,:) + result.start(idx,:) + result.RF(idx,:) / 2;
    LB = max(LB, [1, 1]);
    UB = min(UB, inputSize(1:2));
end 
RF = [LB; UB];
disp(RF)
    

