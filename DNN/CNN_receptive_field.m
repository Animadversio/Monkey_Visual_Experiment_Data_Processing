function result = CNN_receptive_field(net)
layerNames = string({net.Layers.Name});
inputSize = net.Layers(1).InputSize;
result = repmat(struct("Name","","feat_tsr_size",[],"RF",[],"start",[],"jump",[]), 1, length(net.Layers));
for i = 1:length(net.Layers)
    L = net.Layers(i);
    switch class(net.Layers(i))
        case 'nnet.cnn.layer.ImageInputLayer'
            result(i).feat_tsr_size = net.Layers(i).InputSize;
            result(i).RF = [1, 1];
            result(i).jump = [1, 1];
            result(i).start = [0.5, 0.5];
        case {'nnet.cnn.layer.Convolution2DLayer', 'nnet.cnn.layer.GroupedConvolution2DLayer'}
            kernel_size = net.Layers(i).FilterSize;
            stride = net.Layers(i).Stride;
            padding = net.Layers(i).PaddingSize(1:2);
            % kernel_size, stride, padding = map(check_same, [kernel_size, stride, padding])
            result(i).jump = result(i-1).jump .* stride;
            result(i).RF = result(i-1).RF + (kernel_size - 1) .* result(i-1).jump;
            result(i).start = result(i-1).start + (floor((kernel_size - 1) ./ 2) - padding) .* result(i-1).jump;
            size_prev = result(i-1).feat_tsr_size(1:2);
            size_cur = ceil((size_prev + 2 .* padding - kernel_size + 1) ./ stride);
            if isa(L, 'nnet.cnn.layer.GroupedConvolution2DLayer')
                result(i).feat_tsr_size = [size_cur, net.Layers(i).NumGroups * net.Layers(i).NumFiltersPerGroup];
            else
                result(i).feat_tsr_size = [size_cur, net.Layers(i).NumFilters];
            end
        case 'nnet.cnn.layer.MaxPooling2DLayer'
            kernel_size = net.Layers(i).PoolSize;
            stride = net.Layers(i).Stride;
            padding = net.Layers(i).PaddingSize(1:2);
            % kernel_size, stride, padding = map(check_same, [kernel_size, stride, padding])
            result(i).jump = result(i-1).jump .* stride;
            result(i).RF = result(i-1).RF + (kernel_size - 1) .* result(i-1).jump;
            result(i).start = result(i-1).start + (floor((kernel_size - 1) ./ 2) - padding) .* result(i-1).jump;
            size_prev = result(i-1).feat_tsr_size(1:2);
            size_cur = ceil((size_prev + 2 .* padding - kernel_size + 1) ./ stride);
            result(i).feat_tsr_size = [size_cur, result(i-1).feat_tsr_size(3)]; % doesn't change the filter number
        case 'nnet.cnn.layer.ReLULayer'
            result(i) = result(i-1);
        case 'nnet.cnn.layer.FullyConnectedLayer'
            result(i) = result(i-1);
            result(i).RF = inputSize(1:2);
            result(i).feat_tsr_size = [1, 1, net.Layers(i).OutputSize];
        otherwise
            result(i) = result(i-1);
    % isa(net.Layers(i),'nnet.cnn.layer.Convolution2DLayer')
    end
    result(i).Name = net.Layers(i).Name;
end
result=struct2table(result);
disp(result)
end
%%
