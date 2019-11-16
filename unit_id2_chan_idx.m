function [first_idx, last_idx] = unit_id2_chan_idx(unit_arr, spikeID)
first_idx = zeros(length(unit_arr),1);
last_idx = zeros(length(unit_arr),1);
for i = 1:numel(unit_arr)
    unit = unit_arr(i);
    first_idx(i) = find(spikeID==unit,1,'first');
    last_idx(i) = find(spikeID==unit,1,'last');
end
end