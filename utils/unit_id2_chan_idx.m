function [first_idx, last_idx] = unit_id2_chan_idx(unit_arr, spikeID, activ_msk)
% Map channel number to the idx in rasters. Currently provide the first and
% the last index for **active channel**. (Umpty channels are masked out)
%   unit_arr provides the channels to find
%   spikeID is meta.spikeID
%   activ_msk is the units that are active, judged by
%       `check_channel_active_label` function. 
first_idx = zeros(length(unit_arr),1);
last_idx = zeros(length(unit_arr),1);
for i = 1:numel(unit_arr)
    unit = unit_arr(i);
    if ~ isempty(find(spikeID==unit & activ_msk, 1)) 
        first_idx(i) = find(spikeID==unit & activ_msk,1,'first');
        last_idx(i) = find(spikeID==unit & activ_msk,1,'last');
    else
        first_idx(i) = nan;
        last_idx(i) = nan;
    end
end
end