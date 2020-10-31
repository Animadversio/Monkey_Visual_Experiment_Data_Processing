function unit_name_arr = generate_unit_labels_new(spikeID, unitID, varargin)
if nargin == 3 % decide if we want 0 padding for the matrix
    fmt = varargin{2};
else
    fmt = '%d';
end
unit_name_arr = {}; % name tag for each unit 
for i = 1:numel(spikeID)
    cur_chan = spikeID(i);
    if unitID(i) == 0
        unit_name_arr{i} = [num2str(cur_chan, fmt), 'U'];
    else
        unit_name_arr{i} = [num2str(cur_chan, fmt), char(64+unitID(i))];
    end
end
unit_name_arr = string(unit_name_arr);
end