function unit_name_arr = generate_unit_labels(spikeID, varargin)
% Generate the unit labels e.g. 17B from the spikeID variable (and save it)

Unit_id = spikeID; % meta.spikeID;
unit_name_arr = {}; % name tag for each unit 
if nargin == 3 % decide if we want 0 padding for the matrix
    fmt = varargin{2};
else
    fmt = '%d';
end
for i = 1:length(Unit_id)
    cur_chan = Unit_id(i);
    if sum(Unit_id == cur_chan) == 1
        unit_name_arr{i} = num2str(cur_chan, fmt);
    else
        cur_chan = Unit_id(i);
        rel_idx = find(find(Unit_id == cur_chan) == i);
        unit_name_arr{i} = [num2str(cur_chan, fmt), char(64+rel_idx)];
    end
end
if nargin > 1 
    savepath = varargin{1};
    if isempty(savepath)
        return
    end
    fid = fopen(fullfile(savepath, "Unit_Label.txt"),'w');
    fprintf(fid,'%s\n', unit_name_arr{:});
    fclose(fid);
end
end 
