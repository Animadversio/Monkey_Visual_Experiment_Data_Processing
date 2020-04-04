% check response and change the unit_name_arr. An advanced version of
% generate_unit_labels.m
function [activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, spikeID, rasters, varargin)
% you can put `savepath` in `varargin`
if nargin == 5 % decide if we want 0 padding for the matrix
    fmt = varargin{2};
else
    fmt = '%d';
end

%EMPTYTHR = 800;
EMPTYFRTHR = 0.25;
empty_msk = mean(rasters(:,:,:),[2,3]) < EMPTYFRTHR; 
empty_unit_tmp = find(empty_msk); % channel index in empty mask.
for l = 1:numel(empty_unit_tmp) % if there is only one channel then it is not unsorted channel
    cur_chan = spikeID(empty_unit_tmp(l));
    if sum(spikeID==cur_chan) == 1 % if there is only one unit it cannot be empty
        empty_msk(empty_unit_tmp(l)) = 0; %it's no longer empty
    end
    if find(find(spikeID == cur_chan)==empty_unit_tmp(l))~=1% the unsorted unit always happens the first in it's channel 
        empty_msk(empty_unit_tmp(l)) = 0; %it's no longer empty
    end
end
activ_msk = ~ empty_msk;
empty_labels = unit_name_arr(empty_msk);
if ~ all(contains(empty_labels,"A")) 
    % Heurist rule: empty unit should be the first in its channel, and there is usually an active unit present there. 
    % May be better ways of handling this is process them by hand
    fprintf("Empty Channel original labels :")
    fprintf("%s ",string(empty_labels))
    fprintf("\nPlease Check\n")
    keyboard
end 
% TODO: we may conbine the generate unit_labels here! 
unit_name_arr = repmat("",length(spikeID),1); % name tag for each unit 
unit_num_arr = zeros(length(spikeID),1);
for i = 1:length(spikeID) % point of this loop is to figure out the unit number! 
    cur_chan = spikeID(i); 
    if empty_msk(i) % this channel is empty! (effectively)
        unit_name_arr{i} = [num2str(cur_chan, fmt), 'U']; % Unsorted unit
        unit_num_arr(i) = 0; % means this unit is unsorted unit. 
    else % if this channel is not empty
        unit_num_arr(i) = find(find(spikeID == cur_chan & activ_msk) == i); % in all the active units in this channel find the numbering of this.
        if sum(spikeID(activ_msk) == cur_chan) == 1 % if this is the only unit active (or only unit present)
            unit_name_arr{i} = num2str(cur_chan, fmt);
        else % if multiple unit active, then attach unit label to it.
            unit_name_arr{i} = [num2str(cur_chan, fmt), char(64+unit_num_arr(i))];
        end
    end
end

assert(length(unique(unit_name_arr))==length(unit_name_arr),"Some problem with unit label, Examine the label!")
if nargin >= 4 % If there is varargin, then we parse it as savepath, then write the Unit_Label.txt There
    savepath = varargin{1};
    if isempty(savepath)
        return
    end
    fid = fopen(fullfile(savepath, "Unit_Label.txt"),'w');
    fprintf(fid,'%s\n', unit_name_arr{:});
    fclose(fid);
end
end