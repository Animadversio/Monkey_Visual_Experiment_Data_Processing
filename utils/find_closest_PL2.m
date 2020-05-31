function closest_idx = find_closest_PL2(ManifPL2, RFchoices)
choices_idx = 1:length(RFchoices);
target = regexp(ManifPL2,"-(?<date>\d\d\d\d\d\d\d\d)-(?<PLnum>\d\d\d)","names");
RFparse = regexp(RFchoices,"-(?<date>\d\d\d\d\d\d\d\d)-(?<PLnum>\d\d\d)","names");
targetdate = datetime(target.date,'InputFormat','ddMMyyyy');
targetPLnum = str2double(target.PLnum);
dates = cellfun(@(s)datetime(s.date,'InputFormat','ddMMyyyy'),RFparse);
PLnums = cellfun(@(s)str2double(s.PLnum),RFparse);
sameday_mask = (dates == targetdate);
sameday_idx = choices_idx(sameday_mask);
if length(sameday_idx)==0
    fprintf("No RF experiment on the same day!\n choose closest experiment")
    % choose the closest day's all experiment
    gapdur = targetdate - dates;
    insertIdx = sum(gapdur>0);
    %find((gapdur(1:end-1)>0) .* (gapdur(2:end)<0)); % find the position that target date could split the dates in 2 half
    closest_idx = insertIdx;
%     [nearestexp_gap, nearest_idx] = min(abs(targetdate - dates));
%     nearest_idx = find(dates == dates(nearest_idx));
    
else
    [mindist, relidx] = min(abs(PLnums(sameday_mask)-targetPLnum));
    closest_idx = sameday_idx(relidx);
    closest_idx = find(RFchoices(closest_idx) == RFchoices);
end
if length(closest_idx) > 1
    fprintf("Multiple entries with same name are closest, return them all!\n")
end
end
