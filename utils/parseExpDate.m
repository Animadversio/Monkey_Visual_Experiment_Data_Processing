function date = parseExpDate(ephysfn, bhv2fn)
% util function to parse out experiment date from ephys file or bhv file 
    if length(ephysfn) ~= 0 
        parts = split(ephysfn,'-');
        date_str = parts{2};
        ephys_id = str2num(parts{3});
        date = datetime(date_str,'InputFormat','ddMMyyyy');
    elseif length(bhv2fn) ~= 0 
        parts = split(bhv2fn,'_');
        date_str = parts{1};
        date = datetime(date_str,'InputFormat','yyMMdd');
    else
        error("msg")
    end
end