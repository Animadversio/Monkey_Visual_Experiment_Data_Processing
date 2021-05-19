function [chanmsk, targ_area] = parse_mode2mask(mode, chan_arr, unit_arr)
% Parse the score_mode string into a channel mask of which channels are
% used in the exp.
% chanmsk = parse_mode2mask("corr_V4IT");
MAXUNUM = 4;
if nargin ==1, chan_arr=[1:64]'; unit_arr = 0:MAXUNUM; end
chanMat = repmat(chan_arr,1,numel(unit_arr));
unitMat = repmat(unit_arr,numel(chan_arr),1);
targ_area = "";
if ~contains(mode,["V1","V4","IT"]) || contains(mode,"All")
    chanmsk = ones(size(chanMat),'logical');
    targ_area = "V1V4IT";
else
    chanmsk = zeros(size(chanMat),'logical');
    if contains(mode,"IT")
        chanmsk = chanmsk | ((chan_arr<=32));
        targ_area=targ_area+"IT";
    end
    if contains(mode,"V4")
        chanmsk = chanmsk | ((chan_arr>=49));
        targ_area=targ_area+"V4";
    end
    if contains(mode,"V1")
        chanmsk = chanmsk | ((chan_arr<49) & (chan_arr>32));
        targ_area=targ_area+"V1";
    end
end
end