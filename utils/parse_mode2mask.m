function [chanmsk, targ_area] = parse_mode2mask(mode, array_layout, chan_arr, unit_arr)
% Parse the score_mode string into a channel mask of which channels are
% used in the exp.
% chanmsk = parse_mode2mask("corr_V4IT","Alfa");
% chanmsk = parse_mode2mask("corr_V4IT","Beto_new");
MAXUNUM = 4;
if nargin <=1, array_layout="Alfa"; end
if nargin <=2, chan_arr=[1:64]'; unit_arr = 0:MAXUNUM; end
if array_layout == "Alfa"
    ITmsk = (chan_arr<=32);
    V4msk = (chan_arr>=49);
    V1msk = (chan_arr<49) & (chan_arr>32);
elseif array_layout == "Beto_new"
    ITmsk = (chan_arr>=33);
    V4msk = (chan_arr<=16);
    V1msk = (chan_arr<=32) & (chan_arr>16);
elseif array_layout == "Beto"
    fprintf("Warning: when was the experiment? the new array or old one? (Default is Beto_new for new array)\n")
    ITmsk = (chan_arr>=33);
    V4msk = (chan_arr<=16);
    V1msk = (chan_arr<=32) & (chan_arr>16);
    keyboard
end
chanMat = repmat(chan_arr,1,numel(unit_arr));
unitMat = repmat(unit_arr,numel(chan_arr),1);
targ_area = "";
if ~contains(mode,["V1","V4","IT"]) || contains(mode,"All")
    chanmsk = ones(size(chanMat),'logical');
    targ_area = "V1V4IT";
else
    chanmsk = zeros(size(chanMat),'logical');
    if contains(mode,"IT")
        chanmsk = chanmsk | (ITmsk);
        targ_area=targ_area+"IT";
    end
    if contains(mode,"V4")
        chanmsk = chanmsk | (V4msk);
        targ_area=targ_area+"V4";
    end
    if contains(mode,"V1")
        chanmsk = chanmsk | (V1msk);
        targ_area=targ_area+"V1";
    end
end
end