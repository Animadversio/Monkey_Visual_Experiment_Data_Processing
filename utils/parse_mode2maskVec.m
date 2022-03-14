function [chanmsk, targ_area] = parse_mode2maskVec(mode, array_layout, spikeID)
% Parse the score_mode string into a channel mask of which channels are
% used in the exp.
% chanmsk = parse_mode2mask("corr_V4IT","Alfa");
% chanmsk = parse_mode2mask("corr_V4IT","Beto_new");
MAXUNUM = 4;
if nargin <=1, array_layout="Alfa"; end
if nargin <=2, spikeID=[1:64]'; end
if array_layout == "Alfa"
    ITmsk = (spikeID<=32);
    V4msk = (spikeID>=49);
    V1msk = (spikeID<49) & (spikeID>32);
elseif array_layout == "Beto_new"
    ITmsk = (spikeID>=33);
    V4msk = (spikeID<=16);
    V1msk = (spikeID<=32) & (spikeID>16);
elseif array_layout == "Beto"
    fprintf("Warning: when was the experiment? the new array or old one? (Default is Beto_new for new array)\n")
    ITmsk = (spikeID>=33);
    V4msk = (spikeID<=16);
    V1msk = (spikeID<=32) & (spikeID>16);
    keyboard
end

targ_area = "";
if ~contains(mode,["V1","V4","IT"]) || contains(mode,"All")
    chanmsk = ones(size(spikeID),'logical');
    targ_area = "V1V4IT";
else
    chanmsk = zeros(size(spikeID),'logical');
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