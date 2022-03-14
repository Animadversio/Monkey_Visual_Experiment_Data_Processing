function [V1msk, V4msk, ITmsk] = get_areamask(spikeID, array_layout)
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
end