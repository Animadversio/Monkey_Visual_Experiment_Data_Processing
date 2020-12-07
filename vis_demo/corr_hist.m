figure; hold on
this_c = c(ismember(meta.spikeID,1:32),ismember(meta.spikeID,1:32) ) ; histogram(this_c(this_c<1),20)
idx = 33:48; this_c = c(ismember(meta.spikeID,idx),ismember(meta.spikeID,idx) ) ; histogram(this_c(this_c<1),20)
idx = 49:64; this_c = c(ismember(meta.spikeID,idx),ismember(meta.spikeID,idx) ) ; histogram(this_c(this_c<1),20)
idx = 33:48; this_c = c(ismember(meta.spikeID,1:32),ismember(meta.spikeID,32:64) ) ; histogram(this_c(this_c<1),20)