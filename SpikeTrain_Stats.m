%% Spike Train analysis
load('E:\Monkey_Data\Beto-13082020-004_formatted.mat')
load('E:\Monkey_Data\Beto-13082020-004.mat')
%% Visualize ISI for each channel
for iCh = 1:66
spkT = find(spikeChans(iCh,:));
ISI = spkT(2:end)-spkT(1:end-1);
figure(6);histogram(ISI,0:300,'EdgeColor','none')
title(meta.spikeID(iCh))
pause
end