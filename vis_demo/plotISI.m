for iCh =1:80
SpkT = find(spikeChans(iCh,:));
ISIs = SpkT(2:end) - SpkT(1:end-1);
figure(1);histogram(ISIs,0:200)
disp(iCh)
pause
end