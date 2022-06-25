%% Align RFs for CosineEvolution and RFmapping 
function Rchans = findAlignUnits(srcchans,CStat,RFStat)
% This utils function find the unit indices that match the srcchans in
% source experiments in he target expeirments. 
% 
% srcchans: a list of indices 
% CStat: Stats for the source exp (like an evolution or cosine evolution)
% RFStat: Stats for the target exp (like an RF mapping exp)
if islogical(srcchans), 
    srcchans=find(srcchans); 
end
Rchans = [];
for iCh = reshape(srcchans,1,[])
    Cchan = CStat.meta.spikeID(iCh);
    Cunit = CStat.meta.unitID(iCh);
    Rch = find(RFStat.meta.spikeID==Cchan &...
               RFStat.meta.unitID==Cunit );
    if isempty(Rch)
        Rch = nan;
    end
    Rchans = [Rchans; Rch];
end
end
