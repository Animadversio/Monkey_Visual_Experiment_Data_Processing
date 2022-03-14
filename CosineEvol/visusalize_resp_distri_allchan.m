function visusalize_resp_distri_allchan(SelS_col,K)
if nargin == 1, K =25; end
for iTr = 1:numel(SelS_col)
for iCh = 1:numel(SelS_col(iTr).units.spikeID)
chan = SelS_col(iTr).units.spikeID(iCh);
unit = SelS_col(iTr).units.unit_num_arr(iCh);
label = SelS_col(iTr).units.unit_name_arr(iCh);
figh1 = vis_selectivity_static(SelS_col(iTr),true,chan,unit,'top',K, 1);
figh2 = vis_selectivity_static(SelS_col(iTr),true,chan,unit,'bottom',K, 2);
figh3 = vis_selectivity_static(SelS_col(iTr),false,chan,unit,'top',K, 3);
saveallform(SelS_col(iTr).meta.figdir,compose("resp_dist_%s-top",label),figh1,["png"])%,"pdf"
saveallform(SelS_col(iTr).meta.figdir,compose("resp_dist_%s-bottom",label),figh2,["png"])%,"pdf"
saveallform(SelS_col(iTr).meta.figdir,compose("resp_dist_%s-top_unsort",label),figh3,["png"])%,"pdf"
end
end
end