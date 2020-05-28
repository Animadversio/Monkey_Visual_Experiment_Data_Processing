ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
% Animal="Beto";Expi=45;
% ExpType = "Evol"
for layername = ["conv3_1", "conv4_3", "conv5_3"] %"conv5_3";
outfn = fullfile(ccmat_dir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername));
if ExpType=="Evol"
load(outfn,'cc_tsr','MFeat','StdFeat');
elseif ExpType=="Manif"
load(outfn,'corr_tsr','MFeat','StdFeat');
cc_tsr = corr_tsr;
end
variableInfo = who('-file', outfn);
if ismember('wdw_vect', variableInfo), load(outfn,'wdw_vect'); end
savedir = fullfile(ccmat_dir,compose("%s_%s_Exp%d",Animal,ExpType,Expi));
%%
cc_tsr_sum.L1 = squeeze(mean(abs(cc_tsr),3));
cc_tsr_sum.L2 = squeeze(sqrt(mean(cc_tsr.^2,3)));
cc_tsr_sum.max = squeeze(max(cc_tsr,[],3));
%% 
% plot_tsr = cc_tsr_L1;Meanstr = "L1";
for mean_meth = ["L1"]%,"L2","max"]
plot_tsr = getfield(cc_tsr_sum, string(mean_meth));
Meanstr = string(mean_meth);
v = VideoWriter(fullfile(savedir,compose("%s_%s_Exp%d_%s_cc_%s.mp4",Animal,ExpType,Expi,layername,Meanstr)));
v.FrameRate = 2;open(v);
CMIN = prctile(plot_tsr, [ 2], 'all');
CMAX = prctile(plot_tsr, [98], 'all'); 
figure(18);set(18,'Position',[0   435   552   543])
IMS = imagesc(plot_tsr(:,:,1));axis image
caxis([CMIN,CMAX]);colorbar;
for fi = 1:size(plot_tsr,3)
    wdw = [1, 20] + 10 * (fi - 1);
    IMS.CData = plot_tsr(:,:,fi);
    IMS.Parent.Title.String = sprintf("Exp %d Pref chan %d\n %s CorreCoef in of VGG16 %s feature\n with [%d,%d] ms firing rate", ...
            Expi, EStats(Expi).units.pref_chan, Meanstr, layername, wdw(1), wdw(2));
    pause(0.2)
    Fs = getframe(18);
    writeVideo(v,Fs);
    drawnow;
end
close(v);
end

end
