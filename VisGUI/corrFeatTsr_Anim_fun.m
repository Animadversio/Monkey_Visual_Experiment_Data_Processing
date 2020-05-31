%% 
function corrFeatTsr_Anim_fun(EStats,ExpType,Animal,Expi,option)
if nargin==4
    option = struct("save",true);
    option.result_dir = "E:\OneDrive - Washington University in St. Louis\Evol_Manif_Movies";
    option.cent = []; % can be used to specify the centor of exploration. 
end

layerlist = ["conv1_2","conv2_2","conv3_1", "conv4_3", "conv5_3"];
sum_method = ["L1"];%,"L2","max"];

ccmat_dir = "E:\OneDrive - Washington University in St. Louis\CNNFeatCorr";
for layername = layerlist%["conv3_1", "conv4_3", "conv5_3"] %"conv5_3";
outfn = fullfile(ccmat_dir, compose("%s_%s_Exp%d_%s.mat",Animal,ExpType,Expi,layername));
variableInfo = who('-file', outfn);
if ismember('wdw_vect', variableInfo), load(outfn,'wdw_vect'); end
if ExpType=="Evol"
load(outfn,'cc_tsr', 'MFeat', 'StdFeat');
elseif ExpType=="Manif"
if ismember('corr_tsr', variableInfo)
    load(outfn,'corr_tsr', 'MFeat', 'StdFeat');
    cc_tsr = corr_tsr;
else
    load(outfn,'cc_tsr', 'MFeat', 'StdFeat');
end
end
savedir = fullfile(ccmat_dir,compose("%s_%s_Exp%d",Animal,ExpType,Expi));
%%
cc_tsr_sum.L1 = squeeze(mean(abs(cc_tsr),3));
cc_tsr_sum.L2 = squeeze(sqrt(mean(cc_tsr.^2,3)));
cc_tsr_sum.max = squeeze(max(cc_tsr,[],3));
%% 
% plot_tsr = cc_tsr_L1;Meanstr = "L1";
for mean_meth = sum_method%,"L2","max"]
plot_tsr = getfield(cc_tsr_sum, string(mean_meth));
Meanstr = string(mean_meth);
if option.save
v = VideoWriter(fullfile(savedir,compose("%s_%s_Exp%d_%s_cc_%s.avi",Animal,ExpType,Expi,layername,Meanstr)));
v.FrameRate = 2;open(v);
fprintf("Writing video to %s\n",fullfile(savedir,compose("%s_%s_Exp%d_%s_cc_%s.avi",Animal,ExpType,Expi,layername,Meanstr)))
end
CMIN = prctile(plot_tsr, [ 2], 'all');
CMAX = prctile(plot_tsr, [98], 'all'); 
figure(18);set(18,'Position',[0   435   552   543])
IMS = imagesc(plot_tsr(:,:,1));axis image
caxis([CMIN,CMAX]);colorbar;
for fi = 1:size(plot_tsr,3)
    wdw = wdw_vect(fi,:); %[1, 20] + 10 * (fi - 1);
    IMS.CData = plot_tsr(:,:,fi);
    IMS.Parent.Title.String = sprintf("%s %s Exp %d Pref chan %d\n %s CorreCoef in of VGG16 %s feature\n with [%d,%d] ms firing rate", ...
            Animal, ExpType, Expi, EStats(Expi).units.pref_chan, Meanstr, layername, wdw(1), wdw(2));
    pause(0.2)
    Fs = getframe(18);
    if option.save, writeVideo(v,Fs); end
    drawnow;
end
if option.save, close(v); end
end
end

end