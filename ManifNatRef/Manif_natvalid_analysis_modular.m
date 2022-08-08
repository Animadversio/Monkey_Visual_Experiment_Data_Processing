set(groot,'defaultAxesTickLabelInterpreter','none');  
%% Manif_natvalid_analysis_modular
saveroot = "E:\OneDrive - Harvard University\Manifold_Invariance";
Animal = "Beto";Set_Path;
%"220118", "220119", "220225", "220228", "220302","220307", "220309", 
%"220516", "220524", "220531", "220607", "220613"
currows = find(contains(ExpRecord.expControlFN,["220613"])); 
% currows = 437:558;%find(contains(ExpRecord.expControlFN,["220311"])); 
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(currows, Animal, false);
bhvfns = ExpRecord.expControlFN(currows);
%%
rfidx = contains(bhvfns,"rfMapper");
RFS_col = RF_Calc_Stats_fun(meta_new(rfidx), rasters_new(rfidx), Trials_new(rfidx));
%  Process RF data and get the masks saved to disk. 
for iRF = 2:numel(RFS_col)
    RFStat = RFS_col(iRF);
    maskS = RFStats_indiv_chan_gen_mask(RFStat);
    expdir = fullfile(saveroot, compose("%s-%s-RF",datestr(RFStat.meta.datetime,"yyyy-mm-dd"),RFStat.Animal));
    mkdir(expdir)
    save(fullfile(expdir,'RFStat.mat'),'RFStat')
    save(fullfile(expdir,'maskStat.mat'),'maskS')
end
%%

%%
BFEStats = Evol_BigGAN_FC6_Collect_Stats_fun(meta_new(3), rasters_new(3), Trials_new(3));
Evol_BigGAN_FC6_Animation_fun(BFEStats)
%%
selidx = contains(bhvfns,"selectivity_basic") & ~cellfun(@isempty,meta_new');
SelS_col = selectivity_Collect_Stats_fun(meta_new(selidx), rasters_new(selidx), Trials_new(selidx));
% plot the response distribution for each channel. 
visusalize_resp_distri_allchan(SelS_col);
%% Invariance Manifold analysis
visualize_invariance_map_allchan(SelS_col(1))

%% Dynamic visualization of the selectivity experiments 
calc_invariance_fun(SelS_col(1))
%%
calc_sparseness_fun(SelS_col(1))
%%
Fstats = [];
for iCh = 1:numel(ReprStat.units.spikeID)
Fstat = anova_cells(cellfun(@(A)A(iCh,:),ReprStat.resp.trial_col,'uni',0));
Fstats = [Fstats; Fstat];
end
%% find the idx mat 
[objnames, objcats] = object_names();
sfxs = ["left", "med", "right", "small", "large", "background"];
remap_idxmat = sort_invariance_sel_exp(SelS_col(1));

respmat_col = arrayfun(@(remapid)S.resp.meanMat(:, remapid), remap_idxmat,'uni',0);
respmat_tsr = cell2mat(cellfun(@(resp)reshape(resp,1,1,[]),respmat_col,'uni',0));
respmat_col = arrayfun(@(remapid)S.resp.semMat(:, remapid), remap_idxmat,'uni',0);
respsem_tsr = cell2mat(cellfun(@(resp)reshape(resp,1,1,[]),respmat_col,'uni',0));

%%
addrix = [1,2,3,4,2,5,6,2];
figh = 11;figure(figh);set(figh,'pos',[680   402   780   560])
figh2 = 12;figure(figh2);set(figh2,'pos',[680   402   630   590])
objN = numel(objcats); 
clrs = brewermap(objN,'Spectral');
for iCh = 1:numel(S.units.spikeID)
figure(figh);clf;
imagesc(respmat_tsr(addrix,:,iCh))
axis equal tight
yticklabels(sfxs(addrix))
xticks(1:numel(objnames))
xticklabels(objcats)
cb = colorbar();
title(compose("%s Object Response Invariance",S.units.unit_name_arr(iCh)))
hline([3.5,6.5],'r')
saveas(figh, fullfile(S.meta.figdir,...
    compose("tfm_tuning_map_%s.png",S.units.unit_name_arr(iCh))))
%%
figure(figh2);clf;
% rowidx2plot = [1,2,3,4,2,5,6,2];
bsl_ch = S.resp.bslmean(iCh);
bslsem_ch = S.resp.bslsem(iCh);
respmat_ch = [respmat_tsr(:,:,iCh)];% nan(1,objN)];
respsem_ch = [respsem_tsr(:,:,iCh)];% nan(1,objN)];
xarr = [1,2,3, nan, 4,5,6, nan, 7,8];
axes('ColorOrder',clrs,'NextPlot','replacechildren')
for obji = 1:objN
lineprop = {'Color', [clrs(obji,:),0.7], 'LineWidth', 1.5};
shadedErrorBar([1,2,3], respmat_ch([1,2,3],obji), respsem_ch([1,2,3],obji), ...
    'LineProp',lineprop)
lineprop = cat(2, lineprop, {'HandleVis', 'off'});
shadedErrorBar([4,5,6], respmat_ch([4,2,5],obji), respsem_ch([4,2,5],obji), ...
    'LineProp',lineprop)
shadedErrorBar([7,8], respmat_ch([6,2],obji), respsem_ch([6,2],obji), ...
    'LineProp',lineprop)
end
xlim([0.5,8.5])
xticks(1:numel(addrix))
xticklabels(sfxs(addrix))
hline(bsl_ch,'k')
hline(bsl_ch+bslsem_ch*[-1,1],'-.k')
xlabel("Transformation");ylabel("Firing rate")
legend(objcats,'interpreter','none','location','bestoutside')
title(compose("%s Object Response Invariance",S.units.unit_name_arr(iCh)))
saveas(figh2, fullfile(S.meta.figdir,...
    compose("tfm_tuning_curve_%s.png",S.units.unit_name_arr(iCh))))
end
% respmat_2plot = respmat_ch(rowidx2plot,:);
% respsem_2plot = respsem_ch(rowidx2plot,:);
% [respmat_ch([1,2,3],:);
%                  nan(1,size(respmat_ch,2));
%                  respmat_ch([4,2,5],:);
%                  nan(1,size(respmat_ch,2));
%                  respmat_ch([6,2],:);];
% plot(xarr, respmat_2plot,'LineWidth',2)
% for obji = 1:objN
% shadedErrorBar(xarr, respmat_2plot(:,obji), respsem_2plot(:,obji), ...
%     'LineProp',{'Color', [clrs(obji,:),0.7]})
% end
function remap_idxmat = sort_invariance_sel_exp(S)
imgname_uniq = S.stim.imgname_uniq;
sfxs = ["_left", "_med", "_right", "_small", "_lrg", "_background"];
objnames = strip(string(ls("N:\Stimuli\Invariance\Project_Manifold\*.jpg")));
objparts = split(objnames,'.');
objnames = objparts(:, 1)';
remap_idxmat = zeros(numel(sfxs), numel(objnames));
for i = 1:numel(sfxs)
	for j = 1:numel(objnames)
		sfx = sfxs(i);
		objname = objnames(j);
		imgnm = string(objname) + string(sfx); 
		remap_idxmat(i, j) = find(strcmp(imgname_uniq, imgnm)); 
	end
end
end
function [objnames, objcats, mapper] = object_names()
objnames = strip(string(ls("N:\Stimuli\Invariance\Project_Manifold\*.jpg")));
objparts = split(objnames,'.');
objnames = objparts(:, 1)';
objcats = ["birdcage", "squirrel", "monkeyL", "monkeyM", "gear", "guitar", "fruits", "pancake", "tree", "magiccube"];
mapper = containers.Map(["bing_birdcage_0001_seg",
"n02357911_47_seg",
"n02487347_3641_seg",
"n02487547_1709_seg",
"n03430551_637_seg",
"n03716887_63_seg",
"n07753592_1991_seg",
"n07880968_399_seg",
"n13912260_18694_seg",
"n13914608_726_seg"], objcats);
end