
Animal = "Beto";Set_Path;
%"220118", "220119", "220225", "220228", "220302","220307", "220309",
%"220311", "220404", "220406"
currows = find(contains(ExpRecord.expControlFN,["220406"])); 
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(currows, Animal, false);
bhvfns = ExpRecord.expControlFN(currows);
saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_Cosine";

%%
% Alfa-23032021-006 and 210323_Alfa_generate_BigGAN_cosine(3) files have different numbers of words! Check inputs
%% Process RFMapping Experiments 
rfidx = contains(bhvfns,"rfMapper");
RFS_col = RF_Calc_Stats_fun(meta_new(rfidx), rasters_new(rfidx), Trials_new(rfidx));
%  Process RF data and get the masks saved to disk. 
for iRF = 1:numel(RFS_col)
    RFStat = RFS_col(iRF);
    maskS = RFStats_indiv_chan_gen_mask(RFStat);
    expdir = fullfile(saveroot, compose("%s-%s-RF",datestr(RFStat.meta.datetime,"yyyy-mm-dd"),RFStat.Animal));
    mkdir(expdir)
    save(fullfile(expdir,'RFStat.mat'),'RFStat')
    save(fullfile(expdir,'maskStat.mat'),'maskS')
end

%% Process the Selectivity Representation Ecoding Experiments 
selidx = contains(bhvfns,"selectivity_basic") & ~cellfun(@isempty,meta_new');
SelS_col = selectivity_Collect_Stats_fun(meta_new(selidx), rasters_new(selidx), Trials_new(selidx));
% for i = 1:numel(SelS_col)
%     seldir = fullfile(saveroot,SelS_col(i).meta.fdrnm);
%     SelS_col(i).meta.figdir = seldir; ReprStat = SelS_col(i);
%     mkdir(seldir);save(fullfile(seldir,"ReprStat.mat"),'ReprStat')
%     fprintf("Selectivity Exp stats saved to %s\n",fullfile(seldir,"ReprStat.mat"))
% end
%% Visualize response distribution of all channels 
visusalize_resp_distri_allchan(SelS_col);

%% Extract and Visualize Cosine Experiments. 
evoidx = contains(bhvfns,"generate_BigGAN_cosine") & ~cellfun(@isempty,meta_new');
CStats = Evol_Cosine_Collect_Stats_fun(meta_new(evoidx), rasters_new(evoidx), Trials_new(evoidx));
%%
CStats(10).meta % "2022-02-28-Beto-01-cosine_V1V4"
%% Visualizations
visualize_TargetImg(CStats(:))
visualize_Cosine_PopEvol(CStats(:),9);
visualize_Cosine_score_traj(CStats(:),10);
visualize_PCCosine_imageEvol(CStats(:),7,8)
calc_Cosine_RFmask_fun(CStats(:))
visusalize_resp_distri_allchan(SelS_col);
%% 
animate_Cosine_Evol_summary(CStats(:),15)

%%
visualize_Cosine_score_traj(CStats(10),10);

%%
msktype = "rect_wt_thr";
visualize_PCCosine_imageEvol_wMask(CStat, alpha_msk_flip, msktype, figmtg)
%%


function visualize_TargetImg(CStats)
for i = 1:numel(CStats)
targnm = CStats(i).targ.target_cfg{1}{2};
repr_dir = fileparts(CStats(i).targ.repr_path); % parent dir of repr path. 
targfullpath = map2fullpath({targnm},repr_dir);
targfullpath = targfullpath{1};
img = imread(targfullpath);
imwrite(img, fullfile(CStats(i).meta.figdir,"targetimg.png"))
copyfile(targfullpath,CStats(i).meta.figdir)
end
end