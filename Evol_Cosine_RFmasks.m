Set_Path;
RFdir = "O:\RFstats";
refdir = "N:\Stimuli\2020-CosineEvol\RefCollection";
%%

expdir="O:\Evol_Cosine\2021-03-02-Alfa-03-MSE_ITV4";
%%
load(fullfile(expdir, 'expStat.mat'),'expStat')
load(fullfile(expdir, 'ExpImMatchStat.mat'),'ExpImMatch')
datastr = gen_rfdatestr(expStat.meta);
load(fullfile(RFdir, compose('Alfa_%s_RFStat.mat',datastr)))
% figh = RF_contour_plot(RFStat, "edge", targmask);
%%
imgsize = expStat.evol.imgsize;
imgpos = expStat.evol.imgpos;
baseMask = expStat.targ.baseMask; % F mask and get rid of 0 units.
chan_arr = (1:64)';
areamsk = area2chanmsk(expStat.evol.targ_area);
targmask = baseMask & areamsk; 
iCh_mat = zeros(size(baseMask));
targmask_vec = zeros(size(RFStat.unit.chan_num_arr), 'logical');
for iCh = 1:numel(RFStat.unit.chan_num_arr)
    chan = RFStat.unit.chan_num_arr(iCh);
    unit = RFStat.unit.unit_num_arr(iCh);
    iCh_mat(chan,unit+1)=iCh;
    if targmask(chan,unit+1), targmask_vec(iCh) = true; end
end
targmsk_chans = iCh_mat(targmask)';
if any(targmsk_chans==0)
   fprintf("Warning!!\n") 
end
targmsk_chans(targmsk_chans==0)=[];
%%
figh = RF_contour_plot(RFStat, "edge", targmask_vec);
saveallform(expdir,"targetPopul_RF_contours",figh)
%% set the x, y axes aligned with the image size
axs = findobj( get(figh, 'Children'), '-depth', 1, 'type', 'axes');
for ax = axs
    xlim(ax, [-0.5,0.5]*imgsize + imgpos(1))
    ylim(ax, [-0.5,0.5]*imgsize + imgpos(2))
end
saveallform(expdir,"targetPopul_RF_contours_imgcoord",figh)

function chanmsk = area2chanmsk(area_str)
% spikeID, the mapping from iCh to the actual channel number.
chanmsk = zeros(64,1,'logical');
chan_arr = (1:64)';
if contains(area_str,"V1")
chanmsk = chanmsk | ((chan_arr <= 48) & (chan_arr >= 33));
end
if contains(area_str,"V4")
chanmsk = chanmsk |  (chan_arr >= 49);
end
if contains(area_str,"IT")
chanmsk = chanmsk |  (chan_arr <= 32);
end
end

function datestr = gen_rfdatestr(meta)
prts = split(meta.expControlFN,"_");
datestr = ['20',prts{1}];
end