%% RF on cortex
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa"; Set_Path;
load(fullfile(mat_dir, Animal+"_Manif_RFMaps.mat"), 'MaskStats')
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
figIT = figure(9);clf;hold on
figV1 = figure(10);clf;hold on
figV4 = figure(11);clf;hold on
[ax_arrA,tIT,tV1,tV4] = Cortex_Channel_Tile_Layout_All(figIT, figV1, figV4);
% figITB = figure(12);clf;hold on
% figV1B = figure(13);clf;hold on
% figV4B = figure(14);clf;hold on
% [ax_arrB,tITB,tV1B,tV4B] = Cortex_Channel_Tile_Layout_All(figITB, figV1B, figV4B);
%%
xlin = linspace(-10,10,201);
ylin = linspace(-10,10,201);
%%
for Expi = 34:numel(MaskStats)
[first_idx, last_idx] = unit_id2_chan_idx(1:64, MaskStats(Expi).meta.spikeID, MaskStats(Expi).unit.activ_msk);
for chan = 1:64
    iCh = first_idx(chan);
    if isnan(iCh), set(ax_arrA{chan},'visible', 'off'); continue; end
    if ~MaskStats(Expi).signf_chans(iCh), set(ax_arrA{chan},'visible', 'off');end
    set(ax_arrA{chan},'visible', 'on');
    imagesc(ax_arrA{chan}, xlin, ylin, double(MaskStats(Expi).interpmasks(:,:,iCh)));
    axis image;xlim(ax_arrA{chan},[-8,8]);ylim(ax_arrA{chan},[-8,8])
end
pause
end