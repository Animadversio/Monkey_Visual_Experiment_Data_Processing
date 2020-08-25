%% Manif_Animation
%% Really compelling visualization
system("subst S: E:\Network_Data_Sync")
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
mat_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'))
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
%% Set up figure stage for the 3 arrays 
% (no need to redo this for different channels. Just use the layout in ax_arr,tIT,tV1,tV4)
tic
% savepath = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Pop_PSTH_Anim";
% savepath = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Pop_PSTH_Anim";
savepath = "S:\Pop_PSTH_Anim";
mkdir(savepath)
figIT = figure(2);
figV1 = figure(3);
figV4 = figure(4);
[ax_arr,tIT,tV1,tV4] = Cortex_Channel_Tile_Layout_All(figIT, figV1, figV4);
toc
%
tic
% paint it once with fake data, so that the layout if correct. 
imgsc_list = {};
for arr_chan = 1:64 
    imgsc_list{arr_chan} = plot_heatmap(ax_arr{arr_chan}, nan(11), [0,1],"nan");
end
toc
%%
savepath = "S:\Pop_PSTH_Anim";
for Expi = 1:length(meta_new)
%load(fullfile("S:\Data-Ephys-MAT",Stats(Expi).meta.ephysFN+"_formatted.mat"))
rasters = rasters_new{Expi};
fprintf("Processing Exp %d...\n",Expi)
%
si = 1; % here we do subspace 1 first! There is subspace 2 and 3
% manif_psth_all = cellfun(@(idx)rasters(:, :, idx), Stats(Expi).manif.idx_grid{si}, "UniformOutput", false);
manif_psth_avg = cellfun(@(idx)mean(rasters(:, :, idx),3), Stats(Expi).manif.idx_grid{si}, "UniformOutput", false);
manif_psth_avg_tsr = cell2mat(reshape(manif_psth_avg,1,1,11,11)); % reshape and transform, so size 86, 200, 11, 11, key data! 
% load in manifold data and pre-compute the data for each frame. PreCompute CLIM for each channel
Wlen=20; shift_step=2; % setting for averaging the frames
start_arr = 0:shift_step:200-Wlen; % array of the start time for the moving window.  
act_map_tsr = zeros(size(rasters,1),11,11,length(start_arr)); 
fi = 1;
for start = start_arr
wdw = start+1:start+Wlen;
act_map_tsr(:,:,:,fi) = squeeze(mean(manif_psth_avg_tsr(:,wdw,:,:),2)); fi=fi+1;
end
% Color limit for each channel through percentile, prctile(act_map_tsr,[2,98],[2,3,4]);
CMIN_arr = prctile(act_map_tsr,2.5,[2,3,4]);
CMAX_arr = prctile(act_map_tsr, 98,[2,3,4]);

% split the channels into group a and group b to plot
[chan_idxA, chan_idxB] = unit_id2_chan_idx(1:64, Stats(Expi).units.spikeID, Stats(Expi).units.activ_msk);

% Set up saving destination, title, label for group A
Exp_label = sprintf("%s Exp %d pref chan %d PC23 space", Animal, Expi, Stats(Expi).units.pref_chan);
if false
ITgif = fullfile(savepath,compose('%s_Exp%d_manif_IT_calign_A.gif',Animal,Expi));
V4gif = fullfile(savepath,compose('%s_Exp%d_manif_V4_calign_A.gif',Animal,Expi));
V1gif = fullfile(savepath,compose('%s_Exp%d_manif_V1_calign_A.gif',Animal,Expi));
% Set up title and color axis.
tic
for arr_chan = 1:64 % array channel! not number in the resulting array
    ch_j = chan_idxA(arr_chan);
    if isnan(ch_j), continue; end
    title_str = Stats(Expi).units.unit_name_arr(ch_j);
    set(get(ax_arr{arr_chan},"Title"),"string",title_str) % set axes title for this recording
    caxis(ax_arr{arr_chan}, [CMIN_arr(ch_j),CMAX_arr(ch_j)]) % set axes color limit for this recording
end
toc % 3 second to set it up
% Put the color data in for each frame and record them in gif
for fi = 1:length(start_arr)
wdw = start_arr(fi)+1:start_arr(fi)+Wlen;
% Add the time to the title
set(get(tIT,"title"),"string",compose(Exp_label+" IT array\nWindow: [%d,%d] ms",wdw(1),wdw(end)))
set(get(tV1,"title"),"string",compose(Exp_label+" V1V2 array\nWindow: [%d,%d] ms",wdw(1),wdw(end)))
set(get(tV4,"title"),"string",compose(Exp_label+" V4 array\nWindow: [%d,%d] ms",wdw(1),wdw(end)))
for arr_chan = 1:64 % update Color data for each channel.
    ch_j = chan_idxA(arr_chan);
    if isnan(ch_j), continue; end
    imgsc_list{arr_chan}.CData = squeeze(act_map_tsr(ch_j,:,:,fi));
    %set(get(ax_arr{arr_chan},"Children"),"CData",squeeze(act_map_tsr(ch_j,:,:,fi)))
end
drawnow
write2gif(figIT,ITgif,fi); 
write2gif(figV4,V4gif,fi); 
write2gif(figV1,V1gif,fi); 
end
toc
end
if Expi <= 39
% Set up destination for group B
ITgif = fullfile(savepath,compose('%s_Exp%d_manif_IT_calign_B.gif',Animal,Expi));
V4gif = fullfile(savepath,compose('%s_Exp%d_manif_V4_calign_B.gif',Animal,Expi));
V1gif = fullfile(savepath,compose('%s_Exp%d_manif_V1_calign_B.gif',Animal,Expi));
% Set up title and color axis.
tic
for arr_chan = 1:64 % array channel! not number in the resulting array
    ch_j = chan_idxB(arr_chan);
    if isnan(ch_j), continue; end
    title_str = Stats(Expi).units.unit_name_arr(ch_j);
    set(get(ax_arr{arr_chan},"Title"),"string",title_str) % set axes title for this recording
    caxis(ax_arr{arr_chan}, [CMIN_arr(ch_j),CMAX_arr(ch_j)]) % set axes color limit for this recording
end
toc % 3 second to set it up
% Put the color data in for each frame and record them in gif
for fi = 1:length(start_arr)
wdw = start_arr(fi)+1:start_arr(fi)+Wlen;
% Add the time to the title
set(get(tIT,"title"),"string",compose(Exp_label+" IT array\nWindow: [%d,%d] ms",wdw(1),wdw(end)))
set(get(tV1,"title"),"string",compose(Exp_label+" V1V2 array\nWindow: [%d,%d] ms",wdw(1),wdw(end)))
set(get(tV4,"title"),"string",compose(Exp_label+" V4 array\nWindow: [%d,%d] ms",wdw(1),wdw(end)))
for arr_chan = 1:64 % update Color data for each channel.
    ch_j = chan_idxB(arr_chan);
    if isnan(ch_j), continue; end
    imgsc_list{arr_chan}.CData = squeeze(act_map_tsr(ch_j,:,:,fi));
    %set(get(ax_arr{arr_chan},"Children"),"CData",squeeze(act_map_tsr(ch_j,:,:,fi)))
end
drawnow
write2gif(figIT,ITgif,fi); 
write2gif(figV4,V4gif,fi); 
write2gif(figV1,V1gif,fi);
end
toc
end
end


function write2gif(h,gifname,fi,Delay)
% do a frame shot on a figure and append it to the gif collection.
if nargin == 3
    Delay = 0.05;
end
frame = getframe(h);
% frame = getframe(); % This version without border! 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
if fi == 1 
  imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',Delay); 
else 
  imwrite(imind,cm,gifname,'gif', 'WriteMode','append','DelayTime',Delay); 
end 
end
function imgsc = plot_heatmap(ax, score_mat_avg, CLIM, title_str)
    % Given a averaged score matrix, plot a heatmap to certain axis! 
    set(0,"CurrentFigure",ancestor(ax,'figure'));
    set(ancestor(ax,'figure'),"CurrentAxes",ax); % set 
    imgsc = imagesc(squeeze(score_mat_avg), CLIM); % sum(score_mat,3)./cnt_mat -90:18:90, -90:18:90, 
    colormap('parula')
    title(title_str)
    %ylabel("PC 2 degree");xlabel("PC 3 degree")
    axis off image;shading flat;colorbar
end