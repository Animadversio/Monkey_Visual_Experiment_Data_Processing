%% Manif_Population Dynamics in another fashion. (show the array together)
system("subst S: E:\Network_Data_Sync")
addpath .\utils
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir,Animal+"_ManifPopDynamics.mat"));
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"));
load(fullfile(mat_dir,Animal+"_Evol_stats.mat"));
%%
Extract_Channel_Organization;
%%
Expi = 36;
[chan_idxA, chan_idxB] = unit_id2_chan_idx(1:32, Stats(Expi).units.spikeID, Stats(Expi).units.activ_msk);
% plot the channel layout 
h = figure;clf;hold on
scatter(IT_chan_XY(:,1), IT_chan_XY(:,2),64,ManifDyn(Expi).psth_tsr(chan_idxA,1,1,1),'filled')
axis equal
scatter(V4_chan_XY(:,1), 1600 + V4_chan_XY(:,2),64,ManifDyn(Expi).psth_tsr(chan_idxA(1:16),1,1,1),'filled')
scatter(2200 + V1_chan_XY(:,1), 1600 + V1_chan_XY(:,2),64,ManifDyn(Expi).psth_tsr(chan_idxA(17:32),1,1,1),'filled')
set(gca,'color',[0.5,0.5,0.5])
%%
Ysft = 1700; V1Xsft = 2200;
chan64_XY = [[IT_chan_XY(:,1);       V4_chan_XY(:,1);V1Xsft + V1_chan_XY(:,1)], ...
             [IT_chan_XY(:,2);Ysft + V4_chan_XY(:,2);  Ysft + V1_chan_XY(:,2)]];
%%
[chan_idxA, chan_idxB] = unit_id2_chan_idx(1:64, Stats(Expi).units.spikeID, Stats(Expi).units.activ_msk);
h = figure(11);clf;hold on
text(chan64_XY(:,1)-50, chan64_XY(:,2)+100, arrayfun(@string,1:64))% 
sct = scatter(chan64_XY(:,1), chan64_XY(:,2),81,ManifDyn(Expi).psth_tsr(chan_idxA,1,1,1),'filled');
axis equal off
set(gca,'color',[0.,0.,0.])
set(gca,'position',[0.05,0.05,0.9,0.9])
for fi = 1:200
sct.CData = ManifDyn(Expi).psth_tsr(chan_idxA,fi,1,1);
drawnow
end
%%
h = figure(12);clf;hold on
text(chan64_XY(:,1)-50, chan64_XY(:,2)+100, arrayfun(@string,1:64))% 
sct = scatter(chan_all_XY(:,1), chan_all_XY(:,2),81,ManifDyn(Expi).psth_tsr(:,1,1,1),'filled');
axis equal off
set(gca,'color',[0.,0.,0.])
set(gca,'position',[0.05,0.05,0.9,0.9])
for fi = 1:200
sct.CData = ManifDyn(Expi).psth_tsr(:,fi,1,1);
drawnow
end
%% Calculate all XY coordinates  
chan_all_XY = calcXY_all(chan64_XY, Stats(Expi).units.spikeID, Stats(Expi).units.unit_num_arr);
%% the normalization scheme should be really careful. Normalize across channel is not informative.
Wlen=20; shift_step=2; % setting for averaging the frames
start_arr = 0:shift_step:200-Wlen; % array of the start time for the moving window.  
act_map_tsr = zeros(length(Stats(Expi).units.spikeID),length(start_arr),11,11); 
fi = 1;
for start = start_arr
wdw = start+1:start+Wlen;
act_map_tsr(:,fi,:,:) = squeeze(mean(ManifDyn(Expi).psth_tsr(:,wdw,:,:),2)); fi=fi+1;
end
%%
zscore_tsr = zscore(act_map_tsr,0,[2,3,4]);%ManifDyn(Expi).psth_tsr
CMIN = prctile(zscore_tsr,2.5,'all');
CMAX = prctile(zscore_tsr,97.5,'all');
%%
h = figure(12);clf;hold on
text(chan64_XY(:,1)-50, chan64_XY(:,2)+100, arrayfun(@string,1:64))% 
sct = scatter(chan_all_XY(:,1), chan_all_XY(:,2),81,zscore_tsr(:,1,1,1),'filled',"MarkerFaceAlpha",0.7);
title(compose("[%d,%d] ms",start_arr(fi),start_arr(fi)+Wlen))
caxis([CMIN, CMAX]);axis equal off
set(gca,'color',[0.,0.,0.])
set(gca,'position',[0.05,0.05,0.9,0.9])
for fi = 1:2:size(zscore_tsr,2)
sct.CData = zscore_tsr(:,fi,1,1);
set(get(gca,"title"),"string",compose("[%d,%d] ms",start_arr(fi),start_arr(fi)+Wlen))
drawnow
pause(0.06)
%title(compose("[%d] ms",fi))
% title(compose("%s Manif Exp %d\n Pref Chan %s\nps: [%d,%d] ms",Animal,Expi,Unit_str,fi))
end
%%
figure(3)
sct_arr = repmat(matlab.graphics.chart.primitive.Scatter,0,0);
ax_arr = repmat(matlab.graphics.axis.Axes,0,0);
for i=1:11
    for j=1:11
    ax_arr(i,j) = subplottight(11,11,j+(i-1)*11,0.05,0.08);
    sct_arr(i,j) = scatter(chan_all_XY(:,1),chan_all_XY(:,2),49,zscore_tsr(:,1,i,j),'filled',"MarkerFaceAlpha",0.7);
    caxis([CMIN, CMAX]);axis equal off
    set(ax_arr(i,j),'color',[0.,0.,0.])
%     set(ax,'position',[0.05,0.05,0.9,0.9])
    end
end
%%
pause(0.02)
for fi = 1:1:size(zscore_tsr,2)
for i=1:11
    for j=1:11
    sct_arr(i,j).CData = zscore_tsr(:,fi,i,j);    
    end
end
drawnow
end
%%
ch_msk = 48 < Stats(Expi).units.spikeID & Stats(Expi).units.spikeID <= 64;
figure(5);clf;set(5,'position',[426          67        1108         914])
sct_arr = cell(11);%repmat(struct(),0,0);
ax_arr = cell(11);%repmat(struct(),0,0);
for i=1:11
    for j=1:11
    ax_arr{i,j} = subplottight(11,11,j+(i-1)*11,0.09,0.08,0.0);
    sct_arr{i,j} = scatter(chan_all_XY(ch_msk,1),chan_all_XY(ch_msk,2),49,zscore_tsr(ch_msk,1,i,j),'filled',"MarkerFaceAlpha",0.7);
    caxis([CMIN, CMAX]);axis equal off
    set(ax_arr{i,j},'color',[0.,0.,0.])
%     set(ax,'position',[0.05,0.05,0.9,0.9])
    end
end
% T = suptitle(compose("[%d,%d] ms",start_arr(fi),start_arr(fi)+Wlen));
%%
drawnow
pause(0.02)
for fi = 1:1:size(zscore_tsr,2)
for i=1:11
    for j=1:11
    sct_arr{i,j}.CData = zscore_tsr(ch_msk,fi,i,j);    
    end
end
drawnow
end
%%
savepath = "C:\Users\binxu\OneDrive - Washington University in St. Louis\PopDyn_PSTH_Anim";
mkdir(savepath)
%%
Unit_str = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
% Unit_str = Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui));
ch_msk = 48 < Stats(Expi).units.spikeID & Stats(Expi).units.spikeID <= 64;
figure(8);clf;set(8, 'position',[544          71        1035         907])
tile = tiledlayout(11,11,'TileSpacing','Compact','Padding','compact');
for i=1:11
    for j=1:11
    ax_arr{i,j} = nexttile(j+(i-1)*11);
    sct_arr{i,j} = scatter(chan_all_XY(ch_msk,1),chan_all_XY(ch_msk,2),49,zscore_tsr(ch_msk,1,i,j),'filled',"MarkerFaceAlpha",0.7);
    caxis([CMIN, CMAX]);axis tight off
%     set(ax_arr{i,j},'color',[0.,0.,0.])
    end
end
title(compose("%s Manif Exp %d\n Pref Chan %s ps: [%d,%d] ms",Animal,Expi,Unit_str,start_arr(fi),start_arr(fi)+Wlen))
% title(tile, compose("[%d,%d] ms",start_arr(fi),start_arr(fi)+Wlen))
%% IT image
set(8,'position',[206          64        1572         927]) % setting for IT image
for Expi = 1:46
tic
Unit_str = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
ITgif = fullfile(savepath,compose('%s_Exp%d_manif_IT_z.gif',Animal,Expi));
V4gif = fullfile(savepath,compose('%s_Exp%d_manif_V4_z.gif',Animal,Expi));
V1gif = fullfile(savepath,compose('%s_Exp%d_manif_V1_z.gif',Animal,Expi));
chan_all_XY = calcXY_all(chan64_XY, Stats(Expi).units.spikeID, Stats(Expi).units.unit_num_arr);

Wlen=20; shift_step=2; % setting for averaging the frames
start_arr = 0:shift_step:200-Wlen; % array of the start time for the moving window.  
act_map_tsr = zeros(length(Stats(Expi).units.spikeID),length(start_arr),11,11); 
fi = 1;
for start = start_arr
wdw = start+1:start+Wlen;
act_map_tsr(:,fi,:,:) = squeeze(mean(ManifDyn(Expi).psth_tsr(:,wdw,:,:),2)); fi=fi+1;
end
zscore_tsr = zscore(act_map_tsr,0,[2,3,4]);%ManifDyn(Expi).psth_tsr
CMIN = prctile(zscore_tsr,2.5,'all');
CMAX = prctile(zscore_tsr,97.5,'all');
toc
IT_msk = 0 < Stats(Expi).units.spikeID & Stats(Expi).units.spikeID <= 32;
for i=1:11
    for j=1:11
    sct_arr{i,j}.XData = chan_all_XY(IT_msk,1);
    sct_arr{i,j}.YData = chan_all_XY(IT_msk,2);
    sct_arr{i,j}.CData = zscore_tsr(IT_msk,1,i,j);
    % axis(ax_arr{i,j},'equal','tight','off')% this line really takes time.
    end
end
toc
% approx 1 min to reset all axes. 
tic
for fi = 1:1:size(zscore_tsr,2)
    for i=1:11
        for j=1:11
        sct_arr{i,j}.CData = zscore_tsr(IT_msk,fi,i,j);    
        end
    end
    % set(tile.Title,"string",compose("[%d,%d] ms",start_arr(fi),start_arr(fi)+Wlen))
    set(tile.Title,"string",compose("%s Manif Exp %d IT\n Pref Chan %s ps: [%d,%d] ms",...
        Animal,Expi,Unit_str,start_arr(fi),start_arr(fi)+Wlen))
    drawnow
    write2gif(8,ITgif,fi,0.10); 
    toc
end
end

%% V1, V4 Image
set(8,'position',[529          27        1249         964]) % setting for V4,V1
for Expi = 1:46
fprintf("Processing Exp %d\n",Expi)
tic
Unit_str = EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id);
ITgif = fullfile(savepath,compose('%s_Exp%d_manif_IT_z.gif',Animal,Expi));
V4gif = fullfile(savepath,compose('%s_Exp%d_manif_V4_z.gif',Animal,Expi));
V1gif = fullfile(savepath,compose('%s_Exp%d_manif_V1_z.gif',Animal,Expi));
chan_all_XY = calcXY_all(chan64_XY, Stats(Expi).units.spikeID, Stats(Expi).units.unit_num_arr); % plot location for all the units
% Compute moving average and zscore activation in each unit 
Wlen=20; shift_step=2; % setting for averaging the frames
start_arr = 0:shift_step:200-Wlen; % array of the start time for the moving window.  
act_map_tsr = zeros(length(Stats(Expi).units.spikeID),length(start_arr),11,11); 
fi = 1;
for start = start_arr
wdw = start+1:start+Wlen;
act_map_tsr(:,fi,:,:) = squeeze(mean(ManifDyn(Expi).psth_tsr(:,wdw,:,:),2)); fi=fi+1;
end
zscore_tsr = zscore(act_map_tsr,0,[2,3,4]);%ManifDyn(Expi).psth_tsr
CMIN = prctile(zscore_tsr,2.5,'all');% compute CMin and CMax for color normalization
CMAX = prctile(zscore_tsr,97.5,'all');
toc % Finish computing

V4_msk = 48 < Stats(Expi).units.spikeID & Stats(Expi).units.spikeID <= 64;
for i=1:11
    for j=1:11
    sct_arr{i,j}.XData = chan_all_XY(V4_msk,1);
    sct_arr{i,j}.YData = chan_all_XY(V4_msk,2);
    sct_arr{i,j}.CData = zscore_tsr(V4_msk,1,i,j);
    %axis(ax_arr{i,j},'equal','tight','off') % this line really takes
    %time. and is useless. Change the data and the axis limit will change
    %automatically
    end
end
toc % finish setting up stage
% approx 1 min to reset all axes. about 0.22 sec to set the data instead
tic
for fi = 1:1:size(zscore_tsr,2)
    for i=1:11
        for j=1:11
        sct_arr{i,j}.CData = zscore_tsr(V4_msk,fi,i,j);    
        end
    end
    % set(tile.Title,"string",compose("[%d,%d] ms",start_arr(fi),start_arr(fi)+Wlen))
    set(tile.Title,"string",compose("%s Manif Exp %d V4\n Pref Chan %s ps: [%d,%d] ms",...
        Animal,Expi,Unit_str,start_arr(fi),start_arr(fi)+Wlen))
    drawnow
    write2gif(8,V4gif,fi,0.10); % each frame takes around 1.2 sec
end
toc % ploting the map of V4 array 

tic
V1_msk = 32 < Stats(Expi).units.spikeID & Stats(Expi).units.spikeID <= 48;
for i=1:11
    for j=1:11
    sct_arr{i,j}.XData = chan_all_XY(V1_msk,1);
    sct_arr{i,j}.YData = chan_all_XY(V1_msk,2);
    sct_arr{i,j}.CData = zscore_tsr(V1_msk,1,i,j);
    end
end
toc % finish setting up stage
% approx 1 min to reset all axes. about 0.22 sec to set the data instead
tic
for fi = 1:1:size(zscore_tsr,2)
    for i=1:11
        for j=1:11
        sct_arr{i,j}.CData = zscore_tsr(V1_msk,fi,i,j);    
        end
    end
    % set(tile.Title,"string",compose("[%d,%d] ms",start_arr(fi),start_arr(fi)+Wlen))
    set(tile.Title,"string",compose("%s Manif Exp %d V1\n Pref Chan %s ps: [%d,%d] ms",...
        Animal,Expi,Unit_str,start_arr(fi),start_arr(fi)+Wlen))
    drawnow
    write2gif(8,V1gif,fi,0.10);  
end
toc % ploting the map of V1 array 
end
%%
function chan_all_XY = calcXY_all(chan64_XY, spikeID, unit_num_arr)
intv = 400;
offset = [[0,0];[120, 0];[0, -120];[120, -120];[240, -120];];
chan_all_XY = nan(length(spikeID),2);
for ch_i = 1:length(spikeID)
unit = unit_num_arr(ch_i);
if unit > 0
    chan_all_XY(ch_i,:) = chan64_XY(spikeID(ch_i),:) + offset(unit, :);
end
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