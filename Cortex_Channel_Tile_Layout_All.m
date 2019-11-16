function [ax_arr,tIT,tV1,tV4] = Cortex_Channel_Tile_Layout_All(figIT, figV1, figV4)
Extract_Channel_Organization;
ax_arr = {};
if isempty(figIT), figIT = figure(); end
figure(figIT)
set(figIT,'position',[  1          41        2560         963]);%set(figIT,'Visible','off')
tIT = tiledlayout(8,20,'TileSpacing','Compact');
for i = 1:32
    init_idx = calc_tile_init_idx(IT_chan_XY(i, 1), IT_chan_XY(i, 2), 8, 20);
    ax = nexttile(init_idx, [2, 2]);
    axis image
    %title(num2str(i))
    ax_arr{end+1} = ax;
end

if isempty(figV1), figV1 = figure(); end
figure(figV1)
set(figV1,'position',[ 50          33        1220         963]);%set(figV1,'Visible','off')
tV1 = tiledlayout(8,10,'TileSpacing','Compact');
for i = 1:16
    init_idx = calc_tile_init_idx(V1_chan_XY(i, 1), V1_chan_XY(i, 2), 8, 10);
    ax = nexttile(init_idx, [2, 2]);
    axis image
    %title(num2str(i))
    ax_arr{end+1} = ax;
end

if isempty(figV4), figV4 = figure(); end
figure(figV4)
set(figV4,'position',[ 50          33        1220         963]);%set(figV4,'Visible','off')
tV4 = tiledlayout(8,10,'TileSpacing','Compact');
for i = 1:16
    init_idx = calc_tile_init_idx(V4_chan_XY(i, 1), V4_chan_XY(i, 2), 8, 10);
    ax = nexttile(init_idx, [2, 2]);
    axis image
    %title(num2str(i))
    ax_arr{end+1} = ax;
end
end
function idx = calc_tile_init_idx(X, Y, rown, coln)
    Yi = int32(Y / 400 / sin(pi/3));
    Xi = int32(X / 200) ;
    idx = coln * Yi * 2 + Xi + 1 ;
end