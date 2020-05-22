%% Clustering PSTH dynamics (the spatial temporal tensor) for each channel

%% Really compelling visualization
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'))
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
%%
savepath = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Manif_PSTH_clster";
mkdir(savepath)
%% Only those tuned channel
lab_arr = [];
tc_patterns = zeros(46,prod([200    11    11]));
for Expi = 1:length(EStats)
ui = EStats(Expi).units.pref_chan_id; % note this pref_chan id can change in Evolution and Manifold experiments, so the pref_chan_id in evolution and Manifold exp may not be the same id.
unitlab = string(Expi)+" "+EStats(Expi).units.unit_name_arr(ui); % note this pref_chan id can change in Evolution and Manifold experiments
tc = ManifDyn(Expi).psth_tsr(ui,1:end,1:end,1:end);
tc_patterns(Expi, :) = tc(:);
end
tc_patterns = zscore(tc_patterns,0,2);
%%
U = UMAP();
embed = U.fit_transform(tc_patterns);
figure()
scatter(embed(:,1),embed(:,2))
%%
%% All those channels, and do moving average to make it better looking. 
lab_arr = [];
tc_patterns = [];%zeros(0,prod([200    11    11]));
for Expi = 1:length(EStats)
% ui = EStats(Expi).units.pref_chan_id;
unitlab = compose("E%02d ",Expi)+Stats(Expi).units.unit_name_arr;
lab_arr = [lab_arr; unitlab];
tc = movmean(ManifDyn(Expi).psth_tsr,20,2,'Endpoints','discard'); % smoothed 
tc = tc(:,1:2:end,:,:);
tc_patterns = [tc_patterns ; reshape(tc,size(tc,1),[])];
end
tc_patterns = zscore(tc_patterns,0,2); % Do zscore to normalize. 
%%
[~, umap_corr]=run_umap(tc_patterns,'n_components',2,'metric','correlation');
[~, umap_eucl]=run_umap(tc_patterns,'n_components',2,'metric','euclidean');
%% We need a good coloring to understand this distribution. So compute different labeling here.
Exp_arr = []; ui_arr = []; is_driver_arr = []; uinc_arr = []; ch_arr = []; driverch_arr = [];
for Expi = 1:length(EStats)
% ui = EStats(Expi).units.pref_chan_id;
unitlab = compose("E%02d ",Expi)+Stats(Expi).units.unit_name_arr;
is_driver = zeros(length(unitlab), 1); is_driver(Stats(Expi).units.pref_chan_id) = 1;
is_driver_arr = [is_driver_arr; is_driver];
lab_arr = [lab_arr; unitlab];
Exp_arr = [Exp_arr; Expi * ones(length(unitlab),1)];
ui_arr = [ui_arr; [1:length(unitlab)]'];
uinc_arr = [uinc_arr; Stats(Expi).units.unit_num_arr];
ch_arr = [ch_arr; Stats(Expi).units.spikeID];
driverch_arr = [driverch_arr; Stats(Expi).units.pref_chan * ones(length(unitlab),1)];
end
area_lab = zeros(length(ch_arr),1);
area_lab(ch_arr <= 48 & ch_arr >= 33) = -1; % V1 label is -1
area_lab(ch_arr <= 32) = 1; % IT label is 1
%%
zmedian = median(tc_patterns,2);
%%
embed_coord = umap.embedding;
seed = [22.2, 4.127];
window = [0.8, 0.3];
cluster1 = find(extract_idx_from_wdw(embed_coord, [22.2, 4.127], [0.8, 0.3]));
seed = [9.214, 0.2804];
window = [0.8, 0.5];
cluster2 = find(extract_idx_from_wdw(embed_coord, seed, window));
seed = [5.75, -1.709];
window = [0.8, 0.5];
cluster3 = find(extract_idx_from_wdw(embed_coord, seed, window));
%%
embed_coord = umap_eucl.embedding;
% C = find(extract_idx_from_wdw(embed_coord, [5.524, -11.26], [0.3,0.5]));
% C = find(extract_idx_from_wdw(embed_coord, [4.549, -9.895], [0.3,0.5]));
% C = find(extract_idx_from_wdw(embed_coord, [2.331, -11.42], [0.3,0.5]));
% C = find(extract_idx_from_wdw(embed_coord, [8.87, -14.69], [0.3,0.5])); % no pattern at all? 
% C = find(extract_idx_from_wdw(embed_coord, [8.629, -15.15], [0.3,0.5])); % no pattern at all? 
C = find(extract_idx_from_wdw(embed_coord, [20.37, 9.129], [0.1,0.1])); % no pattern at all? 
C = find(extract_idx_from_wdw(embed_coord, [23.29, 7.994], [0.1,0.1])); % no pattern at all? 
rowi = randsample(C,1); %C(1);
% rowi = randsample(size(embed_coord,1),1);
play_gif_for_ch(ManifDyn,Exp_arr(rowi),ui_arr(rowi))
% U = UMAP();
% embed_all = U.fit_transform(tc_patterns);
% figure()
% scatter(embed_all(:,1),embed_all(:,2))
%%
%%
% Exp_arr = []; ui_arr = []; uinc_arr = []; ch_arr = []; driverch_arr = [];
embed_coord = umap_eucl.embedding;
figure(11)
sct = scatter(embed_coord(:,1), embed_coord(:,2), 32, Exp_arr,'filled','MarkerFaceAlpha',0.3);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("Chan",ch_arr);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("unit in ch",uinc_arr);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("Expi",Exp_arr);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("driver ch",driverch_arr);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("zmed",zmedian);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("rowi",1:length(ch_arr));
colorbar()
embed_coord = umap_corr.embedding;
figure(17);
sct2 = scatter(embed_coord(:,1), embed_coord(:,2), 32, Exp_arr,'filled','MarkerFaceAlpha',0.3);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("Chan",ch_arr);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("unit in ch",uinc_arr);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("Expi",Exp_arr);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("driver ch",driverch_arr);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("zmed",zmedian);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("rowi",1:length(ch_arr));
colorbar()
%% Color the cluster based on area of recording.
figure(18)
sct = scatter(umap_eucl.embedding(:,1), umap_eucl.embedding(:,2), ...
    25 * (1+is_driver_arr).^3, area_lab,'filled','MarkerFaceAlpha',0.5);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("Chan",ch_arr);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("unit in ch",uinc_arr);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("Expi",Exp_arr);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("driver ch",driverch_arr);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("zmed",zmedian);
sct.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("rowi",1:length(ch_arr));
title("V1: Blue, -1; V4: Cyan 0; IT: Yellow 1")
colorbar()
caxis([-1,1.5])

figure(19);
sct2 = scatter(umap_corr.embedding(:,1), umap_corr.embedding(:,2), ...
    25 * (1+is_driver_arr).^3, area_lab,'filled','MarkerFaceAlpha',0.5);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("Chan",ch_arr);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("unit in ch",uinc_arr);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("Expi",Exp_arr);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("driver ch",driverch_arr);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("zmed",zmedian);
sct2.DataTipTemplate.DataTipRows(end + 1) = dataTipTextRow("rowi",1:length(ch_arr));
title("V1: Blue, -1; V4: Cyan 0; IT: Yellow 1")
colorbar()
caxis([-1,1.5])
%%
figure()
plot(median(tc_patterns,2))
%%
umap_eucl.embedding(rowi,:)
umap_corr.embedding(rowi,:)
%%
% rowi = randsample(size(embed_coord,1),1);
rowi = 3587;%randsample(C,1); %C(1);
play_gif_for_ch(ManifDyn,Animal,Exp_arr(rowi),ui_arr(rowi),lab_arr(rowi), false);
%%
function h=play_gif_for_ch(ManifDyn,Animal,Expi,ch_i,title_str,savegif)
if nargin == 5
    savegif = false;
end
WdwL = 20; stride = 2;
tc = movmean(ManifDyn(Expi).psth_tsr(ch_i,:,:,:), WdwL, 2,'Endpoints','discard'); % smoothed 
tc = tc(:,1:stride:end,:,:);
CLIM = prctile(tc,[2.5,98],'all'); 
h = figure();
imgsc = imagesc(-90:18:90, -90:18:90, squeeze(tc(1,1,:,:))); % sum(score_mat,3)./cnt_mat
axis image;ylabel("PC 2 degree");xlabel("PC 3 degree")
caxis(CLIM); colorbar()
title(compose(title_str+"\n%s Manif Exp %d Unit %d\n [%d,%d] ms",Animal,Expi,ch_i, 1,1+WdwL))%Animal,

if savegif
savepath = "C:\Users\binxu\OneDrive - Washington University in St. Louis\PSTH_anim\OtherChanAnim";
gifname = fullfile(savepath, compose("%s_manif_%s.gif",Animal,title_str));%Animal
end
for fi = 1:size(tc,2)
    imgsc.CData = squeeze(tc(1,fi,:,:));
    set(get(gca,"Title"),"String",compose(title_str+"\n%s Manif Exp %d Unit %d\n [%d,%d] ms",...
        Animal, Expi,ch_i,stride*(fi-1)+1,stride*(fi-1)+1+WdwL))
    drawnow
    pause(0.02)
    if savegif
    frame = getframe(h);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    if fi == 1 
      imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',0.05); 
    else 
      imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',0.05); 
    end 
    end
end
end

function idx = extract_idx_from_wdw(coord, seed, window)
idx = all((coord < seed + window) & (coord > seed - window),2);
end