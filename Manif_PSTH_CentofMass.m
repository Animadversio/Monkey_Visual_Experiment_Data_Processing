%% 
global ManifDyn Stats EStats
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'))
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
load(fullfile(mat_dir, Animal+'_ManifPopDynamics.mat'),'ManifDyn')
savepath = "E:\OneDrive - Washington University in St. Louis\PSTH_anim";
mkdir(savepath)
%% Dynamics of `prototype`


% What is the weight for center of mass. 
% - The activation. 
% 

PSTHDynViewer(ManifDyn(1).psth_tsr(1,:,:,:), sleep)

% for Expi = 20:numel(Stats)
% for si = 1:length(Stats(Expi).manif.psth) % space idx
% for ui = 1:length(Stats(Expi).units.pref_chan_id) % units idx in pref chan
%%
Wlen=20; shift_step=2;
h=figure(1);set(1,'position',[680   436   552   542]);
% index for SU and Ha within each struct
Expi = 3; ui = 1;
Unit_str = Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui));
set(0,"CurrentFigure",1)
% gifname = fullfile(savepath, compose('%s_Exp%d_manif%s_calign.gif',Animal,Expi,Unit_str));
start_arr =  0:shift_step:200-Wlen;
act_map_col = zeros(11,11,length(start_arr));
fi = 1;
for start = start_arr
wdw = start+1:start+Wlen;
act_map = cellfun(@(psth) mean(psth(ui,wdw,:),[2,3]),Stats(Expi).manif.psth{1},'UniformOutput',true);
act_map_col(:,:,fi) = act_map;fi=fi+1;
end
CMAX = prctile(act_map_col(:),98);
CMIN = prctile(act_map_col(:),2.5);
fi = 1;
for start = 0:shift_step:200-Wlen
wdw = start+1:start+Wlen;
act_map = act_map_col(:, :, fi); fi=fi+1;
imagesc(-90:18:90, -90:18:90, act_map) % sum(score_mat,3)./cnt_mat
caxis([CMIN CMAX]);
axis image;ylabel("PC 2 degree");xlabel("PC 3 degree")
title(compose("%s Manif Exp %d\n Pref Chan %s\nps: [%d,%d] ms",Animal,Expi,Unit_str,wdw(1),wdw(end)))
colorbar()
drawnow
% frame = getframe(h);% frame = getframe(); % This version without border! 
% im = frame2im(frame); 
% [imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
% if start == 0 
%   imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',0.05); 
% else 
%   imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',0.05); 
% end 
end
%%
Expi = 11;
ui = Stats(Expi).units.pref_chan_id(1);
PSTHDynViewer(ManifDyn(Expi).psth_tsr(ui,:,:,:))


function DynViewer(Expi, ui, pref)
global ManifDyn Stats
PSTHDynViewer(ManifDyn(Expi).psth_tsr(ui,:,:,:), sleep)
end
function PSTHDynViewer(psth_map, sleep)
% PSTH map is a 200 by 11 by 11 array
if nargin==1, sleep=0.01; end
if ndims(psth_map)==4, psth_map=squeeze(psth_map); end
Wlen=20; shift_step=2;
% gifname = fullfile(savepath, compose('%s_Exp%d_manif%s_calign.gif',Animal,Expi,Unit_str));
start_arr =  0:shift_step:200-Wlen;
act_map_col = cell2mat(arrayfun(@(strt)mean(psth_map(strt+1:strt+Wlen,:,:), [1]),start_arr','Uni',false)); % length(start_arr), 11, 11
act_map_col = permute(act_map_col,[2,3,1]); % 11, 11, length(start_arr)
CMAX = prctile(act_map_col(:),98);
CMIN = prctile(act_map_col(:),2.5);
imsc = imagesc(-90:18:90, -90:18:90, act_map_col(:,:,1));

h = figure(1);set(1,'position',[680   436   552   542]);
set(0,"CurrentFigure",h)
caxis([CMIN CMAX]); colorbar();
axis image; ylabel("PC 2 degree");xlabel("PC 3 degree")
title(compose("ps: [%d,%d] ms",1,20))
fi = 1;
for start = 0:shift_step:200-Wlen
wdw = start+1:start+Wlen;
imsc.CData = act_map_col(:, :, fi); fi=fi+1;
title(compose("ps: [%d,%d] ms",wdw(1),wdw(end))) % ,Animal,Expi,Unit_str,
drawnow; pause(sleep);
end

end
