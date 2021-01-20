%% Manif_Animation  
%% Really compelling visualization
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'))
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
load(fullfile(mat_dir, Animal+'_ManifPopDynamics.mat'),'ManifDyn')
%%
savepath = "E:\OneDrive - Washington University in St. Louis\PSTH_anim";
mkdir(savepath)
%%
Wlen=20; shift_step=2;
% Expi = 10;
% ui=2; 
h=figure(1);set(1,'position',[680   436   552   542]);
for Expi = 1:numel(Stats)
 % index for SU and Ha within each struct
% for si = 1:length(Stats(Expi).manif.psth) % space idx
for ui = 1:length(Stats(Expi).units.pref_chan_id)
Unit_str = Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui));
set(0,"CurrentFigure",1)
% M=repmat(getframe(gca),0,0);
gifname = fullfile(savepath,compose('%s_Exp%d_manif%s.gif',Animal,Expi,Unit_str));
for start = 0:shift_step:200-Wlen
wdw = start+1:start+Wlen;
act_map = cellfun(@(psth) mean(psth(ui,wdw,:),[2,3]),Stats(Expi).manif.psth{1},'UniformOutput',true);
imagesc(-90:18:90, -90:18:90, act_map) % sum(score_mat,3)./cnt_mat
axis image
ylabel("PC 2 degree");xlabel("PC 3 degree")
title(compose("%s Manif Exp %d\n Pref Chan %s\nps: [%d,%d] ms",Animal,Expi,Unit_str,wdw(1),wdw(end)))
colorbar()
drawnow
frame = getframe(h);
% frame = getframe(); % This version without border! 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
if start == 0 
  imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',0.05); 
else 
  imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',0.05); 
end 
% M=[M,getframe()];
end
% end
end
end
%% Use the same caxis
Wlen=20; shift_step=2;
%Expi = 10;
%ui=2; 
h=figure(1);set(1,'position',[680   436   552   542]);
for Expi = 20:numel(Stats)
 % index for SU and Ha within each struct
% for si = 1:length(Stats(Expi).manif.psth) % space idx
for ui = 1:length(Stats(Expi).units.pref_chan_id)
Unit_str = Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui));
set(0,"CurrentFigure",1)
% M=repmat(getframe(gca),0,0);
gifname = fullfile(savepath,compose('%s_Exp%d_manif%s_calign.gif',Animal,Expi,Unit_str));
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
axis image
ylabel("PC 2 degree");xlabel("PC 3 degree")
title(compose("%s Manif Exp %d\n Pref Chan %s\nps: [%d,%d] ms",Animal,Expi,Unit_str,wdw(1),wdw(end)))
colorbar()
drawnow
frame = getframe(h);
% frame = getframe(); % This version without border! 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
if start == 0 
  imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',0.05); 
else 
  imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',0.05); 
end 
% M=[M,getframe()];
end
% end
end

end
%%
for Expi = 1:numel(Stats)
 % index for SU and Ha within each struct
% for si = 1:length(Stats(Expi).manif.psth) % space idx
for ui = 1:length(Stats(Expi).units.pref_chan_id)
Unit_str = Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui));
% M=repmat(getframe(gca),0,0);
movefile(fullfile(savepath,compose('Exp%d_manif%s_calign.gif',Expi,Unit_str)),...
        fullfile(savepath,compose('%s_Exp%d_manif%s_calign.gif',Animal,Expi,Unit_str)));
movefile(fullfile(savepath,compose('Exp%d_manif%s.gif',Expi,Unit_str)),...
        fullfile(savepath,compose('%s_Exp%d_manif%s.gif',Animal,Expi,Unit_str)));
end
end