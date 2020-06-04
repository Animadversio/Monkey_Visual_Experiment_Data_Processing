% Loading and visualizing the Manifold experiment (image change and psth change w.r.t image)
function Manif_Stats_Anim_fun(Stats, EStats, Animal, Expi, option)
global G
if nargin==4
    option = struct("save",true);
    option.result_dir = "E:\OneDrive - Washington University in St. Louis\Evol_Manif_Movies";
    option.result_dir = "S:\Evol_Manif_Movies";
    option.cent = []; % can be used to specify the centor of exploration. 
end
mvWdL = 10;
fprintf("Processing %s Manifold Exp %d\n",Animal,Expi)
Window = 50:200;
ui=1;si=1;
imgnm_grid = cellfun(@(idx)string(unique(Stats(Expi).imageName(idx))),Stats(Expi).manif.idx_grid{si});
imgfullfn = arrayfun(@(fn) string(ls(fullfile(Stats(Expi).meta.stimuli, fn+"*"))),imgnm_grid);
imgpath_grid = arrayfun(@(fn) fullfile(Stats(Expi).meta.stimuli, fn),imgfullfn);

manif_psth_avg = cellfun(@(psth)mean(movmean(psth(ui,:,:),mvWdL,2),3),Stats(Expi).manif.psth{si},'UniformOutput',false);
manif_psth_avg = cell2mat(reshape(manif_psth_avg,[1,1,11,11]));
manif_psth_sem = cellfun(@(psth)std(movmean(psth(ui,:,:),mvWdL,2),0,3)/sqrt(size(psth,3)),Stats(Expi).manif.psth{si},'UniformOutput',false);
manif_psth_sem = cell2mat(reshape(manif_psth_sem,[1,1,11,11]));
%% Load the source basis data for Manifold Experiments from python
% basis_path = fullfile(Stats(Expi).meta.stimuli,"PC_vector_data.npz");
basis_path = fullfile(EStats(Expi).meta.stimuli,"PC_imgs","PC_vector_data.npz");
f = py.numpy.load(basis_path);
PC_Vec = f.get('PC_vecs').double;
sphere_norm = f.get('sphere_norm').double;
f.close();
% Get the mat containing all the codes of the last generation. To see
% whether we should inverse PC1 
matfns = string(ls(fullfile(EStats(Expi).meta.stimuli,"*.mat")));
code_tmp = load(fullfile(EStats(Expi).meta.stimuli,matfns(end)));
proj_coord = mean(code_tmp.codes,1) * PC_Vec';
if proj_coord(1)>0
    basis = PC_Vec(1:3,:);
else
    fprintf("The evolution direction is inverse to the PC1 direction of PCA. Inverse PC1 as basis\n")
    basis = [-1,1,1]' .* PC_Vec(1:3,:);% Note the final PC may need to reverse! not always the same dir!
end
clear code_tmp matfns
%%
% set up trajectory for the focal point
scoremap = squeeze(mean(manif_psth_avg(:,Window,:,:),[1,2]));
[~,peakidx]=max(scoremap,[],'all','linear');
[peaki,peakj] = ind2sub(size(scoremap),peakidx);
cent = ([peaki,peakj] - 6) * pi/10;
vect = [0, 1];
coord_col{1} = cent + ([0:0.1:4,4:-0.1:-4,-4:0.1:0])'*vect;
vect = [1, 0];
coord_col{2} = cent + ([0:0.1:4,4:-0.1:-4,-4:0.1:0])'*vect;
vect = [1, 1];
coord_col{3} = cent + ([0:0.1:4,4:-0.1:-4,-4:0.1:0])'*vect;
vect = [1, -1];
coord_col{4} = cent + ([0:0.1:4,4:-0.1:-4,-4:0.1:0])'*vect;

result_dir = option.result_dir; % "E:\OneDrive - Washington University in St. Louis\Evol_Manif_Movies";
savepath = result_dir; 
if option.save 
v = VideoWriter(fullfile(savepath,compose('%s_Manif_Exp%02d_Avg_PSTH.mov',Animal,Expi)));
v.FrameRate = 4;
open(v);
end

h3=figure(3);set(3,'Position',[505   235   760   742]);clf;
ax1 = subplot('position',[0.1300    0.5601    0.3347    0.3234]);%2,2,1);
imshow(nan(256,256,3))
ax3 = subplot('position',[0.5474    0.5647    0.3316    0.3208]);%2,2,2-=[
imagesc(linspace(-pi/2,pi/2,11),linspace(-pi/2,pi/2,11),scoremap);axis image;colorbar;hold on
ylabel("PC2");xlabel("PC3")
focalpoint = plot(cent(2),cent(1),"Color",'r','Marker','o');
ax2 = subplot('position',[0.0869    0.0809    0.8552    0.3787]);cla%212
hold on % background data % bgpsth = plot(psth2plot','Color',[0.7,0.7,0.7]);
sEB = shadedErrorBar([],squeeze(manif_psth_avg(1,:,peaki,peakj)),manif_psth_sem(1,:,peaki,peakj),...
        'lineprops',{'Color','r','LineWidth',1.5},'transparent',0,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
ylabel("PSTH (Hz)");xlabel("time(ms)")
ylim([0,max(manif_psth_avg+manif_psth_sem,[],'all')])
title("Evoked PSTH")
ST = suptitle(compose("%s Manif Exp %02d pref chan %s", Animal, Expi, ...
    Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui))));
for k = 1:numel(coord_col)
coord = coord_col{k};
validmask = coord(:,1)<=pi/2 & coord(:,1)>=-pi/2 & coord(:,2)<=pi/2 & coord(:,2)>=-pi/2;
coordfin = coord(validmask,:);
% Collect the psths that we will encounter
matidxs = unique(round(clip(6 + coord / pi * 10, [1,11])),"rows");
psth2plot = cell2mat(arrayfun(@(i,j)manif_psth_avg(1,:,i,j),matidxs(:,1),matidxs(:,2),'UniformOutput', false));
set(h3,"CurrentAxes",ax2);
bgpsth = plot(psth2plot','Color',[0.7,0.7,0.7]);
for fi = 1:size(coordfin,1)
    theta=coordfin(fi,1); phi=coordfin(fi,2);
    % interpolate data 
    codecur = [cos(theta)*cos(phi), sin(theta)*cos(phi), sin(phi)] * basis * sphere_norm;
    imgcur = G.visualize(codecur);
    % actually this is not nearest neighbor in the spherical sense
    rowi = round(clip(6 + theta / pi * 10, [1,11]));
    colj = round(clip(6 + phi   / pi * 10, [1,11]));
    psthcur = manif_psth_avg(:,:,rowi,colj);
    psthsemcur = manif_psth_sem(:,:,rowi,colj);
    % visualize 
    focalpoint.YData = theta;
    focalpoint.XData = phi;
    ax3.Title.String = compose("Theta=%.1f Phi=%.1f", theta / pi * 180, phi / pi * 180);
    set(0,"CurrentFigure",h3)
    set(h3,"CurrentAxes",ax1);%subplot(211);
    imshow(imgcur); 
    set(h3,"CurrentAxes",ax2);%subplot(212);%cla(gca,'reset')%hold on
%     sEB = shadedErrorBar([],psthcur,psthsemcur,...
%         'transparent',1,'patchSaturation',0.15);%'lineprops',{'Color',[color_seq(blocki, :),0.85]},
    ax2.Title.String = compose("PSTH, Evoked Rate %.1f", mean(psthcur(Window)));
    uE = psthcur + psthsemcur; lE = psthcur - psthsemcur;
    sEB.mainLine.YData = psthcur;
    sEB.edge(1).YData = lE;
    sEB.edge(2).YData = uE;
    sEB.patch.Vertices(:,2) = [lE,fliplr(uE)];
    drawnow;
    Fcur = getframe(h3);
    if option.save, writeVideo(v,Fcur); end
end
if k~=numel(coord_col), delete(bgpsth); end
end
if option.save, close(v); end
end