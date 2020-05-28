%% Evolution Experiment Visualization from a global Embedding
addpath DNN
addpath utils
system("subst S: E:\Network_Data_Sync")
system("subst N: E:\Network_Data_Sync")
% addpath("C:\Users\binxu\OneDrive - Washington University in St. Louis\Matlab_add_on\npy-matlab-master\npy-matlab")
pe = pyenv('Version','C:\Users\binxu\.conda\envs\caffe36\python.exe'); % Note the python env could not be changed in a matlab session
py.importlib.import_module('numpy');
%% Loading up Data and G
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir,"Alfa_ManifPopDynamics.mat"));
load(fullfile(mat_dir,"Alfa_Manif_stats.mat"));
load(fullfile(mat_dir,"Alfa_Evol_stats.mat"));
load(fullfile(mat_dir,"Alfa_Manif_RFMaps.mat"));
%
global G 
G = FC6Generator("matlabGANfc6.mat");
%%
% code = randn(1,4096);
% imgs = G.visualize(code);
%% Load the Basis vectors
Expi = 3;
%% Load the source basis data for Manifold Experiments from python
basis_path = fullfile(Stats(Expi).meta.stimuli,"PC_vector_data.npz");
f = py.numpy.load(basis_path);
PC_Vec = f.get('PC_vecs').double;
sphere_norm = f.get('sphere_norm').double;
f.close();
% Get the mat containing all the codes of the last generation. 
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
%% Prepare neural response data 
global psth_avg_tsr psth_std_tsr
si=1; ui=1; iCh=Stats(Expi).units.pref_chan_id(ui);
meanpsth = cellfun(@(psth) mean(psth(ui,:,:),3), Stats(Expi).manif.psth{si}, "UniformOutput", false);
stdpsth = cellfun(@(psth)  std(psth(ui,:,:),0,3)/sqrt(size(psth,3)), Stats(Expi).manif.psth{si}, "UniformOutput", false);
psth_avg_tsr = cell2mat(reshape(meanpsth,1,1,11,11)); % reshape and transform, so size 86, 200, 11, 11, key data! 
psth_std_tsr = cell2mat(reshape(stdpsth,1,1,11,11)); % reshape and transform, so size 86, 200, 11, 11, key data! 
scoremap = squeeze(mean(psth_avg_tsr(1,50:200,:,:),[1,2]));
%% RF Masks
% Assume the Manifold and Evolution Exp happens at the same location
imgsize = EStats(Expi).evol.imgsize;
imgpos = EStats(Expi).evol.imgpos;
x_ext = imgpos(1) + [- imgsize / 2, imgsize / 2];
y_ext = imgpos(2) + [- imgsize / 2, imgsize / 2];
x_grid = linspace(x_ext(1), x_ext(2), 256);
y_grid = linspace(y_ext(1), y_ext(2), 256);
[XX, YY] = meshgrid(x_grid,y_grid);
% universal map grid. 
ntick = 201;
visualField = [-10 10]; 
coli = linspace(visualField(1),visualField(2),ntick);
rowi = linspace(visualField(1),visualField(2),ntick);
[mapgridX,mapgridY]  = meshgrid(coli,rowi); 
% Search for the closest exp % Find the RF map exp on the same day
Mapi = 3;
% Get the corresponding iCh
target_ui = Stats(Expi).units.unit_num_arr(iCh);
target_chi = Stats(Expi).units.spikeID(iCh);
RF_iCh = find((MaskStats(Mapi).unit.unit_num_arr == target_ui) & (MaskStats(Mapi).unit.chan_num_arr == target_chi));
% Get the mask
interpmask = MaskStats(Mapi).interpmasks(:,:,RF_iCh);
convmask = MaskStats(Mapi).convmasks(:,:,RF_iCh);
mask_conv = griddata(mapgridX(:),mapgridY(:),double(convmask(:)),XX, YY); 
mask_interp = griddata(mapgridX(:),mapgridY(:),double(interpmask(:)),XX, YY); 

mask = mask_interp > max(mask_interp,[],'all') * 0.65;
% %% Dummy RF mask
% [XX, YY] = meshgrid(1:256,1:256);
% D = sqrt((XX-127).^2 + (YY-127).^2);
% mask = exp(-(D-40).^2/100);
% mask = max((mask ./ max(mask,[],'all')), D<40);
%% True visualizatin GUI 
figh = figure(3);clf;figh.Position = [28         524        1446         454];
data = struct('basis', basis, 'norm', sphere_norm, 'mask', mask, ...
            'Addmask', 0, "Theta", 0, "Phi", 0, "curcode", basis(1,:), "scoremap", scoremap);
subplot('Position',[0.0913    0.1100    0.2144    0.7446]); hold on
data.tunemap = imagesc(-90:18:90,-90:18:90,scoremap);data.cbar = colorbar();
data.focalpoint = plot(data.Phi,data.Theta,"Color",'r','Marker','o');
axis image;
subplot("Position",[ 0.4108    0.1091    0.2134    0.7530])
data.imax = imshow(G.visualize(sphere_norm*basis(1,:)));
subplot("Position",[0.6694    0.1806    0.3036    0.6784])
data.psthplot = plot(squeeze(psth_avg_tsr(1,:,6,6))); 
ylim([min(psth_avg_tsr,[],'all'),max(psth_avg_tsr,[],'all')])
data.timeL = line([1,1],ylim(),'Color','red','Visible',0);
data.suptitle = suptitle(compose("%s Manif Exp %02d pref chan %s", Animal, Expi, EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id)));

TheSLD = uicontrol(figh, 'Style', "slider", 'Position', [600 30 250 10], "String", "Theta", ...
            'Value', 0, "Callback", @(src, evt) ThetaSlider_Callback(src, evt));
TheSLD.Min = -pi/2; TheSLD.Max = pi/2;TheSLD.SliderStep = [0.005 0.05];
PhiSLD = uicontrol(figh, 'Style', "slider", 'Position', [600 5 250 10], "String", "Phi", ...
            'Value', 0, "Callback", @(src, evt) PhiSlider_Callback(src, evt));
PhiSLD.Min = -pi/2; PhiSLD.Max = pi/2;PhiSLD.SliderStep = [0.005 0.05];
AxisSLD = uicontrol(figh, 'Style', "slider", 'Position', [150 30 250 15], "String", "Dev on Axis", ...
            'Value', 0, "Callback", @(src, evt) AxisSlider_Callback(src, evt));
AxisSLD.Min = -2; AxisSLD.Max = 2;AxisSLD.SliderStep = [0.005 0.05];
TimeSLD = uicontrol(figh, 'Style', "slider", 'Position', [150 5 250 15], "String", "Time after onset", ...
            'Value', 1, "Callback", @(src, evt) TimeSlider_Callback(src, evt));
TimeSLD.Min = 1; TimeSLD.Max = 200;TimeSLD.Value=1;TimeSLD.SliderStep = [1 2]/199;
TG = uicontrol(figh, 'Style', 'togglebutton', "String", "RF Mask",...
    "Callback", @(src, evt) Mask_Callback(src, evt));
CT = uicontrol(figh, 'Style', 'pushbutton', 'Position', [20 40 80 30], "String", "Set Image",...
    "Callback", @(src, evt) SetImage_Callback(src, evt));
SA = uicontrol(figh, 'Style', 'pushbutton', 'Position', [20 70 80 30], "String", "Set Axis",...
    "Callback", @(src, evt) SetAxis_Callback(src, evt)); % Set 
PT = uicontrol(figh, 'Style', 'pushbutton', 'Position', [20 100 80 30], "String", "Play Tuning",...
    "Callback", @(src, evt) PlayTuning_Callback(src, evt));
PD = uicontrol(figh, 'Style', 'pushbutton', 'Position', [20 160 80 30], "String", "Play Dynamics",...
    "Callback", @(src, evt) PlayDynamics_Callback(src, evt));
TD = uicontrol(figh, 'Style', 'togglebutton', 'Position', [20 130 80 30], "String", "Toggle Time",...
    "Callback", @(src, evt) toggleDynamics_Callback(src, evt));
data.thetaSLD = TheSLD;data.phiSLD = PhiSLD;data.axisSLD = AxisSLD;data.timeSLD = TimeSLD;
data.toggleMask=TG;data.setImage=CT;data.setAxis=SA;data.playDynamics=PD;data.toggleDyn=TD;
guidata(figh, data);
%%
%% Code for bilinear interpolation in nd array
interpi = 6 + data.Theta / pi * 10;
interpj = 6 + data.Phi / pi * 10;
i_grid = [floor(interpi), ceil(interpi)];
j_grid = [floor(interpj), ceil(interpj)];
if interpi == floor(interpi)
    i_W = [1, 0];
else
    i_W = [i_grid(2) - interpi, interpi - i_grid(1)];
end
if interpj == floor(interpj)
    j_W = [1, 0];
else
    j_W = [j_grid(2) - interpj, interpj - j_grid(1)];
end
sum(psth_avg_tsr(:,:,i_grid,j_grid) .* reshape(i_W'*j_W,[1,1,2,2]),[1,3,4])
%%
function SetImage_Callback(hObj, evt)
[x,y] = ginput(1);
fprintf("%d,%d\n",x,y)
data = guidata(hObj);
data.Theta = y / 180 * pi;
data.Phi = x / 180 * pi;
data.thetaSLD.Value = data.Theta;
data.phiSLD.Value = data.Phi;
guidata(hObj,data);
draw(data.imax)
drawpsth(data.psthplot)
drawfocal(data.focalpoint)
end

function SetAxis_Callback(hObj, evt)
[x,y] = ginput(2);
fprintf("%d,%d\n",x,y)
data = guidata(hObj);
disp(size([x;y]))
data.axis = [x, y] / 180 * pi; % x, y are column vector
data.axis_cent = [x(1), y(1)] / 180 * pi;
data.axis_dir = [x(2) - x(1), y(2) - y(1)];
data.axis_dir = data.axis_dir / norm(data.axis_dir);
% disp(data.axis_cent)
% disp(data.axis_dir)
data.Theta = y(1) / 180 * pi;
data.Phi = x(1) / 180 * pi;
data.thetaSLD.Value = data.Theta;
data.phiSLD.Value = data.Phi;
data.axisSLD.Value = 0;
guidata(hObj,data);
draw(data.imax)
drawAxis(guidata(hObj).focalpoint, x, y)
end

function PlayTuning_Callback(hObj, evt)
data = guidata(hObj);
for devAxis = [0:0.05:data.axisSLD.Max, data.axisSLD.Max:-0.05:data.axisSLD.Min, data.axisSLD.Min:0.05:0]
data.DevAxis = devAxis;
pos = data.axis_cent + data.DevAxis * data.axis_dir;
data.Theta = pos(2);
data.Phi = pos(1);
guidata(hObj,data);
draw(data.imax)
drawpsth(data.psthplot)
drawfocal(data.focalpoint)
pause(0.025)
end
end

function toggleDynamics_Callback(hObj, evt)
data = guidata(hObj);
if hObj.Value
    global psth_avg_tsr psth_std_tsr
    caxis(data.tunemap.Parent, prctile(psth_avg_tsr,[2,98],'all'))
    data.timeL.Visible=true;
else
    data.timeL.Visible=false;
    data.tunemap.CData = data.scoremap;
    caxis(data.tunemap.Parent, prctile(data.scoremap,[2,98],'all'))
    draw(data.imax)
    drawpsth(data.psthplot)
    drawfocal(data.focalpoint)
end
end

function PlayDynamics_Callback(hObj, evt)
data = guidata(hObj);data.toggleDyn.Value=true;
global psth_avg_tsr psth_std_tsr
toggleDynamics_Callback(data.toggleDyn, evt)
for fi = 1:size(psth_avg_tsr,2)
data.tunemap.CData = squeeze(psth_avg_tsr(1,fi,:,:));
data.tunemap.Parent.Title.String = compose("%d ms",fi);
data.timeL.XData = [fi,fi];
pause(0.025)
end
end

function TimeSlider_Callback(hObj, evt)
data = guidata(hObj);
if data.toggleDyn.Value
global psth_avg_tsr psth_std_tsr
% caxis(data.tunemap.Parent, prctile(psth_avg_tsr,[2,98],'all'))
fi = clip(round(hObj.Value),[1,200]);
data.tunemap.CData = squeeze(psth_avg_tsr(1,fi,:,:));
data.tunemap.Parent.Title.String = compose("%d ms",fi);
data.timeL.XData = [fi,fi];
drawfocal(data.focalpoint)
pause(0.02)
end
end

function AxisSlider_Callback(hObj, evt)
% hObj.Min = -pi/2;hObj.Max = pi/2;
data = guidata(hObj);
data.DevAxis = hObj.Value;
pos = data.axis_cent + data.DevAxis * data.axis_dir;
data.Theta = pos(2);
data.Phi = pos(1);
guidata(hObj,data);
draw(data.imax)
drawpsth(data.psthplot)
drawfocal(data.focalpoint)
end

function ThetaSlider_Callback(hObj, evt)
% disp(hObj.Value)  
data = guidata(hObj);
data.Theta = hObj.Value;
guidata(hObj,data);
draw(data.imax)
drawpsth(data.psthplot)
drawfocal(data.focalpoint)
end

function PhiSlider_Callback(hObj, evt)
data = guidata(hObj);
data.Phi = hObj.Value;
guidata(hObj,data);
draw(data.imax)
drawpsth(data.psthplot)
drawfocal(data.focalpoint)
end

function Mask_Callback(hObj, evt)
data = guidata(hObj);
data.Addmask = hObj.Value;
guidata(hObj,data);
% guidata(hObj).Addmask = hObj.Value;
draw(data.imax)
end

function draw(imax)
global G 
data = guidata(imax);
data.curcode = [cos(data.Theta)*cos(data.Phi), sin(data.Theta)*cos(data.Phi), sin(data.Phi)] * data.basis * data.norm;
if ~ data.Addmask
    imax.CData = G.visualize(data.curcode);
else
    imax.CData = uint8( data.mask .* double(G.visualize(data.curcode)));
end
imax.Parent.Title.String = compose("%.1f, (%.1f, %.1f)",norm(data.curcode), data.Theta, data.Phi);
drawnow;
guidata(imax,data);
end

function drawpsth(psthplot)
global psth_avg_tsr psth_std_tsr
data = guidata(psthplot);
% nearest neighbor of interpolation 
% interpi = 6 + data.Theta / pi * 10;
% interpj = 6 + data.Phi / pi * 10;
% psthplot.YData = squeeze(psth_avg_tsr(1,:,round(interpi),round(interpj)));

% custom version of bilinear interpolation of the psth and sem tensor
interpi = clip(6 + data.Theta / pi * 10, [1,11]);
interpj = clip(6 + data.Phi   / pi * 10, [1,11]);
i_grid = [floor(interpi), ceil(interpi)];
j_grid = [floor(interpj), ceil(interpj)];
if interpi == floor(interpi)
    i_W = [1, 0];
else
    i_W = [i_grid(2) - interpi, interpi - i_grid(1)];
end
if interpj == floor(interpj)
    j_W = [1, 0];
else
    j_W = [j_grid(2) - interpj, interpj - j_grid(1)];
end
psth = sum(psth_avg_tsr(:,:,i_grid,j_grid) .* reshape(i_W'*j_W,[1,1,2,2]),[1,3,4]);
err = sum(psth_std_tsr(:,:,i_grid,j_grid) .* reshape(i_W'*j_W,[1,1,2,2]),[1,3,4]);
psthplot.YData = [psth,nan,psth+err,nan,psth-err]; % add nan to make 3 discontinuous lines
psthplot.XData = [1:200,nan,1:200,nan,1:200];
drawnow;
end

function drawAxis(focalpoint, X, Y)
focalpoint.XData = X;
focalpoint.YData = Y;
drawnow;
end

function drawfocal(focalpoint)
focalpoint.XData = guidata(focalpoint).Phi / pi * 180;
focalpoint.YData = guidata(focalpoint).Theta / pi * 180;
drawnow;
end