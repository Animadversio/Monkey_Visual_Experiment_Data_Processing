addpath DNN
addpath utils
system("subst N: E:\Network_Data_Sync")
% addpath("C:\Users\binxu\OneDrive - Washington University in St. Louis\Matlab_add_on\npy-matlab-master\npy-matlab")
pe = pyenv('Version','C:\Users\binxu\.conda\envs\caffe36\python.exe'); % Note the python env could not be changed in a matlab session
%%
system("subst S: E:\Network_Data_Sync")
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir,"Alfa_ManifPopDynamics.mat"));
load(fullfile(mat_dir,"Alfa_Manif_stats.mat"));
load(fullfile(mat_dir,"Alfa_Evol_stats.mat"));
%%
global G 
G = FC6Generator("matlabGANfc6.mat");
%%
code = randn(1,4096);
imgs = G.visualize(code);
%%
% mask = zeros(256);
[XX, YY] = meshgrid(1:256,1:256);
D = sqrt((XX-127).^2 + (YY-127).^2);
mask = exp(-(D-40).^2/100);
mask = max((mask ./ max(mask,[],'all')), D<40);
%%
Expi = 3;
basis_path = fullfile(Stats(Expi).meta.stimuli,"PC_vector_data.npz");
%% Load the source data for Manifold Experiments from python
py.importlib.import_module('numpy');
f = py.numpy.load(basis_path);
% "S:\Stimuli\2019-Manifold\alfa-191119a\backup_11_19_2019_11_58_11\PC_imgs\PC_vector_data.npz"
PC_Vec = f.get('PC_vecs').double;
sphere_norm = f.get('sphere_norm').double;
basis = [-1,1,1]' .* PC_Vec(1:3,:);% Note the final PC may need to reverse! not always the same dir!
% f.close();
%% Prepare data 
global psth_avg_tsr psth_std_tsr
si=1;ui=1;
meanpsth = cellfun(@(psth) mean(psth(ui,:,:),3), Stats(Expi).manif.psth{si}, "UniformOutput", false);
stdpsth = cellfun(@(psth)  std(psth(ui,:,:),0,3), Stats(Expi).manif.psth{si}, "UniformOutput", false);
psth_avg_tsr = cell2mat(reshape(meanpsth,1,1,11,11)); % reshape and transform, so size 86, 200, 11, 11, key data! 
psth_std_tsr = cell2mat(reshape(stdpsth,1,1,11,11)); % reshape and transform, so size 86, 200, 11, 11, key data! 
scoremap = squeeze(mean(psth_avg_tsr(1,50:200,:,:),[1,2]));
%%
figh = figure(3);clf;figh.Position = [28         524        1446         454];
subplot(131);hold on
tunemap = imagesc(-90:18:90,-90:18:90,scoremap);
focalpoint = plot(data.Phi,data.Theta,"Color",'r','Marker','o');
axis image;
subplot(132)
imax = imshow(G.visualize(sphere_norm*basis(1,:)));
subplot(133)
psthplot = plot(squeeze(psth_avg_tsr(1,:,6,6)));
ylim([min(psth_avg_tsr,[],'all'),max(psth_avg_tsr,[],'all')])
data = struct('basis', basis, 'norm', sphere_norm, 'mask', mask, ...
            'Addmask', 0, "Theta", 0, "Phi", 0, "curcode", basis(1,:), ...
            "imax", imax, "psthplot", psthplot, "focalpoint", focalpoint);
guidata(figh, data);
TheSLD = uicontrol(figh, 'Style', "slider", 'Position', [600 30 250 10], "String", "Theta", ...
            'Value', 0, "Callback", @(src, evt) ThetaSlider_Callback(src, evt));
TheSLD.Min = -pi/2; TheSLD.Max = pi/2;TheSLD.SliderStep = [0.005 0.05];
PhiSLD = uicontrol(figh, 'Style', "slider", 'Position', [600 5 250 10], "String", "Phi", ...
            'Value', 0, "Callback", @(src, evt) PhiSlider_Callback(src, evt));
PhiSLD.Min = -pi/2; PhiSLD.Max = pi/2;PhiSLD.SliderStep = [0.005 0.05];
TG = uicontrol(figh, 'Style', 'togglebutton', "String", "RF Mask",...
    "Callback", @(src, evt) Mask_Callback(src, evt));
CT = uicontrol(figh, 'Style', 'pushbutton', 'Position', [20 40 60 30], "String", "Set Image",...
    "Callback", @(src, evt) SetImage_Callback(src, evt));
SA = uicontrol(figh, 'Style', 'pushbutton', 'Position', [20 70 60 30], "String", "Set Axis",...
    "Callback", @(src, evt) SetAxis_Callback(src, evt));
AxisSLD = uicontrol(figh, 'Style', "slider", 'Position', [100 10 250 15], "String", "Dev on Axis", ...
            'Value', 0, "Callback", @(src, evt) AxisSlider_Callback(src, evt));
AxisSLD.SliderStep = [0.005 0.05];
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
disp(data.axis_cent)
disp(data.axis_dir)
data.Theta = y(1) / 180 * pi;
data.Phi = x(1) / 180 * pi;
guidata(hObj,data);
draw(data.imax)
drawAxis(guidata(hObj).focalpoint, x, y)
end
function AxisSlider_Callback(hObj, evt)
hObj.Min = -pi/2;
hObj.Max = pi/2;
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
% interpi = 6 + data.Theta / pi * 10;
% interpj = 6 + data.Phi / pi * 10;
% psthplot.YData = squeeze(psth_avg_tsr(1,:,round(interpi),round(interpj)));
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
psthplot.YData = sum(psth_avg_tsr(:,:,i_grid,j_grid) .* reshape(i_W'*j_W,[1,1,2,2]),[1,3,4]);
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