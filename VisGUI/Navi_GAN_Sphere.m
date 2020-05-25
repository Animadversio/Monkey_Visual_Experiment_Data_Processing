%% GUI for Navigate Manifold Hemisphere through theta phi scroll bar
addpath DNN
global G 
G = FC6Generator("matlabGANfc6.mat");
%%
pe = pyenv('Version','C:\Users\binxu\.conda\envs\caffe36\python.exe'); % Note the python env could not be changed in a matlab session
py.importlib.import_module('numpy');
f = py.numpy.load("S:\Stimuli\2019-Manifold\alfa-191119a\backup_11_19_2019_11_58_11\PC_imgs\PC_vector_data.npz");
PC_Vec = f.get('PC_vecs').double;
sphere_norm = f.get('sphere_norm').double;
basis = PC_Vec(1:3,:);
f.close();
%%
figh = figure(2);
imax = imshow(G.visualize(basis(1,:)));
guidata(figh,struct('basis', basis, 'norm', sphere_norm, 'mask', mask, ...
            'Addmask', 0, "Theta", 0, "Phi", 0, "curcode", basis(1,:)));
TheSLD = uicontrol(figh, 'Style', "slider", 'Position', [100 30 250 10], "String", "Norm", ...
            'Value', 0, "Callback", @(src, evt) ThetaSlider_Callback(src, evt, imax));
TheSLD.Min = -pi; TheSLD.Max = pi;TheSLD.SliderStep = [0.005 0.05];
PhiSLD = uicontrol(figh, 'Style', "slider", 'Position', [100 5 250 10], "String", "Norm", ...
            'Value', 0, "Callback", @(src, evt) PhiSlider_Callback(src, evt, imax));
PhiSLD.Min = -pi/2; PhiSLD.Max = pi/2;PhiSLD.SliderStep = [0.005 0.05];
TG = uicontrol(figh, 'Style', 'togglebutton', "String", "RF Mask","Callback", @(src, evt) Mask_Callback(src, evt, imax));
%%
function ThetaSlider_Callback(hObj, evt, imax)
% disp(hObj.Value)  
data = guidata(hObj);
data.Theta = hObj.Value;
guidata(hObj,data);
draw(imax)
end
function PhiSlider_Callback(hObj, evt, imax)
data = guidata(hObj);
data.Phi = hObj.Value;
guidata(hObj,data);
draw(imax)
end
function Mask_Callback(hObj, evt, imax)
data = guidata(hObj);
data.Addmask = hObj.Value;
guidata(hObj,data);
% guidata(hObj).Addmask = hObj.Value;
draw(imax)
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