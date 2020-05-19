addpath DNN
addpath utils
%%
global G
G = FC6Generator("matlabGANfc6.mat");
%%
code = randn(1,4096);
imgs = G.visualize(code);
%%
figh = figure(1);
imax = imshow(G.visualize(code));
SLD = uicontrol(figh, 'Style', "slider", 'Position', [100 20 250 20], "String", "Norm", ...
            'Value', 0, "Callback", @(src, evt) Slider_Callback(src, evt, imax, code));
SLD.Min = -5; SLD.Max = 5; SLD.SliderStep = [0.005 0.05];
TG = uicontrol(figh, 'Style', 'togglebutton', "String", "RF Mask");
%% 
function Slider_Callback(hObj, eventdata, imax, code)
global G
disp(hObj.Value)    
imax.CData = G.visualize(code * hObj.Value);
imax.Parent.Title.String = compose("%.1f",hObj.Value*norm(code));
drawnow;
end