global G curcode mask Addmask
G = FC6Generator("matlabGANfc6.mat");
%%
code = randn(1,4096);
curcode = code;
imgs = G.visualize(code);
%%
% mask = zeros(256);
[XX, YY] = meshgrid(1:256,1:256);
D = sqrt((XX-127).^2 + (YY-127).^2);
mask = exp(-(D-40).^2/100);
mask = max((mask ./ max(mask,[],'all')), D<40);
Addmask = False;
%%
figh = figure(1);
imax = imshow(G.visualize(code));draw(imax);
SLD = uicontrol(figh, 'Style', "slider", 'Position', [100 20 250 20], "String", "Norm", ...
            'Value', 0, "Callback", @(src, evt) Slider_Callback(src, evt, imax, code));
SLD.Min = -5; SLD.Max = 5; SLD.SliderStep = [0.005 0.05];
TG = uicontrol(figh, 'Style', 'togglebutton', "String", "RF Mask","Callback", @(src, evt) Mask_Callback(src, evt, imax));
%% 
function Slider_Callback(hObj, evt, imax, code)
global curcode
disp(hObj.Value)    
curcode = code * hObj.Value;
draw(imax)
end
function Mask_Callback(hObj, evt, imax)
global Addmask
Addmask = hObj.Value;
draw(imax)
end
function draw(imax)
global G curcode Addmask mask
if ~ Addmask
    imax.CData = G.visualize(curcode);
else
    imax.CData = uint8( mask .* double(G.visualize(curcode)));
end
imax.Parent.Title.String = compose("%.1f",norm(curcode));
drawnow;
end
%%
%%
% figh = figure(2);
% imax = imshow(G.visualize(code));
% guidata(figh,struct('curcode', code, 'mask', mask, 'Addmask', 0));
% SLD = uicontrol(figh, 'Style', "slider", 'Position', [100 20 250 20], "String", "Norm", ...
%             'Value', 0, "Callback", @(src, evt) Slider_Callback(src, evt, imax, code));
% SLD.Min = -5; SLD.Max = 5;SLD.SliderStep = [0.005 0.05];
% TG = uicontrol(figh, 'Style', 'togglebutton', "String", "RF Mask","Callback", @(src, evt) Mask_Callback(src, evt, imax));
% %%
% function Slider_Callback(hObj, evt, imax, code)
% % disp(hObj.Value)  
% data = guidata(hObj);
% data.curcode = code * hObj.Value;
% guidata(hObj,data);
% draw(imax)
% end
% function Mask_Callback(hObj, evt, imax)
% data = guidata(hObj);
% data.Addmask = hObj.Value;
% guidata(hObj,data);
% % guidata(hObj).Addmask = hObj.Value;
% draw(imax)
% end
% function draw(imax)
% global G 
% if ~ guidata(imax).Addmask
%     imax.CData = G.visualize(guidata(imax).curcode);
% else
%     imax.CData = uint8( guidata(imax).mask .* double(G.visualize(guidata(imax).curcode)));
% end
% imax.Parent.Title.String = compose("%.1f",norm(guidata(imax).curcode));
% drawnow;
% end