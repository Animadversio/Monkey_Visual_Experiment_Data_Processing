%% Corr_interpretation (Kinda obsolete.... use python interface)
exportONNXNetwork(net,"N:\vgg16mat.onnx")
%% Python interface building
% compare the vgg16 in matlab vs that in python

pe = pyenv('Version','C:\ProgramData\Anaconda3\envs\tf-torch\python.exe'); % Set up the python executable
py.importlib.import_module('torch');
py.importlib.import_module('torchvision');
py.importlib.import_module('torchvision.models'); % Load some deep learning framework
pyvgg = py.torchvision.models.vgg16(pyargs("pretrained",1)); % Load a pretrained neural net 
% out = SNet.forward(py.torch.rand(py.int(1),py.int(3),py.int(224),py.int(224))); % process sth in python torch framework
% out.data.numpy().double();
%%
matvgg = vgg16;
% matlab doesn't have std rescale, only 0 mean centered
RGBmean = [123.6800  116.7790  103.9390];squeeze(matvgg.Layers(1).Mean(1,1,:));
%%
vggmatW = struct();
for i=1:numel(matvgg.Layers)
%     fprintf(layer.Name) 
    layer = matvgg.Layers(i);
    if any(strcmp('Weights', fieldnames(layer)))
        fprintf("%s:W [%s], B [%s]\n",layer.Name,num2str(size(layer.Weights)),num2str(size(layer.Bias))) 
%         vggmatW.(compose("%s_num",layer.Name)) = i;
        vggmatW.(compose("%s_weight",layer.Name)) = layer.Weights;
        vggmatW.(compose("%s_bias",layer.Name)) = layer.Bias;
    end
end
save("N:\vgg16W.mat","vggmatW","-v6")
%%
%%
outtmp = squeeze(activations(matvgg,ones(224,224,3,1),'fc8'));
corr(pyact,outtmp)
%%
pyact = py.numpy.load("S:\fc8out.npy").single';
%%

%%
catimg = imread("E:\Monkey_Data\Generator_DB_Windows\nets\upconv\Cat.jpg");
cattsr = catimg(39:end-38,39:end-38,:);
outtmp = squeeze(activations(matvgg,cattsr,'fc8'));
corr(pyact,outtmp)
%%
midtmp = squeeze(activations(matvgg,cattsr,'pool5'));
pyactp5 = py.numpy.load("S:\poo5out.npy").single';
pool5_cc = corr(reshape(permute(pyactp5,[3,4,2,1]),[],1), reshape(midtmp,[],1));
fprintf("%.8f",pool5_cc)
%% fc6
fc6tmp = squeeze(activations(matvgg,cattsr,'fc6'));
pyactfc6 = py.numpy.load("S:\fc6out.npy").single';
fc6_cc = corr(reshape(permute(pyactfc6,[3,4,2,1]),[],1), reshape(fc6tmp,[],1));
fc6_L1 = mean(abs(pyactfc6-fc6tmp));
fprintf("%.8f, %.3f\n",fc6_cc,fc6_L1)
%%
pyorgfc6 = py.numpy.load("S:\fc6out_orig.npy").single';
%%
figure(1);clf;hold on;plot(pyorgfc6);plot(fc6tmp);plot(pyactfc6);
title(compose("VGG fc6 output\nCorrCoef %.3f L1 %.3f\nOriginal CorrCoef %.3f L1 %.3f",fc6_cc,fc6_L1,...
    corr(pyorgfc6,fc6tmp),mean(abs(pyorgfc6-fc6tmp))))
legend(["pytorch original VGG","matlab VGG","pytorch VGG"])
%%
figure;hold on;plot(outtmp);plot(pyact)
title(compose("VGG fc8 output\nCorrCoef %.3f",corr(pyact,outtmp)))
legend(["matlab VGG","pytorch VGG"])
%%
figure;hold on;
plot(pyact-outtmp);plot(pyact)
%%

%%
figure;hist(pyact-outtmp)