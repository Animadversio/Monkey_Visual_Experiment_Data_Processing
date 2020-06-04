% compare the vgg16 in matlab vs that in python, to make the network
% weights uniform. 
%% Corr_interpretation (Kinda obsolete.... use python interface)
exportONNXNetwork(net,"N:\vgg16mat.onnx")
%% Python interface building

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
%% export weight of matlab vgg16 into mat file for import in python. 
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
outtmp = squeeze(activations(matvgg,ones(224,224,3,1),'fc8'));
corr(pyact,outtmp)
%%
catimg = imread("E:\Monkey_Data\Generator_DB_Windows\nets\upconv\Cat.jpg");
cattsr = catimg(39:end-38,39:end-38,:);
outtmp = squeeze(activations(matvgg,cattsr,'fc8'));
corr(pyact,outtmp)
%%
figure;hold on;plot(outtmp);plot(pyact)