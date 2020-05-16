%% Using Python interface to use pytorch in matlab
pe = pyenv('Version','C:\ProgramData\Anaconda3\envs\tf-torch\python.exe'); % Set up the python executable
py.importlib.import_module('torch');
py.importlib.import_module('torchvision');
py.importlib.import_module('torchvision.models'); % Load some deep learning framework
SNet = py.torchvision.models.vgg16(pyargs("pretrained",1)); % Load a pretrained neural net 
out = SNet.forward(py.torch.rand(py.int(1),py.int(3),py.int(224),py.int(224))); % process sth in python torch framework
out.data.numpy().double();% Get the data back in matlab data type
% py.torch.Tensor([1,2,3])