pyenv('Version',"C:\Users\binxu\.conda\envs\tf\python.exe")
syspath = py.sys.path(); % add the official stylegan2 repo. 
syspath.append("E:\Github_Projects\BigBiGAN");
py.importlib.import_module('BigBiGAN'); % doesn't work