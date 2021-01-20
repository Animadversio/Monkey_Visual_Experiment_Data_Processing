G = FC6Generator();
D = torchImDist();
%%
py.importlib.import_module("numpy")
data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\HessTune\NullSpace\Texture_Avg_Hess.npz");
% data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy\Nullspace\Pasu_Space_Avg_Hess.npz");
eigvect = data.get("eigvect_avg").double; % each column is an eigen vector. 
eigvals = data.get("eigv_avg").double;
eigvals = eigvals(end:-1:1);
eigvect = eigvect(:,end:-1:1);