function [evc_all, eva_all, evc_cls, eva_cls, evc_nos, eva_nos] = loadHessian(GAN, param)
if nargin<2, param = [];end
if nargin<1, GAN="BigGAN";end
switch getenv('COMPUTERNAME')
    case 'PONCELAB-ML2A'
	HessRoot = fullfile(getenv('HOMEPATH'),'Documents','Python');
	case 'PONCELAB-ML2B'
	HessRoot = fullfile(getenv('HOMEPATH'),'Documents','Python');
	otherwise
	HessRoot = "E:\OneDrive - Washington University in St. Louis\Hessian_summary\BigGAN";
end
if strcmp(GAN, "BigGAN")
Hpath = fullfile(HessRoot, "H_avg_1000cls.npz");
% So the eigen vectors are the columns of the `evc_*` matrix. Sorted in
% descending order, the same as `eva_*`
Hdata = py.numpy.load(Hpath);
eva_all = Hdata.get('eigvals_avg').double;
evc_all = Hdata.get('eigvects_avg').double;
evc_all = evc_all(:,end:-1:1);
eva_all = eva_all(end:-1:1);
eva_cls = Hdata.get('eigvals_clas_avg').double;
evc_cls = Hdata.get('eigvects_clas_avg').double;
evc_cls = evc_cls(:,end:-1:1);
eva_cls = eva_cls(end:-1:1);
eva_nos = Hdata.get('eigvals_nois_avg').double;
evc_nos = Hdata.get('eigvects_nois_avg').double;
evc_nos = evc_nos(:,end:-1:1);
eva_nos = eva_nos(end:-1:1);
evc_nos_ag = [evc_nos;zeros(128)];
evc_cls_ag = [zeros(128);evc_cls];

elseif strcmp(GAN, "FC6GAN")
switch param
	case "evol"
		npyname = "Evolution_Avg_Hess.npz";
	case "pasu"
		npyname = "Pasu_Space_Avg_Hess.npz";
	otherwise	
		npyname = "Evolution_Avg_Hess.npz";
end
Hpath = fullfile(HessRoot, npyname);
Hdata = py.numpy.load(Hpath);
eva_all = Hdata.get('eigv_avg').double;
evc_all = Hdata.get('eigvect_avg').double;
evc_all = evc_all(:,end:-1:1);
eva_all = eva_all(end:-1:1);
evc_cls = [];
eva_cls = [];
evc_nos = [];
eva_nos = [];

end
end