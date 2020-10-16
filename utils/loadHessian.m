function [evc_all, eva_all, evc_cls, eva_cls, evc_nos, eva_nos] = loadHessian(GAN, param)
if nargin<2, param = [];end
if nargin<1, GAN="BigGAN";end
if strcmp(GAN, "BigGAN")
Hdata = py.numpy.load("E:\OneDrive - Washington University in St. Louis\Hessian_summary\BigGAN\H_avg_1000cls.npz");
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
end
end