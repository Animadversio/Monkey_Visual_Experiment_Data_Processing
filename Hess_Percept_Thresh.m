%% This script dedicates to find the perceptible threshold for each eigen axes 
% and use it to threshold projection
%  
%%
G = FC6Generator();
%%
py.importlib.import_module("numpy")
% data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\HessTune\NullSpace\Evolution_Avg_Hess.npz");
data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy\Nullspace\Pasu_Space_Avg_Hess.npz");
eigvect = data.get("eigvect_avg").double; % each column is an eigen vector. 
eigvals = data.get("eigv_avg").double;
eigvals = eigvals(end:-1:1);
eigvect = eigvect(:,end:-1:1);
%%
%%
sel_list = [1,2,5,10,20,50,80,120,160,200,250,300,400,500,600,800,1000,1500,2000,2500,3000,3500];
figure;montage(G.visualize(2*-percp_thrsh(sel_list).*eigvect(:,sel_list)))
%%
percp_thrsh = 0.3*abs(eigvals).^(-1/3);
expi = [1,5,10,20,30,80];
expcodes_cc = code_proj_cc(expi,:);
expcodes_cc((expcodes_cc<percp_thrsh) & (expcodes_cc>-percp_thrsh)) = 0;
%
figure(25);
subplot(121)
montage(G.visualize(expcodes_cc*eigvect'))
subplot(122)
montage(G.visualize(code_proj_cc(expi,:)*eigvect'))
fprintf("Sparsity\n")
fprintf("%d/4096 \n",4096-sum(expcodes_cc==0,2))
%% Threshold the evolved codes and see where are their genes residing in. 
percept_dims = code_proj_cc>percp_thrsh | code_proj_cc<-percp_thrsh;
percept_dims_perm = perm_proj_cc > percp_thrsh | perm_proj_cc < -percp_thrsh;
percept_dims_ctrl = ctrl_proj_cc > percp_thrsh | ctrl_proj_cc < -percp_thrsh;
percept_dims_ctrl_perm = ctrl_perm_proj_cc > percp_thrsh | ctrl_perm_proj_cc < -percp_thrsh;

figure(10);clf;hold on
plot(mean(percept_dims,1))
plot(mean(percept_dims_perm,1))
plot(mean(percept_dims_ctrl,1))
plot(mean(percept_dims_ctrl_perm,1))
%%
