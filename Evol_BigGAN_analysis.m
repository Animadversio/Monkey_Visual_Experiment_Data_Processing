%% Evol BigGAN Analysis
Animal = "Both";Set_Path;
%"200803","200804","200805"
expftr = contains(ExpRecord.expControlFN,["200807","200810"]); %& contains(ExpRecord.Exp_collection,"BigGAN_Hessian");% & contains(ExpRecord.Exp_collection,"BigGAN");
fllist = find(expftr);no_return=false;
[meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(fllist(1:end),Animal,no_return);%find(expftr)%find(expftr)
%%
median(param,1) 
mean(param,1)
%%
% all -1.0000   -0.5000   -2.4374   -4.6619   -3.2552    0.4563
%%
expnm = "E:\Cluster_Backup\BigGAN_invert\cute_cat_rsz_all";
expresult = readtable(fullfile(expnm,"BigGAN_Hess_Adam_optim_BO_tune600.csv"));
%
[sorted, idx]=sort(expresult.score,'ascend');
expresult(idx,:);
median(table2array(expresult(idx(1:10),:)),1)