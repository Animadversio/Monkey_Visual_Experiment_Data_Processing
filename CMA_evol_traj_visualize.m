% 
% abs(Get) inspired by BigGAN version BigGAN_Hess_Frame_Evol_Analysis.
% updated from Hess_Frame_Evol_Test_Summary
% used for visualize the PC structure of evolutionary trajectory
%% 
load("D:\Project_CMA_Monkeys.mat");
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
ctrl_evol_dir = "N:\Data-Ephys-MAT\Project_CMA_ghostEvolutions_matchingNorm";
fnlist = string(ls(ctrl_evol_dir));
%% Load the initial codes
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir,"evo_drift_codes_all.mat"),'evol_codes','drft_codes','evo_meta')
load(fullfile(mat_dir,"evol_ctrl_codes.mat"),'drft_codes_col') % Control code by stop at the same mean norm
load(fullfile(mat_dir,"evol_ctrl_gen_codes.mat"),'drft_codes_gen_col') % Control code by stop at the same generation number
Data = load("N:\Code\Generator_DB_Windows\init_population\init_code.mat");
initcodes = Data.codes;
drft_codes_nrm = cell2mat(drft_codes_col); % stop at the same mean norm
drft_codes_gen = cell2mat(drft_codes_gen_col); % stop at the same generation. 
%% Collect data into form 
evol_fincodes = [];
PC_coef_col = {};
MonkName = ["Alfa","Beto"];
tic
for Mi = 1:2
   Animal = MonkName(Mi);
   for Expi = 1:numel(MonkeysCMA{Mi}.Stats)
      generations = MonkeysCMA{Mi}.Stats{Expi}.genes.gen;
      nGen = max(generations);
      meancodes_all = arrayfun(@(iGen)mean(MonkeysCMA{Mi}.Stats{Expi}.genes.all(iGen==generations,:),1),1:nGen,'uni',0);
      meancodes_all = cell2mat(meancodes_all');
      [PCvec, projcoef, latent] = pca(meancodes_all,'num',10);
      meanfincode = mean(MonkeysCMA{Mi}.Stats{Expi}.genes.all(nGen==generations,:),1);
      evol_fincodes = [evol_fincodes; meanfincode];
      PC_coef_col{end+1} = projcoef;
      fprintf("%s-Exp%03d %.2f\n",Animal,Expi,toc)
   end
end
%%
figdir = "E:\OneDrive - Harvard University\GECCO2022\Figures\Sinusoidal_invivo";
%% Plot the projected curves 
figure(2);clf
T = tiledlayout(4,1,'TileSp','Compact','Padding','compact');
for expi = 1:numel(PC_coef_col)
    for PCi = 1:4
    nexttile(T,PCi);hold on
    sgn = -1 * sign(PC_coef_col{expi}(1,PCi));
    plot(sgn*PC_coef_col{expi}(:,PCi),'Color',[0,0,0,.2])
    end
end
for PCi = 1:4
    nexttile(T,PCi);
    ylabel(compose("PC %d",PCi));
    xlim([0,130]);line([1,130],[0,0],'color','r','linestyle','-.')
end
for PCi = 1:3, nexttile(T,PCi);xticklabels([]); end
nexttile(T,4);xlabel("Generations")
title(T,compose("Top PC Projection Coefficients of Evol Trajectory\n Monkey A,B  N=%d",numel(PC_coef_col)))
%%
saveallform(figdir,"evol_PCA_projcoefs_CMA_all",2)


%%
figure;
eigid = 1:numel(expvar);
N = numel(eigid);
expvar_thry = 3./(1-cos(pi*eigid/N))/(N^2-1);
expvar_thry2 = 6./pi^2./eigid.^2;
loglog(1:numel(expvar),expvar/100);hold on
loglog(1:numel(expvar),expvar_thry)
loglog(1:numel(expvar),expvar_thry2)





%% Structure in Hessian Eigenframe
py.importlib.import_module("numpy")
% data = py.numpy.load("N:\Code\Hessian\Texture_Avg_Hess.npz");
data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\HessTune\NullSpace\Texture_Avg_Hess.npz");
% data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy\Nullspace\Pasu_Space_Avg_Hess.npz");
eigvect = data.get("eigvect_avg").double; % each column is an eigen vector. 
eigvals = data.get("eigv_avg").double;
eigvals = eigvals(end:-1:1);
eigvect = eigvect(:,end:-1:1);
%%
evol_fincodes_norm = evol_fincodes ./ vecnorm(evol_fincodes,2,2);
[pop_PCvec, pop_projcoef, ~,~,pop_expvar] = pca(evol_fincodes);
[pop_PCvec, pop_projcoef, ~,~,pop_expvar] = pca(evol_fincodes_norm);
%%
initprojs = initcodes * eigvect;
projcoefs = evol_fincodes * eigvect;
devicoefs = (evol_fincodes-mean(initcodes,1)) * eigvect;
normprojcoefs = evol_fincodes_norm * eigvect;
%%
cutoff = 600;
eigenrange = 1:700;%800:4000;%1:cutoff;
cc_vec         = corr(abs(devicoefs(:,eigenrange)'),eigvals(eigenrange)');
cc_vec_drftnrm = corr(abs(drftnrm_devicoefs(:,eigenrange)'),eigvals(eigenrange)');
cc_vec_drftgen = corr(abs(drftgen_devicoefs(:,eigenrange)'),eigvals(eigenrange)');
cc_vec_log         = corr(abs(devicoefs(:,eigenrange)'),log10(eigvals(eigenrange)'));
cc_vec_log_drftnrm = corr(abs(drftnrm_devicoefs(:,eigenrange)'),eigvals(eigenrange)');
cc_vec_log_drftgen = corr(abs(drftgen_devicoefs(:,eigenrange)'),log10(eigvals(eigenrange)'));

ttest2corr_print(cc_vec,cc_vec_drft)
ttest2corr_print(cc_vec_log,cc_vec_log_drftgen)
%%
figure;hold on
histogram(cc_vec)
% histogram(cc_vec_drftnrm)
histogram(cc_vec_drftgen)

figure;hold on
histogram(cc_vec_log)
% histogram(cc_vec_drftnrm)
histogram(cc_vec_log_drftgen)
%%
%%
figure()
cdfplot(projcoefs(:,end));hold on
cdfplot(projcoefs(:,end-2000))
%%
figure;hold on
scatter(log10(abs(eigvals)),abs(mean(initprojs,1)),12,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
scatter(log10(abs(eigvals)),abs(mean(devicoefs,1)),12,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
scatter(log10(abs(eigvals)),abs(mean(projcoefs,1)),12,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
%%
figure;hold on
% scatter(log10(abs(eigvals)),mean(abs(initprojs),1),12,...
%     'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
% scatter(log10(abs(eigvals)),mean(abs(devicoefs),1),12,...
%     'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
scatter(log10(abs(eigvals)),mean(abs(projcoefs),1),12,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
scatter(log10(abs(eigvals)),mean(abs(drftgen_projcoefs(1:5:end,:)),1),12,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
scatter(log10(abs(eigvals)),mean(abs(drftnrm_projcoefs),1),12,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)

%%
figure;hold on 
for iexp  = 1:264
scatter(log10(abs(eigvals)),normprojcoefs(iexp,:),6,'k',...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.01)
end
%% KS test for each dimension independently ,,, not a good idea. 
Pcol = [];
KScol = [];
for idim = 1:4096
    [H,P,KSSTAT] = kstest2(reshape(drftgen_projcoefs(:,idim),1,[]),...
                           reshape(evol_fincodes(:,idim),1,[]));
    Pcol = [Pcol,P];
    KScol = [KScol,KSSTAT];
end
%%
figure;
hold on;
histogram(vecnorm(evol_fincodes,2,2))
histogram(vecnorm(drftnrm_projcoefs(1:5:end,:),2,2))
histogram(vecnorm(drftgen_projcoefs,2,2))
%%

%%
figure;hold on
scatter(log10(abs(eigvals)),mean(abs(projcoefs),1),12,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
% scatter(log10(abs(eigvals)),mean(abs(initprojs),1),12,...
%     'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
scatter(log10(abs(eigvals)),mean(abs(devicoefs),1),12,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
%%
figure;hold on
scatter((abs(eigvals)),mean(abs(projcoefs),1),...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
scatter((abs(eigvals)),mean(abs(initprojs),1),...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
scatter((abs(eigvals)),mean(abs(devicoefs),1),...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
%%
cutoff = sum(eigvals > 1E-4);
[rval,pval] = corr(eigvals(1:cutoff)', mean(abs(initprojs(:,1:cutoff)),1)');
fprintf("init gen vector corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
[rval,pval] = corr(eigvals(1:cutoff)', mean(abs(projcoefs(:,1:cutoff)),1)');
fprintf("Evol vector corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
[rval,pval] = corr(eigvals(1:cutoff)', mean(abs(devicoefs(:,1:cutoff)),1)');
fprintf("Deviation from init gen mean corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
%%
%%
drftnrm_projcoefs = drft_codes_nrm * eigvect;
drftnrm_devicoefs = (drft_codes_nrm-mean(initcodes,1)) * eigvect;
drftgen_projcoefs = drft_codes_gen * eigvect;
drftgen_devicoefs = (drft_codes_gen-mean(initcodes,1)) * eigvect;

[rval,pval] = corr((eigvals(1:cutoff))', mean(abs(drftnrm_projcoefs(:,1:cutoff)),1)');
fprintf("Drift evolution, same norm, gen mean corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
[rval,pval] = corr((eigvals(1:cutoff))', mean(abs(drftnrm_devicoefs(:,1:cutoff)),1)');
fprintf("Drift evolution, same norm Deviation from init gen mean corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
[rval,pval] = corr((eigvals(1:cutoff))', mean(abs(drftgen_projcoefs(:,1:cutoff)),1)');
fprintf("Drift evolution, same gen corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
[rval,pval] = corr((eigvals(1:cutoff))', mean(abs(drftgen_devicoefs(:,1:cutoff)),1)');
fprintf("Drift evolution, same norm Deviation from init gen mean corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
%%
[rval,pval] = corr(eigvals(1:cutoff)', abs(mean(drftnrm_projcoefs(:,1:cutoff),1))');
fprintf("Drift evolution, same norm, gen mean corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
[rval,pval] = corr(eigvals(1:cutoff)', abs(mean(drftnrm_devicoefs(:,1:cutoff),1))');
fprintf("Drift evolution, same norm Deviation from init gen mean corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
[rval,pval] = corr(eigvals(1:cutoff)', abs(mean(drftgen_projcoefs(:,1:cutoff),1))');
fprintf("Drift evolution, same gen corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
[rval,pval] = corr(eigvals(1:cutoff)', abs(mean(drftgen_devicoefs(:,1:cutoff),1))');
fprintf("Drift evolution, same norm Deviation from init gen mean corr=%.3f P=%.1e (H eig cutoff=%d)\n",rval,pval,cutoff)
%%
[H,P,KSSTAT] = kstest2(drftnrm_projcoefs(:,end-1),drftnrm_projcoefs(:,end-2))
%%
[H,P,KSSTAT] = kstest2(reshape(drftgen_projcoefs(:,cutoff:end),1,[]),reshape(evol_fincodes(:,cutoff:end),1,[]))
%%

figure('pos',[680   498   430   400]);hold on
% scatter(log10(abs(eigvals)),mean(abs(initprojs),1),12,...
%     'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
% scatter(log10(abs(eigvals)),abs(mean(devicoefs,1)),12,...
%     'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
% scatter(log10(abs(eigvals)),abs(mean(drftnrm_devicoefs(1:5:end,:),1)),12,...
%     'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
% scatter(log10(abs(eigvals)),abs(mean(drftgen_devicoefs(1:5:end,:),1)),12,...
%     'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15)
scatter(log10(abs(eigvals)),mean(abs(projcoefs),1),6,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.10,...
    'Disp',"Evolution (N=264)")
scatter(log10(abs(eigvals)),mean(abs(drftgen_projcoefs(1:5:end,:)),1),6,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.10,...
    'Disp',"Drift Evolution (match gen) (N=264)")
scatter(log10(abs(eigvals)),mean(abs(drftnrm_projcoefs(1:5:end,:)),1),6,...
    'filled','MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.10,...
    'Disp',"Drift Evolution (match norm) (N=264)")
legend('location','best')
xlabel("log10(eigval)")
ylabel("mean abs amplitude")
ylim([2.5,7.5]);xlim([-12,-1])
vline(log10(eigvals(801)),'-.k')
%%
mdl = fitlm(log10(abs(eigvals(1:800))),mean(abs(projcoefs(:,1:800)),1));
plot(mdl)
mdl2 = fitlm(log10(abs(eigvals(1:800))),mean(abs(drftgen_projcoefs(1:5:end,1:800)),1));
plot(mdl2)
mdl3 = fitlm(log10(abs(eigvals(1:800))),mean(abs(drftnrm_projcoefs(1:5:end,1:800)),1));
plot(mdl3)
legend('location','best')
xlabel("log10(eigval)")
ylabel("mean abs amplitude")
ylim([2.5,7.5]);xlim([-12,-1])
figdir = "E:\OneDrive - Harvard University\GECCO2022\Figures\EigenSpaceAlignment";
saveallform(figdir,"meanabs_log_scatter_ctrl_cmp_invivo_with_regression",1)


%%
figdir = "E:\OneDrive - Harvard University\GECCO2022\Figures\EigenSpaceAlignment";
saveallform(figdir,"meanabs_log_scatter_ctrl_cmp_invivo",4)
%%
cutoff = 800;
[H,P,KSSTAT] = kstest2(reshape(projcoefs(:,1:cutoff),1,[]),reshape(projcoefs(:,cutoff+1:end),1,[]))
%%
[H,P,KSSTAT] = kstest2(reshape(drftgen_projcoefs(:,1:cutoff),1,[]),reshape(drftgen_projcoefs(:,cutoff+1:end),1,[]))
%%
[rval,pval] = corr(log10(eigvals(1:cutoff))', mean(abs(projcoefs(:,1:cutoff)),1)')
[rval,pval] = corr(log10(eigvals(1:cutoff))', mean(abs(drftgen_projcoefs(:,1:cutoff)),1)')
[rval,pval] = corr(log10(eigvals(1:cutoff))', mean(abs(drftnrm_projcoefs(:,1:cutoff)),1)')
%%

%%
test_corr_cutoff(eigvals,{initprojs,projcoefs,devicoefs,...
    drftnrm_projcoefs,drftnrm_devicoefs,...
    drftgen_projcoefs,drftgen_devicoefs},...
    ["init","evol","evol - mean init",...
    "drftevol match norm", "drftevol - mean init",...
    "drftevol match gen", "drftevol - mean init"],...
    cutoff,"meanabs")
test_corr_cutoff(eigvals,{initprojs,projcoefs,devicoefs,...
    drftnrm_projcoefs,drftnrm_devicoefs,...
    drftgen_projcoefs,drftgen_devicoefs},...
    ["init","evol","evol - mean init",...
    "drftevol match norm", "drftevol - mean init",...
    "drftevol match gen", "drftevol - mean init"],...
    cutoff,"absmean")
%%
test_corr_cutoff(eigvals,{initprojs,projcoefs,devicoefs,...
    drftnrm_projcoefs,drftnrm_devicoefs,...
    drftgen_projcoefs,drftgen_devicoefs},...
    ["init","evol","evol - mean init",...
    "drftevol match norm", "drftevol - mean init",...
    "drftevol match gen", "drftevol - mean init"],...
    cutoff,"meanabs",'log')
test_corr_cutoff(eigvals,{initprojs,projcoefs,devicoefs,...
    drftnrm_projcoefs,drftnrm_devicoefs,...
    drftgen_projcoefs,drftgen_devicoefs},...
    ["init","evol","evol - mean init",...
    "drftevol match norm", "drftevol - mean init",...
    "drftevol match gen", "drftevol - mean init"],...
    cutoff,"absmean",'log')
function test_corr_cutoff(eigvals,projcoef_col,label_col,cutoff,oper,eigop)
if nargin <= 5, eigop = "";end
if nargin <= 4, oper = "meanabs";end
if oper=="absmean"
    fprintf("abs of mean proj coef\n")
elseif oper=="meanabs"
    fprintf("mean of abs proj coef\n")
end
if eigop == "log"
    fprintf("do log transform to eigenvalues\n")
end
for mi = 1:numel(projcoef_col)
    label = label_col(mi);
    projcoef = projcoef_col{mi};
    if oper=="absmean"
        coefs = abs(mean(projcoef(:,1:cutoff),1));
    elseif oper=="meanabs"
        coefs = mean(abs(projcoef(:,1:cutoff)),1);
    end
    if eigop=="log"
        eigcoefs = log10(eigvals(1:cutoff));
    else
        eigcoefs = eigvals(1:cutoff);
    end
    [rval,pval] = corr(eigcoefs', coefs');
    fprintf("%s mean corr=%.3f P=%.1e (H eig cutoff=%d)\n",label,rval,pval,cutoff)
end
end

