% Analyzing the Evolution Trajectory in the Frame of Hessian eigenvectors.
% 
% Get inspired by BigGAN version BigGAN_Hess_Frame_Evol_Analysis.
%% 
load("D:\Project_CMA_Monkeys.mat");
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
ctrl_evol_dir = "N:\Data-Ephys-MAT\Project_CMA_ghostEvolutions_matchingNorm";
fnlist = string(ls(ctrl_evol_dir));
%% Load the initial codes
Data = load("N:\Code\Generator_DB_Windows\init_population\init_code.mat");
initcodes = Data.codes;
%% Collect the evolved codes from Carlos' collection
evol_codes = [];
MonkName = ["Alfa","Beto"];
for Mi = 1:2
   Animal = MonkName(Mi);
   for Expi = 1:numel(MonkeysCMA{Mi}.Stats)
      nGen = max(MonkeysCMA{Mi}.Stats{Expi}.genes.gen);
      meancode = mean(MonkeysCMA{Mi}.Stats{Expi}.genes.all(nGen==MonkeysCMA{Mi}.Stats{Expi}.genes.gen,:),1);
      evol_codes = [evol_codes; meancode];
   end
end
%
evo_meta = repmat(struct(),size(evol_codes,1),1);
csr=1;
for Mi = 1:2
   Animal = MonkName(Mi);
   for Expi = 1:numel(MonkeysCMA{Mi}.Stats)
       evo_meta(csr).nGen = max(MonkeysCMA{Mi}.Stats{Expi}.genes.gen);
       evo_meta(csr).Animal = Animal;
       evo_meta(csr).Expi = Expi;
       evo_meta(csr).prefchan = MonkeysCMA{Mi}.Stats{Expi}.prefChan;
       csr = csr+1;
   end
end
%% Collect the controled drift evolution codes
drft_codes = [];
for Mi = 1:2
   Animal = MonkName(Mi);
   for Expi = 1:numel(MonkeysCMA{Mi}.Stats)
       chan = MonkeysCMA{Mi}.Stats{Expi}.prefChan;
       matfn = compose("%s-exp-%d-chan-%02d-001.mat",Animal,Expi,chan);
%        matpart = compose("%s-exp-%d-",Animal,Expi);
%        matfn = fnlist(find(contains(fnlist,matpart)));
%        fprintf("%s %s\n",matfn,matfn1)
       assert(exist(fullfile(ctrl_evol_dir,matfn),'File')>0)
       D = load(fullfile(ctrl_evol_dir, matfn));
       drft_codes = [drft_codes; D.genes(end,:)]; 
   end
end
%%
save(fullfile(mat_dir,"evo_drift_codes_all.mat"),'evol_codes','drft_codes','evo_meta')
%%

load(fullfile(mat_dir,"evol_ctrl_codes.mat"),'drft_codes_col')
%% Load up the Hessian for GAN
py.importlib.import_module("numpy")
% data = py.numpy.load("N:\Code\Hessian\Texture_Avg_Hess.npz");
data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\HessTune\NullSpace\Texture_Avg_Hess.npz");
% data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy\Nullspace\Pasu_Space_Avg_Hess.npz");
eigvect = data.get("eigvect_avg").double; % each column is an eigen vector. 
eigvals = data.get("eigv_avg").double;
eigvals = eigvals(end:-1:1);
eigvect = eigvect(:,end:-1:1);
%% simulated drift data. 
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir,"evo_drift_codes_all.mat"),'evol_codes','drft_codes','evo_meta')
load(fullfile(mat_dir,"evol_ctrl_codes.mat"),'drft_codes_col') % Control code by stop at the same mean norm
load(fullfile(mat_dir,"evol_ctrl_gen_codes.mat"),'drft_codes_gen_col') % Control code by stop at the same generation number
Data = load("N:\Code\Generator_DB_Windows\init_population\init_code.mat");
initcodes = Data.codes;
drft_codes_nrm = cell2mat(drft_codes_col); % stop at the same mean norm
drft_codes_gen = cell2mat(drft_codes_gen_col); % stop at the same generation. 
%%
evol_proj = evol_codes * eigvect;
drft_proj = drft_codes * eigvect;
init_proj = initcodes * eigvect;
drft_proj_gen = drft_codes_gen * eigvect; % stop at the same generation. 
drft_proj_nrm = drft_codes_nrm * eigvect; % stop at the same mean norm. 
%%
init_devi_shfl_codes = [];
init_mean = mean(initcodes);
for ci = 1:size(evol_codes, 1)
    devi_code = evol_codes(ci, :) - init_mean;
end

%% Deviation shuffling from initial code.
tic
cutoff = 800;
init_mean = mean(initcodes); % mean of initial generation
cc_idx_dist = []; cc_eva_dist = []; cc_logeva_dist = []; cc_spear_dist = [];
cc_idx_800_dist = []; cc_eva_800_dist = []; cc_logeva_800_dist = []; cc_spear_800_dist = [];
for triali = 1:1000
fprintf("%d ", triali)
evol_devi_codes = evol_codes - init_mean;
perm_evol_mean = []; % one round of shuffling
for i = 1:size(evol_codes,1)
perm_evol_mean = [perm_evol_mean; evol_devi_codes(i,randperm(4096))];
end
perm_evol_mean = perm_evol_mean + init_mean;
perm_evol_devi_proj = perm_evol_mean * eigvect;

cc_idx_shfl = corr([1:numel(eigvals)]', mean(abs(perm_evol_devi_proj(:,:)),1)');
cc_eva_shfl = corr(eigvals', mean(abs(perm_evol_devi_proj(:,:)),1)');
cc_logeva_shfl = corr(real(log10(eigvals))', mean(abs(perm_evol_devi_proj(:,:)),1)');
cc_spear_shfl = corr(eigvals', mean(abs(perm_evol_devi_proj(:,:)),1)','Type','Spearman');
cc_idx_dist = [cc_idx_dist, cc_idx_shfl];
cc_eva_dist = [cc_eva_dist, cc_eva_shfl];
cc_logeva_dist = [cc_logeva_dist, cc_logeva_shfl];
cc_spear_dist = [cc_spear_dist, cc_spear_shfl];

cc_idx_shfl800 = corr([1:cutoff]', mean(abs(perm_evol_devi_proj(:,1:cutoff)),1)');
cc_eva_shfl800 = corr(eigvals(1:cutoff)', mean(abs(perm_evol_devi_proj(:,1:cutoff)),1)');
cc_logeva_shfl800 = corr(real(log10(eigvals(1:cutoff)))', mean(abs(perm_evol_devi_proj(:,1:cutoff)),1)');
cc_spear_shfl800 = corr(eigvals(1:cutoff)', mean(abs(perm_evol_devi_proj(:,1:cutoff)),1)','Type','Spearman');
cc_idx_800_dist = [cc_idx_800_dist, cc_idx_shfl800];
cc_eva_800_dist = [cc_eva_800_dist, cc_eva_shfl800];
cc_logeva_800_dist = [cc_logeva_800_dist, cc_logeva_shfl800];
cc_spear_800_dist = [cc_spear_800_dist, cc_spear_shfl800];
end
toc
%%
figdir = "E:\OneDrive - Washington University in St. Louis\HessEvolStruct";
save(fullfile(figdir,"FC6_evol_devi_perm_corr.mat"), ...
  "cc_idx_dist", "cc_eva_dist", "cc_logeva_dist", "cc_spear_dist", ...
  "cc_idx_800_dist", "cc_eva_800_dist", "cc_logeva_800_dist", "cc_spear_800_dist")

%% Compute the correlations 
% Real Evolved codes 
[cc_idx, cc_eva, cc_logeva, cc_spear] = all_eigproj_corr(eigvals, evol_proj);
% Initial Generation codes
[cc_idx_init, cc_eva_init, cc_logeva_init, cc_spear_init] = all_eigproj_corr(eigvals, init_proj);
% Drift evolution, terminate at same gen
[cc_idx_drft_gen, cc_eva_drft_gen, cc_logeva_drft_gen, cc_spear_drft_gen] = all_eigproj_corr(eigvals, drft_proj_gen);
% Drift evolution, terminate at same Norm
[cc_idx_drft_nrm, cc_eva_drft_nrm, cc_logeva_drft_nrm, cc_spear_drft_nrm] = all_eigproj_corr(eigvals, drft_proj_nrm);
% Real Evolved codes, cutoff in spectrum
[cc_idx_800, cc_eva_800, cc_logeva_800, cc_spear_800] = all_eigproj_corr(eigvals(1:cutoff), evol_proj(:,1:cutoff));
% Initial Generation codes, cutoff in spectrum
[cc_idx_800_init, cc_eva_800_init, cc_logeva_800_init, cc_spear_800_init] = all_eigproj_corr(eigvals(1:cutoff), init_proj(:,1:cutoff));
% Drift evolution, terminate at same gen, cutoff in spectrum
[cc_idx_800_drft_gen, cc_eva_800_drft_gen, cc_logeva_800_drft_gen, cc_spear_800_drft_gen] = all_eigproj_corr(eigvals(1:cutoff), drft_proj_gen(:,1:cutoff));
% Drift evolution, terminate at same Norm, cutoff in spectrum
[cc_idx_800_drft_nrm, cc_eva_800_drft_nrm, cc_logeva_800_drft_nrm, cc_spear_800_drft_nrm] = all_eigproj_corr(eigvals(1:cutoff), drft_proj_nrm(:,1:cutoff));

%% Plot the Shuffling test of the real evolved correlation
figdir = "E:\OneDrive - Washington University in St. Louis\HessEvolStruct";
H = figure('position',[300,400, 1100, 400]);
T = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
nexttile(1)
histogram(cc_idx_dist, 30)
line([cc_idx, cc_idx], ylim(), 'Color', 'red')
line([cc_idx_init, cc_idx_init], ylim(), 'Color', 'k')
line([cc_idx_drft_gen, cc_idx_drft_gen], ylim(), 'Color', 'magenta')
line([cc_idx_drft_nrm, cc_idx_drft_nrm], ylim(), 'Color', 'green')
legend(["Permutation Control","Real Value","Initial Gen","Drift same gen#","Drift same norm"])
xlabel("Corr with eigen rank")
ylabel("Density")
title(compose("cc=%.3f (p=%.4f)", cc_idx, 1-sum(cc_idx_dist>cc_idx)/numel(cc_idx_dist)))
box off;
nexttile(2)
histogram(cc_eva_dist, 30)
line([cc_eva, cc_eva], ylim(), 'Color', 'red')
line([cc_eva_init, cc_eva_init], ylim(), 'Color', 'k')
line([cc_eva_drft_gen, cc_eva_drft_gen], ylim(), 'Color', 'magenta')
line([cc_eva_drft_nrm, cc_eva_drft_nrm], ylim(), 'Color', 'green')
xlabel("Corr with eigval")
title(compose("cc=%.3f (p=%.4f)", cc_eva, sum(cc_eva_dist>cc_eva)/numel(cc_eva_dist)))
box off;
nexttile(3)
histogram(cc_logeva_dist, 30)
line([cc_logeva, cc_logeva], ylim(), 'Color', 'red')
line([cc_logeva_init, cc_logeva_init], ylim(), 'Color', 'k')
line([cc_logeva_drft_gen, cc_logeva_drft_gen], ylim(), 'Color', 'magenta')
line([cc_logeva_drft_nrm, cc_logeva_drft_nrm], ylim(), 'Color', 'green')
xlabel("Corr with Log(eigval)")
title(compose("cc=%.3f (p=%.4f)", cc_logeva, sum(cc_logeva_dist>cc_logeva)/numel(cc_logeva_dist)))
box off;
nexttile(4)
histogram(cc_spear_dist, 30)
line([cc_spear, cc_spear], ylim(), 'Color', 'red')
line([cc_spear_init, cc_spear_init], ylim(), 'Color', 'k')
line([cc_spear_drft_gen, cc_spear_drft_gen], ylim(), 'Color', 'magenta')
line([cc_spear_drft_nrm, cc_spear_drft_nrm], ylim(), 'Color', 'green')
xlabel("Spearman Corr")
title(compose("cc=%.3f (p=%.4f)", cc_spear, sum(cc_spear_dist>cc_spear)/numel(cc_spear_dist)))
box off;
title(T,["Correlation of Projection Amplitude of Evolved Code and Hessian Eigen Structures (FC6, N=264)",...
    "Shuffle the Deviation vector of evolved code from initial code"],'FontSize',16)
savefig(H, fullfile(figdir,"FC6_proj_eig_corr_devi_permtest.fig"))
saveas(H, fullfile(figdir,"FC6_proj_eig_corr_devi_permtest.png"))
saveas(H, fullfile(figdir,"FC6_proj_eig_corr_devi_permtest.pdf"))
%% Plot the Shuffling test of the real evolved correlation,  with cutoff
figdir = "E:\OneDrive - Washington University in St. Louis\HessEvolStruct";
H2 = figure('position',[300,400, 1100, 400]);
T = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
nexttile(1)
histogram(cc_idx_800_dist, 30)
line([cc_idx_800, cc_idx_800], ylim(), 'Color', 'red')
line([cc_idx_800_init, cc_idx_800_init], ylim(), 'Color', 'k')
line([cc_idx_800_drft_gen, cc_idx_800_drft_gen], ylim(), 'Color', 'magenta')
line([cc_idx_800_drft_nrm, cc_idx_800_drft_nrm], ylim(), 'Color', 'green')
legend(["Permutation Control","Real Value","Initial Gen","Drift same gen#","Drift same norm"])
xlabel("Corr with eigen rank")
ylabel("Density")
title(compose("cc=%.3f (p=%.4f)", cc_idx_800, 1-sum(cc_idx_800_dist>cc_idx_800)/numel(cc_idx_800_dist)))
box off;
nexttile(2)
histogram(cc_eva_800_dist, 30)
line([cc_eva_800, cc_eva_800], ylim(), 'Color', 'red')
line([cc_eva_800_init, cc_eva_800_init], ylim(), 'Color', 'k')
line([cc_eva_800_drft_gen, cc_eva_800_drft_gen], ylim(), 'Color', 'magenta')
line([cc_eva_800_drft_nrm, cc_eva_800_drft_nrm], ylim(), 'Color', 'green')
xlabel("Corr with eigval")
title(compose("cc=%.3f (p=%.4f)", cc_eva_800, sum(cc_eva_800_dist>cc_eva_800)/numel(cc_eva_800_dist)))
box off;
nexttile(3)
histogram(cc_logeva_800_dist, 30)
line(real([cc_logeva_800, cc_logeva_800]), ylim(), 'Color', 'red')
line(real([cc_logeva_800_init, cc_logeva_800_init]), ylim(), 'Color', 'k')
line([cc_logeva_800_drft_gen, cc_logeva_800_drft_gen], ylim(), 'Color', 'magenta')
line([cc_logeva_800_drft_nrm, cc_logeva_800_drft_nrm], ylim(), 'Color', 'green')
xlabel("Corr with Log(eigval)")
title(compose("cc=%.3f (p=%.4f)", cc_logeva_800, sum(cc_logeva_800_dist>cc_logeva_800)/numel(cc_logeva_800_dist)))
box off;
nexttile(4)
histogram(cc_spear_800_dist, 30)
line([cc_spear_800, cc_spear_800], ylim(), 'Color', 'red')
line([cc_spear_800_init, cc_spear_800_init], ylim(), 'Color', 'k')
line([cc_spear_800_drft_gen, cc_spear_800_drft_gen], ylim(), 'Color', 'magenta')
line([cc_spear_800_drft_nrm, cc_spear_800_drft_nrm], ylim(), 'Color', 'green')
xlabel("Spearman Corr")
title(compose("cc=%.3f (p=%.4f)", cc_spear_800, sum(cc_spear_800_dist>cc_spear_800)/numel(cc_spear_800_dist)))
box off;
title(T,["Correlation of Projection Amplitude of Evolved Code and Hessian Eigen Structures (Top 800) (FC6, N=264)",...
    "Shuffle the Deviation vector of evolved code from initial code"],'FontSize',16)
savefig(H2, fullfile(figdir,"FC6_proj_eig800_corr_devi_permtest.fig"))
saveas(H2, fullfile(figdir,"FC6_proj_eig800_corr_devi_permtest.png"))
saveas(H2, fullfile(figdir,"FC6_proj_eig800_corr_devi_permtest.pdf"))
%%
H4=figure;clf;hold on;set(H4,'pos',[1000,612,420,380])
sampN = size(evol_proj,1);
scatter(eigvals,mean(abs(evol_proj)),10,'k')
plotRegressLine(eigvals,mean(abs(evol_proj)),'k')
scatter(eigvals,mean(abs(init_proj)),10,'r')
% plotRegressLine(eigvals,mean(abs(init_proj)),'r')
scatter(eigvals,mean(abs(drft_proj_gen)),10,'m')
plotRegressLine(eigvals,mean(abs(drft_proj_gen)),'m')
scatter(eigvals,mean(abs(drft_proj_nrm)),10,'g')
plotRegressLine(eigvals,mean(abs(drft_proj_nrm)),'g')
xlabel("eigen value");ylabel("Projection Amplitude");
title(["Correlation of Projection Amplitude and Eigenvalue","(FC6 N="+num2str(sampN)+")"])
legend(["Evolved","Initial Gen","Drift same gen#","Drift same norm"])
savefig(H4, fullfile(figdir,"FC6_proj_eig_scatter.fig"))
saveas(H4, fullfile(figdir,"FC6_proj_eig_scatter.png"))
saveas(H4, fullfile(figdir,"FC6_proj_eig_scatter.pdf"))
%%
H5=figure;clf;hold on;set(H5,'pos',[1000,612,420,380])
scatter(eigvals(1:cutoff),mean(abs(evol_proj(:,1:cutoff))),10,'k')
plotRegressLine(eigvals(1:cutoff),mean(abs(evol_proj(:,1:cutoff))),'k')
scatter(eigvals(1:cutoff),mean(abs(init_proj(:,1:cutoff))),10,'r')
% plotRegressLine(eigvals(1:cutoff),mean(abs(init_proj(:,1:cutoff))),'r')
scatter(eigvals(1:cutoff),mean(abs(drft_proj_gen(:,1:cutoff))),10,'m')
plotRegressLine(eigvals(1:cutoff),mean(abs(drft_proj_gen(:,1:cutoff))),'m')
scatter(eigvals(1:cutoff),mean(abs(drft_proj_nrm(:,1:cutoff))),10,'g')
plotRegressLine(eigvals(1:cutoff),mean(abs(drft_proj_nrm(:,1:cutoff))),'g')
xlabel("eigen value");ylabel("Projection Amplitude");
title(["Correlation of Projection Amplitude and Eigenvalue","(Cutoff at 800) (FC6 N="+num2str(sampN)+")"])
legend(["Evolved","Initial Gen","Drift same gen#","Drift same norm"])
savefig(H5, fullfile(figdir,"FC6_proj_eig800_scatter.fig"))
saveas(H5, fullfile(figdir,"FC6_proj_eig800_scatter.png"))
saveas(H5, fullfile(figdir,"FC6_proj_eig800_scatter.pdf"))
%%
H6=figure;clf;hold on
scatter(real(log10(eigvals)),mean(abs(drft_proj_gen)),10,'g')
plotRegressLine(real(log10(eigvals)),mean(abs(drft_proj_gen)),'g')
scatter(real(log10(eigvals)),mean(abs(drft_proj_nrm)),10,'magenta')
plotRegressLine(real(log10(eigvals)),mean(abs(drft_proj_nrm)),'magenta')
scatter(real(log10(eigvals)),mean(abs(init_proj)),10,'r')
plotRegressLine(real(log10(eigvals)),mean(abs(init_proj)),'r')
scatter(real(log10(eigvals)),mean(abs(evol_proj)),10,'k')
plotRegressLine(real(log10(eigvals)),mean(abs(evol_proj)),'k')


%%
figure;clf;hold on
scatter(real(eigvals),mean(abs(drft_proj_gen)),10,'g')
scatter(real(eigvals),mean(abs(drft_proj_nrm)),10,'magenta')
scatter(real(eigvals),mean(abs(init_proj)),10,'r')
scatter(real(eigvals),mean(abs(evol_proj)),10,'k')
%%
% [ctrl_rksm_P, ~, ctrl_rksm_Stat] = arrayfun(@(id)ranksum(code_proj_cc(:,id),ctrl_proj_cc(:,id)),1:4096,'uni',1);
% [~, ctrl_krsm_P, ctrl_krsm_Stat] = arrayfun(@(id)kstest2(code_proj_cc(:,id),ctrl_proj_cc(:,id),'alpha',0.01),1:4096,'uni',1);
[~, ctrl_ttst_P, ~, ctrl_ttst_Stat] = arrayfun(@(id)ttest(evol_proj(:,id),drft_proj(:,id)),1:4096,'uni',0);
ctrl_ttst_Stat = cell2mat(ctrl_ttst_Stat); ctrl_ttst_P = cell2mat(ctrl_ttst_P); 
%%
[~, ctrl_ttst_abs_P, ~, ctrl_ttst_abs_Stat] = arrayfun(@(id)ttest(abs(evol_proj(:,id)),abs(drft_proj(:,id))),1:4096,'uni',0);
ctrl_ttst_abs_Stat = cell2mat(ctrl_ttst_abs_Stat); ctrl_ttst_abs_P = cell2mat(ctrl_ttst_abs_P); 
%%
[ctrl_sgrk_P, ~, ctrl_sgrk_Stat] = arrayfun(@(id)signrank(evol_proj(:,id),drft_proj(:,id)),1:4096,'uni',1);
[ctrl_sgrk_abs_P, ~, ctrl_sgrk_abs_Stat] = arrayfun(@(id)signrank(abs(evol_proj(:,id)),abs(drft_proj(:,id)),'tail','both'),1:4096,'uni',1);
%%
corr_P = mafdr(ctrl_sgrk_abs_P,'BHFDR',true); % nothing fall out of it
find(corr_P<0.1)
%% To match the number. 
rep_evol_proj = repelem(evol_proj,5,1);
[ctrl2_sgrk_abs_P, ~, ctrl2_sgrk_abs_Stat] = arrayfun(@(id)signrank(rep_evol_proj(:,id),drft_proj_gen(:,id),'tail','both'),1:4096,'uni',1);
corr_P = mafdr(ctrl2_sgrk_abs_P,'BHFDR',true); % nothing fall out of it
find(corr_P<0.001)
%%
[~, ctrl_ttst_abs_P, ~, ctrl_ttst_abs_Stat] = arrayfun(@(id)ttest(rep_evol_proj(:,id),drft_proj_nrm(:,id)),1:4096,'uni',0);
ctrl_ttst_abs_Stat = cell2mat(ctrl_ttst_abs_Stat); ctrl_ttst_abs_P = cell2mat(ctrl_ttst_abs_P); 
corr_P = mafdr(ctrl_ttst_abs_P,'BHFDR',true); % nothing fall out of it
find(corr_P<0.00001)
%%
figure;plot(arrayfun(@(S)S.tstat,ctrl_ttst_abs_Stat))
%%
thresh_covobs = covobs;
thresh_covobs(covobs<5)=nan;
figure;
imagesc(thresh_covobs)
%%
ccobs = corrcoef(obs);
ccobs_nodiag = ccobs + diag(nan(1,4096));
thresh_ccobs = ccobs_nodiag;
thresh_ccobs(ccobs_nodiag<0.1)=nan;
%%
figure;plot(sort(ccobs_nodiag(:)))
%%
figure;imagesc(thresh_ccobs);axis image

function plotRegressLine(X, Y, varargin) % Util function to add a linear reg line to a scatter
if nargin==2, varargin={};end
[r,m,b] = regression(X, Y);
xmin = min(X); xmax = max(X);
p = plot([xmin,xmax],[xmin,xmax].*m+b,varargin{:});
set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
function [cc_idx, cc_eva, cc_logeva, cc_spear] = all_eigproj_corr(eigvals, evol_proj)
cc_idx = corr([1:numel(eigvals)]', mean(abs(evol_proj(:,:)),1)');
cc_eva = corr(eigvals', mean(abs(evol_proj(:,:)),1)');
cc_logeva = corr(log10(eigvals)', mean(abs(evol_proj(:,:)),1)');
cc_spear = corr(eigvals', mean(abs(evol_proj(:,:)),1)','Type','Spearman');
end
