% Analyzing the Evolution Trajectory in the Frame of Hessian eigenvectors.
% 
%% 
load("D:\Project_CMA_Monkeys.mat");
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
ctrl_evol_dir = "N:\Data-Ephys-MAT\Project_CMA_ghostEvolutions_matchingNorm";
fnlist = string(ls(ctrl_evol_dir));
%% Load the initial codes
Data = load("N:\Code\Generator_DB_Windows\init_population\init_code.mat");
initcodes = Data.codes;
%% Collect the evolved codes
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
%%
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
load(fullfile(mat_dir,"evol_ctrl_codes.mat"),'drft_codes_col')
load(fullfile(mat_dir,"evol_ctrl_gen_codes.mat"),'drft_codes_gen_col')
drft_codes_nrm = cell2mat(drft_codes_col); % stop at the same mean norm
drft_codes_gen = cell2mat(drft_codes_gen_col); % stop at the same generation. 
%%
evol_proj = evol_codes * eigvect;
drft_proj = drft_codes * eigvect;
init_proj = initcodes * eigvect;
drft_proj_gen = drft_codes_gen * eigvect;
drft_proj_nrm = drft_codes_nrm * eigvect;
%%
figure;clf;hold on
scatter(1:4096,mean(drft_proj_gen),10,'g')
scatter(1:4096,mean(drft_proj_nrm),10,'magenta')
scatter(1:4096,mean(init_proj),10,'r')
scatter(1:4096,mean(evol_proj),10,'k')
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