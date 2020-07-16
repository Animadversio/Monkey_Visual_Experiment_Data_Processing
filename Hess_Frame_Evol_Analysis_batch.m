%% Hess_Frame_Evol_Analysis_batch
Animal="Alfa"; Set_Path;
evol_ids = find(contains(ExpRecord.expControlFN,"generate_integrated") & ExpRecord.Exp_collection=="ReducDimen_Evol") ;
%%
G = FC6Generator();
%% load the Hessian Matrices
py.importlib.import_module("numpy")
% data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\HessTune\NullSpace\Evolution_Avg_Hess.npz");
data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy\Nullspace\Pasu_Space_Avg_Hess.npz");
eigvect = data.get("eigvect_avg").double; % each column is an eigen vector. 
eigvals = data.get("eigv_avg").double;
eigvals = eigvals(end:-1:1);
eigvect = eigvect(:,end:-1:1);
%% load the evolution codes
for evoli = evol_ids(31:end)'
disp(ExpRecord.comments(evoli))
stimuli_path = ExpRecord.stimuli{evoli}; %ExpRecord.stimuli{241}; % Manif Expi 11
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1); % each row is an evolved code
% [codes_all_rd, img_ids_rd, code_geni_rd] = load_codes_all(stimuli_path, 2); % each row is an evolved code
%
Hproj_coef = codes_all * eigvect;
Hproj_coef_rd = codes_all_rd * eigvect;
%
geni = max(code_geni);
[~,P,~,STATS] = ttest2(abs(mean(Hproj_coef(code_geni==geni,1:100),1)), abs(mean(Hproj_coef(code_geni==geni,300:400),1)));
fprintf("1-100 vs 300-400, t=%.3f (p=%.3f)\n",STATS.tstat,P)
[~,P,~,STATS] = ttest2(abs(mean(Hproj_coef(code_geni==geni,1:100),1)), abs(mean(Hproj_coef(code_geni==geni,3996:4096),1)));
fprintf("1-100 vs 3996-4096, t=%.3f (p=%.3f)\n",STATS.tstat,P)
[~,P,~,STATS] = ttest2(abs(mean(Hproj_coef(code_geni==geni,300:400),1)), abs(mean(Hproj_coef(code_geni==geni,3996:4096),1)));
fprintf("300-400 vs 3996-4096, t=%.3f (p=%.3f)\n",STATS.tstat,P)
end