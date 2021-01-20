%% Evol OU process
%  this code tries to find the interesting genes that grows in evolution.
Animal="Beto"; Set_Path;
ftrrows = find(contains(ExpRecord.expControlFN,"generate_") & ExpRecord.Exp_collection=="Manifold");
%%
py.importlib.import_module("numpy")
data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\HessTune\NullSpace\Texture_Avg_Hess.npz");
% data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy\Nullspace\Pasu_Space_Avg_Hess.npz");
eigvect = data.get("eigvect_avg").double; % each column is an eigen vector. 
eigvals = data.get("eigv_avg").double;
eigvals = eigvals(end:-1:1);
eigvect = eigvect(:,end:-1:1);
%%
stimuli_path = ExpRecord.stimuli{22}; %ExpRecord.stimuli{241}; % Manif Expi 11
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1); % each row is an evolved code
%%
avg_codes = arrayfun(@(geni)mean(codes_all(code_geni==geni,:)), 1:max(code_geni),'Uni',0);
avg_codes = cell2mat(avg_codes');
sigma_arr = arrayfun(@(geni)mean(std(codes_all(code_geni==geni,:),1,1)), 1:max(code_geni),'Uni',1);
%%
avg_proj = avg_codes * eigvect;
%%
steps = avg_codes(2:end,:)-avg_codes(1:end-1,:);
devia = avg_codes(1:end,:)-avg_codes(1,:);
%%
steps_proj = avg_proj(2:end,:)-avg_proj(1:end-1,:);%%
devia_proj = avg_proj(:,:)-avg_proj(1,:);
%%
[sgrk_P,~,sgrk_STATS] = arrayfun(@(eigi)signrank(steps_proj(:,eigi)),1:4096); % not significant, flat p distribution
%%
[~,tt_P,~,tt_STATS] = arrayfun(@(eigi)ttest(steps_proj(:,eigi)),1:4096,'Uni',0); % not significant, flat p distribution
tt_P = cell2mat(tt_P);tt_STATS=cell2mat(tt_STATS);
%%
[~,tt_P_o,~,tt_STATS_o] = arrayfun(@(eigi)ttest(steps(:,eigi)),1:4096,'Uni',0); % not significant, flat p distribution
tt_P_o = cell2mat(tt_P_o);tt_STATS_o=cell2mat(tt_STATS_o);
%%
figure;histogram(tt_P,50)
%%
G = FC6Generator();
%%
figure;imshow(G.visualize(avg_codes(end,:)))
%%
figure;imshow(G.visualize(50*mean(steps_proj,1)*eigvect'))
%%
sparse_proj = avg_proj(end,:);
sparse_proj(sgrk_P>0.1)= 0;
figure;imshow(G.visualize(sparse_proj*eigvect'))
%%
[reg_slp,slp_int,~,~,regstats] = arrayfun(@(eigi)regress(devia_proj(:,eigi), [1:size(devia_proj,1)]'),1:4096,'uni',0);
slp_int = cell2mat(slp_int'); regstats = cell2mat(regstats');reg_slp = cell2mat(reg_slp);
% Note the slope has a strong stron correlation with the projection values
% in last generation.
%%
fprintf("The t stat of the steps is super correlated with final Projection values (corr%.3f)\n",corr(arrayfun(@(T)T.tstat,tt_STATS)',avg_proj(end,:)'))
fprintf("Same for non-parametric: The z value of `signrank` of the steps is super correlated with final Projection values (corr%.3f)\n",corr(arrayfun(@(T)T.signedrank,sgrk_STATS)',avg_proj(end,:)'))
fprintf("This Correlation with `t` is still high without the Hessian eigenframe projection (corr%.3f)\n", corr(arrayfun(@(T)T.tstat,tt_STATS_o)',avg_codes(end,:)'))
fprintf("Correlation between slope and final gen projection is high (corr%.3f)\n",corr(reg_slp',avg_proj(end,:)'))
fprintf("Correlation between Rsquare and absolute value of final gen projection (corr%.3f)\n",corr(regstats(:,1),abs(avg_proj(end,:))'))
%% Create statistical null distribution from the data generating process
Data = load("N:\Code\Generator_DB_Windows\init_population\init_code.mat");
initcodes = Data.codes;
%%
ctrlcodes = [];
tic
for triali = 1:1000
codes = initcodes;
for geni = 2:max(code_geni)
% ctrl_codes_all = [ctrl_codes_all; codes];
scores = randn(1,size(codes,1));
newcodes = rankweight(size(codes,1), 20)*codes + sigma_arr(geni) * randn(40, 4096);
codes = newcodes;
end
ctrlcodes = [ctrlcodes; mean(newcodes,1)];
if mod(triali,10)==0, fprintf("%d/1000\n",triali);end
end
ctrlproj = ctrlcodes * eigvect;
toc;% 161 sec for 1000 iterations
%%
% ctrlcodes_col = [ctrlcodes_col; mean(newcodes,1)];
% end
% save(fullfile(summarydir,"ctrl_code_col.mat"),'ctrlcodes_col')
empir_P = arrayfun(@(eigi)min(sum(avg_proj(end,eigi)<ctrlproj(:,eigi)), sum(avg_proj(end,eigi)>ctrlproj(:,eigi)))/500,1:4096);
figure;plot(sort(empir_P))
mafdr(empir_P,'BHFDR',1)
