%% Hessian Frame Evolution Analysis
%  Project the evolution trajectory onto the Hessian frames 
Animal="Beto"; Set_Path;
find(contains(ExpRecord.expControlFN,"generate_integrated") & ExpRecord.Exp_collection=="ReducDimen_Evol") 
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
stimuli_path = ExpRecord.stimuli{182}; %ExpRecord.stimuli{241}; % Manif Expi 11
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1); % each row is an evolved code
[codes_all_rd, img_ids_rd, code_geni_rd] = load_codes_all(stimuli_path, 2); % each row is an evolved code
%% Project onto Hessian eigen frame
Hproj_coef = codes_all * eigvect;
Hproj_coef_rd = codes_all_rd * eigvect;
%% Within generation std is kind of same across spectrum.
figure(3);
for geni = min(code_geni):max(code_geni)
std_vect = std(Hproj_coef(code_geni==geni,:),1,1);
plot(std_vect);ylim([1.5,6])
title(compose("%d std %.3f stdstd %.3f (across dim)", geni, mean(std_vect), std(std_vect)))
pause(0.3)
end
%%
figure(1);clf;
for geni = min(code_geni):max(code_geni)
% scatter(repmat(1:4096,sum(code_geni==geni),1)', Hproj_coef(code_geni==geni,:)')
codes_gen = mean(codes_all(code_geni==geni,:),1);
codes_gen_rd = mean(codes_all_rd(code_geni_rd==geni,:),1);
bsize = sum(code_geni==geni);
subtightplot(1,3,1)
cla; hold on 
for ci = find(code_geni==geni) % plot each codes separately 
scatter(1:4096, Hproj_coef(ci,:), 10)
end
scatter(1:4096, mean(Hproj_coef(code_geni==geni,:),1), 45,'k','filled')
xlim([1,100]);ylim([-30,35])

subtightplot(1,3,2)
cla; hold on 
for ci = find(code_geni==geni) % plot each codes separately 
scatter(1:4096, Hproj_coef(ci,:), 10)
end
scatter(1:4096, mean(Hproj_coef(code_geni==geni,:),1), 45,'k','filled')
xlim([300,400]);ylim([-30,35])
subtightplot(1,3,3)
cla; hold on 
for ci = find(code_geni==geni) % plot each codes separately 
scatter(1:4096, Hproj_coef(ci,:), 10)
end
scatter(1:4096, mean(Hproj_coef(code_geni==geni,:),1), 45,'k','filled')
xlim([4000,4096]);ylim([-30,35])
suptitle(num2str(geni))
[~,P,~,STATS] = ttest2(abs(mean(Hproj_coef(code_geni==geni,1:100),1)), abs(mean(Hproj_coef(code_geni==geni,300:400),1)));
fprintf("1-100 vs 300-400, t=%.3f (p=%.3f)\n",STATS.tstat,P)
[~,P,~,STATS] = ttest2(abs(mean(Hproj_coef(code_geni==geni,1:100),1)), abs(mean(Hproj_coef(code_geni==geni,3996:4096),1)));
fprintf("1-100 vs 3996-4096, t=%.3f (p=%.3f)\n",STATS.tstat,P)
[~,P,~,STATS] = ttest2(abs(mean(Hproj_coef(code_geni==geni,300:400),1)), abs(mean(Hproj_coef(code_geni==geni,3996:4096),1)));
fprintf("300-400 vs 3996-4096, t=%.3f (p=%.3f)\n",STATS.tstat,P)
pause
end
%%
figure(2);clf;
for geni = min(code_geni):max(code_geni)
% scatter(repmat(1:4096,sum(code_geni==geni),1)', Hproj_coef(code_geni==geni,:)')
codes_gen = mean(codes_all(code_geni==geni,:),1);
codes_gen_rd = mean(codes_all_rd(code_geni_rd==geni,:),1);
bsize = sum(code_geni==geni);
subtightplot(1,3,1)
cla; hold on 
for ci = find(code_geni==geni) % plot each codes separately 
scatter(1:4096, codes_all(ci,:), 10)
end
scatter(1:4096, mean(codes_all(code_geni==geni,:),1), 45,'k','filled')
xlim([1,100]);ylim([-30,35])

subtightplot(1,3,2)
cla; hold on 
for ci = find(code_geni==geni) % plot each codes separately 
scatter(1:4096, codes_all(ci,:), 10)
end
scatter(1:4096, mean(codes_all(code_geni==geni,:),1), 45,'k','filled')
xlim([200,300]);ylim([-30,35])
subtightplot(1,3,3)
cla; hold on 
for ci = find(code_geni==geni) % plot each codes separately 
scatter(1:4096, codes_all(ci,:), 10)
end
scatter(1:4096, mean(codes_all(code_geni==geni,:),1), 45,'k','filled')
xlim([4000,4096]);ylim([-30,35])
suptitle(num2str(geni))
[~,P,~,STATS] = ttest2(abs(mean(codes_all(code_geni==geni,1:100),1)), abs(mean(codes_all(code_geni==geni,300:400),1)))
pause
end
%%
figure(14);clf;
for geni = min(code_geni):max(code_geni)
% scatter(repmat(1:4096,sum(code_geni==geni),1)', Hproj_coef(code_geni==geni,:)')
codes_gen = codes_all(code_geni==geni,:);
codes_gen_rd = codes_all_rd(code_geni_rd==geni,:);
bsize = sum(code_geni==geni);
bins = linspace(0,20,41);
subtightplot(1,1,1);cla;hold on 
histogram(abs(codes_gen(1,1:50)),bins,'Normalization','probability')
% histogram(codes_gen(1,51:200),40,'Normalization','probability')
histogram(abs(codes_gen(1,4046:4096)),bins,'Normalization','probability')
% xlim([-20,20])
% xlim([4000,4096]);ylim([-20,25])
title(num2str(geni))
pause
end
%%
figure;clf;
for geni = min(code_geni):max(code_geni)-1
cla; hold on 
% scatter(repmat(1:4096,sum(code_geni==geni),1)', Hproj_coef(code_geni==geni,:)')
codes_gen = mean(codes_all(code_geni==geni,:),1);
% codes_gen_rd = mean(codes_all_rd(code_geni_rd==geni,:),1);
scatter(1:4096, codes_gen * eigvect / norm(codes_gen),10)
% scatter(1:4096, codes_gen_rd * eigvect / norm(codes_gen_rd),10)
% scatter(1:4096, codes_gen(1,:) * eigvect / norm(codes_gen), 10)
% scatter(1:4096, mean(Hproj_coef(code_geni==geni,:),1),10)
% scatter(1:4096, codes_gen(1,randperm(4096)) * eigvect / norm(codes_gen), 10);
% ylim([-25,30])
title(num2str(geni))
pause
end
%%
figure(23);clf;
cutoff = 4000;
for geni = min(code_geni):max(code_geni)
cla;hold on 
% scatter(repmat(1:4096,sum(code_geni==geni),1)', Hproj_coef(code_geni==geni,:)')
codes_gen = mean(codes_all(code_geni==geni,:),1);
proj_code = codes_gen(1,:) * eigvect(:,1:cutoff);
scatter(1:cutoff, proj_code, 20,'filled')
for ci = 1:size(proj_code,1)%find(code_geni==geni) % plot each codes separately 
scatter(1:cutoff, proj_code(ci,:), 10)
end
% scatter(1:4096, mean(Hproj_coef(code_geni==geni,:),1),10)
proj_code_shf = codes_gen(1,randperm(4096)) * eigvect(:,1:cutoff);
% scatter(1:cutoff, proj_code_shf / norm(proj_code_shf), 10);
% ylim([-25,30])
xlim([0,4096])
title(num2str(geni))
pause
end
%% Just correlate code generation with its coef you will get nothing. 
cc_eigdir = corr(code_geni', Hproj_coef);
figure(2);clf;hold on
scatter(1:4096,cc_eigdir)
plot(movmean(cc_eigdir,10))
%% 
figure;
imagesc(Hproj_coef)
%% Project out later dimensions and see the projected codes
%% Many dimensions doesn't matter. 
cutoff_arr = [1,10,50,100,200,400,800,1000,1500,2000,3000,4096];
cutoff_codes = [];
for cutoff = cutoff_arr
cutoff_codes = [cutoff_codes; Hproj_coef(end, 1:cutoff) * eigvect(:,1:cutoff)'];
end
% cutoff_codes = renormalize(cutoff_codes, 300);
cutoff_imgs = G.visualize(cutoff_codes);
figure;montage(cutoff_imgs,'Size',[3,4])
%%
bandwidth = 200;
cutoff_arr = [1,10,20,50,75,100,200,400,800,1000,1500,2000,3000,3600];
cutoff_codes = [];
for cutoff = cutoff_arr
cutoff_codes = [cutoff_codes; Hproj_coef(end, cutoff:cutoff+bandwidth) * eigvect(:,cutoff:cutoff+bandwidth)'];
end
% cutoff_codes = renormalize(cutoff_codes, 300);
cutoff_imgs = G.visualize(cutoff_codes);
figure;montage(cutoff_imgs)
%%
G = FC6Generator();