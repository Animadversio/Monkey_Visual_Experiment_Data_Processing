%% Hessian Frame Evolution Analysis
%  Project the evolution trajectory onto the Hessian frames 
Animal="Beto"; Set_Path;
ftrrows = find(contains(ExpRecord.expControlFN,"generate_") & ExpRecord.Exp_collection=="Manifold") 
%%
G = FC6Generator();
%% load the Hessian Matrices
py.importlib.import_module("numpy")
data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\HessTune\NullSpace\Evolution_Avg_Hess.npz");
% data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy\Nullspace\Pasu_Space_Avg_Hess.npz");
eigvect = data.get("eigvect_avg").double; % each column is an eigen vector. 
eigvals = data.get("eigv_avg").double;
eigvals = eigvals(end:-1:1);
eigvect = eigvect(:,end:-1:1);
%% load the evolution codes
stimuli_path = ExpRecord.stimuli{10}; %ExpRecord.stimuli{241}; % Manif Expi 11
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1); % each row is an evolved code
[codes_all_rd, img_ids_rd, code_geni_rd] = load_codes_all(stimuli_path, 2); % each row is an evolved code
%% Project codes onto Hessian eigen frame
Hproj_coef = codes_all * eigvect;
% Hproj_coef_rd = codes_all_rd * eigvect;
%% Within generation std is kind of same across spectrum.
figure(3);
for geni = min(code_geni):max(code_geni)
std_vect = std(Hproj_coef(code_geni==geni,:),1,1);
plot(std_vect);ylim([1.5,6])
title(compose("%d std %.3f stdstd %.3f (across dim)", geni, mean(std_vect), std(std_vect)))
pause(0.3)
end
%% Well made Demo of code evolution in the Hessian Eigen Frame
figure(1);clf;
for geni = min(code_geni):max(code_geni)
% scatter(repmat(1:4096,sum(code_geni==geni),1)', Hproj_coef(code_geni==geni,:)')
codes_gen = mean(codes_all(code_geni==geni,:),1);
% codes_gen_rd = mean(codes_all_rd(code_geni_rd==geni,:),1);
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
%% Find a proper test for the notion we want 
% ttest of the amplitude is quite hard.
%%
Xmat = [code_geni;ones(size(code_geni))];
coefmat = inv(Xmat*Xmat')*Xmat*Hproj_coef;
residue = (Hproj_coef - Xmat'*coefmat);
varres = mean(residue.^2,1);
vartot = var(Hproj_coef,0,1);

sum(varres./vartot>0.95) % 1108 channels
mask = (varres ./ vartot);
%%
figure(2);clf;
for geni = min(code_geni):max(code_geni)
% scatter(repmat(1:4096,sum(code_geni==geni),1)', Hproj_coef(code_geni==geni,:)')
codes_gen = mean(codes_all(code_geni==geni,:),1);
% codes_gen_rd = mean(codes_all_rd(code_geni_rd==geni,:),1);
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
cutoff_codes = [cutoff_codes; Hproj_coef(end-2, 1:cutoff) * eigvect(:,1:cutoff)'];
end
% cutoff_codes = renormalize(cutoff_codes, 300);
cutoff_imgs = G.visualize(cutoff_codes);
figure(16);montage(cutoff_imgs,'Size',[3,4])
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
threshs = [0, 0.1,1,2,4,6,8,10,13];
cutoff_codes = [];
eigid = 1:4096;
for cutoff = threshs
cutoff_code = Hproj_coef(end, :);
cutoff_code((cutoff_code<cutoff & cutoff_code>-cutoff) | eigid>400) = 0;
cutoff_codes = [cutoff_codes; cutoff_code];
fprintf("%d ",sum(abs(cutoff_code)>0))
end
% cutoff_codes = renormalize(cutoff_codes, 300);
cutoff_imgs = G.visualize(cutoff_codes * eigvect');
figure;montage(cutoff_imgs)
%%
fincodes = mean(Hproj_coef(code_geni==geni, :),1);
cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((fincodes>cutoff | fincodes<-cutoff) & eigid<=eigthr,2)),...
            prctile(fincodes,[90,95,99,99.5,99.9,99.95,99.99]))./eigthr,...[0,0.2,1,2,4,6,8,10,13,15,18]
        [100,200,400,800,1000,2000,4096]','Uni',0))
%%
ngen = max(code_geni);
fincodes = (Hproj_coef(code_geni==ngen, :));%,1);
cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((fincodes>cutoff | fincodes<-cutoff) & eigid>=eigthr-200+1 & eigid<=eigthr,2)),...
            prctile(fincodes(:),[80,90,95,99,99.5,99.9,99.95,99.99],'all')')./200,...[0,0.2,1,2,4,6,8,10,13,15,18]
        [100,200,400,800,1000,2000,4096]','Uni',0))
%%
figure;hold on
plot([300 295 266 228 172 110 69 41 16]./300)
plot([400 394 353 303 229 151 95 54 24 ]./400)
plot([600 591 530 455 336 225 142 79 28 ]./600);
plot([4096 4038 3581 3075 2134 1401 826 469 143]./4096);
%% Threshold and number of axes surpassing threshold in different part of spectra
cutoff = 13;
cutoff_code = fincodes;
cutoff_code((cutoff_code<cutoff & cutoff_code>-cutoff) ) = 0; %| eigid>400
figure;stem(1:4096,cutoff_code,'Marker','none')

%%
G = FC6Generator();
%%
Animal="Beto"; Set_Path;
ftrrows = find(contains(ExpRecord.expControlFN,"generate_") & ExpRecord.Exp_collection=="Manifold");
%%
% 
stimuli_path = ExpRecord.stimuli{ftrrows(11)}; %ExpRecord.stimuli{241}; % Manif Expi 11
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1); % each row is an evolved code
Hproj_coef = codes_all * eigvect;
%%
ngen = max(code_geni);
fincodes = mean(Hproj_coef(code_geni==ngen, :),1);%,1);
prctl_thr = prctile(abs(fincodes(:)),[80,90,95,99,99.5,99.9,99.95],'all')';
pct_dstr_mat1 = cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((fincodes>cutoff | fincodes<-cutoff) & eigid>=eigthr-200+1 & eigid<=eigthr,2)),...
            prctl_thr)./min(eigthr,200),...[0,0.2,1,2,4,6,8,10,13,15,18]
        [25,50,100,200,400,800,1000,2000,3000,4096]','Uni',0))

pct_dstr_mat2 = cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((fincodes>cutoff | fincodes<-cutoff) & eigid<=eigthr,2)),...
            prctl_thr)./eigthr,...[0,0.2,1,2,4,6,8,10,13,15,18]
        [25,50,100,200,400,600,800,1000,2000,3000,4096]','Uni',0))
%%
ctrl_mat_col1 = [];
ctrl_mat_col2 = [];
for shfi = 1:500
shfl_codes = fincodes(randperm(4096));
ctrl_mat1 = cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((shfl_codes>cutoff | shfl_codes<-cutoff) & eigid>=eigthr-200+1 & eigid<=eigthr,2)),...
            prctl_thr)./min(eigthr,200), [25,50,100,200,400,800,1000,2000,3000,4096]','Uni',0));
ctrl_mat2 = cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((shfl_codes>cutoff | shfl_codes<-cutoff) & eigid<=eigthr,2)),...
            prctl_thr)./eigthr, [25,50,100,200,400,600,800,1000,2000,3000,4096]','Uni',0));
ctrl_mat_col1=cat(3,ctrl_mat_col1,ctrl_mat1);
ctrl_mat_col2=cat(3,ctrl_mat_col2,ctrl_mat2);
end
%%
pct_dstr_pool1 = [];
pct_dstr_pool2 = [];
for Expi = 1:45
stimuli_path = ExpRecord.stimuli{ftrrows(Expi)}; %ExpRecord.stimuli{241}; % Manif Expi 11
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1); % each row is an evolved code
Hproj_coef = codes_all * eigvect;
ngen = max(code_geni);
fincodes = mean(Hproj_coef(code_geni==ngen, :),1);%,1);
prctl_thr = prctile(abs(fincodes(:)),[80,90,95,99,99.5,99.9,99.95],'all')';
pct_dstr_mat1 = cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((fincodes>cutoff | fincodes<-cutoff) & eigid>=eigthr-200+1 & eigid<=eigthr,2)),...
            prctl_thr)./min(eigthr,200),...[0,0.2,1,2,4,6,8,10,13,15,18]
        [25,50,100,200,400,600,800,1000,2000,3000,4096]','Uni',0));

pct_dstr_mat2 = cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((fincodes>cutoff | fincodes<-cutoff) & eigid<=eigthr,2)),...
            prctl_thr)./eigthr,...[0,0.2,1,2,4,6,8,10,13,15,18]
        [25,50,100,200,400,600,800,1000,2000,3000,4096]','Uni',0));
pct_dstr_pool1 = cat(3,pct_dstr_pool1,pct_dstr_mat1);
pct_dstr_pool2 = cat(3,pct_dstr_pool2,pct_dstr_mat2);
end
%% Extreme values (meaningful genes) occurs more in the top Hessian spectra > bottom Hessian spectra
summarydir = "E:\OneDrive - Washington University in St. Louis\HessEvolStruct";
%%
figure(7);
subtightplot(1,2,1,0.06,0.08,0.08)
plot(mean(pct_dstr_pool2,3)'./[0.2,0.1,0.05,0.01,0.005,0.001,0.0005]')
xticks(1:7);xticklabels([0.2,0.1,0.05,0.01,0.005,0.001,0.0005])
legend(compose("Top %d d",[25,50,100,200,400,600,800,1000,2000,3000,4096]),'location','best')
xlabel("Percentile of Amplitude in the Projected Vector")
ylabel("Frequency Normalized to Percentage")
subtightplot(1,2,2,0.06,0.08,0.08);
plot(mean(pct_dstr_pool2,3)')
xticks(1:7);xticklabels([0.2,0.1,0.05,0.01,0.005,0.001,0.0005])
legend(compose("Top %d d",[25,50,100,200,400,600,800,1000,2000,3000,4096]),'location','best')
xlabel("Percentile of Amplitude in the Projected Vector")
ylabel("Frequency of Occuring")
suptitle("Distribution of Extreme Projected Coefficients Along the Spectra")
saveas(7,fullfile(summarydir,"Prctile_onTopHess.png"))
savefig(7,fullfile(summarydir,"Prctile_onTopHess.fig"))
%%
figure(8);
subtightplot(1,2,1,0.06,0.08,0.08)
plot(mean(pct_dstr_pool1,3)'./[0.2,0.1,0.05,0.01,0.005,0.001,0.0005]')
xticks(1:7);xticklabels([0.2,0.1,0.05,0.01,0.005,0.001,0.0005])
legend(compose("200d around %d d",[25,50,100,200,400,600,800,1000,2000,3000,4096]),'location','best')
xlabel("Percentile of Amplitude in the Projected Vector")
ylabel("Frequency Normalized to Percentage")
subtightplot(1,2,2,0.06,0.08,0.08);
plot(mean(pct_dstr_pool1,3)')
xticks(1:7);xticklabels([0.2,0.1,0.05,0.01,0.005,0.001,0.0005])
legend(compose("200d around %d d",[25,50,100,200,400,600,800,1000,2000,3000,4096]),'location','best')
xlabel("Percentile of Amplitude in the Projected Vector")
ylabel("Frequency of Occuring")
suptitle("Distribution of Extreme Projected Coefficients Along the Spectra")
saveas(8,fullfile(summarydir,"Prctile_onHessBand.png"))
savefig(8,fullfile(summarydir,"Prctile_onHessBand.fig"))
%% test the significance of the percentage value
pct_dstr_pool1_m = mean(pct_dstr_pool1,3);
pct_dstr_pool2_m = mean(pct_dstr_pool2,3);
cutoff_list = [25,50,100,200,400,600,800,1000,2000,3000,4096];
prct_list = [0.2,0.1,0.05,0.01,0.005,0.001,0.0005];
prob_mat2 = zeros(size(pct_dstr_pool2_m));
prob_mat1 = zeros(size(pct_dstr_pool2_m));
for ci = 1:numel(cutoff_list)
    for pi = 1:numel(prct_list)
        spacen = cutoff_list(ci);
        seln = floor(4096*prct_list(pi));
        prob_mat2(ci,pi) = hygecdf(spacen*pct_dstr_pool2_m(ci,pi),4096,seln,spacen);
%         prob_mat2(ci,pi) = discr_Binomial(spacen*pct_dstr_pool2_m(ci,pi),4096,seln,spacen);
        spacen = min(200,cutoff_list(ci));
%         prob_mat1(ci,pi) = discr_Binomial(spacen*pct_dstr_pool1_m(ci,pi),4096,seln,spacen);
        prob_mat1(ci,pi) = hygecdf(spacen*pct_dstr_pool1_m(ci,pi),4096,seln,spacen);
    end
end
%%
discr_Binomial(1:100,4096,100,400)
%% Comparing code to shuffling 
final_gen_xexp = [];
perm_gen_xexp = [];
for Expi = 1:45
stimuli_path = ExpRecord.stimuli{ftrrows(Expi)}; %ExpRecord.stimuli{241}; % Manif Expi 11
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1,"last"); % each row is an evolved code
mean_gene = mean(codes_all,1);
final_gen_xexp = [final_gen_xexp; mean_gene];
for shufi = 1:10
perm_gen_xexp = [perm_gen_xexp; mean_gene(randperm(4096))];
end
end
%%
code_proj_cc = final_gen_xexp * eigvect;
perm_proj_cc = perm_gen_xexp * eigvect;
code_proj_cc = code_proj_cc./norm_axis(code_proj_cc,2).*300;
perm_proj_cc = perm_proj_cc./norm_axis(perm_proj_cc,2).*300;
%%
figure;hold on 
histogram(code_proj_cc(:,3))
histogram(perm_proj_cc(:,3))
%%
rksm_col = arrayfun(@(id)ranksum(abs(code_proj_cc(:,id)),abs(perm_proj_cc(:,id))),1:4096,'uni',1);
%%
krsm_col = arrayfun(@(id)kstest2(abs(code_proj_cc(:,id)),abs(perm_proj_cc(:,id)),'alpha',0.01),1:4096,'uni',1);
%%
% N = 4096;
% sel = round(4096*0.2);
% n = 30;
% k = 0:30;
% prob = discr_Binomial(k,N,n,sel)
% function prob = discr_Binomial(k,N,n,sel)
% % Hand written version of Hyper geometric distribution.
% if all(k>n | k>sel)
% prob = zeros(size(k));
% return
% end
% 
% prob = beta(n+1,N-n+1)./beta(n-k+1,N-sel-n+k+1)./beta(sel-k+1,k+1)*(N+1)/(sel+1)/(N-sel+1);
% end
% 
