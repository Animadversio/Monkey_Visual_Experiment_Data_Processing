% Plot a sample Gabor orientation tuning experiments with a recent
%   Evolution Dataset. 
% Load up a recent experiment 
Set_Path;
rowid = find(contains(ExpRecord.ephysFN,'Beto-11112021-003'));
[meta_new,rasters_new,~,Trials_new] = loadExperiments(rowid, "Beto", false); %'Beto-11112021-003'
EStats = Evol_Collect_Stats_fun(meta_new, rasters_new, Trials_new);
%%
psth_cat = cat(3,EStats.ref.psth_arr{:});
refimg_bsl_mean = mean(psth_cat(:,1:50,:),'all');
refimg_bsl_sem = sem(mean(psth_cat(:,1:50,:),[2]));
refimg_mean = cellfun(@(P) mean(P(1,51:200,:), 'all'),EStats.ref.psth_arr);
refimg_sem  = cellfun(@(P) sem(mean(P(1,51:200,:), [1,2])),EStats.ref.psth_arr);
refimg_std  = cellfun(@(P) std(mean(P(1,51:200,:), [1,2])),EStats.ref.psth_arr);
%%
img_param = regexp(EStats.ref.imgnm,"gab_ori_(?<ori>[\d.]*)_(?<sf>[\d.]*)_",'names');
img_param = cat(1,img_param{:});
orientations=str2double([img_param.ori]);
spfreqs=str2double([img_param.sf]);
[~, sortidx] = sortrows([spfreqs',orientations']);
%%
scatter(refimg_mean,refimg_sem)
%%
figure(1);clf;
shadedErrorBar([], refimg_mean(sortidx([1:6,1])), refimg_sem(sortidx([1:6,1])), 'LineProp', {'Display', "spatialFreq 0.5"});
shadedErrorBar([], refimg_mean(sortidx([7:12,7])), refimg_sem(sortidx([7:12,7])), 'LineProp', {'Display', "spatialFreq 1.5"});
shadedErrorBar([], repmat(refimg_bsl_mean,1,7), repmat(refimg_bsl_sem,1,7), 'LineProp', {'-.', 'Display', "Baseline firing"});
legend()
% xticklabels(orientations(sortidx([1:6,1])))
xticklabels([0:30:180])
ylabel("Evoked Firing Rate (sp/sec)")
xlabel("Orientation of Gabor Grating")
ylabel("Evoked Firing Rate (sp/sec)")
title("Example V1 neuron")
saveallform("O:\Manuscript_Manifold\Figure1","tuning_curv_example",1)
%%
sort_psth_col = EStats.ref.psth_arr(sortidx([1:6,1]));
rate_col = cellfun(@(P)squeeze(mean(P(1,51:200,:),[1,2])),sort_psth_col,'Unif',0);
ori_col = arrayfun(@(i)i * ones(numel(rate_col{i}),1), 1:7,'Unif',0);
rate_vec = cat(1,rate_col{:});
idx_vec = cat(1,ori_col{:});
%%
sort_psth_col = EStats.ref.psth_arr(sortidx([7:12,7]));
rate_col = cellfun(@(P)squeeze(mean(P(1,51:200,:),[1,2])),sort_psth_col,'Unif',0);
ori_col = arrayfun(@(i)i * ones(numel(rate_col{i}),1), 1:7,'Unif',0);
rate_vec2 = cat(1,rate_col{:});
idx_vec2 = cat(1,ori_col{:});
%%
figure(2);clf;
shadedErrorBar([], refimg_mean(sortidx([1:6,1])), refimg_sem(sortidx([1:6,1])), 'LineProp', {'Display', "spatialFreq 0.5"});
shadedErrorBar([], refimg_mean(sortidx([7:12,7])), refimg_sem(sortidx([7:12,7])), 'LineProp', {'Display', "spatialFreq 1.5"});
hold on;
scatter(idx_vec, rate_vec,9,'blue','LineWidth',2)
scatter(idx_vec2, rate_vec2,9,'magenta','LineWidth',2)
shadedErrorBar([], repmat(refimg_bsl_mean,1,7), repmat(refimg_bsl_sem,1,7), 'LineProp', {'-.k', 'Display', "Baseline firing"});
legend()
% xticklabels(orientations(sortidx([1:6,1])))
xticklabels([0:30:180])
ylabel("Evoked Firing Rate (sp/sec)")
xlabel("Orientation of Gabor Grating")
ylabel("Evoked Firing Rate (sp/sec)")
title("Example V1 neuron")
saveallform("O:\Manuscript_Manifold\Figure1","tuning_curv_example_with_data",2)
%% Fit gp regress on the raw single trial firing 
gprMdl = fitrgp((idx_vec-1)*30, rate_vec);%, 'SigmaLowerBound', 5E-3)
[ypred,ysd,yint] = predict(gprMdl, [0:1:180]','Alpha',0.05);
gprMdl2 = fitrgp((idx_vec2-1)*30, rate_vec2);%, 'SigmaLowerBound', 5E-3)
[ypred2,ysd2,yint2] = predict(gprMdl2, [0:1:180]','Alpha',0.05);
%% Fit gp regress on the bootstrapped mean
rate_col_all = cellfun(@(P)squeeze(mean(P(1,51:200,:),[1,2])),EStats.ref.psth_arr(sortidx),'unif',0);
meanrate_btrps = cellfun(@(rates) bootstrp(100,@mean,rates),rate_col_all,'unif',0);
%%
btrprate_col1 = meanrate_btrps([1:6,1]);
ori_col1 = arrayfun(@(i)30 * (i-1) * ones(numel(btrprate_col1{i}),1), 1:7, 'Unif',0);
btrprate_vec1 = cat(1,btrprate_col1{:});
ori_vec1 = cat(1,ori_col1{:});

btrprate_col2 = meanrate_btrps([7:12,7]);
ori_col2 = arrayfun(@(i)30 * (i-1) * ones(numel(btrprate_col2{i}),1), 1:7, 'Unif',0);
btrprate_vec2 = cat(1,btrprate_col2{:});
ori_vec2 = cat(1,ori_col2{:});

rate_col1 = rate_col_all([1:6,1]);
rate_vec1 = cat(1, rate_col1{:});
idx_col1 = arrayfun(@(i)i * ones(numel(rate_col1{i}),1), 1:7,'Unif',0);
idx_vec1 = cat(1,idx_col1{:});
rate_col2 = rate_col_all([7:12,7]);
rate_vec2 = cat(1, rate_col2{:});
idx_col2 = arrayfun(@(i)i * ones(numel(rate_col2{i}),1), 1:7,'Unif',0);
idx_vec2 = cat(1,idx_col2{:});

gprMdl1 = fitrgp(ori_vec1, btrprate_vec1);%, 'SigmaLowerBound', 5E-3)
[ypred1,ysd1,yint1] = predict(gprMdl1, [0:1:180]','Alpha',0.05);
gprMdl2 = fitrgp(ori_vec2, btrprate_vec2);%, 'SigmaLowerBound', 5E-3)
[ypred2,ysd2,yint2] = predict(gprMdl2, [0:1:180]','Alpha',0.05);
%% Gaussian process model of the tuning curve... 
figure(4);clf;
% shadedErrorBar([], refimg_mean(sortidx([1:6,1])), refimg_sem(sortidx([1:6,1])), 'LineProp', {'Display', "spatialFreq 0.5"});
% shadedErrorBar([], refimg_mean(sortidx([7:12,7])), refimg_sem(sortidx([7:12,7])), 'LineProp', {'Display', "spatialFreq 1.5"});
shadedErrorBar([0:1:180], ypred1, ysd1, 'LineProp', {'Display', "spatialFreq 0.5"});
shadedErrorBar([0:1:180], ypred2, ysd2, 'LineProp', {'Display', "spatialFreq 1.5"});
hold on;
scatter(30*(idx_vec-1), rate_vec1,9,'blue','LineWidth',2,'MarkerEdgeAlpha',0.4)
scatter(30*(idx_vec2-1), rate_vec2,9,'magenta','LineWidth',2,'MarkerEdgeAlpha',0.4)
shadedErrorBar(0:180, repmat(refimg_bsl_mean,1,181), repmat(refimg_bsl_sem,1,181), 'LineProp', {'-.k', 'Display', "Baseline firing"});
legend('location','best')
xlim([0,180]);xticks([0:30:180])
% xticklabels([0:30:180])
ylabel("Evoked Firing Rate (sp/sec)")
xlabel("Orientation of Gabor Grating")
ylabel("Evoked Firing Rate (sp/sec)")
title("Example V1 Orientation Tuning")
saveallform("O:\Manuscript_Manifold\Figure1","tuning_curv_example_gprfit_bstrpmean_wdata",4)
%%
figure(5);clf;
shadedErrorBar([0:1:180], ypred1, ysd1, 'LineProp', {'Display', "spatialFreq 0.5"});
shadedErrorBar([0:1:180], ypred2, ysd2, 'LineProp', {'Display', "spatialFreq 1.5"});
hold on;
shadedErrorBar(0:180, repmat(refimg_bsl_mean,1,181), repmat(refimg_bsl_sem,1,181), 'LineProp', {'-.k', 'Display', "Baseline firing"});
scatter([0:30:180], refimg_mean(sortidx([1:6,1])),  25, 'k');
scatter([0:30:180], refimg_mean(sortidx([7:12,7])), 25, 'k');
legend('location','best')
xlim([0,180]);xticks([0:30:180])
% xticklabels([0:30:180])
ylabel("Evoked Firing Rate (sp/sec)")
xlabel("Orientation of Gabor Grating")
ylabel("Evoked Firing Rate (sp/sec)")
title("Example V1 Orientation Tuning")
saveallform("O:\Manuscript_Manifold\Figure1","tuning_curv_example_gprfit_bstrpmean",5)