%% Process experiment with CMA vs GA optimizers. 
%  loading => Collect Stats => Generate summary figure 
!ExpRecordBackup.bat 
%%
Animal = "Both";Set_Path;
expftr = contains(ExpRecord.expControlFN,"generate") & ...
        contains(ExpRecord.Exp_collection, "CMAGA_cmp");
row_idx = find(expftr);
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(row_idx, Animal, false);
%%
% Data Structures to store score traj and zscore traj 
% exp num-by-thread num, in each cell there is a 1d vector [1, image-num] 
score_traces = cell(length(meta_new), 2); 
zscore_traces = cell(length(meta_new), 2); 
block_traces = cell(length(meta_new), 2); % vector of block num for plotting
ref_traces = cell(length(meta_new), 2);
ref_block_traces = cell(length(meta_new), 2); % vector of block num for plotting
for Triali = 1:length(meta_new)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);
fprintf("Processing Exp %d:\n",Expi)
disp(ExpRecord.comments(exp_rowi))
% assert(Expi_tab == Expi, "Data Expi doesn't match that in exp record! Check your indexing or record.")
%% Sort channel id
% finding spike ID, note for multi-thread optimizer, we will have multiple
% pref_chan for different optimizers 
pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_prefchan = Trials.TrialRecord.User.evoConfiguration{1,4};
assert(pref_chan(1) == pref_chan(2))
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
try
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
catch
unit_name_arr = generate_unit_labels(meta.spikeID);
end
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
pref_chan_id = find(meta.spikeID==pref_chan(1)&unit_num_arr==unit_in_prefchan(1)); % the id in the raster and lfps matrix 

%% Optimizer Names 
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, string(Trials.TrialRecord.User.evoConfiguration{i,end})];
end
if contains(Optim_names(1),"CMAES")&&contains(Optim_names(2),"GA"),
optimMap = [1,2]; % make sure CMA in the 1st entry in summary cell. 
else % make sure GA in the 2nd entry in summary cell. 
optimMap = [2,1];
end
%% Sort the images
imgnm = Trials.imageName;
% seperate the thread natural images and generated images 
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm); % first 2 characters are not digits
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the thread 1 and 2 
row_thread0 = contains(imgnm, compose("thread%03d", 0));
row_thread1 = contains(imgnm, compose("thread%03d", 1));
assert(sum(row_thread0)+sum(row_thread1) == length(imgnm))
thread_msks = {row_thread0, row_thread1}; % store masks in a structure for the ease to iterate
% get the generation number 
block_arr = cell2mat(Trials.block);
% if needed, analyze the image names to see the block and thread
% information. see Evol_Traj_Cmp
%% Compute score for evolution 
block_list = min(block_arr):max(block_arr)-1;% last (incompleted) block is already croped here. 
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2)); % [unit_num, img_nums]
zscores_tsr = zscore(scores_tsr,1,2); % [unit_num, img_nums]
bsl_tsr = squeeze(mean(rasters(:, 1:40, :), 2));
bslmean = mean(bsl_tsr,2);
bslsem = sem(bsl_tsr,2);
meanscore_syn = nan(size(rasters, 1), length(block_list), 2); % [unit_num, gen_nums, threads]
stdscore_syn = nan(size(rasters, 1), length(block_list), 2); 
meanscore_nat = nan(size(rasters, 1), length(block_list), 2);
stdscore_nat = nan(size(rasters, 1), length(block_list), 2);
for threadi = 1:thread_num 
    for blocki = block_list
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        meanscore_syn(:, blocki, threadi) = mean(scores_tsr(:, gen_msk), 2);
        meanscore_nat(:, blocki, threadi) = mean(scores_tsr(:, nat_msk), 2);
        stdscore_syn(:, blocki, threadi)  = std(scores_tsr(:, gen_msk), 1, 2) / sqrt(sum(gen_msk));
        stdscore_nat(:, blocki, threadi)  = std(scores_tsr(:, nat_msk), 1, 2) / sqrt(sum(nat_msk));
    end
    % Store the exp stats into the pools
    optimi=optimMap(threadi); % CMA at 1st GA at 2nd
    score_traces{Triali,optimi} = scores_tsr(pref_chan_id, row_gen & thread_msks{threadi}) - bslmean(pref_chan_id); % [1, gen_imgN]
    zscore_traces{Triali,optimi} = zscores_tsr(pref_chan_id, row_gen & thread_msks{threadi}); % [1, gen_imgN]
    block_traces{Triali,optimi} = block_arr(row_gen & thread_msks{threadi}); % [1, gen_imgN]
    ref_traces{Triali,optimi} = scores_tsr(pref_chan_id, row_nat & thread_msks{threadi}) - bslmean(pref_chan_id); % [1, nat_pres]
    ref_block_traces{Triali,optimi} = block_arr(row_nat & thread_msks{threadi}); % [1, nat_pres]
end
% %% Compute Average PSTH and sem for evolved image
% evol_stim_fr = nan(size(rasters, 1), size(rasters, 2), length(block_list), thread_num);
% evol_stim_sem = nan(size(rasters, 1), size(rasters, 2), length(block_list), thread_num);
% for threadi = 1:thread_num
% for blocki = 1:length(block_list)
%     gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
%     evol_stim_fr(:, :, blocki, threadi) = mean(rasters(:,:, gen_msk),3);
%     evol_stim_sem(:, :, blocki, threadi) = std(rasters(:,:, gen_msk),1,3) / sqrt(sum(gen_msk));
% end
% end
end
%% Note the last exp should be exluded, not enough trial to do anything.

%% A demo of testing the score of last 5 gen
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
scores_thread1 = scores_tsr(pref_chan_id, row_gen & (block_arr > max(block_arr)-6) & thread_msks{1});
scores_thread2 = scores_tsr(pref_chan_id, row_gen & (block_arr > max(block_arr)-6) & thread_msks{2});
[~,Pval,CI] = ttest2(scores_thread1, scores_thread2);

%% Plot the Score Evolution Trajectory of CMA vs GA
%  Note the key is proper normalization to show all score traj together. 
figdir = "O:\Optimizer_Cmp\summary";
Corder = colororder;
%% Scatter the zscored score trajectory for each experiment
figure;hold on;fignm="zscore_scatter_all";
for Expi=1:size(zscore_traces,1)-1
    scatter(block_traces{Expi,1},zscore_traces{Expi,1},'MarkerEdgeColor','none','MarkerFaceColor',Corder(2,:),'MarkerFaceAlpha',0.05)
    scatter(block_traces{Expi,2},zscore_traces{Expi,2},'MarkerEdgeColor','none','MarkerFaceColor',Corder(1,:),'MarkerFaceAlpha',0.05)
end
ylabel("z-scored Activation each session")
xlabel("Generation")
legend(["CMAES","GA-classic"])
saveallform(figdir,fignm);
xlim([0,50]);fignm = fignm+"_Xlim";
saveallform(figdir,fignm);
%% Plot zscored score trajeoctory for each experiment
figure;hold on;fignm="zscore_scoreTraj_all";
for Expi=1:size(zscore_traces,1)-1
    [zscore_C_m,zscore_C_s,blockvec] = sort_scoreblock(block_traces{Expi,1},zscore_traces{Expi,1});
    shadedErrorBar(blockvec,zscore_C_m,zscore_C_s,'lineProps',{'Color',Corder(2,:)})
    [zscore_G_m,zscore_G_s,blockvec] = sort_scoreblock(block_traces{Expi,2},zscore_traces{Expi,2});
    shadedErrorBar(blockvec,zscore_G_m,zscore_G_s,'lineProps',{'Color',Corder(1,:)})
end
ylabel("z-scored Activation each session")
xlabel("Generation")
legend(["CMAES","GA-classic"])
saveallform(figdir,fignm);
xlim([0,50]);fignm = fignm+"_Xlim";
saveallform(figdir,fignm);
%% Plot raw score trajeoctory for each experiment
figure;hold on;fignm="scoreTraj_example";
for Expi=7%1:size(zscore_traces,1)-1
    [score_C_m,score_C_s,blockvec] = sort_scoreblock(block_traces{Expi,1},score_traces{Expi,1});
    [score_G_m,score_G_s,blockvec] = sort_scoreblock(block_traces{Expi,2},score_traces{Expi,2});
    shadedErrorBar(blockvec,score_C_m,score_C_s,'lineProps',{'Color',[Corder(2,:),0.6],'lineWidth',1})
    shadedErrorBar(blockvec,score_G_m,score_G_s,'lineProps',{'Color',[Corder(1,:),0.6],'lineWidth',1})
end
title(compose("Example GA-CMAES comparison Exp%d %s Chan%02d",Expi,meta_new{Expi}.ephysFN,Trials_new{Expi}.TrialRecord.User.prefChan(1)))
ylabel("z-scored Activation each session")
xlabel("Generation")
legend(["CMAES","GA-classic"])
saveallform(figdir,fignm);
%% Plot max normalized score trajeoctory for each experiment
score_C_traj_col = [];
score_G_traj_col = [];
block_traj_col = [];
figure;hold on;fignm="MaxNorm_scoreTraj_all";
for Expi=1:size(score_traces,1)-1
    [score_C_m,score_C_s,blockvec] = sort_scoreblock(block_traces{Expi,1},score_traces{Expi,1});
    [score_G_m,score_G_s,blockvec] = sort_scoreblock(block_traces{Expi,2},score_traces{Expi,2});
    normmin = 0;%(score_C_m(1)+score_G_m(1))/2;
    normmax = max(score_C_m);
    scaling = abs(normmax - normmin);
    score_C_m_norm = (score_C_m-normmin)/scaling;
    score_G_m_norm = (score_G_m-normmin)/scaling;
    score_C_traj_col = [score_C_traj_col,score_C_m_norm];
    score_G_traj_col = [score_G_traj_col,score_G_m_norm];
    block_traj_col = [block_traj_col,blockvec];
    shadedErrorBar(blockvec,(score_C_m-normmin)/scaling,score_C_s/scaling,'lineProps',{'Color',[Corder(2,:),0.6],'lineWidth',1},'patchSaturation',0.1)
    shadedErrorBar(blockvec,(score_G_m-normmin)/scaling,score_G_s/scaling,'lineProps',{'Color',[Corder(1,:),0.6],'lineWidth',1},'patchSaturation',0.1)
end
xlabel("Generation")
ylabel("Activation / Max each session")
legend(["CMAES","GA-classic"])
% saveallform(figdir,fignm);
xlim([0,50]);fignm = fignm+"_Xlim";
% saveallform(figdir,fignm);
%%
figure;hold on;fignm="MaxNorm_scoreTraj_ColAvg";
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_traj_col,score_C_traj_col);
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(block_traj_col,score_G_traj_col);
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title("Averaged Optimization Trajectory")
legend(["CMAES","GA-classic"])
% saveallform(figdir,fignm);
% xlim([0,50]);fignm="MaxNorm_scoreTraj_ColAvg_Xlim";
% saveallform(figdir,fignm);
%% Filter the experiments, get rid of non successful ones. 
% score_traces{Triali,optimi} = scores_tsr(pref_chan_id, row_gen & thread_msks{threadi}); % [1, gen_imgN]
% zscore_traces{Triali,optimi} = zscores_tsr(pref_chan_id, row_gen & thread_msks{threadi}); % [1, gen_imgN]
% block_traces{Triali,optimi} = block_arr(row_gen & thread_msks{threadi}); % [1, gen_imgN]
optimstr = ["CMAES","GA-cl"];
tval_end_arr = []; pval_end_arr = [];
tval_max_arr = []; pval_max_arr = [];
for Triali = 1:size(score_traces, 1)
for optimi=1:2
fprintf(optimstr(optimi)+" ")
blockN = max(block_traces{Triali,optimi});
[tval,pval] = ttest2_print(score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[1,2],2)),...
    score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[blockN-1,blockN],2)),"init12","last12");
tval_end_arr(Triali,optimi) = tval;
pval_end_arr(Triali,optimi) = pval;
[score_m,score_s,blockvec] = sort_scoreblock(block_traces{Triali,optimi},...
                score_traces{Triali,optimi});
[~,maxN] = max(score_m);
[tval,pval] = ttest2_print(score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[1,2],2)),...
    score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[maxN-1,maxN],2)),"init12","last12");
tval_max_arr(Triali,optimi) = tval;
pval_max_arr(Triali,optimi) = pval;
end
end
sucsmsk = any((tval_end_arr<0)&(pval_end_arr<0.001), 2);
sucsmsk_max = any((tval_max_arr<0)&(pval_max_arr<0.001), 2);
%% Collect the mean of each trajectory, Add Extrapolation to match the block number `blockN_extrap` for better mean and sem curve. 
score_C_traj_col = [];
score_G_traj_col = [];
block_traj_col = [];
blockN_extrap = 60;
score_C_traj_extrap_col = [];
score_G_traj_extrap_col = [];
block_traj_extrap_col = [];
for Expi=1:size(score_traces,1)
    [score_C_m,score_C_s,blockvec] = sort_scoreblock(block_traces{Expi,1},score_traces{Expi,1});
    [score_G_m,score_G_s,blockvec] = sort_scoreblock(block_traces{Expi,2},score_traces{Expi,2});
    normmin = 0;%(score_C_m(1)+score_G_m(1))/2;
    normmax = max(score_C_m);
    scaling = abs(normmax - normmin);
    score_C_m_norm = (score_C_m-normmin)/scaling;
    score_G_m_norm = (score_G_m-normmin)/scaling;
    score_C_traj_col{Expi} = score_C_m_norm;
    score_G_traj_col{Expi} = score_G_m_norm;
    block_traj_col{Expi} = blockvec;
    
    if numel(score_G_m_norm)>=2
    extrap_val_C = mean(score_C_m_norm(end-1:end));
    extrap_val_G = mean(score_G_m_norm(end-1:end));
    else
    extrap_val_C = mean(score_C_m_norm(end));
    extrap_val_G = mean(score_G_m_norm(end));
    end
    if numel(blockvec) < blockN_extrap
        score_C_traj_extrap_col{Expi} = [score_C_m_norm,extrap_val_C * ones(1,blockN_extrap-numel(score_C_m_norm))];
        score_G_traj_extrap_col{Expi} = [score_G_m_norm,extrap_val_G * ones(1,blockN_extrap-numel(score_G_m_norm))];
        block_traj_extrap_col{Expi} = 1:blockN_extrap;
    else
        score_C_traj_extrap_col{Expi} = score_C_m_norm;
        score_G_traj_extrap_col{Expi} = score_G_m_norm;
        block_traj_extrap_col{Expi} = blockvec;
    end
end
%% Plot the mean sem curve with extraplotated individual exp.
figure;hold on;fignm="MaxNorm_scoreTraj_sucsonly_ColAvg";
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_traj_col(sucsmsk),score_C_traj_col(sucsmsk));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(block_traj_col(sucsmsk),score_G_traj_col(sucsmsk));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title("Averaged Optimization Trajectory (success only)")
legend(["CMAES","GA-classic"])
xlim([0,50])
saveallform(figdir,fignm)
% TODO represent the extrapolation with dashed line. 
figure;hold on;fignm="MaxNorm_scoreTraj_sucsonly_extrap_ColAvg";
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_traj_extrap_col(sucsmsk),score_C_traj_extrap_col(sucsmsk));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(block_traj_extrap_col(sucsmsk),score_G_traj_extrap_col(sucsmsk));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory (success only, extrap to %d)",blockN_extrap))
legend(["CMAES","GA-classic"])
xlim([0,50])
saveallform(figdir,fignm)
% TODO represent the extrapolation with dashed line. 
figure;hold on;fignm="MaxNorm_scoreTraj_sucsonly_extrap_ColAvg_Indiv";
[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_traj_extrap_col(sucsmsk),score_C_traj_extrap_col(sucsmsk));
[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(block_traj_extrap_col(sucsmsk),score_G_traj_extrap_col(sucsmsk));
shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:),'LineWidth',2})
shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:),'LineWidth',2})
for Expi=1:size(score_traces,1)
    if ~sucsmsk(Expi),continue;end
    plot(block_traj_col{Expi},movmean(score_C_traj_col{Expi},3),'Color',[Corder(2,:),0.5],'LineWidth',1)
    plot(block_traj_col{Expi},movmean(score_G_traj_col{Expi},3),'Color',[Corder(1,:),0.5],'LineWidth',1)
end
xlabel("Generation")
ylabel("Activation / Max each session")
title(compose("Averaged Optimization Trajectory (success only, extrap to %d)",blockN_extrap))
legend(["CMAES","GA-classic"],'Location','Best')
xlim([0,50])
saveallform(figdir,fignm)
% winopen(figdir)
%% Statistics summary scatter plot
initact_mean = []; initact_sem = [];
endact_mean = []; endact_sem = [];
maxact_mean = []; maxact_sem = [];
for Triali = 1:size(score_traces, 1)
for optimi=1:2
blockN = max(block_traces{Triali,optimi});
[score_m,score_s,blockvec] = sort_scoreblock(block_traces{Triali,optimi}, score_traces{Triali,optimi});
[~,maxN] = max(score_m);
initact_mean(Triali,optimi) = mean(score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[1,2],2)));
initact_sem(Triali,optimi) = sem(score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[1,2],2)));

endact_mean(Triali,optimi) = mean(score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[blockN-1,blockN],2))); 
endact_sem(Triali,optimi) = sem(score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[blockN-1,blockN],2)));

maxact_mean(Triali,optimi) = mean(score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[maxN-1,maxN],2))); 
maxact_sem(Triali,optimi) = sem(score_traces{Triali,optimi}(any(block_traces{Triali,optimi}==[maxN-1,maxN],2)));
end
end
figure;set(gcf,'pos', [1000         421         575         550])
errorbar(endact_mean(sucsmsk,1),endact_mean(sucsmsk,2),...
    endact_sem(sucsmsk,2),endact_sem(sucsmsk,2),endact_sem(sucsmsk,1),endact_sem(sucsmsk,1),'o')
xlabel("CMAES");ylabel("GA classic");title(["Last Generation activation compare","(sp/sec, - bsl, sem of 2 blocks)"])
axis equal square;box off
add_diagonal(gca,'k-.');
saveallform(figdir, "LastGen_GACMA_Cmp")
%% Statistical test. 
%% t test the last 5 generations 
CMA_score = zeros(length(meta_new),1);
CMA_err = zeros(length(meta_new),1);
GA_score = zeros(length(meta_new),1);
GA_err = zeros(length(meta_new),1);
CMA_ref_score = zeros(length(meta_new),1);
CMA_ref_err = zeros(length(meta_new),1);
GA_ref_score = zeros(length(meta_new),1);
GA_ref_err = zeros(length(meta_new),1);
for Triali = 1:length(meta_new)
    scores = score_traces{Triali,1}(block_traces{Triali,1} > max(block_traces{Triali,1})-5);
    CMA_score(Triali) = mean(scores);
    CMA_err(Triali) = std(scores)/sqrt(length(scores));
    scores = score_traces{Triali,2}(block_traces{Triali,2} > max(block_traces{Triali,2})-5);
    GA_score(Triali) = mean(scores);
    GA_err(Triali) = std(scores)/sqrt(length(scores));
    
    scores = ref_traces{Triali,1}(ref_block_traces{Triali,1} > max(ref_block_traces{Triali,1})-5);
    CMA_ref_score(Triali) = mean(scores);
    CMA_ref_err(Triali) = std(scores)/sqrt(length(scores));
    scores = ref_traces{Triali,2}(ref_block_traces{Triali,2} > max(ref_block_traces{Triali,2})-5);
    GA_ref_score(Triali) = mean(scores);
    GA_ref_err(Triali) = std(scores)/sqrt(length(scores));
end
%% stats for normalized CMA, GA. 
normact_CMA = cellfun(@(traj)mean(traj(end-1:end)),score_C_traj_col(sucsmsk));
normact_GA  = cellfun(@(traj)mean(traj(end-1:end)),score_G_traj_col(sucsmsk));
ttest2_print(normact_CMA,normact_GA,"CMA","GA",true)
%%
% blockarr = block_traces{Expi,1};
% scorearr = zscore_traces{Expi,1};
% [zscore_m,zscore_s,blockvec] = sort_scoreblock(block_traces{Expi,1},zscore_traces{Expi,1})
% function [score_m,score_s,blockvec] = sort_scoreblock(blockarr,scorearr)
% blockvec = min(blockarr):max(blockarr);
% for blocki = min(blockarr):max(blockarr)
%     score_m(blocki) = mean(scorearr(blockarr==blocki));
%     score_s(blocki) = sem(scorearr(blockarr==blocki));
% end
% end
% function saveallform(figdir,fignm,h)
% if nargin == 2, h=gcf; end
% savefig(h,fullfile(figdir,fignm+".fig"))
% saveas(h,fullfile(figdir,fignm+".pdf"))
% saveas(h,fullfile(figdir,fignm+".png"))
% end