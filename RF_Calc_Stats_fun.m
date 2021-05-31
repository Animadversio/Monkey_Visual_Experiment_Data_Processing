function S_col = RF_Calc_Stats_fun(meta_new, rasters_new, Trials_new)
pixperdeg = 40;
P = struct();
P.plotEachChan = false;
P.collectstat = false;
savedir = "O:\RFstats";
if P.collectstat, S_col = []; end
for iTr = 1:numel(meta_new)
Trials = Trials_new{iTr}; 
rasters = rasters_new{iTr}; 
meta = meta_new{iTr};
if isempty(meta), continue;end 
fprintf("Processing RFMap %d\n",iTr)
S = struct();
if isfield(meta,"unitID")
unit_num_arr = meta.unitID;
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
activ_msk = unit_num_arr~=0;
else
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
end
S.unit.unit_name_arr = unit_name_arr;
S.unit.activ_msk = activ_msk;
S.unit.unit_num_arr = unit_num_arr;
S.unit.chan_num_arr = meta.spikeID;
S.meta = meta;

xy_all = Trials.TrialRecord.User.xy; % basis to sort the trials 
size_all = Trials.width(:,1);
uniqsize_all = sort(unique(size_all));
uniqsize_deg_all = sort(Trials.TrialRecord.User.uniqueSizes_deg); % in deg
nImSize = length(uniqsize_deg_all); % before its always 1, recently, we use multiple image sizes! 
S.stim.xy_all = xy_all;
S.stim.size_all = size_all;
S.stim.size_deg = size_all / pixperdeg;

idx_grid_all = cell(numel(uniqsize_all),1);% collect the maps together. 
x_grid_all = cell(numel(uniqsize_all),1);
y_grid_all = cell(numel(uniqsize_all),1);
for iSz = 1:nImSize % iterate image size 
imgsize = uniqsize_all(iSz); 
imgsize_deg = uniqsize_deg_all(iSz); 
uniqpos = unique(xy_all(size_all==imgsize, :), 'rows');%Trials.TrialRecord.User.uniquePositions;
xpos = unique(uniqpos(:,1));
ypos = unique(uniqpos(:,2));
nPos = size(uniqpos,1);
% This is a key assumption! The position could be factorized into a square grid. 
assert(nPos==length(xpos)*length(xpos)); % if there is multiple sets of xpos this will fail
posgrid = reshape([uniqpos, (1:nPos)'],[length(xpos),length(ypos),3]);
% sort neural responses by getting the trial indices
idx_grid = cell(length(ypos),length(xpos)); % trial index
iPos_grid = cell(length(ypos),length(xpos)); % position index for grid
x_grid = nan(length(ypos),length(xpos));
y_grid = nan(length(ypos),length(xpos));
for iPos = 1:nPos % iterate position
    tPos = uniqpos(iPos,:) ;
    idx = find( xy_all(:,1) == tPos(1) & xy_all(:,2) == tPos(2) & size_all==imgsize) ; % debug this May26. before without adding the size criterion. 
    gridi = find(xpos == tPos(1));
    gridj = find(ypos == tPos(2));
    % posgrid(gridi,gridj,1) == tPos(1)
    % posgrid(gridi,gridj,2) == tPos(2)
    idx_grid{gridj, gridi} = idx;
    pos_grid{gridj, gridi} = tPos;
    iPos_grid{gridj, gridi} = repmat([iPos],1,length(idx));
    x_grid(gridj, gridi) = tPos(1); 
    y_grid(gridj, gridi) = tPos(2); 
end
idx_grid_all{iSz} = idx_grid; 
x_grid_all{iSz} = x_grid;
y_grid_all{iSz} = y_grid;
nrep = cellfun(@length, idx_grid);
%% averaged psth for each location. 
psth_mean = cellfun(@(idx) mean(rasters(:,:,idx),3), idx_grid, 'uni', 0); % Each cell is channel # by 200
psth_sem = cellfun(@(idx) std(rasters(:,:,idx),0,3) / sqrt(numel(idx)), idx_grid, 'uni', 0); % Each cell is channel # by 200
%  Compute the score (psth time averaged)
wdwbsl = 1:50; wdwevk = 51:200;
act_col   = cellfun(@(idx) squeeze(mean(rasters(:,wdwevk,idx),[2])), idx_grid, 'uni', 0); 
score_col = cellfun(@(idx) squeeze(mean(rasters(:,wdwevk,idx),[2])) - squeeze(mean(rasters(:,wdwbsl,idx),[2,3])), idx_grid, 'uni', 0); 
act_mean   = cellfun(@(idx) mean(rasters(:,wdwevk,idx),[2,3]), idx_grid, 'uni', 0); 
act_mean   = cell2mat(cellfun(@(act)reshape(act,1,1,[]), act_mean,'uni',0));
score_mean = cellfun(@(idx) mean(rasters(:,wdwevk,idx),[2,3]) - mean(rasters(:,wdwbsl,idx),[2,3]), idx_grid, 'uni', 0); 
score_mean = cell2mat(cellfun(@(act)reshape(act,1,1,[]),score_mean,'uni',0));
score_std  = cellfun(@(idx) squeeze(std(mean(rasters(:,wdwevk,idx),[2]),1,3)), idx_grid, 'uni', 0); 
score_std  = cell2mat(cellfun(@(act)reshape(act,1,1,[]),score_std,'uni',0));
score_sem  = cellfun(@(idx) squeeze(sem(mean(rasters(:,wdwevk,idx),[2]),3)), idx_grid, 'uni', 0); 
score_sem  = cell2mat(cellfun(@(act)reshape(act,1,1,[]),score_sem,'uni',0));
psthbest_mean = nan(numel(unit_num_arr),200);
psthbest_sem = nan(numel(unit_num_arr),200);
for iCh = 1:numel(unit_num_arr)
[maxscore, idx] = max(score_mean(:,:,iCh),[],'all','linear');
[idx_j,idx_i] = ind2sub(size(score_mean(:,:,iCh)), idx);
psthbest_mean(iCh,:) = psth_mean{idx_j,idx_i}(iCh,:);
psthbest_sem(iCh,:)  = psth_sem{idx_j,idx_i}(iCh,:);
end
S.psth.idx_grid{iSz} = idx_grid; 
S.psth.psth_mean{iSz} = psth_mean; % cell array of chan-200
S.psth.psth_sem{iSz} = psth_sem;
S.psth.psth_best_mean{iSz} = psthbest_mean; % chan N, 200 mat
S.psth.psth_best_sem{iSz} = psthbest_sem; % chan N, 200 mat
S.psth.act_col{iSz} = act_col; % cell array of chan-repN
S.psth.act_mean{iSz} = act_mean; % ypos N, xpos N, chan N tsr
S.psth.score_mean{iSz} = score_mean; % ypos N, xpos N, chan N tsr
S.psth.score_std{iSz} = score_std; % ypos N, xpos N, chan N tsr
S.psth.score_sem{iSz} = score_sem; % ypos N, xpos N, chan N tsr
S.stim.nrep{iSz} = nrep; % matrix of nrep
S.stim.xpos{iSz} = xpos; % col vector 
S.stim.ypos{iSz} = ypos; % col vector 
S.stim.x_grid{iSz} = x_grid; % matrix of x pos 
S.stim.y_grid{iSz} = y_grid; % matrix of y pos. 
end
bsl_mean = mean(rasters(:,wdwbsl,:), [2,3]);
bsl_sem = sem(mean(rasters(:,wdwbsl,:), [2]), 3);
rsp_vec = squeeze(mean(rasters(:,wdwevk,:),[2])); 
S.psth.bsl_mean = bsl_mean; % baseline is the same across sizes (assume). nChan by 1 
S.psth.bsl_sem  = bsl_sem;  % nChan by 1 
S.psth.rsp_vec  = squeeze(mean(rasters(:,wdwevk,:),[2])); % nChan by nTrial 
re_tmp = regexp(meta.ephysFN,"(?<Animal>.*)-(?<datestr>.*)-",'names');
Animal = string(re_tmp.Animal);
date_str = string(re_tmp.datestr);
S.Animal = Animal;
S.meta.datestr = date_str;
S.meta.datetime = datetime(date_str,'InputFormat','ddMMyyyy');

act_col_vec = cellfun(@(A)reshape(A,[],1), S.psth.act_col,'Unif',0);
act_col_vec = cat(1, act_col_vec{:});
allF = zeros(size(rasters,1),1);
allF_P = zeros(size(rasters,1),1);
allT = zeros(size(rasters,1),1);
allT_P = zeros(size(rasters,1),1);
for iCh = 1:size(rasters,1) % Go through all channels see if there is some modulation (w/ ANOVA)
    [~,p,~,TST] = ttest(rsp_vec(iCh,:));
    allT(iCh,1) = TST.tstat;
    allT_P(iCh,1) = p;
    STS = anova_cells(cellfun(@(act)act(iCh,:),act_col_vec,'uni',0));
    allF(iCh,1) = STS.F;
    allF_P(iCh,1) = STS.F_P ; 
end
S.stats.F = allF;
S.stats.F_P = allF_P;
S.stats.T = allT;
S.stats.T_P = allT_P;
%%
% gprMdl_arr = [];
% for iCh = 1:numel(unit_num_arr)
% X = [xy_all,size_all/pixperdeg];
% y = S.psth.rsp_vec(iCh,:)'-bsl_mean(iCh);
% gprMdl = fitrgp(X,y);
% gprMdl_arr = [gprMdl_arr,gprMdl];
% end
%% Synthesized information across sizes
% Save stat to the folder
if P.collectstat, S_col = [S_col,S]; end
RFStat = S;
save(fullfile(savedir, compose("%s_%s_RFStat.mat",Animal,datestr(S.meta.datetime,"yyyymmdd"))), "RFStat");
%%
if P.plotEachChan
figure;plot(psthbest_mean')
for iCh = 1:numel(unit_num_arr)
pixperdeg = 40;
X = [xy_all,size_all/pixperdeg];
y = squeeze(mean(rasters(iCh,wdwevk,:),[2]))-bsl_mean(iCh);
gprMdl = fitrgp(X,y);
%%
Xq = -8:0.2:8; Yq = -8:0.2:8;
[XX,YY] = meshgrid(Xq,Yq);
figure(2);clf;T=tiledlayout(3,3,'tilespac','compact');
set(2,'pos',[1000   193   955   912])
for iSz = 1:3
sz = uniqsize_deg_all(iSz);
nexttile(T,iSz)
imagesc(S.stim.xpos{iSz},S.stim.ypos{iSz},S.psth.score_mean{iSz}(:,:,iCh))
axis image;set(gca,'YDir','normal');xlim([min(Xq),max(Xq)]);ylim([min(Yq),max(Yq)])
title(compose("size %.1f deg",sz))
pred_score = gprMdl.predict([reshape(XX,[],1),reshape(YY,[],1),ones(numel(XX),1)*sz]);
pred_rfmat = reshape(pred_score,size(XX));
nexttile(T,3+iSz)
imagesc(Xq,Yq,pred_rfmat);
title(compose("size %.1f deg",sz))
axis image;set(gca,'YDir','normal')
nexttile(T,6+iSz)
shadedErrorBar(1:200,S.psth.psth_best_mean{iSz}(iCh,:),S.psth.psth_best_sem{iSz}(iCh,:))
title(compose("size %.1f deg",sz))
end
title(T,compose("%s Unit %s",date_str,unit_name_arr(iCh)))
axs = get(T,'Children');
AlignAxisLimits({axs(1),axs(4),axs(7)});
AlignAxisCLimits(axs(2:3:9));
AlignAxisCLimits(axs(3:3:9));
pause
end
end
end 