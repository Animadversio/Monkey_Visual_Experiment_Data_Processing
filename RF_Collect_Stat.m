%% Collect the PSTH and Score in the more raw form
Animal = "Alfa"; Set_Path;
%expftr = (contains(ExpRecord.expControlFN,"200319"));
expftr = (contains(ExpRecord.Exp_collection,"Manifold") &...
            contains(ExpRecord.expControlFN, "rf"));%&...
            %ExpRecord.Expi > 0);
rowis = find(expftr);
[meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(rowis,Animal);
%%
% addpath D:\Lab_Codeshare\project_rfmap
RFStats = repmat(struct(),length(meta_new),1);
%% Go through the exps 
tic
Mapi = 1; % seperate the different size RFMapping exps from the start
for Triali = 1:length(meta_new)
Trials = Trials_new{Triali}; 
rasters = rasters_new{Triali}; 
meta = meta_new{Triali};
%%
fprintf("Processing RFMap %d\n",Triali)
xy = Trials.TrialRecord.User.xy;
size_all = Trials.width(:,1);
uniqsize_all = sort(unique(size_all));
uniqsize_deg_all = sort(Trials.TrialRecord.User.uniqueSizes_deg); % in deg
nImSize = length(uniqsize_deg_all); % before its always 1, recently, we use multiple image sizes! 

for iSz = 1:nImSize
imgsize = uniqsize_all(iSz); 
imgsize_deg = uniqsize_deg_all(iSz); 
uniqpos = unique(xy(size_all==imgsize, :), 'rows');%Trials.TrialRecord.User.uniquePositions;
xpos = unique(uniqpos(:,1));
ypos = unique(uniqpos(:,2));
nPos = size(uniqpos,1);
% This is a key assumption! The position could be factorized into a square grid. 
assert(nPos==length(xpos)*length(xpos)); % if there is multiple sets of xpos this will fail
posgrid = reshape([uniqpos, (1:nPos)'],[length(xpos),length(ypos),3]);
% sort neural responses by getting the trial indices
idx_grid = cell(length(xpos),length(ypos)); % trial index
iPos_grid = cell(length(xpos),length(ypos)); % position index for grid
for iPos = 1:nPos
    tPos = uniqpos(iPos,:) ;
    idx = find( xy(:,1) == tPos(1) & xy(:,2) == tPos(2) ) ;
    gridi = find(xpos == tPos(1));
    gridj = find(ypos == tPos(2));
    % posgrid(gridi,gridj,1) == tPos(1)
    % posgrid(gridi,gridj,2) == tPos(2)
    idx_grid{gridi, gridj} = idx;
    iPos_grid{gridi, gridj} = repmat([iPos],1,length(idx));
end
nrep = cellfun(@length, idx_grid);
psth_mean = cellfun(@(idx) mean(rasters(:,:,idx),3), idx_grid, 'UniformOutput', false);
psth_sem = cellfun(@(idx) std(rasters(:,:,idx),0,3) / sqrt(length(idx)), idx_grid, 'UniformOutput', false);
%% Compute the score (psth time averaged)
Windows.late = 1:49;
Windows.late = 50:200;
score_all = cellfun(@(idx) squeeze(mean(rasters(:,Windows.late,idx),[2])), idx_grid, 'UniformOutput', false);
idx_all_cat = cell2mat(reshape(iPos_grid,1,[]));
score_all_cat = cell2mat(reshape(score_all,1,[]));
score_mean = cellfun(@(idx) mean(rasters(:,Windows.late,idx),[2,3]), idx_grid, 'UniformOutput', false);
score_mean = cell2mat(reshape(score_mean,[1,size(score_mean)]));
score_sem = cellfun(@(idx) std(mean(rasters(:,Windows.late,idx),[2]),0,3)/sqrt(length(idx)), idx_grid, 'UniformOutput', false);
score_sem = cell2mat(reshape(score_sem,[1,size(score_sem)]));
%% Compute ANOVA F, T percentile
resp = squeeze( mean( rasters(:,Windows.late,:) , 2 ) - ...
    mean( rasters(:,Windows.early,:), 2 ) );

allPercentiles = [99.9 99 97.5 95 90 75 50]';
allP = zeros(size(rasters,1),1);
allF = zeros(size(rasters,1),1);
allT = zeros(size(rasters,1),1);
allP_t = zeros(size(rasters,1),1);
prctile_per_chan = zeros(size(rasters,1),length(allPercentiles));
for iCh = 1:size(rasters,1) % Go through all channels see if there is some modulation (w/ ANOVA)
    [H,P,CI,STATS] = ttest(resp(iCh,:));
    allT(iCh,1) = STATS.tstat;
    allP_t(iCh,1) = P;
    [P,ANOVATAB,STATS] = anova1( score_all_cat(iCh,:), idx_all_cat ,'off') ;
    allP(iCh,1) = P ; 
    allF(iCh,1) = nan;
    if ~isempty(ANOVATAB{2,5}) % Get F value
        allF(iCh,1) = ANOVATAB{2,5} ; 
    end
    prctile_per_chan(iCh,:) = prctile(score_all_cat(iCh,:), allPercentiles)';
end
allP_sr = rowfun(@signrank,table(resp)) ;
%% Save stats to the dict 
RFStats(Mapi).stim.offset = Trials.TrialRecord.User.offsets;
RFStats(Mapi).stim.nSize = nImSize;
RFStats(Mapi).stim.imgsize = imgsize;
RFStats(Mapi).stim.imgsize_deg = imgsize_deg;
RFStats(Mapi).stim.imgname = string(unique(Trials.imageName)); 
RFStats(Mapi).stim.uniqpos = uniqpos;
RFStats(Mapi).stim.xpos = unique(uniqpos(:,1));
RFStats(Mapi).stim.ypos = unique(uniqpos(:,2));
RFStats(Mapi).psth.nrep = nrep;
RFStats(Mapi).psth.psth_mean = psth_mean;
RFStats(Mapi).psth.psth_sem = psth_sem;
RFStats(Mapi).psth.score_mean = score_mean;
RFStats(Mapi).psth.score_sem = score_sem;
RFStats(Mapi).stats.anovaP = allP;
RFStats(Mapi).stats.F = allF;
RFStats(Mapi).stats.T = allT;
RFStats(Mapi).stats.ttestP = allP_t;
RFStats(Mapi).stats.signrankP = allP_sr;
RFStats(Mapi).stats.prctile = prctile_per_chan;
RFStats(Mapi).meta = meta;
date_str = string(regexp(meta.ephysFN,"-(.*)-",'tokens'));
RFStats(Mapi).meta.datestr = date_str;
RFStats(Mapi).meta.datetime = datetime(date_str,'InputFormat','ddmmyyyy');
Mapi = Mapi + 1;
end
end
toc
%%
savepath = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Mat_Statistics";
save(compose("D:\\%s_Manif_RFstats.mat", Animal), 'RFStats','-v6')
save(fullfile(savepath, compose("%s_Manif_RFstats.mat", Animal)), 'RFStats','-v6')
%% Carlos' original code
% note seems this function doesn't take the image size into account 
[allRFs,StatsRF] = Project_RfMap_computeRFs(Trials.TrialRecord.User.xy, rasters, meta);
%%

%%
figure;
colorbar()
for chi = 1:size(allRFs,3)
   imagesc(allRFs(:,:,chi));
   drawnow
   pause(0.1)
end
function eq=float_eq(A, B)
eq = abs(A(:) - B(:)) < 1E-7;
return 
end