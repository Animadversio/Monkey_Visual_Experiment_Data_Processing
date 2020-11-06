%% Movie Movie Correlation
Animal = "Alfa";Set_Path;
setMatlabTitle("Movie Movie Correlation")
%%
ftr = find(contains(ExpRecord.expControlFN,"Alfa") & contains(ExpRecord.Exp_collection,"Movie"));
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr(6:7), Animal);
%%
meta_mov = meta_new{1};
Trials_mov = Trials_new{1};

wdw = meta_mov.rasterWindow;
% Get View Time and the Frame Ticks for marking plot
viewTime = Trials_mov.TrialRecord.User.viewTime;
vid = VideoReader(fullfile(meta_mov.stimuli, Trials_mov.imageName{1}+".avi"));
vidLenth = vid.Duration * 1000;
fps = vid.FrameRate; 
frameN = int32(viewTime / (1000/fps) + 1);
frameTick = (1000/fps)*double(0:frameN); 
% Get movie names and sort the trials into movies
movnm = string(unique(Trials_mov.imageName));
[movnm_sorted, sortedProp, sortIdx] = sortMovieNames(movnm);
% Get the trial index for movies and the psth and sem for it 
mov_idx_arr = {}; psth_means = {}; psth_sems = {};
mov_idx_arr{1} = arrayfun(@(mv)find(contains(Trials_new{1}.imageName, mv)),movnm_sorted,'Uni',0);
mov_idx_arr{2} = arrayfun(@(mv)find(contains(Trials_new{2}.imageName, mv)),movnm_sorted,'Uni',0);
psth_means{1} = cellfun(@(idx) mean(rasters_new{1}(:,:,idx),3), mov_idx_arr{1}, 'Uni', 0); 
psth_means{2} = cellfun(@(idx) mean(rasters_new{2}(:,:,idx),3), mov_idx_arr{2}, 'Uni', 0); 
psth_sems{1} = cellfun(@(idx) std(rasters_new{1}(:,:,idx),1,3)/sqrt(numel(idx)), mov_idx_arr{1}, 'Uni', 0); 
psth_sems{2} = cellfun(@(idx) std(rasters_new{2}(:,:,idx),1,3)/sqrt(numel(idx)), mov_idx_arr{2}, 'Uni', 0); 
%% Visualize PSTH for the same channel across days. 
kerL = 21;
iMv = 3;
for iCh = 1:64
% iCh = 34; iU = 1; 
idCh1 = find((meta_new{1}.spikeID == iCh) & (meta_new{1}.unitID == iU));
psth1 = psth_means{1}{iMv}(idCh1,:); 
% iCh = 34; iU = 1; 
idCh2 = find((meta_new{2}.spikeID == iCh) & (meta_new{2}.unitID == iU));
psth2 = psth_means{2}{iMv}(idCh2,:); 
if isempty(idCh1) || isempty(idCh2), continue; end
corr(movmean(psth1(251:2200),kerL,2)',movmean(psth2(251:2200),kerL,2)')
%
figure(1);clf;hold on
shadedLine(wdw, psth1, psth_sems{1}{iMv}(idCh1,:), kerL, 'LineWidth', 1.5)
shadedLine(wdw, psth2, psth_sems{2}{iMv}(idCh2,:), kerL, 'LineWidth', 1.5)
% plot(wdw(1)+1:wdw(2),movmean(psth1,kerL,2))
% plot(wdw(1)+1:wdw(2),movmean(psth2,kerL,2))
xlim(wdw);title(num2str(iCh))
pause
end
%%
pvalmat = []; MvCorrmat_x = [];
for iCh = 1:64 
kerL = 21;
% iCh = 9; iU = 1; 
idCh1 = find((meta_new{1}.spikeID == iCh) & (meta_new{1}.unitID == iU));
% iCh = 9; iU = 1; 
idCh2 = find((meta_new{2}.spikeID == iCh) & (meta_new{2}.unitID == iU));
if isempty(idCh1) || isempty(idCh2), continue; end
for iMv1 = 1:12
    for iMv2 = 1:12
    psth1 = psth_means{1}{iMv1}(idCh1,:);
    psth2 = psth_means{2}{iMv2}(idCh2,:); 
    [MvCorrmat_x(iMv1,iMv2), pvalmat(iMv1,iMv2)] = corr(movmean(psth1(101:2250),kerL,2)',movmean(psth2(101:2250),kerL,2)');
    end
end
[diag(MvCorrmat_x), diag(pvalmat)]
MvCorrmat_nod = MvCorrmat_x + diag(nan(1,12));
ttest2(MvCorrmat_nod(:),diag(MvCorrmat_x))
figure(2);clf;
imagesc(MvCorrmat_x);
colorbar();axis image
pause
end
%%
unit_name_arr = generate_unit_labels_new(meta_new{1}.spikeID,meta_new{1}.unitID);
%% Correlation Matrix for the neural code 
figdir = "E:\OneDrive - Washington University in St. Louis\MovieDynamics\2020-10-27+29-Alfa_Stability";
corrWdw = 101:2250; kerL = 21;
xcorrStats = repmat(struct(),1,numel(meta_new{1}.spikeID));
for chid1 = 1:numel(meta_new{1}.spikeID) % go through channel id list in experiment 1 and find the same channel x unit list in exp 2. 
    chanstr = unit_name_arr(chid1);
    iCh = meta_new{1}.spikeID(chid1); iU = meta_new{1}.unitID(chid1); 
    chid2 = find((meta_new{2}.spikeID == iCh) & (meta_new{2}.unitID == iU));
    if isempty(chid1) || isempty(chid2), continue; end
    MvCorrmat_x = nan(12); pvalmat_x = nan(12);
    for iMv1 = 1:12
        for iMv2 = 1:12
        psth1 = psth_means{1}{iMv1}(chid1,:);
        psth2 = psth_means{2}{iMv2}(chid2,:); 
        [MvCorrmat_x(iMv1,iMv2), pvalmat_x(iMv1,iMv2)] = corr(movmean(psth1(corrWdw),kerL,2)',movmean(psth2(corrWdw),kerL,2)');
        end
    end
    [diag(MvCorrmat_x), diag(pvalmat_x)]
    xcorrStats(chid1).corrmat = MvCorrmat_x;
    xcorrStats(chid1).pvalmat = pvalmat_x;
    MvCorrmat_nod = MvCorrmat_x + diag(nan(1,size(MvCorrmat_x,2)));
    [H,P,~,STAT] = ttest2(diag(MvCorrmat_x),MvCorrmat_nod(:),'tail','right'); % One sided as we hypothesize diag to be larger
    xcorrStats(chid1).diag_t.H = H;
    xcorrStats(chid1).diag_t.P = P;
    xcorrStats(chid1).diag_t.T = STAT.tstat;
    xcorrStats(chid1).kerL = kerL;
    xcorrStats(chid1).corrWdw = corrWdw;
    figure(2);clf; 
    imagesc(MvCorrmat_x);
    colorbar();axis image;set(gca,'pos',[0.2054    0.0722    0.7617    0.8082]) % remove the margin
    set(gca, 'TickLabelInterpreter', 'none')
    yticks(1:12);yticklabels(movnm_sorted)
    ylabel(meta_new{1}.ephysFN); xlabel(meta_new{2}.ephysFN); % First index is movie in exp 1
    title(compose("%s %s  %s Chan %d Unit%d\n Median Corr %.3f(%.1e) Kernel L %d Window [%d,%d]\n Diag T = %.2f(%.1e)", ...
           meta_new{1}.ephysFN, meta_new{2}.ephysFN, Animal, iCh, iU,...
           median(MvCorrmat_x(:)), median(pvalmat_x(:)), ... median of corrmat and p value
           kerL, corrWdw(1)+wdw(1), corrWdw(end)+wdw(1), ... Smooth and window info
           STAT.tstat, P)) % t statistics of diagonal versus all corr. 
    saveas(2,fullfile(figdir,compose("movie_code_xcorrmat_%s.png", chanstr)))
    savefig(2,fullfile(figdir,compose("movie_code_xcorrmat_%s.fig", chanstr)))
    pause
end
%%
save(fullfile(figdir,"PSTHcorrStats.mat"),'xcorrStats') % ,'-append'
%%
P_col = [];
for S = xcorrStats
if ~isempty(S.diag_t), P_col(end+1) = S.diag_t.P; end
end
fprintf("%d/%d unit are has a significantly larger diagonal correlation by p thresh %.4f\n",sum(P_col<0.05),numel(P_col),0.05)
fprintf("%d/%d unit are has a significantly larger diagonal correlation by p thresh %.4f\n",sum(P_col<0.01),numel(P_col),0.01)
fprintf("%d/%d unit are has a significantly larger diagonal correlation by p thresh %.4f\n",sum(P_col<0.001),numel(P_col),0.001)
%% Within day consistency estimated by bootstrapping. 
kerL = 21;
selfcorrStats1 = repmat(struct(),1,0);
selfcorrStats2 = repmat(struct(),1,0);
for chid1 = 1:numel(meta_new{1}.spikeID) % go through channel id list in experiment 1 and find the same channel x unit list in exp 2. 
chanstr = unit_name_arr(chid1);
iCh = meta_new{1}.spikeID(chid1); iU = meta_new{1}.unitID(chid1); 
chid2 = find((meta_new{2}.spikeID == iCh) & (meta_new{2}.unitID == iU));
if isempty(chid1) || isempty(chid2), continue; end
% Compute the bootstrapping consistency
MvCorrmat_1_col = {};
MvCorrmat_2_col = {};
for i = 1:50 % bootstrap reps
% 2 bootstrapped mean of raster see the correlation
psth1_all_rnd1 = cell2mat(cellfun(@(idx) mean(rasters_new{1}(chid1,:,randsample(idx,numel(idx),1)),3), ...
                mov_idx_arr{1}, 'Uni', 0));
psth1_all_rnd2 = cell2mat(cellfun(@(idx) mean(rasters_new{1}(chid1,:,randsample(idx,numel(idx),1)),3), ...
                mov_idx_arr{1}, 'Uni', 0));
psth2_all_rnd1 = cell2mat(cellfun(@(idx) mean(rasters_new{2}(chid2,:,randsample(idx,numel(idx),1)),3), ...
                mov_idx_arr{2}, 'Uni', 0));
psth2_all_rnd2 = cell2mat(cellfun(@(idx) mean(rasters_new{2}(chid2,:,randsample(idx,numel(idx),1)),3), ...
                mov_idx_arr{2}, 'Uni', 0));
[MvCorrmat_1, pvalmat_1] = corr(movmean(psth1_all_rnd1(:,corrWdw),kerL,2)',...
                                movmean(psth1_all_rnd2(:,corrWdw),kerL,2)');
[MvCorrmat_2, pvalmat_2] = corr(movmean(psth2_all_rnd1(:,corrWdw),kerL,2)',...
                                movmean(psth2_all_rnd2(:,corrWdw),kerL,2)');
MvCorrmat_1_col{end+1} = MvCorrmat_1;
MvCorrmat_2_col{end+1} = MvCorrmat_2;
end
MvCorrmat_1_avg = mean(cat(3,MvCorrmat_1_col{:}),3);
MvCorrmat_2_avg = mean(cat(3,MvCorrmat_2_col{:}),3);
selfcorrStats1(chid1).corrmat = MvCorrmat_1_avg;
selfcorrStats2(chid2).corrmat = MvCorrmat_2_avg;
% figure(3);clf; set(3,'pos',[50         400        1120         540])
% subplot(121)
% imagesc(MvCorrmat_1_avg);
% colorbar();axis image;set(gca, 'position', [0.183   0.073  0.34  0.76])
% set(gca, 'TickLabelInterpreter', 'none')
% yticks(1:12);yticklabels(movnm_sorted)
% xlabel("bootstrap 1");ylabel("bootstrap 2");
% title(compose("%s\nMedian Corr %.3f(%.1e)",meta_new{1}.ephysFN,...
%                 median(MvCorrmat_1_avg(:)), median(pvalmat_1(:)))); 
% subplot(122)
% imagesc(MvCorrmat_2_avg);
% colorbar();axis image;set(gca, 'position', [0.602   0.073  0.34  0.76])
% xlabel("bootstrap 1");ylabel("bootstrap 2");yticklabels([])
% title(compose("%s\nMedian Corr %.3f(%.1e)",meta_new{2}.ephysFN,...
%                 median(MvCorrmat_2_avg(:)), median(pvalmat_2(:)))); % First index is movie in exp 1
% suptitle(compose("%s Chan %d Unit%d Trial Bootstrapping Mean Consistency\n  Kernel L %d Window [%d,%d]", ...
%        Animal, iCh, iU,...
%        kerL, corrWdw(1)+wdw(1), corrWdw(end)+wdw(1))) ... %Smooth and window info
%        
% saveas(3,fullfile(figdir,compose("movie_code_self_corrmat_%s.png", chanstr)))
% savefig(3,fullfile(figdir,compose("movie_code_self_xcorrmat_%s.fig", chanstr)))
end
meta12 = meta_new;
save(fullfile(figdir,"PSTHcorrStats.mat"),'selfcorrStats1','selfcorrStats2','meta12','-append') % 
%% For behavioral analysis see that script Trial_Behav analysis

%% 
function shadedLine(wdw, psthmov, errmov, kerL, varargin)
if isempty(kerL), kerL = 21; end
psthmov = movmean(psthmov,kerL);
errmov = movmean(errmov,kerL);
plot(wdw(1)+1:wdw(2), psthmov, varargin{:})
patch([wdw(1)+1:wdw(2),fliplr(wdw(1)+1:wdw(2))], [psthmov-errmov,fliplr(psthmov+errmov)], 'k','FaceAlpha',0.15,'EdgeColor','none');
end 

function [sortedMovnm, sortedRows, sortIdx] = sortMovieNames(movnm)
% Sort the movie names in **lexicoGraphical order**, by space and by eigen idx. 
% This assumes the names to have the structure like "class_eig17_shortshort"
eigi_cell = cellfun(@(eigi){str2double(eigi{1})},regexp(movnm,"_eig(\d*)_",'tokens'));
space_cell = cellfun(@(eigi){eigi{1}{1}},regexp(movnm,"(.*)_eig(\d*)_",'tokens'));
[sortedRows, sortIdx] = sortrows([space_cell,eigi_cell]);
sortedMovnm = movnm(sortIdx);
end