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
% Get View Time and the Frame Ticks for marking
viewTime = Trials_mov.TrialRecord.User.viewTime;
vid = VideoReader(fullfile(meta_mov.stimuli, Trials_mov.imageName{1}+".avi"));
vidLenth = vid.Duration * 1000;
fps = vid.FrameRate; 
frameN = int32(viewTime / (1000/fps) + 1);
frameTick = (1000/fps)*double(0:frameN); 
% Get movie names and sort the trials into movies
movnm = string(unique(Trials_mov.imageName));
[movnm_sorted, sortedProp, sortIdx] = sortMovieNames(movnm);
%%
mov_idx_arr = {}; psth_means = {}; psth_sems = {};
mov_idx_arr{1} = arrayfun(@(mv)find(contains(Trials_new{1}.imageName, mv)),movnm_sorted,'Uni',0);
mov_idx_arr{2} = arrayfun(@(mv)find(contains(Trials_new{2}.imageName, mv)),movnm_sorted,'Uni',0);
psth_means{1} = cellfun(@(idx) mean(rasters_new{1}(:,:,idx),3), mov_idx_arr{1}, 'Uni', 0); 
psth_means{2} = cellfun(@(idx) mean(rasters_new{2}(:,:,idx),3), mov_idx_arr{2}, 'Uni', 0); 
psth_sems{1} = cellfun(@(idx) std(rasters_new{1}(:,:,idx),1,3)/sqrt(numel(idx)), mov_idx_arr{1}, 'Uni', 0); 
psth_sems{2} = cellfun(@(idx) std(rasters_new{2}(:,:,idx),1,3)/sqrt(numel(idx)), mov_idx_arr{2}, 'Uni', 0); 
%%
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
pvalmat = []; MvCorrmat = [];
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
    [MvCorrmat(iMv1,iMv2), pvalmat(iMv1,iMv2)] = corr(movmean(psth1(101:2250),kerL,2)',movmean(psth2(101:2250),kerL,2)');
    end
end
[diag(MvCorrmat), diag(pvalmat)]
MvCorrmat_nod = MvCorrmat + diag(nan(1,12));
ttest2(MvCorrmat_nod(:),diag(MvCorrmat))
figure(2);clf;
imagesc(MvCorrmat);
colorbar();axis image
pause
end
%%
corrWdw = 101:2250; kerL = 21; 
for chid1 = 1:numel(meta_new{1}.spikeID) % go through channel id list in experiment 1 and find the same channel x unit list in exp 2. 
    iCh = meta_new{1}.spikeID(chid1); iU = meta_new{1}.unitID(chid1); 
    chid2 = find((meta_new{2}.spikeID == iCh) & (meta_new{2}.unitID == iU));
    if isempty(chid1) || isempty(chid2), continue; end
    for iMv1 = 1:12
        for iMv2 = 1:12
        psth1 = psth_means{1}{iMv1}(chid1,:);
        psth2 = psth_means{2}{iMv2}(chid2,:); 
        [MvCorrmat(iMv1,iMv2), pvalmat(iMv1,iMv2)] = corr(movmean(psth1(corrWdw),kerL,2)',movmean(psth2(corrWdw),kerL,2)');
        end
    end
    [diag(MvCorrmat), diag(pvalmat)]
    MvCorrmat_nod = MvCorrmat + diag(nan(1,size(MvCorrmat,2)));
    ttest2(MvCorrmat_nod(:),diag(MvCorrmat),'tail','right')
    figure(2);clf; 
    imagesc(MvCorrmat);
    colorbar();axis image
    set(gca, 'TickLabelInterpreter', 'none')
    yticks(1:12);yticklabels(movnm_sorted)
    title(compose("%s %s\n%s Chan %d U%d", meta_new{1}.ephysFN, meta_new{2}.ephysFN, Animal, iCh, iU))
    pause
end
%%
MvCorrmat_nod = MvCorrmat + diag(nan(1,12));
ttest2(MvCorrmat_nod(:),diag(MvCorrmat),'tail','right')
%%
errornums = arrayfun(@(B)B.TrialError, Trials_new{1}.B);
movtr_mask = arrayfun(@(B)numel(B.TaskObject.Attribute)>1, Trials_new{1}.B);
block_arr = arrayfun(@(B)B.Block, Trials_new{1}.B);
for iblc = 1:max(block_arr)

end
%%
blcidx = arrayfun(@(iblc)find(errornums == 0 & block_arr == iblc & movtr_mask), 1:max(block_arr), 'uni', 0);

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