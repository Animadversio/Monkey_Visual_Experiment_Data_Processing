% Generate image sequence montage, RF montage and etc  for the paper and grant
RFdir = "O:\RFstats";
refdir = "N:\Stimuli\2020-CosineEvol\RefCollection";
expdir="O:\Evol_Cosine\2021-03-02-Alfa-03-MSE_ITV4";
load(fullfile(expdir, 'expStat.mat'),'expStat')
load(fullfile(expdir, 'ExpImMatchStat.mat'),'ExpImMatch')
expday = datetime(expStat.meta.expControlFN(1:6),'InputFormat','yyMMdd');
load(fullfile(RFdir,compose("Alfa_%s_RFStat.mat",datestr(expday,"yyyymmdd"))),'RFStat')
%% Image Evolution Traj
rank_statnm = "MSE_All";
grid_size = [3,13];%[4,9]
imgsize = [128,128];
scorevec = expStat.scores.(rank_statnm+"_vec");
% rank_statnm = "squ_dist";
% scorevec = -ExpImMatch.squ.fulldistvec;
expStat.evol.row_gen;
expStat.imageName;
block_num = max(expStat.evol.block_arr);
reprImgIdx = [];
reprImgs = {};
for blocki = 1:block_num-1
    [sortScore, sortIdx] = sort(scorevec(expStat.evol.gen_idx_seq{blocki}),'descend');
    % [maxScore_nat, maxIdx_nat] = sort(scorevec(expStat.evol.gen_idx_seq{blocki}),'descend');
    % score_mean = mean(imgscore_vec(gen_idx_seq{blocki}));
    bestGlobIdx = expStat.evol.gen_idx_seq{blocki}(sortIdx(1));
    reprImgIdx(blocki) = bestGlobIdx;
    reprImgs{blocki} = imread(fullfile(expStat.meta.stimuli, expStat.imageName(bestGlobIdx)+".bmp"));
end
%%
h=figure(1);
tiledimg = imtile(reprImgs,'BorderSize',[2,2],'GridSize',grid_size,'ThumbnailSize',imgsize);
imshow(tiledimg)
title(compose("%s Exemplar Image (ranked by %s)",expStat.meta.fdrnm,rank_statnm),'interpreter','none')
saveas(h,fullfile(expdir,compose("ReprImg_rank_%s.png",rank_statnm)))
imwrite(tiledimg,fullfile(expdir,compose("ReprImg_rank_%s_tile.png",rank_statnm)))
%%
tiledimg = imtile(reprImgs(1:4:end),'BorderSize',[2,2],'GridSize',[1,10],'ThumbnailSize',imgsize);
imwrite(tiledimg,fullfile(expdir,compose("ReprImg_rank_%s_tile_sparse.png",rank_statnm)))
%%
tiledimg = imtile(reprImgs(1:4:end),'BorderSize',[2,2],'GridSize',[2,6],'ThumbnailSize',[192,192]);
imwrite(tiledimg,fullfile(expdir,compose("ReprImg_rank_%s_tile_sparse2.png",rank_statnm)))
%% Last block top half 
rank_statnm = "MSE_V4";
grid_size = [3,6];%[4,9]
imgsize = [192,192];
scorevec = expStat.scores.(rank_statnm+"_vec");
blocki = block_num-1;
[sortScore, sortIdx] = sort(scorevec(expStat.evol.gen_idx_seq{blocki}),'descend');
% [maxScore_nat, maxIdx_nat] = sort(scorevec(expStat.evol.gen_idx_seq{blocki}),'descend');
% score_mean = mean(imgscore_vec(gen_idx_seq{blocki}));
bestHalfGlobIdx = expStat.evol.gen_idx_seq{blocki}(sortIdx(1:16));
finalreprImg = arrayfun(@(idx)imread(fullfile(expStat.meta.stimuli, expStat.imageName(idx)+".bmp")),...
    bestHalfGlobIdx,'uni',0);
tiledimg = imtile(finalreprImg,'BorderSize',[2,2],'GridSize',grid_size,'ThumbnailSize',imgsize);
figure;imshow(tiledimg)
imwrite(tiledimg,fullfile(expdir,compose("FinalImg_rank_%s_tile.png",rank_statnm)))
% reprImgs{blocki} = imread(fullfile(expStat.meta.stimuli, expStat.imageName(bestGlobIdx)+".bmp"));
%% Show the color map 
clrseq = brewermap(block_num-1,'Spectral');
brewermap_view(block_num-1,'Spectral')

%% Find masks for real images
imgsize = expStat.evol.imgsize;
imgpos = expStat.evol.imgpos;
baseMask = expStat.targ.baseMask;
chan_arr = (1:64)';
targmask = baseMask&(chan_arr<=32|chan_arr>=49);
iCh_mat = zeros(size(baseMask));
for iCh = 1:numel(RFStat.unit.chan_num_arr)
    chan = RFStat.unit.chan_num_arr(iCh);
    unit = RFStat.unit.unit_num_arr(iCh);
    iCh_mat(chan,unit+1)=iCh;
end
targmsk_chans = iCh_mat(targmask)';
if any(targmsk_chans==0)
   fprintf("Warning!!\n") 
end
targmsk_chans(targmsk_chans==0)=[];
%%
%%
uniqsize_deg = unique(RFStat.stim.size_deg)';
activ_chans = find(RFStat.stats.F_P<0.001)';
area_colormap = containers.Map({'V1','V4','IT'},{[1,0,0],[0,1,0],[0,0,1]});
Xq = -8:0.2:8; Yq = -8:0.2:8;
[XX,YY] = meshgrid(Xq,Yq);
Nsize = numel(uniqsize_deg);
figure(2);clf;T=tiledlayout(1,Nsize,'tilespac','tight','padding','compact');
set(2,'pos',[1000   193   306*Nsize  320   ])
for iCh = targmsk_chans%activ_chans%
chan_num = RFStat.unit.chan_num_arr(iCh);
area = area_map(chan_num);
clr = area_colormap(area);
pixperdeg = 40;
X = [RFStat.stim.xy_all,RFStat.stim.size_all/pixperdeg];
y = RFStat.psth.rsp_vec(iCh,:)-RFStat.psth.bsl_mean(iCh);
gprMdl = fitrgp(X,y);
for iSz = 1:Nsize
sz = uniqsize_deg(iSz);
% imagesc(S.stim.xpos{iSz},S.stim.ypos{iSz},S.psth.score_mean{iSz}(:,:,iCh))
% axis image;set(gca,'YDir','normal');xlim([min(Xq),max(Xq)]);ylim([min(Yq),max(Yq)])
% title(compose("size %.1f deg",sz))
pred_score = gprMdl.predict([reshape(XX,[],1),reshape(YY,[],1),ones(numel(XX),1)*sz]);
pred_rfmat = reshape(pred_score,size(XX));
nexttile(T,iSz);hold on
% contour(Xq,Yq,pred_rfmat,[0.5 0.5]*max(pred_rfmat,[],'all'),'linecolor',[clr],'linewidth',0.5);

[C,hcf] = contourf(Xq,Yq,pred_rfmat,[0.5 0.5]*max(pred_rfmat,[],'all'),'facecolor',[clr]);
allH = allchild(hcf);
valueToHide = 0.5*max(pred_rfmat,[],'all');
patchValues = cell2mat(get(allH,'UserData'));
patchesToHide = patchValues == valueToHide;
set(allH(patchesToHide),'FaceColor',[clr],'FaceAlpha',0.05);
end
% pause7
end
% Set the axis limit to the size of image. 
for iSz = 1:Nsize
    nexttile(T,iSz);hold on;axis image
    xlim(imgpos(1)+[-0.5, +0.5]*imgsize)
    ylim(imgpos(2)+[-0.5, +0.5]*imgsize)
    title(compose("Feature size %.1f deg",uniqsize_deg(iSz)))
end
saveallform(expdir,"Obj_PopulRF_maps_fill")
% saveallform(expdir,"All_PopulRF_maps")