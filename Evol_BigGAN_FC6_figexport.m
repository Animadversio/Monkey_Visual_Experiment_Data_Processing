%% Compare Evolved image in BigGAN and FC6 in batch
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
for Animal = ["Alfa"] %, "Beto"
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
end
%%
P = struct();
P.plot_img = true;
P.plot_traj = true;
rootdir = "E:\OneDrive - Harvard University\Evol_BigGAN_FC6";
sumdir = "E:\OneDrive - Harvard University\Evol_BigGAN_FC6\summary_figs";
% N:\Stimuli\2020-BigGAN\2020-07-23-Beto-01\2020-07-23-15-59-38
%%
failures = struct('iter', {}, 'str', {});
%%
for Expi = 7:numel(BFEStats)
try
BFES = BFEStats(Expi);
stimpath = BFES.meta.stimuli;
prefchan = BFES.evol.pref_chan(1);
area  = area_map(prefchan);
expdir = fullfile(rootdir,compose("%s_Exp%02d_Ch%02d",Animal,Expi,prefchan));
mkdir(expdir)
expstr = compose("%s Exp%02d PrefChan %02d (%s)\nThr1: %s x %s Thr2: %s x %s\n%s",Animal,Expi,prefchan,area,...
                 BFES.evol.space_cfg{1}{1},BFES.evol.optim_names(1),...
                 BFES.evol.space_cfg{2}{1},BFES.evol.optim_names(2), BFES.meta.ephysFN);
expstr = strrep(expstr,"_","-");
fprintf(expstr+"\n")
%% Get the responses
bsl_col = cellfun(@(P)squeeze(mean(P(:,1:45,:),[1,2])), BFES.evol.psth, 'uni',0);
act_col = cellfun(@(P)squeeze(mean(P(:,51:200,:),[1,2])), BFES.evol.psth, 'uni',0);
refbsl_col = cellfun(@(P)squeeze(mean(P(:,1:45,:),[1,2])), BFES.ref.psth, 'uni',0);
refact_col = cellfun(@(P)squeeze(mean(P(:,51:200,:),[1,2])), BFES.ref.psth, 'uni',0);
bsl_vec = cellfun(@mean, bsl_col);
sembsl_vec = cellfun(@sem, bsl_col);

meanact_vec = cellfun(@mean, act_col);
semact_vec = cellfun(@sem, act_col);
stdact_vec = cellfun(@std, act_col);

refact_vec = cellfun(@mean, refact_col);
refsem_vec = cellfun(@sem, refact_col);
refstd_vec = cellfun(@std, refact_col);

nGen = size(meanact_vec,2);
if P.plot_traj
%% Evol Trajectory mean
Cord = colororder;
figure(1);clf;hold off;set(1,'pos',[680   400   460   480])
for iThr = 1:2
    % plot(meanact_vec(iThr,1:end-1),'linestyle','-','color',Cord(iThr,:),'Linewidth',1.5);
    % shadedErrorBar([],refact_vec(iThr,1:end-1),refsem_vec(iThr,1:end-1),...
    %     'lineProps',{'Color',Cord(iThr+2,:),'Linewidth',1.,'linestyle','-.'},'patchSaturation',0.7)
    plot(refact_vec(iThr,1:end-1),'Color',Cord(iThr+2,:),'Linewidth',1,'linestyle','-')
    hold on 
    shadedErrorBar([],meanact_vec(iThr,1:end-1),stdact_vec(iThr,1:end-1),...
        'lineProps',{'Color',Cord(iThr,:),'Linewidth',2},'patchSaturation',0.3)
    hold on 
    plot(bsl_vec(iThr,1:end-1),'linestyle',':','color',Cord(iThr,:),'Linewidth',1.)
end
xlabel("Generations")
ylabel("Firing Rate (event/sec)")
title(compose("%s\nEvolution Trajectory",expstr))
% legend(["FC6","FC6 baseline","BigGAN","BigGAN baseline"],'location','best')
legend(["FC6","FC6 ref","FC6 baseline","BigGAN","BigGAN ref","BigGAN baseline"],'location','best')
xlim([0,nGen+1])
%% Evol Trajectory std
% Cord = colororder;
% figure(3);clf;hold off;set(3,'pos',[680   558   460   420])
% for iThr = 1:2
% plot(stdact_vec(iThr,1:end-1),'linestyle','-','color',Cord(iThr,:),'Linewidth',1.5);
% hold on
% end
% xlabel("Generations")
% ylabel("Firing Rate (event/sec)")
% title(compose("%s\nEvolution Trajectory",expstr))
% legend(["FC6","BigGAN"],'location','best')
%% Single trial trajectory
[genvec1, actvec1] = actcell2traj(act_col(1,1:end-1));
[genvec2, actvec2] = actcell2traj(act_col(2,1:end-1));
figure(2);clf;hold off;set(2,'pos',[1140   400    460   480])
for iThr = 1:2
    [genvec, actvec] = actcell2traj(act_col(iThr,1:end-1));
    scatter(genvec,actvec,'markerfacecolor',Cord(iThr,:),...
        'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
    hold on
    plot(bsl_vec(iThr,1:end-1),'color',Cord(iThr,:),...
        'linestyle',':','Linewidth',2)
end
legend(["FC6","FC6 baseline","BigGAN","BigGAN baseline"],'location','best')
xlabel("Generations")
ylabel("Firing Rate (event/sec)")
title(compose("%s\nEvolution Trajectory (single trial)",expstr))
xlim([0,nGen+1])

saveallform([expdir], compose("%s_Exp%02d_EvolTraj",Animal,Expi), 1)
saveallform([sumdir], compose("%s_Exp%02d_EvolTraj",Animal,Expi), 1, ["png"])
saveallform([expdir], compose("%s_Exp%02d_EvolTraj_sgtr",Animal,Expi), 2)
saveallform([sumdir], compose("%s_Exp%02d_EvolTraj_sgtr",Animal,Expi), 2, ["png"])
end


%% Preferred images across Evolution
if P.plot_img
% BFES.imageName(BFES.evol.idx_seq{2,3})
bestFCimg_col = {};
bestBGimg_col = {};
meanFCimg_col = {};
meanBGimg_col = {};
WtmeanFCimg_col = {};
WtmeanBGimg_col = {};
for igen = 1:size(BFES.evol.idx_seq,2)-1
%     [FCim_gen, FCimnames] = loadEvolImages(stimpath, 0, igen);
%     [BGim_gen, BGimnames] = loadEvolImages(stimpath, 1, igen);
    tic
    FCimnms = BFES.imageName(BFES.evol.idx_seq{1,igen});
    BGimnms = BFES.imageName(BFES.evol.idx_seq{2,igen});
    FCimfps = fullfile(stimpath, strcat(FCimnms,".bmp"));
    BGimfps = fullfile(stimpath, strcat(BGimnms,".bmp"));
    FCim_gen = cellfun(@(fp) imread(fp), FCimfps, "uni", 0);
    BGim_gen = cellfun(@(fp) imread(fp), BGimfps, "uni", 0);
    toc
    %%
    FCact = act_col{1,igen};
    BGact = act_col{2,igen};
    % meanFCimg = mean(cat(4, FCim_gen{:}) / 255.0, 4); % Buggy! this will qunatize pixel to 0, 1
    % meanBGimg = mean(cat(4, BGim_gen{:}) / 255.0, 4); % Buggy! this will qunatize pixel to 0, 1
    meanFCimg = mean(double(cat(4, FCim_gen{:})) / 255.0, 4);
    meanBGimg = mean(double(cat(4, BGim_gen{:})) / 255.0, 4);
    WtmeanFCimg = sum(double(cat(4, FCim_gen{:})) / 255.0 ...
                       .* reshape(FCact,1,1,1,[]), 4) / sum(FCact);
    WtmeanBGimg = sum(double(cat(4, BGim_gen{:})) / 255.0 ...
                       .* reshape(BGact,1,1,1,[]), 4) / sum(BGact);
    % uniformly averaged images
    meanFCimg_col{end+1} = meanFCimg;
    meanBGimg_col{end+1} = meanBGimg;
    % activation averaged images
    WtmeanFCimg_col{end+1} = WtmeanFCimg;
    WtmeanBGimg_col{end+1} = WtmeanBGimg;
    % get the max activation image, very noisy
    [~,maxidx] = max(FCact);
    bestFCimg_col{end+1} = FCim_gen{maxidx};
    [~,maxidx] = max(BGact);
    bestBGimg_col{end+1} = BGim_gen{maxidx};
    % figure(3);
    % subplot(121);montage(FCim_gen)
    % subplot(122);montage(BGim_gen)
    % figure(4);
    % subplot(121);imshow(meanFCimg)
    % subplot(122);imshow(meanBGimg)
    % pause
end
%%
% from the first to last - 1 generation, the mean activation, 2% and 98% to avoid the max min unstable
CLIM = prctile(meanact_vec(:,1:end-1),[2,98],'all')'; 
frame_FCimg_col = score_frame_image_arr(WtmeanFCimg_col, meanact_vec(1,1:end-1), CLIM);%, cmap, LineWidth)
frame_BGimg_col = score_frame_image_arr(WtmeanBGimg_col, meanact_vec(2,1:end-1), CLIM);%, cmap, LineWidth)
frame_FCbestimg_col = score_frame_image_arr(bestFCimg_col, meanact_vec(1,1:end-1), CLIM);%, cmap, LineWidth)
frame_BGbestimg_col = score_frame_image_arr(bestBGimg_col, meanact_vec(2,1:end-1), CLIM);%, cmap, LineWidth)
%%
FCbestmtg = imtile(bestFCimg_col,8);
BGbestmtg = imtile(bestBGimg_col,8);
frame_FCbestmtg = imtile(frame_FCbestimg_col,8);
frame_BGbestmtg = imtile(frame_BGbestimg_col,8);
FCmtg = imtile(WtmeanFCimg_col,8);
BGmtg = imtile(WtmeanBGimg_col,8);
frame_FCmtg = imtile(frame_FCimg_col,8);
frame_BGmtg = imtile(frame_BGimg_col,8);
figure(5);set(5,'pos',[ 81         100       1720         750])
T = tiledlayout(1,2,'tilesp','compact','pad','compact');
nexttile(1);imshow(FCmtg) % montage(meanFCimg_col)
nexttile(2);imshow(BGmtg) % montage(meanBGimg_col)
title(T,compose("Image Evolution Color lim [%.2f,%.2f]\n%s",CLIM(1),CLIM(2),expstr))
figure(6);set(6,'pos',[ 81         100       1720         750])
T = tiledlayout(1,2,'tilesp','compact','pad','compact');
nexttile(T,1);imshow(frame_FCmtg)
nexttile(T,2);imshow(frame_BGmtg)
title(T,compose("Image Evolution Color lim [%.2f,%.2f]\n%s",CLIM(1),CLIM(2),expstr))
figure(7);set(7,'pos',[ 81         100       1720         750])
T = tiledlayout(1,2,'tilesp','compact','pad','compact');
nexttile(T,1);imshow(frame_FCbestmtg)
nexttile(T,2);imshow(frame_BGbestmtg)
title(T,compose("Image Evolution Color lim [%.2f,%.2f]\n%s",CLIM(1),CLIM(2),expstr))
%%
saveallform([sumdir,expdir], compose("%s_Exp%02d_FC6BGImageEvol_cmp",Animal,Expi), 5, ["png","pdf"])
saveallform([sumdir,expdir], compose("%s_Exp%02d_FC6BGImageEvol_framed_cmp",Animal,Expi), 6, ["png","pdf"])
saveallform([sumdir,expdir], compose("%s_Exp%02d_FC6BGImageEvol_bestframed_cmp",Animal,Expi), 7, ["png","pdf"])

imwrite(FCmtg, fullfile(expdir,compose("%s_Exp%02d_FC6ImageEvol.png",Animal,Expi)))
imwrite(BGmtg, fullfile(expdir,compose("%s_Exp%02d_BGImageEvol.png",Animal,Expi)))
imwrite(frame_FCmtg, fullfile(expdir,compose("%s_Exp%02d_FC6ImageEvol_framed.png",Animal,Expi)))
imwrite(frame_BGmtg, fullfile(expdir,compose("%s_Exp%02d_BGImageEvol_framed.png",Animal,Expi)))
imwrite(FCbestmtg, fullfile(expdir,compose("%s_Exp%02d_FC6ImageEvol_best.png",Animal,Expi)))
imwrite(BGbestmtg, fullfile(expdir,compose("%s_Exp%02d_BGImageEvol_best.png",Animal,Expi)))
imwrite(frame_FCbestmtg, fullfile(expdir,compose("%s_Exp%02d_FC6ImageEvol_best_framed.png",Animal,Expi)))
imwrite(frame_BGbestmtg, fullfile(expdir,compose("%s_Exp%02d_BGImageEvol_best_framed.png",Animal,Expi)))
end
%% 
catch ME
     failures(end + 1).iter = Expi;
     failures(end).str  = getReport(ME);
     disp(ME)
     save("D:\BigGAN_run_failure.mat",'failures')
     fileID = fopen('D:\BigGAN_run_failure.txt','a+');
     nbytes = fprintf(fileID,"Exp %d\n",Expi);
     nbytes = fprintf(fileID,getReport(ME));
end
end


%%

%% test data loading time
tic 
[FCim_gen, FCimnames] = loadEvolImages(stimpath, 0, 1);
[BGim_gen, BGimnames] = loadEvolImages(stimpath, 1, 1);
toc


%%
function [genvec, actvec] = actcell2traj(actcol)
actvec = cat(1,actcol{:});
gencol = arrayfun(@(geni) geni * ones(numel(actcol{geni}),1), 1:numel(actcol),'uni',0);
genvec = cat(1,gencol{:});
end

%%
function ImageEvol_pair_cmp(FCimg_col, BGimg_col, actmat, prefix, suffix, expdir, rootdir)
CLIM = prctile(actmat(:,:),[2,98],'all')'; 
frame_FCimg_col = score_frame_image_arr(FCimg_col, actmat(1,:), CLIM);%, cmap, LineWidth)
frame_BGimg_col = score_frame_image_arr(BGimg_col, actmat(2,:), CLIM);%, cmap, LineWidth)
FCmtg = imtile(FCimg_col,8);
BGmtg = imtile(BGimg_col,8);
frame_FCmtg = imtile(frame_FCimg_col,8);
frame_BGmtg = imtile(frame_BGimg_col,8);
imwrite(FCmtg, fullfile(expdir,compose("%s_FC6ImageEvol%s.png",prefix,suffix)))
imwrite(BGmtg, fullfile(expdir,compose("%s_BGImageEvol%s.png",prefix,suffix)))

end