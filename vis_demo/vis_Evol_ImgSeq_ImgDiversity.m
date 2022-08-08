% for the NatCompSci paper

Set_Path;
mat_dir = "O:\Mat_Statistics"; 
%%
Animal = "Alfa"; 
load(fullfile(mat_dir, Animal+"_Evol_stats.mat"))
load(fullfile(mat_dir, Animal+"_Manif_stats.mat"))
G = FC6Generator();

%% Evolution stimuli
Expi = 33;
stimdir = EStats(Expi).meta.stimuli;

%% Calc mean scores.
actcol = cellfun(@(P)squeeze(mean(P(1,51:200,:),[1,2])),EStats(Expi).evol.psth,'Unif',0);
bslcol = cellfun(@(P)squeeze(mean(P(1,1:51,:),[1,2])),EStats(Expi).evol.psth,'Unif',0);
actmean = cellfun(@mean, actcol);
actsem = cellfun(@sem, actcol);
bslmean = mean(cat(1,bslcol{:}));
bslsem = sem(cat(1,bslcol{:}));
%% Load up all codes 
[codes_all, img_ids, generations] = load_codes_all(stimdir,1);
%% Get mean codes 
code_mean = arrayfun(@(geni)mean(codes_all(generations==geni,:),1),[min(generations):max(generations)]','uni',0);
code_mean = cat(1,code_mean{:});
genN = size(code_mean,1);
%% Visualize the mean images 
code_mean_imgs = G.visualize(code_mean(:,:));
%% Export the Evolution Image Trajectory
outdir = "E:\OneDrive - Harvard University\FactorModel_NatCompSci\Figures\EvolImageSeq";
img_frame = score_frame_image_arr(code_mean_imgs(:,:,:,[1:1:genN-1]),...
    actmean([1:1:genN-1]),[15.53,90]);
mtg = imtile(img_frame, 'GridSize',[4,6],'Thumb',[296,296]);
imwrite(mtg,fullfile(outdir,compose("%s_Exp%d_EvolImgTraj2.png",Animal,Expi)))
%% 
figure;
montage(img_frame,'Size',[3,7],'Thumb',[296,296])
%% Visualize the mean images 
code_last_gen = codes_all(generations==max(generations)-1,:);
act_last_gen = actcol{end-1};
last_gen_imgs = G.visualize(code_last_gen(:,:));
%% 
img_fram_lastgen = score_frame_image_arr(last_gen_imgs,act_last_gen,[15.53,90]); % ,78.41
mtg2 = imtile(img_fram_lastgen, 'GridSize',[4,10],'Thumb',[296,296]);
imwrite(mtg2,fullfile(outdir,compose("%s_Exp%d_final_gen.png",Animal,Expi)))
% figure;
% montage(img_fram_lastgen,'Size',[5,8],'Thumb',[296,296])
%%
actvec = cat(1,actcol{:});
generation_col = arrayfun(@(i) i*ones(numel(actcol{i}),1),1:numel(actcol),'Uni',0);
generation_vec = cat(1,generation_col{:});
%%
figure('pos',[200,200,400,400]);
scatter(generation_vec, actvec,'filled','markerfacealpha',0.2)
shadedErrorBar([],actmean,actsem,'lineprop',{'color','red'})
shadedErrorBar([1,numel(actmean)],...
    [bslmean,bslmean],[bslsem,bslsem])
xlabel("Generation")
ylabel("firing rate (spike/sec)")
title(compose("%s Exp%02d Evolution Chan %02d",Animal,Expi,Stats(Expi).units.pref_chan))
legend(["single trial firing rate","generation mean firing rate","baseline"],'location','best')
xlim([0,26])
saveallform(outdir,compose("%s_Exp%d_EvolTraj",Animal,Expi))

%%
figure('pos',[ 680   536   488   442]);
heatmap(linspace(15.53,90.00,50))
colorbar;
colormap("parula")
saveallform(outdir,compose("%s_Exp%d_colorbartmp",Animal,Expi))