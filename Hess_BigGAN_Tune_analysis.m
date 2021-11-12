%% Adapted from Hess_Tune_Analysis
%  Generate figures for Selectivity experiments for tuning along Hessian
%  eigen-axes. Compatible for different versions of the naming convention. 
%  Tested. 
Animal = "Both";Set_Path;
expftr = contains(ExpRecord.expControlFN,["210312_Alfa_selectivity_basic"]);%& contains(ExpRecord.expControlFN, "selectivity");; %& ,"200812"contains(ExpRecord.Exp_collection,"BigGAN_Hessian");% & contains(ExpRecord.Exp_collection,"BigGAN");
% expftr = contains(ExpRecord.Exp_collection,"BigGAN_Hessian") & contains(ExpRecord.expControlFN, "selectivity");
fllist = find(expftr);no_return=false;
[meta_new,rasters_new,~,Trials_new] = loadExperiments(fllist(1:end),Animal,no_return);
%% Load the experiments
Animal = "Both";Set_Path;
% expftr = contains(ExpRecord.expControlFN,["200813"]); %& ,"200812"contains(ExpRecord.Exp_collection,"BigGAN_Hessian");% & contains(ExpRecord.Exp_collection,"BigGAN");
expftr = contains(ExpRecord.Exp_collection,"BigGAN_Hessian") & contains(ExpRecord.expControlFN, "selectivity");
fllist = find(expftr);no_return=false;
[meta_new,rasters_new,~,Trials_new] = loadExperiments(fllist(end-5:end),Animal,no_return);
% figdir = "E:\OneDrive - Washington University in St. Louis\HessBigGANTune\Beto_Exp04";
% mkdir(figdir)
%% Prepare Image Metric
D = torchImDist();
%%
saveroot = "E:\OneDrive - Washington University in St. Louis\HessBigGANTune";
P=struct();
P.plotTuneCurve=false;% TO IMplement.
for Triali = numel(meta_new)-3:numel(meta_new)-2
%%
meta = meta_new{Triali};
rasters = rasters_new{Triali};
% lfps = lfps_new{Triali};
Trials = Trials_new{Triali};
if contains(meta.ephysFN,"Alfa"),Animal = "Alfa";elseif contains(meta.ephysFN,"Beto"),Animal = "Beto";end
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
Expi = ExpRecord.Expi(exp_rowi);
fprintf("Processing  Exp %d:\n",Expi)
fprintf([ExpRecord.comments{exp_rowi},'\n'])
if contains(meta.ephysFN,'Alfa-13012021-003'), fprintf("This is a incomplete exp!"); end
% Check the Expi match number
% unit_name_arr = generate_unit_labels(meta.spikeID);
% [activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% newer version generate unit labels
unit_num_arr = meta.unitID;
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
imgname_uniq = unique(Trials.imageName); 
prefchan = Trials.TrialRecord.User.prefChan;
prefchan_ids = find(meta.spikeID == prefchan & meta.unitID>0);

% Old naming convention
% figdir = fullfile("E:\OneDrive - Washington University in St. Louis\HessBigGANTune",compose("%s_Exp%02d",Animal,Expi));
% if meta.ephysFN == "Alfa-10082020-004", figdir = figdir + "_post"; end
% mkdir(figdir);fprintf("writing to %s\n",figdir)

stimparts = split(meta.stimuli,"\");
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
fdrnm = compose("%s-Ch%02d", stimparts{end}, prefchan);
% fdrnm = compose("%s-%s-Chan%02d-1",datestr(expday,'yyyy-mm-dd'), Animal, pref_chan(1));
figdir = fullfile(saveroot, fdrnm);
if exist(figdir,'dir'),warning("%s figure directory exist! Beware",figdir);keyboard;end
mkdir(figdir)

%% Load the images in the class space and noise space
% Identify the naming convention used in this Exp. (Different version of py code...)
[noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq);
% Extract parameters of images in noise space from img name
namepart_uniq = regexp(imgname_uniq, noise_pattern, 'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
eig_id_arr_nos = unique(cellfun(@(U)str2num(U.eig_id),namepart_uniq));
dist_arr_nos = unique(cellfun(@(U)str2num(U.dist),namepart_uniq));
if ~any(dist_arr_nos==0), dist_arr_nos = unique([dist_arr_nos;0]); end
% collect image name and trial indices in cell
imgnm_arr_nos = strings(length(eig_id_arr_nos), length(dist_arr_nos));
idx_arr_nos = cell(length(eig_id_arr_nos), length(dist_arr_nos));
for i = 1:length(eig_id_arr_nos)
    eig_id = eig_id_arr_nos(i);
    for j = 1:length(dist_arr_nos)
    dist = dist_arr_nos(j);
    imgnm = compose(noise_imgnm,eig_id,dist); % e.g. "noise_eig%d_lin%.1f"
    idx_arr_nos{i,j} = find(contains(Trials.imageName, imgnm));
    imgnm_arr_nos(i,j) = imgnm;
    if numel(idx_arr_nos{i,j}) == 0 && dist== 0 
    % Recently, we deleted images like "noise_eig2_lin0.0" so should use "noise_eig1_lin0.0" instead
    imgnm_arr_nos(i,j) = compose(noise_imgnm,eig_id_arr_nos(1),0);
    end
    end
end
img_noise = arrayfun(@(imgnm)imread(fullfile(meta.stimuli, imgnm+".jpg")),imgnm_arr_nos,"Un",false);
% Extract parameters of images in class space from img name
namepart_uniq = regexp(imgname_uniq,class_pattern,'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
eig_id_arr_cls = unique(cellfun(@(U)str2num(U.eig_id),namepart_uniq));
dist_arr_cls = unique(cellfun(@(U)str2num(U.dist),namepart_uniq));
% collect image name and trial indices in cell
imgnm_arr_cls = strings(length(eig_id_arr_cls), length(dist_arr_cls));
idx_arr_cls = cell(length(eig_id_arr_cls), length(dist_arr_cls));
for i = 1:length(eig_id_arr_cls)
    eig_id = eig_id_arr_cls(i);
    for j = 1:length(dist_arr_cls)
    dist = dist_arr_cls(j);
    imgnm = compose(class_imgnm,eig_id,dist); % "class_eig%d_lin%.1f"
    idx_arr_cls{i,j} = find(contains(Trials.imageName, imgnm));
    imgnm_arr_cls(i,j) = imgnm;
    if numel(idx_arr_cls{i,j}) == 0 && dist== 0 
    % Recently, we deleted images like "noise_eig2_lin0.0" so should use "noise_eig1_lin0.0" instead
    imgnm_arr_cls(i,j) = compose(class_imgnm,eig_id_arr_cls(1),0);
    end
    end
end
img_class = arrayfun(@(imgnm)imread(fullfile(meta.stimuli, imgnm+".jpg")),imgnm_arr_cls,"Un",false);
%  Plot the image tile in the experiment
cls_tile = imtile(img_class','GridSize',size(img_class));
imwrite(cls_tile, fullfile(figdir, compose("Class_Images_Tile.jpg")))
nos_tile = imtile(img_noise','GridSize',size(img_noise));
imwrite(nos_tile, fullfile(figdir, compose("Noise_Images_Tile.jpg")))

nrow_cls = length(eig_id_arr_cls);
ncol_cls = length(dist_arr_cls);
nrow_nos = length(eig_id_arr_nos);
ncol_nos = length(dist_arr_nos);
%%  Compute Image Dissimilarity and make heatmap
imdist_class = zeros(size(img_class));
cent_i = (ncol_cls + 1) / 2;
for rowi = 1:size(img_class,1)
imrow = cat(4,img_class{rowi,:});
distrow = D.distance(imrow, imrow(:,:,:,cent_i));
imdist_class(rowi, :) = distrow;
end
imdist_noise = zeros(size(img_noise));
cent_i = (ncol_nos + 1) / 2;
for rowi = 1:size(img_noise,1)
imrow = cat(4,img_noise{rowi,:});
distrow = D.distance(imrow, imrow(:,:,:,cent_i));
imdist_noise(rowi, :) = distrow;
end
% Visualize image dissimilarity with heatmap
figure(1);set(1,'pos',[1000  215  560  763])
imagesc(imdist_class)
colorbar();axis image
title(["Deviation from center image by SqueezeNet ImDist","Class space"])
ylabel("eigen id");xlabel("Linear Distance")
xticks(1:length(dist_arr_cls));xticklabels(dist_arr_cls)
yticks(1:length(eig_id_arr_cls));yticklabels(eig_id_arr_cls)
saveas(1,fullfile(figdir, "class_space_ImDist.jpg"))
figure(2);set(2,'pos',[1000  215  560  763])
imagesc(imdist_noise)
colorbar();axis image
title(["Deviation from center image by SqueezeNet ImDist","Noise space"])
ylabel("eigen id");xlabel("Linear Distance")
xticks(1:length(dist_arr_nos));xticklabels(dist_arr_nos)
yticks(1:length(eig_id_arr_nos));yticklabels(eig_id_arr_nos)
saveas(2,fullfile(figdir, "noise_space_ImDist.jpg"))
%% PSTH collections for different classes of images

%% Plot the neural response to image tiles AND compute statistics
figure(12);set(12,'pos',[325    67   717   850]);%set(12,'position',[315   506   560   444])
figure(13);set(13,'position',[961    64   676   850])%set(13,'position',[315   159   560   608])
ExpLabel = compose("%s Pilot Exp %d pref chan %d", Animal, Expi, prefchan);
for iCh = prefchan_ids'%1:size(rasters,1) %prefchan_ids'%%
resp_mat_nois = cellfun(@(idx)mean(rasters(iCh,51:200,idx),'all'), idx_arr_nos);
resp_mat_cls = cellfun(@(idx)mean(rasters(iCh,51:200,idx),'all'), idx_arr_cls);
% Form the psth cell array and compute F, T statistics for rows and cols
psth_col_nois = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_nos, 'uni',0);
psth_col_clas = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_cls, 'uni',0);
stats_nois = calc_tune_stats(psth_col_nois);
stats_clas = calc_tune_stats(psth_col_clas);
stats_row_nois = arrayfun(@(i)calc_tune_stats(psth_col_nois(i,:)), 1:nrow_nos);
stats_col_nois = arrayfun(@(i)calc_tune_stats(psth_col_nois(:,i)), 1:ncol_nos);
stats_row_clas = arrayfun(@(i)calc_tune_stats(psth_col_clas(i,:)), 1:nrow_cls);
stats_col_clas = arrayfun(@(i)calc_tune_stats(psth_col_clas(:,i)), 1:ncol_cls);
% Plot Score Heatmap in Noise space
set(0,'CurrentFigure',12);clf;
imagesc(resp_mat_nois)
colorbar();axis image
ylabel("eigen id")
xlabel("Linear Distance")
xticks(1:ncol_nos);xticklabels(dist_arr_nos)
yticks(1:nrow_nos);yticklabels(eig_id_arr_nos)
title(compose("%s Ch %s\nBigGAN Noise Space Tuning Along Hessian EigenVectors\nF=%.1f(p=%.1e) T=%.1f(p=%.1e)",...
    ExpLabel, unit_name_arr(iCh),stats_nois.F,stats_nois.F_P,stats_nois.T,stats_nois.t_P))
rowFstr = arrayfun(@(S)compose("F=%.1f(%.e)",S.F,S.F_P),stats_row_nois)';
colFstr = arrayfun(@(S)compose("F=%.1f\n(%.e)",S.F,S.F_P),stats_col_nois)';
text(ncol_nos + 1+ones(nrow_nos,1), 1:nrow_nos, rowFstr')
text(1:ncol_nos, nrow_nos+2*ones(ncol_nos,1), colFstr,'HorizontalAlignment','center')
saveas(12, fullfile(figdir, compose("Noise_TuneMap_%s.jpg",unit_name_arr(iCh))))
% Plot Score Heatmap in Class space
set(0,'CurrentFigure',13);clf;
imagesc(resp_mat_cls)
colorbar();axis image
ylabel("eigen id")
xlabel("Linear Distance")
xticks(1:ncol_cls);xticklabels(dist_arr_cls)
yticks(1:nrow_cls);yticklabels(eig_id_arr_cls)
title(compose("%s Ch %s\nBigGAN Class Space Tuning Along Hessian EigenVectors\nF=%.1f(p=%.1e) T=%.1f(p=%.1e)",...
    ExpLabel, unit_name_arr(iCh),stats_clas.F,stats_clas.F_P,stats_clas.T,stats_clas.t_P))
rowFstr = arrayfun(@(S)compose("F=%.1f(%.e)",S.F,S.F_P),stats_row_clas)';
colFstr = arrayfun(@(S)compose("F=%.1f\n(%.e)",S.F,S.F_P),stats_col_clas)';
text(ncol_cls + 1+ones(nrow_cls,1), 1:nrow_cls, rowFstr')
text(1:ncol_cls, nrow_cls+2*ones(ncol_cls,1), colFstr,'HorizontalAlignment','center')
saveas(13, fullfile(figdir, compose("Class_TuneMap_%s.jpg",unit_name_arr(iCh))))
% Plot Score framed image tile in class and noise space
frame_img_class = score_frame_image_arr(img_class, resp_mat_cls, prctile(resp_mat_cls,[0,100],'all')'+[0,0.001], parula, 10);
frame_img_noise = score_frame_image_arr(img_noise, resp_mat_nois, prctile(resp_mat_nois,[0,100],'all')'+[0,0.001], parula, 10);
cls_scr_tile = imtile(frame_img_class','GridSize',[length(eig_id_arr_cls),length(dist_arr_cls)]);
imwrite(cls_scr_tile, fullfile(figdir, compose("Class_TuneTile_%s.jpg",unit_name_arr(iCh))))
nos_scr_tile = imtile(frame_img_noise','GridSize',[length(eig_id_arr_nos),length(dist_arr_nos)]);
imwrite(nos_scr_tile, fullfile(figdir, compose("Noise_TuneTile_%s.jpg",unit_name_arr(iCh))))
end % End Loop of Channels
end % End Loop of Exp


%% Purely Calculate the T and F stats for channels
for iCh = prefchan_ids %1:size(rasters,1)
resp_mat_nois = cellfun(@(idx)mean(rasters(iCh,51:200,idx),'all'), idx_arr_nos);
resp_mat_cls = cellfun(@(idx)mean(rasters(iCh,51:200,idx),'all'), idx_arr_cls);
psth_col_nois = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_nos, 'Un',0);
psth_col_clas = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_cls, 'Un',0);

stats_nois = calc_tune_stats(psth_col_nois);
stats_row_nois = arrayfun(@(i)calc_tune_stats(psth_col_nois(i,:)), 1:nrow_nos);
stats_col_nois = arrayfun(@(i)calc_tune_stats(psth_col_nois(:,i)), 1:ncol_nos);
arrayfun(@(S)S.F_P,stats_row_nois)'

stats_clas = calc_tune_stats(psth_col_clas);
stats_row_clas = arrayfun(@(i)calc_tune_stats(psth_col_clas(i,:)), 1:nrow_cls);
stats_col_clas = arrayfun(@(i)calc_tune_stats(psth_col_clas(:,i)), 1:ncol_cls);
arrayfun(@(S)S.F_P,stats_row_clas)'
end
%%
%% Add statistics to side annotation of the response heatmap
figure(5);clf;hold on
% subplot('position',[0.05,0.2,0.72,0.7])
imagesc(resp_mat_cls)
colorbar();axis image;set(gca,"YDir",'reverse')
ylabel("eigen id")
xlabel("Linear Distance")
xticks(1:ncol_cls);xticklabels(dist_arr_cls)
yticks(1:nrow_cls);yticklabels(eig_id_arr_cls)
title(compose("%s %s Pilot Exp\nBigGAN Class Space Tuning Along Hessian EigenVectors\nF=%.1f(p=%.1e) T=%.1f(p=%.1e)",...
    Animal, unit_name_arr(iCh),stats_clas.F,stats_clas.F_P,stats_clas.T,stats_clas.t_P))
rowFstr = arrayfun(@(S)compose("F=%.1f(%.e)",S.F,S.F_P),stats_row_clas)';
colFstr = arrayfun(@(S)compose("F=%.1f\n(%.e)",S.F,S.F_P),stats_col_clas)';
text(ncol_cls + 1+ones(nrow_cls,1), 1:nrow_cls, rowFstr')
text(1:ncol_cls, nrow_cls+2*ones(ncol_cls,1), colFstr,'HorizontalAlignment','center')
% subplot('position',[0.80,0.2,0.15,0.7])
% imagesc(arrayfun(@(S)S.F,stats_row_clas)')
% axis image;axis off
% colorbar()
% subplot('position',[0.05,0.05,0.72,0.12])
% imagesc(arrayfun(@(S)S.F,stats_col_clas))
% colorbar()
% axis image;axis off
%% Plot individual trials see the deviation 
cnt_nois = cellfun(@length, idx_arr_nos);
cnt_clas = cellfun(@length, idx_arr_cls);
cnt_col_nois = arrayfun(@(n,d) d*ones(1,n), cnt_nois, repmat(dist_arr_nos',length(eig_id_arr_nos),1), 'Un', 0);
cnt_col_clas = arrayfun(@(n,d) d*ones(1,n), cnt_clas, repmat(dist_arr_cls',length(eig_id_arr_cls),1), 'Un', 0);
scores_mat = mean(rasters(:,51:200,:),2);
for iCh = prefchan_ids%1:size(rasters,1) % prefchan_ids
resp_mat_nois = cellfun(@(idx)mean(scores_mat(iCh,idx),'all'), idx_arr_nos);
resp_mat_cls = cellfun(@(idx)mean(scores_mat(iCh,idx),'all'), idx_arr_cls);
resp_sem_nois = cellfun(@(idx)std(scores_mat(iCh,idx))/sqrt(length(idx)), idx_arr_nos);
resp_sem_cls = cellfun(@(idx)std(scores_mat(iCh,idx))/sqrt(length(idx)), idx_arr_cls);
resp_col_nois = cellfun(@(idx)scores_mat(iCh,idx), idx_arr_nos, 'uni', 0);
resp_col_cls = cellfun(@(idx)scores_mat(iCh,idx), idx_arr_cls, 'uni', 0);
figure(3);clf;hold on
plot(repmat(dist_arr_nos',nrow_nos,1)',resp_mat_nois')
xlabel("Linear Distance")
% for ri = 1:length(eig_id_arr_nos)
% scatter(cat(2,cnt_col_nois{ri,:}), cat(2,resp_col_nois{ri,:}))
% end
xticks(dist_arr_cls);xticklabels(dist_arr_cls)
end
%%

figure(7);
% for ri = 1:nrow_nos
% for ci = 1:ncol_nos
% ax1 = subtightplot(size(idx_arr_nos,1),size(idx_arr_nos,2),ci+ncol_nos*(ri-1));
% if ci==1, ylabel(eig_id_arr_nos(ri)); end
% if ri==nrow_nos, xlabel(dist_arr_nos(ci)); end
% end
% end
TL = tiledlayout(nrow_nos,ncol_nos,'TileSpacing','compact','Padding','compact');
xlabel(TL,"Im Distance");ylabel(TL,"Eigen Vector id")
%
for iCh = 2:size(rasters,1)%iCh = prefchan_ids;prefchan_ids%
tic
psth_col_nois = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_nos, 'uni',0);
psth_col_clas = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_cls, 'uni',0);
stats_nois = calc_tune_stats(psth_col_nois);
stats_clas = calc_tune_stats(psth_col_clas);
title(TL,compose("%s Ch %s\nBigGAN Noise Space Tuning Along Hessian EigenVectors\n Raw raster plot F=%.1f(p=%.1e) T=%.1f(p=%.1e)",...\n
    ExpLabel, unit_name_arr(iCh), stats_nois.F, stats_nois.F_P, stats_nois.T, stats_nois.t_P))
for ri = 1:nrow_nos
for ci = 1:ncol_nos
ax1 = nexttile(ci+ncol_nos*(ri-1));
% ax1 = subtightplot(size(idx_arr_nos,1),size(idx_arr_nos,2),ci+ncol_nos*(ri-1));
% Obsolete plotting routines
% rasterplot(find(rasters(iCh,:,idx_arr_nos{ri,ci})),4,200,ax1)
% imagesc(squeeze(rasters(iCh,:,idx_arr_nos{ri,ci}))');
% colormap(flipud(gray));box off
% plt_col{ri,ci} = plot([0],[0],'k');
[spkT,spktr] = find(squeeze(rasters(iCh,:,idx_arr_nos{ri,ci})));
plot([spkT,spkT]',[spktr-0.7,spktr]',"color",'k')
xlim([0,200]);ylim([0,length(idx_arr_nos{ri,ci})]);yticks([]);xticks([]);box off
if ci==1, ylabel(eig_id_arr_nos(ri)); end
if ri==nrow_nos, xlabel(dist_arr_nos(ci)); end
end
end
% keyboard
saveas(7, fullfile(figdir, compose("Noise_RasterMap_%s.jpg",unit_name_arr(iCh))))
toc
end
%%
figure(8); set(8,'pos',[666    42   904   954])
TL = tiledlayout(nrow_nos,ncol_nos,'TileSpacing','compact','Padding','compact');
xlabel(TL,"Im Distance");ylabel(TL,"Eigen Vector id")
for iCh = 1:size(rasters,1)%iCh = prefchan_ids;prefchan_ids%
tic
psth_col_nois = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_nos, 'uni',0);
psth_col_clas = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_cls, 'uni',0);
stats_nois = calc_tune_stats(psth_col_nois);
stats_clas = calc_tune_stats(psth_col_clas);
title(TL,compose("%s Ch %s\nBigGAN Class Space Tuning Along Hessian EigenVectors\n Raw raster plot F=%.1f(p=%.1e) T=%.1f(p=%.1e)",...\n
    ExpLabel, unit_name_arr(iCh), stats_clas.F, stats_clas.F_P, stats_clas.T, stats_clas.t_P))
for ri = 1:nrow_nos
for ci = 1:ncol_nos
ax1 = nexttile(ci+ncol_nos*(ri-1));
[spkT,spktr] = find(squeeze(rasters(iCh,:,idx_arr_cls{ri,ci})));
plot([spkT,spkT]',[spktr-0.7,spktr]',"color",'k')
xlim([0,200]);ylim([0,length(idx_arr_cls{ri,ci})]);yticks([]);xticks([]);box off
if ci==1, ylabel(eig_id_arr_cls(ri)); end
if ri==nrow_nos, xlabel(dist_arr_cls(ci)); end
end
end
saveas(8, fullfile(figdir, compose("Class_RasterMap_%s.jpg",unit_name_arr(iCh))))
toc
end
%%
tic
for ri = 1:nrow_nos
for ci = 1:ncol_nos
[spkT,spktr] = find(squeeze(rasters(iCh,:,idx_arr_nos{ri,ci})));
plt_col{ri,ci}.XData = [spkT,spkT]';
plt_col{ri,ci}.YData = [spktr,spktr+.7]';
% plot([spkT,spkT]',[spktr,spktr+.7]',"color",'k')
end
end
toc
%%
figure;
[spkT,spktr] = find(squeeze(rasters(prefchan_ids,:,idx_arr_nos{ri,ci})));
plot([spkT,spkT]',[spktr,spktr+.7]',"color",'k')
% stem(spkT,spktr)
% imagesc(squeeze(rasters(prefchan_ids,:,idx_arr_nos{ri,ci}))');colormap(flipud(gray))
xlim([0,200]);yticks([]);xticks([])
%% Test the functions
[noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq);
namepart_uniq = regexp(imgname_uniq,noise_pattern,'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
%% Use Deconv Lucy to get spike from trials
ker = normpdf(-6:6,0,2);
tic
rasters_tmp = reshape(permute(rasters,[3,1,2]),[],200);
rasters_tr3 = deconvlucy(rasters_tmp(:,:), ker, 100);
toc
%%
rasters_tr3 = reshape(rasters_tr3,[],size(rasters,1),200);
%%
[noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention( unique(Trials.imageName))
function [noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq)
% noise part
namepart_uniq = regexp(imgname_uniq,"noise_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)",'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
if ~isempty(namepart_uniq)
    decim = cellfun(@(nm)split(nm.dist,"."), namepart_uniq, 'uni', false);
    float_L = unique(cellfun(@(dcm)length(dcm{2}),decim));
    if float_L == 1
        fprintf("Old Image Naming Convention in Noise Space\n")
        noise_pattern = "noise_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
        noise_imgnm = "noise_eig%d_lin%.1f";
    elseif float_L == 2
        fprintf("Latest Line Search Image Naming Convention in Noise Space\n")
        noise_pattern = "noise_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
        noise_imgnm = "noise_eig%d_lin%.2f";
    end
else
    fprintf("New Image Naming Convention in Noise Space\n")
    namepart_uniq = regexp(imgname_uniq,"noise_eig(?<eig_id>\d*)_exp(?<expon>[.\d]*)_lin(?<dist>[-.\d]*)",'names');
    namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
    if ~isempty(namepart_uniq)
        expon_str = string(unique(cellfun(@(U)U.expon,namepart_uniq,"Uni",false)));
        noise_pattern = "noise_eig(?<eig_id>\d*)_exp(?<expon>[.\d]*)_lin(?<dist>[-.\d]*)";
        noise_imgnm = "noise_eig%d_exp"+expon_str+"_lin%.1f";
    else
        error("Name pattern not recognized.")
    end
end
% class part
namepart_uniq = regexp(imgname_uniq,"class_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)",'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
if ~isempty(namepart_uniq)
    decim = cellfun(@(nm)split(nm.dist,"."), namepart_uniq, 'uni', false);
    float_L = unique(cellfun(@(dcm)length(dcm{2}),decim));
    if float_L == 1
        fprintf("Old Image Naming Convention in Class Space\n")
        class_pattern = "class_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
        class_imgnm = "class_eig%d_lin%.1f";
    elseif float_L == 2
        fprintf("Latest Line Search Image Naming Convention in Class Space\n")
        class_pattern = "class_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
        class_imgnm = "class_eig%d_lin%.2f";
    end
else
    fprintf("New Image Naming Convention in Class Space\n")
    namepart_uniq = regexp(imgname_uniq,"class_eig(?<eig_id>\d*)_exp(?<expon>[.\d]*)_lin(?<dist>[-.\d]*)",'names');
    namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
    if ~isempty(namepart_uniq)
        expon_str = string(unique(cellfun(@(U)U.expon,namepart_uniq,"Uni",false)));
        class_pattern = "class_eig(?<eig_id>\d*)_exp(?<expon>[.\d]*)_lin(?<dist>[-.\d]*)";
        class_imgnm = "class_eig%d_exp"+expon_str+"_lin%.1f";
    else
        error("Name pattern not recognized.")
    end
end
end