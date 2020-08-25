%% Visualize Evol BigGAN Geometry
%  
Animal="Beto";Set_Path;
rowidx = find(contains(ExpRecord.expControlFN,"BigGAN"));
ExpRecord(find(contains(ExpRecord.expControlFN,"BigGAN")),:)
ExpRecord(end-10:end,:)
%[meta_new,rasters_new,~,Trials_new] = loadExperiments(rowidx, Animal);
%% Load the Trials and Spikes
rowidx = find(contains(ExpRecord.expControlFN,"BigGAN") & contains(ExpRecord.ephysFN,"Beto-28072020-003"));
[meta_new,rasters_new,~,Trials_new] = loadExperiments(rowidx, Animal);
%% Set up GAN and optimizer
G = torchBigGAN();
dummyOptim = CMAES_simple(128,[],struct());
%%
Triali = 1;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
% lfps = lfps_new{Triali};
Trials = Trials_new{Triali};
%% Load some setting variables
bhvfn = meta.expControlFN;
space_opts = Trials.TrialRecord.User.space_opts;
optim_names = Trials.TrialRecord.User.evoConfiguration(:,5);
% save(fullfile(TrialRecord.User.newPicsHome, "space_opts.mat"), "space_opts");
backupdir = meta.stimuli;%Trials.TrialRecord.User.newPicsHome;
threadn = size(Trials.TrialRecord.User.evoConfiguration,1);
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
G=G.select_space(Trials.TrialRecord.User.space_cfg{2}{1},Trials.TrialRecord.User.space_cfg{2}{1});
%% Visualize Evolution Generation by Generation
iThread = 2;
codefnlist = string(ls(fullfile(backupdir,compose("block*_thread%03d_code.mat",iThread-1))));
% assume codefnlist start from block1
fig = figure(14);set(fig,'position',[410       100        1013         969])
ax = subplot(1,1,1);
iGen = 1;
for iGen = 1:length(codefnlist)-1
data = load(fullfile(backupdir,codefnlist(iGen)));
imgids = string(data.ids);
codes_cur = data.codes;
scores_cur = Trials.TrialRecord.User.scores_record{iGen+1, iThread};
[score_sorted, sortidx] = sort(scores_cur, 'descend');
imgs_cur = arrayfun(@(nm)imread(fullfile(backupdir, nm)),imgids(sortidx)','Uni',0);
montage(imgs_cur,'Parent',ax)
if strcmp(space_opts(2).name,'BigGAN') % if in whole space, show norm in 2 subspaces separately and score
noise_norms = norm_axis(codes_cur(:,1:128),2);
class_norms = norm_axis(codes_cur(:,129:end),2);
title(compose("%s Thread %d Space: %s Optim: %s\nGen %d max %.1f min %.1f mean%.1f\n norm: noise %.1f(%.1f) class %.2f(%.2f)", ...
    strrep(bhvfn,'_',' '), iThread, space_opts(iThread).name, optim_names{iThread},...    
    iGen+1, max(scores_cur), min(scores_cur), mean(scores_cur), ...
    mean(noise_norms),std(noise_norms), mean(class_norms),std(class_norms)))
else % else just show code norm and score
code_norms = norm_axis(codes_cur,2);
title(compose("%s Thread %d Space: %s Optim: %s\nGen %d max %.1f min %.1f mean%.1f\n norm: %.2f(%.2f)", ...
    strrep(bhvfn,'_',' '), iThread, space_opts(iThread).name, optim_names{iThread},...
    iGen+1, max(scores_cur), min(scores_cur), mean(scores_cur), ...
    mean(code_norms),std(code_norms)))
end
pause;
end
%% Best and Worst images in evoking neuron
iThread = 2;
codefnlist = string(ls(fullfile(backupdir,compose("block*_thread%03d_code.mat",iThread-1))));
imgids_all = [];
codes_all = [];
scores_all = [];
for iGen = 1:length(codefnlist)-1
data = load(fullfile(backupdir, codefnlist(iGen)));
imgids_all = [imgids_all; string(data.ids)'];
codes_all = [codes_all; data.codes];
scores_all = [scores_all; Trials.TrialRecord.User.scores_record{iGen+1, iThread}];
end
[score_all_sorted, sortidx_all] = sort(scores_all, 'descend');
%%
topidx = sortidx_all(1:16);
botidx = sortidx_all(end-15:end);
imgs_top = arrayfun(@(nm)imread(fullfile(backupdir, nm)),imgids_all(topidx)','Uni',0);
imgs_bot = arrayfun(@(nm)imread(fullfile(backupdir, nm)),imgids_all(botidx)','Uni',0);
%
figure(19);
ax1 = subplot(1,2,1);
montage(imgs_top,'Parent',ax1,'ThumbnailSize',[256,256])
title(compose("max %.1f min %.1f mean%.1f",...\n norm: noise %.1f(%.1f) class %.2f(%.2f)", ...
    max(scores_all(topidx)), min(scores_all(topidx)), mean(scores_all(topidx))))
ax2 = subplot(1,2,2);
montage(imgs_bot,'Parent',ax2,'ThumbnailSize',[256,256])
title(compose("max %.1f min %.1f mean%.1f",...\n norm: noise %.1f(%.1f) class %.2f(%.2f)", ...
    max(scores_all(botidx)), min(scores_all(botidx)), mean(scores_all(botidx))))
suptitle(compose("%s\n Thread %d Space: %s Optim: %s", strrep(bhvfn,'_',' '), ...
    iThread, strrep(space_opts(iThread).name,'_',' '),strrep(optim_names{iThread},'_',' ')))
%%
img = G.visualize(data.codes);figure;imshow(imtile(img))%codes_all(end-3,:)
%% See the Exploration effects. 
iThread = 2;
codefnlist = string(ls(fullfile(backupdir,compose("block*_thread%03d_code.mat",iThread-1))));

%%
fig = figure(15);
set(fig,'position',[1          41        1920         963])
ax1 = subplot('position',[0.02,0.50,0.32,0.45]);colorbar(ax1);%subplot(1,4,1);
ax2 = subplot('position',[0.02,0.02,0.32,0.45]);%subplot(1,4,2);
ax3 = subplot('position',[0.35,0.02,0.32,0.95]);%subplot(1,4,3);
ax4 = subplot('position',[0.68,0.02,0.32,0.95]);%subplot(1,4,4);
select_n = dummyOptim.mu;
iGen = 1;
for iGen = 1:length(codefnlist)-1
data = load(fullfile(backupdir,codefnlist(iGen)));
imgids = string(data.ids)';
codes_cur = data.codes;
data_next = load(fullfile(backupdir,codefnlist(iGen+1)));
imgids_next = string(data_next.ids)';
scores_cur = Trials.TrialRecord.User.scores_record{iGen+1, iThread};
[score_sorted, sortidx] = sort(scores_cur, 'descend');
recomb_codes = codes_cur(sortidx(1:select_n),:);
newmean = dummyOptim.weights * recomb_codes;
newmeanimg = G.visualize(newmean);
imgs_select = arrayfun(@(nm)imread(fullfile(backupdir, nm)),imgids(sortidx(1:select_n))','Uni',0);
imgs_next = arrayfun(@(nm)imread(fullfile(backupdir, nm)),imgids_next','Uni',0);
img_frame = score_frame_image_arr(imgs_select,score_sorted(1:select_n)');
%
sample_codes = newmean + [0.03 * randn(25, 128), 0.06 * randn(25, 128)];
img_samples = G.visualize(sample_codes);

imshow(imtile(img_frame),'Parent',ax1)
colorbar(ax1);caxis(ax1,score_sorted([select_n,1]))
title(ax1,"Gen "+iGen+" selected images")
imshow(newmeanimg,'Parent',ax2)
title(ax2,"Visualize Weighted mean code")
imshow(imtile(img_samples),'Parent',ax3)
title(ax3,"Potential samples when exploring with smaller step in noise space")
imshow(imtile(imgs_next),'Parent',ax4)
title(ax4,"Real samples / next gen "+(iGen+1))
pause
end