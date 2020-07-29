%% Run this after experiment using torchBigGAN. 
%  This will empty the GPU resources for other programs. Crucial for
%  Hessian computation which requires space on GPU
py.torch.cuda.empty_cache()
%% Experiment Visualization
bhvfn = "200728_Beto_generate_BigGAN.bhv2";
[data,MLConfig,TrialRecord] = mlread(bhvfn);
space_opts = TrialRecord.User.space_opts;
optim_names = TrialRecord.User.evoConfiguration(:,5);
% save(fullfile(TrialRecord.User.newPicsHome, "space_opts.mat"), "space_opts");
backupdir = TrialRecord.User.newPicsHome;
%% Visualize the experiments using scores and images
% imgext = ".bmp";
%%  Examine images sorted by scores Generation by Generation
iThread = 2;
codefnlist = string(ls(fullfile(backupdir,compose("block*_thread%03d_code.mat",iThread-1))));
fig = figure(14);set(fig,'position',[410       -1085        1013         969])
ax = subplot(1,1,1);
iGen = 1;
for iGen = 1:length(codefnlist)-1
data = load(fullfile(backupdir,codefnlist(iGen)));
imgids = string(data.ids);
codes_cur = data.codes;
scores_cur = TrialRecord.User.scores_record{iGen+1, iThread};
[score_sorted, sortidx] = sort(scores_cur, 'descend');
imgs_cur = arrayfun(@(nm)imread(fullfile(backupdir, nm)),imgids(sortidx)','Uni',0);
montage(imgs_cur,'Parent',ax)
if strcmp(space_opts(2).name,'BigGAN')
noise_norms = norm_axis(codes_cur(:,1:128),2);
class_norms = norm_axis(codes_cur(:,129:end),2);
title(compose("%s Thread %d Space: %s Optim: %s\nGen %d max %.1f min %.1f mean%.1f\n norm: noise %.1f(%.1f) class %.2f(%.2f)", ...
    strrep(bhvfn,'_',' '), iThread, space_opts(iThread).name, optim_names{iThread},...    
    iGen+1, max(scores_cur), min(scores_cur), mean(scores_cur), ...
    mean(noise_norms),std(noise_norms), mean(class_norms),std(class_norms)))
else
code_norms = norm_axis(codes_cur,2);
title(compose("%s Thread %d Space: %s Optim: %s\nGen %d max %.1f min %.1f mean%.1f\n norm: %.2f(%.2f)", ...
    strrep(bhvfn,'_',' '), iThread, space_opts(iThread).name, optim_names{iThread},...
    iGen+1, max(scores_cur), min(scores_cur), mean(scores_cur), ...
    mean(code_norms),std(code_norms)))
end
pause;
end

%% Visualize Favorite images in all generations
iThread = 2;
codefnlist = string(ls(fullfile(backupdir,compose("block*_thread%03d_code.mat",iThread-1))));
imgids_all = [];
codes_all = [];
scores_all = [];
for iGen = 1:length(codefnlist)-1
data = load(fullfile(backupdir, codefnlist(iGen)));
imgids_all = [imgids_all; string(data.ids)'];
codes_all = [codes_all; codes_cur];
scores_all = [scores_all; TrialRecord.User.scores_record{iGen+1, iThread}];
end
[score_all_sorted, sortidx_all] = sort(scores_all, 'descend');
topidx = sortidx_all(1:16);
botidx = sortidx_all(end-15:end);
imgs_top = arrayfun(@(nm)imread(fullfile(backupdir, nm)),imgids_all(topidx)','Uni',0);
imgs_bot = arrayfun(@(nm)imread(fullfile(backupdir, nm)),imgids_all(botidx)','Uni',0);
%
figure(19);
ax1 = subplot(1,2,1);
montage(imgs_top,'Parent',ax1)
title(compose("max %.1f min %.1f mean%.1f",...\n norm: noise %.1f(%.1f) class %.2f(%.2f)", ...
    max(scores_all(topidx)), min(scores_all(topidx)), mean(scores_all(topidx))))
ax2 = subplot(1,2,2);
montage(imgs_bot,'Parent',ax2)
title(compose("max %.1f min %.1f mean%.1f",...\n norm: noise %.1f(%.1f) class %.2f(%.2f)", ...
    max(scores_all(botidx)), min(scores_all(botidx)), mean(scores_all(botidx))))
suptitle(compose("%s\n Thread %d Space: %s Optim: %s", strrep(bhvfn,'_',' '), ...
    iThread, strrep(space_opts(iThread).name, '_', ' '),optim_names{iThread}))