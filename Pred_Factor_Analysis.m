FC6Evol_Decomp
% Analyze the new Evolution decomposition experiments. 
% 

%%
Animal = "Alfa";Set_Path;
ftrrows = find(...contains(ExpRecord.expControlFN,["generate_BigGAN_cosine"])&...
               contains(ExpRecord.Exp_collection,["FC6Evol_Decomp"])&...
               ~isnan(ExpRecord.Expi)...
               );
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows, Animal, false);
%%
set(0,'defaultTextInterpreter','none');
saveroot = "O:\corrFeatVis_FactorPredict";
P = struct();
P.plotTuneCurve = false;% TO IMplement.
P.savefig = false;
for Triali = 12%numel(meta_new)-3:numel(meta_new)-2
%
meta = meta_new{Triali};
rasters = rasters_new{Triali};
% lfps = lfps_new{Triali};
Trials = Trials_new{Triali};
if contains(meta.ephysFN,"Alfa"),Animal = "Alfa";elseif contains(meta.ephysFN,"Beto"),Animal = "Beto";end
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
Expi = ExpRecord.Expi(exp_rowi);
fprintf("Processing  Exp %d Notes:\n",Expi)
fprintf([ExpRecord.comments{exp_rowi},'\n'])

unit_num_arr = meta.unitID;
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
imgname_uniq = unique(Trials.imageName); 
prefchan = Trials.TrialRecord.User.prefChan;
prefchan_ids = find(meta.spikeID == prefchan & meta.unitID>0);
if strcmp(meta.ephysFN, 'Alfa-22042021-007'), prefchan_ids = prefchan_ids(2);end

stimparts = split(meta.stimuli,"\");
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
fdrnm = compose("%s-Ch%02d", stimparts{end}, prefchan);
figdir = fullfile(saveroot, fdrnm);
mkdir(figdir)
fprintf("Saved to exp dir\n%s\n",figdir)
%% Sort out the images 
img_idxarr = cellfun(@(imgnm)find(strcmp(Trials.imageName,imgnm)),imgname_uniq,'uni',0); % trial id each image
%%
evkwdw = [51:200];
bslwdw = [1:50];
for iCh = prefchan_ids
% score, evoked per trial
score_trials  = squeeze(mean(rasters(iCh, evkwdw, :),[2]) - mean(rasters(iCh, bslwdw, :),[2]));
evoke_trials  = squeeze(mean(rasters(iCh, evkwdw, :),[2]));
evkbsl_trials = squeeze(mean(rasters(iCh, evkwdw, :),[2])) - mean(rasters(iCh, bslwdw, :),'all');

% sort into image names
psth_arr = cell2mat(cellfun(@(idx)mean(rasters(iCh,:,idx),3),img_idxarr,'Un',0));
psth_arr_sem = cell2mat(cellfun(@(idx)sem(rasters(iCh,:,idx),3),img_idxarr,'Un',0));
score_col = cellfun(@(idx)squeeze(mean(rasters(iCh,51:200,idx),2)),img_idxarr,'Un',0);
score_bsl_col = cellfun(@(idx)squeeze(mean(rasters(iCh,evkwdw,idx),2)) ...
                        - mean(rasters(iCh,bslwdw,idx),'all'),img_idxarr,'Un',0);
% summary for each image
img_score_mean = cellfun(@mean, score_bsl_col);
img_score_sem = cellfun(@sem, score_bsl_col);
figure; errorbar([],img_score_mean,img_score_sem); box off
end

%%
%% Group up images to present
[groupdict, group_uniqidxs, grouplabs] = uniqimg_grouping(imgname_uniq, "type");
group_idxarr = cellfun(@(idxuniq) find(contains(Trials.imageName, imgname_uniq(idxuniq))), group_uniqidxs, 'uni', 0);
group_score_col = cellfun(@(idx) evkbsl_trials(idx), group_idxarr,'uni',0);
group_mean = cellfun(@mean, group_score_col);
group_sem = cellfun(@sem, group_score_col);
figure;
X = categorical(grouplabs);
X = reordercats(X,grouplabs);
B = bar(X,group_mean);
set(gca,'TickLabelInterpreter','none')
box off;B.DataTipTemplate.Interpreter = 'none';
%%
Cord = colororder;
figure;hold on
E = errorbar([],group_mean,group_sem,'o','LineWidth',1.5);
row = dataTipTextRow('groups',grouplabs);
E.DataTipTemplate.DataTipRows(end+1) = row;
E.DataTipTemplate.Interpreter = 'none';
xticks(1:numel(group_mean)); xticklabels(grouplabs); xtickangle(30);
set(gca,'TickLabelInterpreter','none');box off
ylabel("Active - Bsl");title("Image Group Summary")
for gi = 1:length(groupdict)
    SPminor = 0.5 / numel(group_uniqidxs{gi});
    xarr = gi + [1:numel(group_uniqidxs{gi})]*SPminor;
    ER = errorbar(xarr, img_score_mean(group_uniqidxs{gi}), img_score_sem(group_uniqidxs{gi}),...
        'o','color',[Cord(2,:),0.4],'LineWidth',0.1,'MarkerSize',3);
    row = dataTipTextRow('imnm',imgname_uniq(group_uniqidxs{gi}));
    ER.DataTipTemplate.DataTipRows(end+1) = row;
    ER.DataTipTemplate.Interpreter = 'none';
end
xlim([0.5,length(groupdict)+1])
%%

end
%%
[groups, group_uniqidxs, grouplabs] = uniqimg_grouping(imgname_uniq, "type");%uniqnm_idxs, grouplabs
function [groups, uniqnm_idxs, grouplabs] = uniqimg_grouping(imgname_uniq, mode)
% group the uniq image names into groups
% groups: containers.Map, dict
% uniqnm_idxs: cell array of indices array
% grouplabs: cell array of chars 
if ~sum(contains(imgname_uniq,"_shuffle"))
    mode = "type";
else
    mode = "type_shfl";    
end
if strcmp(mode, "type")
    imgnm_re = regexp(imgname_uniq,"fac(\d)_",'tokens');
    imgnm_re(cellfun(@isempty,imgnm_re)) = [];
    factor_id = cellfun(@(C)str2num(C{1}{1}),imgnm_re);
    cNum = max(factor_id)+1;
    groups = containers.Map();
    groups("evol_best") = find(contains(imgname_uniq,"evol_best_"));
    groups("fact_tsr") = find(contains(imgname_uniq,"facttsr_"));
    groups("full_tsr") = find(contains(imgname_uniq,"tsr_") & ~contains(imgname_uniq,"facttsr_"));
    for ci = min(factor_id):max(factor_id)
    groups(compose("fact%d_cntpnt",ci)) = find(contains(imgname_uniq,compose("fac%d_cntpnt_",ci)));
    groups(compose("fact%d_full",ci)) = find(contains(imgname_uniq,compose("fac%d_full_",ci)));
    groups(compose("fact%d_wmap",ci)) = find(contains(imgname_uniq,compose("fac%d_map_",ci)));
    end
    grouplabs = groups.keys;
    uniqnm_idxs = groups.values;
    fprintf("\nGroup summary:\n")
    for gi = 1:length(groups)
        fprintf("%s: %d\n",grouplabs{gi},numel(groups(grouplabs{gi})))
    end
elseif strcmp(mode, "type_shfl")
    imgnm_re = regexp(imgname_uniq,"fac(\d)_",'tokens');
    imgnm_re(cellfun(@isempty,imgnm_re)) = [];
    factor_id = cellfun(@(C)str2num(C{1}{1}),imgnm_re);
    cNum = max(factor_id)+1;
    groups = containers.Map();
    groups("evol_best") = find(contains(imgname_uniq,"evol_best_"));
    groups("fact_tsr") = find(contains(imgname_uniq,"facttsr_"));
    groups("full_tsr") = find(contains(imgname_uniq,"tsr_") & ~contains(imgname_uniq,"facttsr_"));
    for ci = min(factor_id):max(factor_id)
    groups(compose("fact%d_cntpnt",ci)) = find(contains(imgname_uniq,compose("fac%d_cntpnt_",ci)) & ~contains(imgname_uniq,compose("fac%d_cntpnt_shuffle",ci)));
    groups(compose("fact%d_full",ci)) = find(contains(imgname_uniq,compose("fac%d_full_",ci)) & ~contains(imgname_uniq,compose("fac%d_full_shuffle",ci)));
    groups(compose("fact%d_wmap",ci)) = find(contains(imgname_uniq,compose("fac%d_map_",ci)) & ~contains(imgname_uniq,compose("fac%d_map_shuffle",ci)));
    groups(compose("fact%d_cntpnt_shfl",ci)) = find(contains(imgname_uniq,compose("fac%d_cntpnt_shuffle",ci)));
    groups(compose("fact%d_full_shfl",ci)) = find(contains(imgname_uniq,compose("fac%d_full_shuffle",ci)));
    groups(compose("fact%d_wmap_shfl",ci)) = find(contains(imgname_uniq,compose("fac%d_map_shuffle",ci)));
    end
    grouplabs = groups.keys;
    uniqnm_idxs = groups.values;
    fprintf("\nGroup summary:\n")
    for gi = 1:length(groups)
        fprintf("%s: %d\n",grouplabs{gi},numel(groups(grouplabs{gi})))
    end
end
end
