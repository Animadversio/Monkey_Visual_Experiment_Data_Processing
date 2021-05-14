% Analyze the new Evolution decomposition experiments. 
% 

%%
Animal = "Alfa";Set_Path;
ftrrows = find(...
               contains(ExpRecord.expControlFN,["select"])&...
               contains(ExpRecord.Exp_collection,["FC6Evol_Decomp"])&...
               ~isnan(ExpRecord.Expi)...
               );
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows, Animal, false);
%%
set(0,'defaultTextInterpreter','none');
saveroot = "O:\corrFeatVis_FactorPredict";
P = struct();
P.plotTuneCurve = false;% TO IMplement.
P.savefig = true;
P.interactive = false;
global ax1 imcanvs groupdict imgname_uniq
for Triali = 1:numel(meta_new)%-3:numel(meta_new)-2
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
impos = Trials.TrialRecord.User.uniquePositions;
imsize_pix = Trials.width(1,1);
imsize_deg = ExpRecord.stim_size(exp_rowi);
% imsize_deg = Trials.width(1,1) / Trials.MLConfig.PixelsPerDegree(1); %
% precise value
unit_num_arr = meta.unitID;
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
imgname_uniq = unique(Trials.imageName); 
prefchan = Trials.TrialRecord.User.prefChan;
prefUi = ExpRecord.pref_unit(exp_rowi);
prefchan_ids = find(meta.spikeID == prefchan & meta.unitID>0);
prefchan_id = prefchan_ids(prefUi);%(ui)
if strcmp(meta.ephysFN, 'Alfa-22042021-007'), prefchan_id = prefchan_ids(2);end
if strcmp(meta.ephysFN, 'Alfa-27042021-007'), prefchan_id = prefchan_ids(1);end
[imgfps, mapper] = map2fullpath(string(imgname_uniq), meta.stimuli); % full path for images to load and show

stimparts = split(meta.stimuli,"\");
expday = datetime(meta.expControlFN(1:6),'InputFormat','yyMMdd');
fdrnm = compose("%s-Ch%02d", stimparts{end}, prefchan);
figdir = fullfile(saveroot, fdrnm);
mkdir(figdir)
fprintf("Saved to exp dir\n%s\n",figdir)
explabel = compose("%s PrefChan %s size %.1f deg pos [%.1f %.1f]",stimparts{end}(1:end-7),unit_name_arr(prefchan_id),...
    imsize_deg,impos(1),impos(2));
%% Sort out the images 
img_idxarr = cellfun(@(imgnm)find(strcmp(Trials.imageName,imgnm)),imgname_uniq,'uni',0); % trial id each image 
% (img name exact match, DONT USE contains, use strcmp)
%% Group up images to present
[groupdict, group_uniqidxs, grouplabs] = uniqimg_grouping(imgname_uniq, "type_shfl");
% group_idxarr = cellfun(@(idxuniq) find(contains(Trials.imageName,
% imgname_uniq(idxuniq))), group_uniqidxs, 'uni', 0); % WRONG OBsolete
group_idxarr = cellfun(@(idxuniq) cell2mat(img_idxarr(idxuniq)), group_uniqidxs, 'uni', 0); % trial index for all images in a group
%%
evkwdw = [51:200];
bslwdw = [1:50];
for iCh = prefchan_id
% score, evoked per trial
score_trials  = squeeze(mean(rasters(iCh, evkwdw, :),[2]) - mean(rasters(iCh, bslwdw, :),[2])); % - single trial bsl score
evoke_trials  = squeeze(mean(rasters(iCh, evkwdw, :),[2])); % evoked firing rate
evkbsl_trials = squeeze(mean(rasters(iCh, evkwdw, :),[2])) - mean(rasters(iCh, bslwdw, :),'all'); % evk - bsl score using overall baseline

% sort rsp into image names
psth_arr = cell2mat(cellfun(@(idx)mean(rasters(iCh,:,idx),3),img_idxarr,'Un',0)); % mean psth for each img (imgN, T=200)
psth_arr_sem = cell2mat(cellfun(@(idx)sem(rasters(iCh,:,idx),3),img_idxarr,'Un',0)); % psth sem for each img (imgN, T=200)
score_col = cellfun(@(idx)squeeze(mean(rasters(iCh,evkwdw,idx),2)),img_idxarr,'Un',0); % cell array of trial scores each img
score_bsl_col = cellfun(@(idx)squeeze(mean(rasters(iCh,evkwdw,idx),2)) ...
                        - mean(rasters(iCh,bslwdw,idx),'all'),img_idxarr,'Un',0);
% score summary for each image
img_score_mean = cellfun(@mean, score_bsl_col); % trial mean for each img (imgN, 1)
img_score_sem = cellfun(@sem, score_bsl_col); % trial sem for each img (imgN, 1)
end
%%
group_score_col = cellfun(@(idx) evkbsl_trials(idx), group_idxarr,'uni',0);
group_mean = cellfun(@mean, group_score_col);
group_sem = cellfun(@sem, group_score_col);

imcanvs = figure(1);clf;ax1 = subplot(1,1,1);set(1,'Name','Stimuli Viewer');hold on
Cord = colororder;

figure(2); set(2,'pos',[705   425   975   435])
X = categorical(grouplabs);
X = reordercats(X,grouplabs);
B = bar(X,group_mean);
ylabel("Evoke - Bsl (sp/s)");xlabel("Image Groups");
title(compose("%s\nImage Group Summary",explabel))
set(gca,'TickLabelInterpreter','none')
box off; B.DataTipTemplate.Interpreter = 'none';
if P.savefig
savefig(2,fullfile(figdir,"group_summary_prfchan.fig"))
saveas(2,fullfile(figdir,"group_summary_prfchan.png"))
end
%
figure(3);clf; set(3,'pos',[1000    150   1325   825])
subplot(211);hold on; 
E = errorbar([],group_mean,group_sem,'o','LineWidth',1.5);
row = dataTipTextRow('groups',grouplabs);
E.DataTipTemplate.DataTipRows(end+1) = row;
E.DataTipTemplate.Interpreter = 'none';
xticks(1:numel(group_mean)); xticklabels(grouplabs); xtickangle(30);
set(gca,'TickLabelInterpreter','none');box off
ylabel("Evoke - Bsl");suptitle(compose("%s\nImage Group Summary",explabel))%xlabel("Image Groups");
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
if P.interactive
    set(imcanvs.Children, 'UserData', struct("imgfps",string(imgfps), "showed_ids",[], "img_col",[], 'groupdict', groupdict))
    set(E,'ButtonDownFcn',@selectgroup, 'HitTest','on')
end
% figure(4); set(4,'pos',[1353  547  1075  375])
subplot(212)
ERR = errorbar([],img_score_mean,img_score_sem); box off
ylabel("Evoke - Bsl (sp/s)");xlabel("Image Id");
row = dataTipTextRow('imnm',imgname_uniq);
ERR.DataTipTemplate.DataTipRows(end+1) = row;
ERR.DataTipTemplate.Interpreter = 'none';
if P.interactive
    set(ERR,'ButtonDownFcn',@selectpoint, 'HitTest','on')
end
if P.savefig
savefig(3,fullfile(figdir,"group_img_scatter_prfchan.fig"))
saveas(3,fullfile(figdir,"group_img_scatter_prfchan.png"))
end
%% Best image summary 

%% Saving zone 
% save(fullfile(figdir,compose(".mat")),)

end
%%

function [groups, uniqnm_idxs, grouplabs] = uniqimg_grouping(imgname_uniq, mode)
% group the uniq image names into groups
% groups: containers.Map, dict
% uniqnm_idxs: cell array of indices array
% grouplabs: cell array of chars 
% Example:
%  [groups, group_uniqidxs, grouplabs] = uniqimg_grouping(imgname_uniq, "type");
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
        if numel(uniqnm_idxs{gi})==0, remove(groups,grouplabs{gi}); end
        fprintf("%s: %d\n",grouplabs{gi},numel(uniqnm_idxs{gi}))
    end
    grouplabs = groups.keys;
    uniqnm_idxs = groups.values;
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
    groups(compose("fact%d_wmap",ci)) = find(contains(imgname_uniq,compose("fac%d_map_",ci))...
                                          & ~contains(imgname_uniq,compose("fac%d_map_shuffle",ci))...
                                          & ~contains(imgname_uniq,compose("fac%d_map_map_patchshuffle",ci))...
                                          & ~contains(imgname_uniq,compose("fac%d_map_maponly_patchshuffle",ci)));
    groups(compose("fact%d_wmap_patchshflmap",ci)) = find(contains(imgname_uniq,compose("fac%d_map_maponly_patchshuffle",ci)));
    groups(compose("fact%d_cntpnt_shfl",ci)) = find(contains(imgname_uniq,compose("fac%d_cntpnt_shuffle",ci)));
    groups(compose("fact%d_full_shfl",ci)) = find(contains(imgname_uniq,compose("fac%d_full_shuffle",ci)));
    groups(compose("fact%d_wmap_shfl",ci)) = find(contains(imgname_uniq,compose("fac%d_map_shuffle",ci)));
    groups(compose("fact%d_wmap_patchshfl",ci)) = find(contains(imgname_uniq,compose("fac%d_map_map_patchshuffle",ci)));
    end
    grouplabs = groups.keys;
    uniqnm_idxs = groups.values;
    fprintf("\nGroup summary:\n")
    for gi = 1:length(groups)
        if numel(uniqnm_idxs{gi})==0, remove(groups,grouplabs{gi}); end
        fprintf("%s: %d\n",grouplabs{gi},numel(uniqnm_idxs{gi}))
    end
    grouplabs = groups.keys;
    uniqnm_idxs = groups.values;
end
end

function [coordSelected, minIdx] = selectpoint(hObj, event)
% call back function
global ax1 imgname_uniq
% disp(event)
% disp(hObj)
X = hObj.XData; 
Y = hObj.YData; 
coordinates = [reshape(X,[],1),reshape(Y,[],1)];
pt = event.IntersectionPoint; 
x = pt(:,1);y = pt(:,2);
dist = pdist2([x,y],coordinates);            %distance between your selection and all points
[~, minIdx] = min(dist);            % index of minimum distance to points
coordSelected = coordinates(minIdx,:);
selected_id = int32(coordSelected(1));
fprintf("Id %03d score %.1f, name %s\n",selected_id,coordSelected(2),imgname_uniq{selected_id})
dt = datatip(hObj,x,y, 'DeleteFcn', @(dtob,~)update_show_imgs(ax1, selected_id));
update_show_imgs(ax1, selected_id)
end

function [coordSelected, minIdx] = selectgroup(hObj, event)
% call back function
global imcanvs ax1 imgname_uniq groupdict
% disp(event)
% disp(hObj)
groupnms = groupdict.keys;
group_uniqids = groupdict.values;
X = hObj.XData; 
Y = hObj.YData; 
coordinates = [reshape(X,[],1),reshape(Y,[],1)];
pt = event.IntersectionPoint; 
x = pt(:,1);y = pt(:,2);
dist = pdist2([x,y],coordinates);            %distance between your selection and all points
[~, minIdx] = min(dist);            % index of minimum distance to points
coordSelected = coordinates(minIdx,:);
selected_id = int32(coordSelected(1));
fprintf("Id %03d score %.1f, name %s\n",selected_id,coordSelected(2),groupnms{selected_id})
dt = datatip(hObj,x,y,'DeleteFcn', @(dtob,~)update_show_imgs(imcanvs.Children, group_uniqids{selected_id}));
update_show_imgs(imcanvs.Children, group_uniqids{selected_id})
end

function update_show_imgs(ax1, imgid)
% disp(imgid)
for imgi = reshape(imgid,1,[])
    if any(ax1.UserData.showed_ids==imgi) 
        % if the clicked img id is in the collection, remove it
        idx = find(ax1.UserData.showed_ids==imgi);
        ax1.UserData.showed_ids(idx) = [];
        ax1.UserData.img_col(idx) = [];
        continue
    else
        % if the clicked img id isn't in the collection, Append it to the end
        ax1.UserData.showed_ids(end+1) = imgi;
        img = imread(ax1.UserData.imgfps(imgi));
        ax1.UserData.img_col{end+1} = img;
    end
end
% axes(ax1);
cla(ax1)
if numel(ax1.UserData.img_col) > 0 
imshow(imtile(ax1.UserData.img_col),'Parent',ax1)
% montage(ax1.UserData.img_col,'Parent',ax1)
end
end

function [imgfps, mapper] = map2fullpath(picnm_arr, imgdir)
% map the imgnm without suffix to the full path, for all images in a folder.
imgnm_wsfx = deblank(string([ls([imgdir+"\*.png"]);ls([imgdir+"\*.jpg"]);ls([imgdir+"\*.jpeg"])]));
[~,imgnms,sfxs] = arrayfun(@fileparts, imgnm_wsfx); % older version matlab
mapper = containers.Map(imgnms,fullfile(imgdir, imgnm_wsfx));
imgfps = cellfun(@(I)mapper(I),picnm_arr,'uni',0); % image full path of tht input picnms 
end