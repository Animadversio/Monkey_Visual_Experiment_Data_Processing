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
P.plotTuneImage = false;true;% Plot best image in each class of tuning
P.savefig = false;true;
P.interactive = false;
P.vis = false;
P.savestat = true;
P.collectstat = true;
global ax1 imcanvs groupdict imgname_uniq
if P.collectstat, S_col = []; end
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

if P.interactive && P.vis
imcanvs = figure(1);clf;ax1 = subplot(1,1,1);set(1,'Name','Stimuli Viewer');hold on
Cord = colororder;
end
%% Best image (or 2) for each category summary 
[maxscore_ingroup,maxid_ingroup] = cellfun(@(idx)max(img_score_mean(idx)),group_uniqidxs);
bestimgnms = arrayfun(@(gi)imgname_uniq{group_uniqidxs{gi}(maxid_ingroup(gi))},1:numel(group_uniqidxs),'uni',0);
bestfns = cellfun(@(imgnm)mapper(imgnm),bestimgnms,'uni',0);
bestimgs = cellfun(@(imgnm)imread(mapper(imgnm)),bestimgnms,'uni',0);

if P.vis
figure(2); set(2,'pos',[705   425   975   435])
X = categorical(grouplabs);
X = reordercats(X,grouplabs);
B = bar(X,group_mean);
ylabel("Evoke - Bsl (sp/s)");xlabel("Image Groups");
title(compose("%s\nImage Group Summary",explabel))
set(gca,'TickLabelInterpreter','none')
box off; B.DataTipTemplate.Interpreter = 'none';
if P.savefig
    saveallform(figdir, "group_summary_prfchan", 2);
end

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
    saveallform(figdir, "group_img_scatter_prfchan", 3);
end

%%
if P.plotTuneImage
figure(4); clf; set(4,'pos',[300   125   1500   720])
T=tiledlayout(2,1,'Padding','compact','TileSpac','compact');
nexttile(T,1);hold on
X = categorical(grouplabs);
X = reordercats(X,grouplabs);
B = bar(X,group_mean);
for gi = 1:length(groupdict)
%     SPminor = 0.5 / numel(group_uniqidxs{gi});
%     xarr = gi + [1:numel(group_uniqidxs{gi})]*SPminor;xarr
    ER = errorbar(repmat(X(gi),numel(group_uniqidxs{gi}),1), img_score_mean(group_uniqidxs{gi}), img_score_sem(group_uniqidxs{gi}),...
        'o','color',[1,0,0,0.7],'LineWidth',0.1,'MarkerSize',3);
    ER.DataTipTemplate.Interpreter = 'none';
end
ylabel("Evoke - Bsl (sp/s)");xlabel("Image Groups");
title(T,compose("%s\nImage Group Summary",explabel))
set(gca,'TickLabelInterpreter','none')
box off; B.DataTipTemplate.Interpreter = 'none';
nexttile(T,2)
% imshow(imtile(bestimgs,'GridSize',[1,numel(bestimgs)]))
if numel(bestimgs)>13, imshow(tile_zigzag(bestimgs, 2));
else, imshow(tile_zigzag(bestimgs, 1));end
if P.savefig
    saveallform(figdir, "group_summary_bestImg_prfchan", 4);
end
end
end
%%
if P.savestat
    S = struct();
    S.Animal = Animal;
    S.meta = meta;
    S.meta.expday = expday;
    S.meta.fdrnm = fdrnm;
    S.meta.figdir = figdir;
    S.meta.explabel = explabel;
    S.imageName = Trials.imageName;
    S.unit.unit_num_arr = unit_num_arr;
    S.unit.unit_name_arr = unit_name_arr;
    S.unit.prefchan = prefchan;
    S.unit.prefunit = prefUi;
    S.unit.prefchan_id = prefchan_id;

    S.stim.imgfn_mapper = mapper;
    S.stim.imgfps = imgfps;
    S.stim.impos = impos;
    S.stim.imsize_pix = imsize_pix;
    S.stim.imsize_deg = imsize_deg;
    S.stim.imgname_uniq = imgname_uniq;
    S.stim.img_idxarr = img_idxarr;
    S.stim.groupdict = groupdict;
    S.stim.group_uniqidxs = group_uniqidxs;
    S.stim.grouplabs = grouplabs;
    S.stim.group_idxarr = group_idxarr;
    
    S.prefresp.psth_arr = psth_arr;
    S.prefresp.img_evoke_col = score_col;
    S.prefresp.img_score_col = score_bsl_col;
    S.prefresp.group_score_col = group_score_col;
    S.prefresp.img_mean = img_score_mean;
    S.prefresp.img_sem = img_score_sem;
    S.prefresp.group_mean = group_mean;
    S.prefresp.group_sem = group_sem;
    S.prefresp.bestimgnms = bestimgnms;
    save(fullfile(figdir, "ExpStat.mat"), "S");
    if P.collectstat, S_col = [S_col,S]; end
end
%% Saving zone 

end
%%
% figure(5);imshow(tile_zigzag(bestimgs, 3))
function img_mtg = tile_zigzag(imgcol,rown)
imgpix = 256;
brdrpix=2;
imgN = numel(imgcol);
imgrow_mtg = {};
for rowi = 1:rown
imgrow_mtg{rowi}=imtile(imgcol(rowi:rown:end),'GridSize',[1,numel(imgcol(rowi:rown:end))],'BorderSize', brdrpix);
end
padimgrow_mtg = {};
rowlen = 1 + (imgN - 1) * (1 / rown);
rowpix = ceil(rowlen*(imgpix+2*brdrpix));
for rowi = 1:rown
    prepad = round((rowi-1)/rown*(imgpix+2*brdrpix));
    postpad = rowpix - prepad - size(imgrow_mtg{rowi}, 2);
    pre_padarr = zeros(imgpix+2*brdrpix,prepad,3,class(imgrow_mtg{rowi}));
    post_padarr = zeros(imgpix+2*brdrpix,postpad,3,class(imgrow_mtg{rowi}));
    padimgrow_mtg{rowi}=cat(2,pre_padarr,imgrow_mtg{rowi},post_padarr);
end
img_mtg = cat(1,padimgrow_mtg{:});
end

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
    groups("fact_tsr") = find(contains(imgname_uniq,"facttsr_")...
                             & ~contains(imgname_uniq,"facttsr_map_patchshuffle")...
                             & ~contains(imgname_uniq,"facttsr_maponly_patchshuffle"));
    groups("fact_tsr_mapfeat_patchshfl") = find(contains(imgname_uniq,"facttsr_map_patchshuffle"));
    groups("fact_tsr_maponly_patchshfl") = find(contains(imgname_uniq,"facttsr_maponly_patchshuffle"));
    groups("full_tsr") = find(contains(imgname_uniq,"tsr_") & ~contains(imgname_uniq,"facttsr_"));
    for ci = min(factor_id):max(factor_id)
    groups(compose("fact%d_cntpnt",ci)) = find(contains(imgname_uniq,compose("fac%d_cntpnt_",ci)) & ~contains(imgname_uniq,compose("fac%d_cntpnt_shuffle",ci)));
    groups(compose("fact%d_full",ci)) = find(contains(imgname_uniq,compose("fac%d_full_",ci)) & ~contains(imgname_uniq,compose("fac%d_full_shuffle",ci)));
    groups(compose("fact%d_wmap",ci)) = find(contains(imgname_uniq,compose("fac%d_map_",ci))...
                                          & ~contains(imgname_uniq,compose("fac%d_map_shuffle",ci))...
                                          & ~contains(imgname_uniq,compose("fac%d_map_map_patchshuffle",ci))...
                                          & ~contains(imgname_uniq,compose("fac%d_map_maponly_patchshuffle",ci)));
    groups(compose("fact%d_wmap_maponly_patchshfl",ci)) = find(contains(imgname_uniq,compose("fac%d_map_maponly_patchshuffle",ci)));
    groups(compose("fact%d_cntpnt_feat_shfl",ci)) = find(contains(imgname_uniq,compose("fac%d_cntpnt_shuffle",ci)));
    groups(compose("fact%d_full_feat_shfl",ci)) = find(contains(imgname_uniq,compose("fac%d_full_shuffle",ci)));
    groups(compose("fact%d_wmap_feat_shfl",ci)) = find(contains(imgname_uniq,compose("fac%d_map_shuffle",ci)));
    groups(compose("fact%d_wmap_featmap_patchshfl",ci)) = find(contains(imgname_uniq,compose("fac%d_map_map_patchshuffle",ci)));
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