%% vis_demo_Dual_Evolution_Factorization
%  This new analysis script is written in more modular fashion. 
%  Each experiment in a bundle is processed first to be a struct and then
%  visualization and statistics are created from it. 
%  [In development]
Animal="Beto"; Set_Path;
rootdir = "O:\CorrFactor_natvalid";
%%
rowlist = find(contains(ExpRecord.ephysFN,"04012022"));
[meta_new,rasters_new,~,Trials_new] = loadExperiments(rowlist, Animal, false);
%% Evolution Experiments 
S = Evol_BigGAN_FC6_Collect_Stats_fun(meta_new(2), rasters_new(2), Trials_new(2));
fdrnm = S.meta.fdrnm;
expdir = fullfile(rootdir,fdrnm); 
mkdir(expdir)
%% Visualize Evolution Experiments 
vidpath = Evol_BigGAN_FC6_Animation_fun(S);
copyfile(vidpath, expdir)
%% Selectivity Experiments 
%% Format the selectivity Exepriemnts data
Snat = selectivity_Collect_Stats_fun(meta_new(3:5), rasters_new(3:5), Trials_new(3:5));
% %% The rig version Obsolete
% D = vis_GUI_selective_prefchan("N:\Data-Behavior (BHV2)\211220_Beto_selectivity_basic.bhv2", ...
%                 "N:\Stimuli\2021-EvolDecomp\2021-12-20-Beto-01-natvalid");
%% Save the formated mat structs 
savenm = [Snat(1).meta.ephysFN,Snat(3).meta.ephysFN(end-3:end)];
save(fullfile(expdir,compose("%s_selStats.mat",savenm)),'Snat')
savenm = [S(1).meta.ephysFN];
save(fullfile(expdir,compose("%s_evolStats.mat",savenm)),'S')
fprintf("Stats Structures of Evolution and Validation (Selectivity) saved to %s\n", expdir)

%% Visualize and Print the Selectivity Experimental Response. 
unitidx = S.units.pref_chan_id(1,1);
chan = S.units.pref_chan(1);
unit = S.units.unit_in_pref_chan(1);%2
unitidx = find((S.units.spikeID==chan)&(S.units.unit_num_arr==unit));
unitlab = S.units.unit_name_arr(unitidx);
for Si = 1:numel(Snat)
stimpath = Snat(Si).meta.stimuli;
if contains(stimpath, "natvalid")
    sellabel = "natvalid";
elseif contains(stimpath, "-decomp-BigGAN")
    sellabel = "BigGANFactPred";
elseif contains(stimpath, "-decomp")
    sellabel = "FC6FactPred";
else 
    error("un identified validation stimuli space (or stimuli set naming convention changed.)")
end
[figh,~,~] = vis_selectivity_static(Snat(Si),true,chan,unit,'top',30);
saveallform(expdir, compose("sel_topK_%s_sorted_%s", sellabel, unitlab), figh, ["png","pdf"])
[figh,~,~] = vis_selectivity_static(Snat(Si),false,chan,unit,'top',30);
saveallform(expdir, compose("sel_topK_%s_%s", sellabel, unitlab), figh, ["png","pdf"])
[figh,~,~] = vis_selectivity_static(Snat(Si),true,chan,unit,'bottom',30);
saveallform(expdir, compose("sel_bottomK_%s_sorted_%s", sellabel, unitlab), figh, ["png","pdf"])
end

%% Interactive Exploratory Analysis of selectivity experiment
%% Natural selectivity experiments
vis_selectivity_interact(Snat(1),true)
%% BigGAN factorization experiments
vis_selectivity_interact(Snat(2),true)
%% FC6 factorization experiments
vis_selectivity_interact(Snat(3),true)
vis_selectivity_interact(Snat(2),true,2,2)
%% The behavior of prefchan in the three validation experiments 
vis_selectivity_interact(Snat(1),true)
vis_selectivity_interact(Snat(2),true)
vis_selectivity_interact(Snat(3),true)
%%
vis_selectivity_static(Snat(1),true,2,1,'bottom',16)
vis_selectivity_static(Snat(2),false,2,1,'bottom',16)
vis_selectivity_static(Snat(3),true,2,1,'bottom',16)


function [h1,ax1,ax2] = vis_selectivity_interact(Snat,sorted,chan,unit)
% thin wrapper around the interactive function to produce arguments for it.
% Parameters:
%   sorted: Boolean, sort the responses in ascending order (default false)
%   chan: Int, channel number to plot 
%   unit: Int, unit number in that channel. 
if nargin < 2, sorted = false; end
if nargin < 3, do_prefchan =true; 
else, do_prefchan=false; end 
if nargin < 4, unit = 1;end
if do_prefchan
    unit_idx = Snat.units.pref_chan_id;
else
    unit_idx = find((Snat.units.spikeID == chan) & ...
                    (Snat.units.unit_num_arr == unit));
end
rspvec = Snat.resp.meanMat(unit_idx, :);
semvec = Snat.resp.semMat(unit_idx, :);
bsl = Snat.resp.bslmean(unit_idx);
imgnm_uniq = Snat.stim.imgname_uniq;
imgfps = Snat.stim.imgfps;
if sorted
    [rsp_sort, sortidx] = sort(rspvec);
    rspvec = rspvec(sortidx);
    semvec = semvec(sortidx);
    imgnm_uniq = imgnm_uniq(sortidx);
    imgfps = imgfps(sortidx);
end
[h1,ax1,ax2]=vis_selectivity_interact_core(imgnm_uniq, imgfps,...
                    rspvec,semvec,bsl);
titlestr = compose("%s chan %s",Snat.meta.ephysFN,Snat.units.unit_name_arr(unit_idx));
sgtitle(h1,titlestr)
end

%% Static visualization 
function [h1,ax1,ax2] = vis_selectivity_static(Snat,sorted,chan,unit,order,K)
% thin wrapper around the interactive function to produce arguments for it.
% Parameters:
%   sorted: Boolean, sort the responses in ascending order (default false)
%   chan: Int, channel number to plot
%   unit: Int, unit number in that channel. 
%   order: "top" (default) or "bottom". 
%   K: Top K images to show in 2nd panel.. 
if nargin < 2, sorted = false; end
if nargin < 3, do_prefchan =true; 
else, do_prefchan=false; end 
if nargin < 4, unit = 1;end
if nargin < 5, order = "top";end
if nargin < 6, K = 16; end
if do_prefchan
    unit_idx = Snat.units.pref_chan_id;
else
    unit_idx = find((Snat.units.spikeID == chan) & ...
                    (Snat.units.unit_num_arr == unit));
end
rspvec = Snat.resp.meanMat(unit_idx, :);
semvec = Snat.resp.semMat(unit_idx, :);
bsl = Snat.resp.bslmean(unit_idx);
imgnm_uniq = Snat.stim.imgname_uniq;
imgfps = Snat.stim.imgfps;
imgN = numel(rspvec);
[rsp_sort, sortidx] = sort(rspvec);
if sorted
    rspvec = rspvec(sortidx);
    semvec = semvec(sortidx);
    imgnm_uniq = imgnm_uniq(sortidx);
    imgfps = imgfps(sortidx);
end
[h1,ax1,ax2]=vis_selectivity_interact_core(imgnm_uniq, imgfps,...
                    rspvec,semvec,bsl,false);

if sorted
    if strcmp(order,"top")
    idx2show = imgN:-1:imgN-K+1;
    elseif strcmp(order,"bottom")
    idx2show = 1:K;
    end
else
    if strcmp(order,"top")
    idx2show = sortidx(imgN:-1:imgN-K+1);
    elseif strcmp(order,"bottom")
    idx2show = sortidx(1:K);
    end
end
axes(ax1);hold on
SCT = scatter(idx2show,rspvec(idx2show),'ro');
SCT.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('img',imgnm_uniq(idx2show));
SCT.DataTipTemplate.Interpreter = 'none';
% SCT.HitTest = 'off';
show_imgs(ax2, idx2show);
titlestr = compose("%s chan %s",Snat.meta.ephysFN,Snat.units.unit_name_arr(unit_idx));
sgtitle(h1,titlestr)
end

function [h1,ax1,ax2]=vis_selectivity_interact_core(pics_unique,imgfps,rspvec,semvec,bsl,react)
% Example: 
%    [h1,ax1,ax2]=vis_selectivity_interact_core(Snat.stim.imgname_uniq,Snat.stim.imgfps,...
%                     Snat.resp.meanMat_pref,Snat.resp.semMat_pref,...
%                     Snat.resp.bslmean(Snat.units.pref_chan_id));
% Parameter:
%    react: if false, the interactivity will be disabled. 
if nargin < 6, react=true; end
h1 = figure;set(h1,'pos',[188     89    1365       655])
figure(h1)
dcm = datacursormode(h1);
ax2 = subplot(122);set(ax2,'pos',[0.54    0.075    0.43    0.85]);hold on
set(ax2, 'UserData', struct("rspvec",rspvec,"semvec",semvec,...
    "imgnm_uniq",string(pics_unique),"imgfps",string(imgfps),...
    "showed_ids",[], "img_col",[]))
ax1 = subplot(121);set(ax1,'pos',[0.05    0.10    0.43    0.85])
barH = shadedErrorBar([],rspvec,semvec);
line([0,numel(rspvec)],[bsl, bsl],'linestyle','-.','linewidth',1.5)
xlim([0,numel(rspvec)+1])
xlabel("image id"); ylabel("Activation: rsp - bsl")
if react
set(barH.mainLine,'ButtonDownFcn',@(h,evt)selectpoint(h,evt,ax2),...
   'HitTest','on')
end
barH.patch.HitTest = 'off';
barH.mainLine.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('img',pics_unique);
barH.mainLine.DataTipTemplate.Interpreter = 'none';
dcm.UpdateFcn()
end

function [coordSelected, minIdx] = selectpoint(hObj, event, ax2)
X = hObj.XData; 
Y = hObj.YData; 
coordinates = [reshape(X,[],1),reshape(Y,[],1)];
pt = event.IntersectionPoint; 
x = pt(:,1);y = pt(:,2);
dist = pdist2([x,y],coordinates);            %distance between your selection and all points
[~, minIdx] = min(dist);            % index of minimum distance to points
coordSelected = coordinates(minIdx,:);
selected_id = int32(coordSelected(1));
% fprintf("Id %02d score %.1f, name %s\n",selected_id,coordSelected(2),ax2.UserData.imgnm_uniq(selected_id))
fprintf("Id %02d score %.1f+-%.1f, name %s\n",selected_id,...
    ax2.UserData.rspvec(selected_id),ax2.UserData.semvec(selected_id),ax2.UserData.imgnm_uniq(selected_id))
dt = datatip(hObj,x,y,'DeleteFcn',@(dtob,~)show_imgs(ax2, selected_id));
show_imgs(ax2, selected_id)
end

function show_imgs(ax2, imgid)
% disp(imgid)
for imgi = reshape(imgid,1,[])
    if any(ax2.UserData.showed_ids==imgi) 
        idx = find(ax2.UserData.showed_ids==imgi);
        ax2.UserData.showed_ids(idx) = [];
        ax2.UserData.img_col(idx) = [];
        continue
    else
        ax2.UserData.showed_ids(end+1) = imgi;
        img = imread(ax2.UserData.imgfps(imgi));
        ax2.UserData.img_col{end+1} = img;
    end
end
axes(ax2);cla(ax2)
montage(ax2.UserData.img_col)
% imshow(img)
% title(ax2.UserData.imgfps(imgi),'interp','none')
end



