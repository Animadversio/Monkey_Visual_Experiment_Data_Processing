%% Evolution Select Axis
system("subst S: E:\Network_Data_Sync")
system("subst N: E:\Network_Data_Sync")
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir,"Alfa_ManifPopDynamics.mat"));
load(fullfile(mat_dir,"Alfa_Manif_stats.mat"));
load(fullfile(mat_dir,"Alfa_Evol_stats.mat"));
%% 
Expi = 3;geni = 1;
global evol_psth
evol_psth = EStats(Expi).evol.psth;
psth_all = squeeze(cell2mat(reshape(evol_psth,1,1,[])))'; % big array
%% Substitute the stimuli path
assert(exist(strrep(EStats(Expi).meta.stimuli,"2019-06-Evolutions","2019-Manifold"),'dir')>0)
EStats(Expi).meta.stimuli = strrep(EStats(Expi).meta.stimuli,"2019-06-Evolutions","2019-Manifold");
%% Collect the generated codes in the same format as the exp 
global code_col idx_col
codematfns = sort(string(ls(EStats(Expi).meta.stimuli+"\*.mat"))); % these are codes generated in matlab
code_col = cell(length(codematfns),1);
idx_col = cell(length(codematfns),1); 
for geni = 1:length(codematfns)
tmpdata = load(fullfile(EStats(Expi).meta.stimuli, codematfns(geni)));
code_col{geni} = double(tmpdata.codes);
idx_col{geni} = string(tmpdata.ids)';
end
codes_all = cat(1,code_col{:});
imgidx_all = cat(1,idx_col{:}); % These are code index in generating images
%%
% psth_sort_gen = cell(length(codematfns),1); 
% for geni = 1:length(codematfns)
% psth_sort_gen{geni} = nan(length(idx_col{geni}), 200);
% for ii = 1:length(idx_col{geni})
% idx = find(contains(imgname_all, idx_col{geni}(ii)));
% psth_sort_gen{geni}(ii, :) = psth_all(idx,:);
% 
% end
% end
%% imagename for all generated images
imgname_gen = cellfun(@(idx)string(EStats(Expi).imageName(idx)),EStats(Expi).evol.idx_seq,'UniformOutput',false); % name of image files 
imgname_all = cat(1, imgname_gen{:});
%% Find the code for each psth.
global code_sort_gen
imgname_blocksp = split(imgname_all,"_gen_");
code_sort_gen = cell(length(imgname_gen), 1);
% code_gen = zeros(size(psth_all,1), 1);
for geni = 1:length(imgname_gen)
code_sort_gen{geni} = zeros(length(imgname_gen{geni}), 4096);
for ii = 1:length(imgname_gen{geni})
tmpsp = split(imgname_gen{geni}(ii),"_gen_");
idx = find(contains(imgidx_all, tmpsp(2)));
code_sort_gen{geni}(ii, :) = codes_all(idx, :);
end
end
code_sort = cat(1, code_sort_gen{:});
% imgname_blocksp = split(imgname_all,"_gen_");
% code_sort = zeros(size(psth_all,1), 4096);
% code_gen = zeros(size(psth_all,1), 1);
% for ii = 1:size(psth_all,1)
% idx = find(contains(imgidx_all, imgname_blocksp(ii,2)));
% code_sort(ii, :) = codes_all(idx, :);
% geni = regexp(imgname_blocksp(ii,1),"block(\d.*)_thread",'tokens'); % block n from the image name not necessarily match block in trial
% code_gen(ii) = str2num(geni{1});
% end

%% Dimension Reduce the Code and visualize it with the image and PSTH. 
global code_red code_umap
[code_red,code_umap] = run_umap(code_sort, 'metric', "correlation");
%% 
[psth_red,psth_umap] = run_umap(psth_all, 'metric', "euclidean"); % correlation is not a good metric for code similarity here
%%
[reduction,umap2] = run_umap(codes_all, 'metric', "euclidean");
%%
global G 
G = FC6Generator("matlabGANfc6.mat");
%%
evofig = figure(4);clf;evofig.Position = [28         524        1446         454];
evoldata = struct('mask', mask, 'addMask', 0, "curcode", basis(1,:), "scoremap", scoremap, ...
    'psth',psth_all(1,:),'meanPSTH', 1, "rank", 1);
% You want multi-views of your data. 
% Neural view (PSTH), Image View, Code Space View
subplot("Position",[ 0.0913    0.1564    0.2579    0.7313]);hold on 
evoldata.DR = scatter(code_red(:,1), code_red(:,2), 50, 1:size(code_red,1));
evoldata.csr = plot(code_red(1,1),code_red(2,2),"Color",'r','Marker','o');
colorbar()
subplot("Position",[ 0.4108    0.1091    0.2134    0.7530])
evoldata.imax = imshow(G.visualize(sphere_norm*basis(1,:)));
subplot("Position",[0.6694    0.0969    0.3036    0.8304])
evoldata.psthplot = plot(squeeze(evol_psth{1}(1,:,1))); 
ylim([min(psth_all,[],"all"),max(psth_all,[],"all")]) % not perfect for showing avg response
% Slider controlling the generation
evoldata.genSLD = uicontrol(evofig, 'Style', "slider", 'Position', [600 30 250 20], "String", "generation Slider", ...
            'Value', 1, "Min", 1, "Max", length(evol_psth), "SliderStep", [1, 2]./(length(evol_psth)-1),...
            "Callback", @(src, evt) genSlider_Callback(src, evt));
evoldata.rankSLD = uicontrol(evofig, 'Style', "slider", 'Position', [600 5 250 20], "String", "Rank Slider", ...
            'Value', 1, "Min", 1, "Max", 40, "SliderStep", [1, 2]./(40-1),...
            "Callback", @(src, evt) rankSlider_Callback(src, evt));
% TheSLD.Min = -pi/2; TheSLD.Max = pi/2;TheSLD.SliderStep = [0.005 0.05];
% Button to play Evolution Animation
evoldata.ST = suptitle(compose("%s Evol (Manif) Exp %02d pref chan %s", ...
    Animal, Expi, EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id)));
evoldata.playEvol = uicontrol(evofig, 'Style', 'pushbutton', 'Position', [20 160 80 30], "String", "Play Evolution",...
    "Callback", @(src, evt) playEvol_Callback(src, evt));
% Button to toggle mask
evoldata.toggleMask = uicontrol(evofig, 'Style', 'togglebutton', 'Position', [20 10 80 30], "String", "RF Mask",...
    "Callback", @(src, evt) Mask_Callback(src, evt));
evoldata.meanToggle = uicontrol(evofig, 'Style', 'checkbox', 'Position', [20 40 80 30], "String", "Block Avg", ...
    "Callback", @(src, evt) meanToggle_Callback(src, evt),'Value', 1);
% LocalPCASLD

guidata(evofig, evoldata);
%%
function playEvol_Callback(hObj, evt)
global evol_psth
data = guidata(hObj);
for geni=data.genSLD.Min:data.genSLD.Max
data.genSLD.Value = geni;
genSlider_Callback(data.genSLD, evt)
pause(0.025);
end
end

function rankSlider_Callback(hObj, evt)
global evol_psth code_sort_gen
data = guidata(hObj);
data.rank = clip(round(hObj.Value),[1,hObj.Max]);

if data.meanPSTH
    data.psth = mean(evol_psth{data.geni}, 3);
    data.curcode = mean(code_sort_gen{data.geni}, 1);
else % individual trial
%     [~,maxi] = max(squeeze(sum(evol_psth{data.geni}(:,50:end,:),2))); 
    [~,maxi] = maxk(squeeze(sum(evol_psth{data.geni}(:,50:end,:),2)),data.rank); % note maxk can handle larger then length rank! Graceful! 
    maxi = maxi(end); % k th largest value
    data.psth = evol_psth{data.geni}(1,:,maxi);
    data.curcode = code_sort_gen{data.geni}(maxi,:); % Note this is wrong since the 2 array are not aligned...
end
guidata(hObj,data)
draw(data.imax)
drawPSTH(data.psthplot)
drawlatent(data.DR)
end

function genSlider_Callback(hObj, evt)
global evol_psth code_sort_gen
data = guidata(hObj);
data.geni = clip(round(hObj.Value),[1,length(evol_psth)]);
data.rankSLD.Max = size(evol_psth{data.geni}, 3); data.rankSLD.SliderStep = [1, 2]/(data.rankSLD.Max - 1);
if data.meanPSTH
    data.psth = mean(evol_psth{data.geni}, 3);
    data.curcode = mean(code_sort_gen{data.geni}, 1);
else % individual trial
%     [~,maxi] = max(squeeze(sum(evol_psth{data.geni}(:,50:end,:),2))); 
    [~,maxi] = maxk(squeeze(sum(evol_psth{data.geni}(:,50:end,:),2)),data.rank); % note maxk can handle larger then length rank! Graceful! 
    maxi = maxi(end); % k th largest value
    data.psth = evol_psth{data.geni}(1,:,maxi);
    data.curcode = code_sort_gen{data.geni}(maxi,:); % Note this is wrong since the 2 array are not aligned...
end
guidata(hObj,data)
draw(data.imax)
drawPSTH(data.psthplot)
drawlatent(data.DR)
end

function Mask_Callback(hObj, evt)
data = guidata(hObj);
data.addMask = hObj.Value;
guidata(hObj,data);
end

function meanToggle_Callback(hObj, evt)
global evol_psth
data = guidata(hObj);
data.meanPSTH = hObj.Value;
if hObj.Value
data.psthplot.Parent.YLim = [0  max(cellfun(@(psth)max(mean(psth,3)),evol_psth'))];
else
data.psthplot.Parent.YLim = [0 max(cellfun(@(psth)max(psth,[],'all'),evol_psth'))];
end
guidata(hObj,data);
drawPSTH(data.psthplot)
drawlatent(data.DR)
end

function draw(imax)
global G 
data = guidata(imax);
if ~ data.addMask
    imax.CData = G.visualize(data.curcode);
else
    imax.CData = uint8( data.mask .* double(G.visualize(data.curcode)));
end
if data.meanPSTH
    imax.Parent.Title.String = compose("Mean code gen%d", data.geni);
else
    imax.Parent.Title.String = compose("Sample code gen%d, rank %d", data.geni, data.rank);
end
% imax.Parent.Title.String = compose("%.1f, (%.1f, %.1f)",norm(data.curcode), data.Theta, data.Phi);
drawnow;
guidata(imax,data);
end
function drawlatent(DR)
global G code_umap
data = guidata(DR);
data.cur_embed = code_umap.transform([data.curcode;data.curcode]);
data.csr.XData = data.cur_embed(1,1); data.csr.YData = data.cur_embed(1,2);
% imax.Parent.Title.String = compose("%.1f, (%.1f, %.1f)",norm(data.curcode), data.Theta, data.Phi);
drawnow;
guidata(DR,data);
end
function drawPSTH(psthplot)
% global evol_psth
data = guidata(psthplot);
% psth = sum(psth_avg_tsr(:,:,i_grid,j_grid) .* reshape(i_W'*j_W,[1,1,2,2]),[1,3,4]);
% err = sum(psth_std_tsr(:,:,i_grid,j_grid) .* reshape(i_W'*j_W,[1,1,2,2]),[1,3,4]);
psthplot.YData = data.psth;
% psthplot.XData = [1:200,nan,1:200,nan,1:200];
drawnow;
end