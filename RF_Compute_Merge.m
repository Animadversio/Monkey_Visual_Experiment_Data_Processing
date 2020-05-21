Animal = "Alfa";

savepath = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(savepath, compose("%s_Manif_RFstats.mat", Animal)), 'RFStats');
%% Make masks from responses, in the new way and the interpolation way.
% Prepare canvas of mask
ntick = 201;
visualField = [-10 10]; 
coli = linspace(visualField(1),visualField(2),ntick);
rowi = linspace(visualField(1),visualField(2),ntick);
[gridX,gridY]  = meshgrid(coli,rowi);
MaskStats = repmat(struct(),length(RFStats),1);
for Mapi = 1:length(RFStats)
uniqpos = RFStats(Mapi).stim.uniqpos;
xpos = RFStats(Mapi).stim.xpos;
ypos = RFStats(Mapi).stim.ypos;
x_vect = uniqpos(:,1); %reshape(posgrid(:,:,1),1,[]);
y_vect = uniqpos(:,2); %reshape(posgrid(:,:,2),1,[]);
nPos = size(uniqpos,1);
posgrid = reshape([uniqpos, (1:nPos)'],[length(xpos),length(ypos),3]);
%% Set up the visual field plot
img_maskstack = zeros([size(gridX), nPos]);
img_radius = RFStats(Mapi).stim.imgsize_deg / 2;
for iPos = 1:nPos
xcent = uniqpos(iPos,1);
ycent = uniqpos(iPos,2);
Dinf = max(abs(gridX - xcent), abs(gridY - ycent)); % L_infty mask matrix to the center of the field
img_maskstack(:, :, iPos) = Dinf < img_radius;
end
% score_mat = zscore(reshape(RFStats(Mapi).psth.score_mean,[],nPos)',1); %
% this is obsolete since the baseline can be negative! 

% Subtract the mean baseline get the raw FR increase, and linearize spatial
% dimension into [spatial, channel] matrix
score_mat = reshape(RFStats(Mapi).psth.score_mean - mean(RFStats(Mapi).psth.baseline_mean,[2,3]), [], nPos)'; 
score_mat_norm = score_mat ./ max(score_mat, [], 1); % Normalize the maximal score to 1. 
% Use these normalized weights as score to weight the masks created from
% images. 
allmasks = einsum(img_maskstack, score_mat_norm, 'ijk,kl->ijl'); % Use 
allmasks_norm = allmasks ./ sum(img_maskstack,3); 
% Get the significant channels from T and F statistics 
signf_chans = RFStats(Mapi).stats.ttestP < 0.001 & RFStats(Mapi).stats.anovaP < 0.001;
% The old way of interpolating a mask.
interpmasks = single(zeros([ntick ntick length(RFStats(Mapi).meta.spikeID)]));
for iCh = 1:length(RFStats(Mapi).meta.spikeID)
resp = squeeze(RFStats(Mapi).psth.score_mean(iCh,:,:)) - mean(RFStats(Mapi).psth.baseline_mean(iCh,:,:),'all');
interp_rsp = griddata(x_vect, y_vect, resp(:), gridX,gridY);
interpmasks(:,:,iCh) = interp_rsp;
end
MaskStats(Mapi).signf_chans = signf_chans;
MaskStats(Mapi).interpmasks = single(interpmasks);
MaskStats(Mapi).convmasks = single(allmasks_norm);
MaskStats(Mapi).meta = RFStats(Mapi).meta;
MaskStats(Mapi).unit = RFStats(Mapi).unit;
end

savepath = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
save(fullfile(savepath, compose("%s_Manif_RFMaps.mat", Animal)), 'MaskStats');

%% Merge the Masks between multiple exps 
rfExpDates = arrayfun(@(S)S.meta.datestr, RFStats); % date for each RF experiment
uniqDates = unique(rfExpDates,'stable'); % uniq with stable to prevent matlab from shuffling it.
MergeMasks = repmat(struct(),length(uniqDates),1); % compute a summary 
for Dayj = 1:length(uniqDates)
cur_date_str = uniqDates(Dayj);
fprintf("Processing Exp Day %s, \n",cur_date_str)
MergeMasks(Dayj).datestr = cur_date_str; 
MergeMasks(Dayj).Mapidx = find(contains(rfExpDates, cur_date_str)); % Which rf Exps did I merged
nMerge = length(MergeMasks(Dayj).Mapidx);
MergeMasks(Dayj).meta = [];
% Collect the processed result and 
Results = [];
for j = 1:nMerge % in these rf Exps that I want to merge
    Mapi = MergeMasks(Dayj).Mapidx(j);
    Result = compute_RFMap(RFStats(Mapi), gridX, gridY);
    Results = [Results; Result];
    MergeMasks(Dayj).meta = [MergeMasks(Dayj).meta; RFStats(Mapi).meta];
    fprintf("Processed Exp %s, \n",RFStats(Mapi).meta.expControlFN)
end
try
    chan_num = arrayfun(@(M) length(M.spikeID), MergeMasks(Dayj).meta);
    assert(all(chan_num == chan_num(1)),compose("Channel number changed in the same exp day %s", cur_date_str))
catch
    fprintf("Channel number changed in the same exp day %s\n", cur_date_str);
    keyboard;
    continue
end
img_maskstack = cat(3,Results.img_maskstack);
score_mat = cat(1,Results.score_mat);
x_vect = cat(1,Results.x_vect);
y_vect = cat(1,Results.y_vect);
resp_mat = arrayfun(@(S) reshape(S.resp, size(S.resp, 1), []), Results, 'UniformOutput', false); % reshape response to linear 
resp_all = cat(2, resp_mat{:});
% convolution method
score_mat_norm = score_mat ./ max(score_mat, [], 1);
allmasks = einsum(img_maskstack, score_mat_norm, 'ijk,kl->ijl'); % Use 
allmasks_norm = allmasks ./ sum(img_maskstack,3); % normalize by the occupancy frequency of that entry in image field.
% interpolation method 
interpmasks = single(zeros([size(gridX), size(resp_all, 1)])); % length(S.meta.spikeID)
for iCh = 1:size(resp_all, 1)
resp = squeeze(resp_all(iCh, :, :));
interp_rsp = griddata(x_vect, y_vect, resp(:), gridX, gridY);
interpmasks(:,:,iCh) = interp_rsp;
end
MergeMasks(Dayj).convmasks = allmasks_norm;
MergeMasks(Dayj).interpmasks = interpmasks;
% MaskStats(Mapi).interpmasks;
end

%%  Visualize and inspect the results
figure;
scatter(x_vect,y_vect,25,resp(:))
hold on 
surf(gridX, gridY, interp_rsp)
figure;
imagesc(coli, rowi, interp_rsp)

%%
figure(3)
for iPos = 1:nPos
    imagesc(img_maskstack(:,:,iPos));
    axis image;set(gca,'YDir','normal')
    pause(0.1)
end

%% Show the masks and norma
figure(2); 
for iCh=1:size(allmasks,3)
if RFStats(Mapi).stats.anovaP(iCh) > 0.001
    fprintf("%d (%d) non modulated \n", RFStats(Mapi).meta.spikeID(iCh), iCh)
    continue
else
    fprintf("%d (%d) Significant \n", RFStats(Mapi).meta.spikeID(iCh), iCh)
end
subplot(131)
imagesc(coli,rowi,allmasks(:,:,iCh));
colorbar()
set(gca,'YDir','normal')
axis image
subplot(132)
imagesc(coli,rowi,allmasks_norm(:,:,iCh)); % this looks pretty good! 
title( compose("Ch %d, F %.1f T %.1f", iCh, RFStats(Mapi).stats.F(iCh), RFStats(Mapi).stats.T(iCh)) )
colorbar()
set(gca,'YDir','normal')
axis image
subplot(133)
resp = squeeze(RFStats(Mapi).psth.score_mean(iCh,:,:)) - mean(RFStats(Mapi).psth.baseline_mean(iCh,:,:),'all');
interp_rsp = griddata(x_vect, y_vect, resp(:), gridX,gridY);
imagesc(coli,rowi,interp_rsp);
axis image
colorbar()
set(gca,'YDir','normal')
pause(.2);
end

%% Loop through psth Animations of Channels (like the Manifold dynamics animation)
for iCh = 1:76
if RFStats(Mapi).stats.anovaP(iCh) > 0.001
    fprintf("%d (%d) non modulated \n", RFStats(Mapi).meta.spikeID(iCh), iCh)
    continue
else
    fprintf("%d (%d) Significant \n", RFStats(Mapi).meta.spikeID(iCh), iCh)
end
% iCh = 76;
psthdyn_mat = cell2mat(cellfun(@(psth) reshape(psth(iCh,:),1,1,[]), RFStats(Mapi).psth.psth_mean, 'UniformOutput', false));
%%
smth_dyn_mat = movmean(psthdyn_mat, 10, 3);
figure(5)
IM = imagesc(xpos,ypos,psthdyn_mat(:,:,1));
axis image
CMIN = prctile(psthdyn_mat,[2.5],'all');
CMAX = prctile(psthdyn_mat,[98],'all');
caxis([CMIN, CMAX]);
colorbar()
for fi = 1:size(smth_dyn_mat,3)
   IM.CData = smth_dyn_mat(:,:,fi);
   IM.Parent.Title.String = compose("Ch %d F %.1f T%.1f\n %d ms",...
       RFStats(Mapi).meta.spikeID(iCh), RFStats(Mapi).stats.F(iCh),RFStats(Mapi).stats.T(iCh),fi);
   pause(0.01)
   drawnow;
end
end

function Result = compute_RFMap(S, gridX, gridY)  % wrap up the RF computation in a function.
uniqpos = S.stim.uniqpos;
xpos = S.stim.xpos;
ypos = S.stim.ypos;
x_vect = uniqpos(:,1); %reshape(posgrid(:,:,1),1,[]);
y_vect = uniqpos(:,2); %reshape(posgrid(:,:,2),1,[]);
nPos = size(uniqpos,1);
posgrid = reshape([uniqpos, (1:nPos)'],[length(xpos),length(ypos),3]);
%% Set up the visual field plot
img_maskstack = single(zeros([size(gridX), nPos]));
img_radius = S.stim.imgsize_deg / 2;
for iPos = 1:nPos
xcent = uniqpos(iPos,1);
ycent = uniqpos(iPos,2);
Dinf = max(abs(gridX - xcent), abs(gridY - ycent)); % L_infty mask matrix to the center of the field
img_maskstack(:, :, iPos) = Dinf < img_radius;
end
% score_mat = zscore(reshape(S.psth.score_mean,[],nPos)',1); %
% this is obsolete since the baseline can be negative! 

% Subtract the **mean** baseline get the raw FR increase, and linearize spatial
% dimension into [spatial, channel] matrix. Use mean baseline instead of
% the paired one, to reduce noise from short window firing rate. 
score_mat = reshape(S.psth.score_mean - mean(S.psth.baseline_mean,[2,3]), [], nPos)'; 
score_mat_norm = score_mat ./ max(score_mat, [], 1); % Normalize the maximal score to 1. 
% Use these normalized weights as score to weight the masks created from
% images. 
allmasks = einsum(img_maskstack, score_mat_norm, 'ijk,kl->ijl'); % Use 
allmasks_norm = allmasks ./ sum(img_maskstack,3); 
% Get the significant channels from T and F statistics 
signf_chans = S.stats.ttestP < 0.001 & S.stats.anovaP < 0.001;
% The old way of interpolating a mask.
interpmasks = single(zeros([size(gridX), length(S.meta.spikeID)]));
resp_all = S.psth.score_mean(:,:,:) - mean(S.psth.baseline_mean(:,:,:),[2,3]);
for iCh = 1:length(S.meta.spikeID)
resp = squeeze(resp_all(iCh, :, :));
interp_rsp = griddata(x_vect, y_vect, resp(:), gridX, gridY);
interpmasks(:,:,iCh) = interp_rsp;
end
Result.x_vect = x_vect;
Result.y_vect = y_vect;
Result.resp = resp_all;
Result.signf_chans = signf_chans;
Result.interpmasks = interpmasks;
Result.img_maskstack = img_maskstack;
Result.convmasks = allmasks_norm;
Result.score_mat = score_mat;
end
