function SS = selectivity_Collect_Stats_fun(meta_new, rasters_new, Trials_new)
% Extract enough statistics from the files to do plotting. 
for iTr = 1:numel(meta_new)
meta = meta_new{iTr};
rasters = rasters_new{iTr};
Trials = Trials_new{iTr};
if contains(meta.ephysFN,["Alfa","ALfa"]), Animal = "Alfa";
elseif contains(meta.ephysFN,["Beto"]), Animal = "Beto";
end
% Meta info 
SS(iTr).Animal = Animal;
% BFEStats(Expi).Expi = Expi; 
SS(iTr).imageName = Trials.imageName;
SS(iTr).meta = meta;

%% units information
if isfield(meta,"unitID")
unit_num_arr = meta.unitID;
unit_name_arr = generate_unit_labels_new(meta.spikeID, meta.unitID);
activ_msk = unit_num_arr~=0;
else
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
end
prefchan = Trials.TrialRecord.User.prefChan;
if isfield(Trials.TrialRecord.User,"prefUnit")
    prefUi = Trials.TrialRecord.User.prefUnit;
    % prefUi = ExpRecord.pref_unit(exp_rowi);
else
    frpintf("pref unit information not found, use unit 1\n")
    prefUi = 1;
end	
prefchan_ids = find(meta.spikeID == prefchan & meta.unitID>0);
pref_chan_id = prefchan_ids(prefUi); 
unit_in_pref_chan = prefUi;
% if strcmp(meta.ephysFN, 'Alfa-22042021-007'), prefchan_id = prefchan_ids(2);end
% if strcmp(meta.ephysFN, 'Alfa-27042021-007'), prefchan_id = prefchan_ids(1);end

SS(iTr).units.unit_name_arr = unit_name_arr;
SS(iTr).units.unit_num_arr = unit_num_arr;
SS(iTr).units.activ_msk = activ_msk;
SS(iTr).units.spikeID = meta.spikeID;
SS(iTr).units.pref_chan = prefchan;
SS(iTr).units.pref_chan_id = pref_chan_id;
SS(iTr).units.unit_in_pref_chan = unit_in_pref_chan;

%% sort images 
fprintf([meta.comments,'\n'])
impos = Trials.TrialRecord.User.uniquePositions;
imsize_pix = Trials.width(1,1);
% imsize_deg = ExpRecord.stim_size(exp_rowi);
imsize_deg = round(imsize_pix/40 * 2) / 2;% usually imsize_deg is set to be multiple to 0.5. 
imgname_uniq = unique(Trials.imageName); 
[imgfps, mapper] = map2fullpath(string(imgname_uniq), meta.stimuli);
img_idxarr = cellfun(@(imgnm)find(strcmp(Trials.imageName,imgnm)),imgname_uniq,'uni',0); % trial id each image 

SS(iTr).stim.imgfn_mapper = mapper;
SS(iTr).stim.imgfps = imgfps;
SS(iTr).stim.impos = impos;
SS(iTr).stim.imsize_pix = imsize_pix;
SS(iTr).stim.imsize_deg = imsize_deg;
SS(iTr).stim.imgname_uniq = imgname_uniq;
SS(iTr).stim.img_idxarr = img_idxarr;

%% format all experiments 

%% sort all responses into images 
evkwdw = [51:200];  bslwdw = [1:40];
bslmean = squeeze(mean(rasters(:, bslwdw, :),[2,3])); %average baseline activity 
bslsem = sem(squeeze(mean(rasters(:, bslwdw, :),[2])),3); % variability between trials in baseline. 
evoke_trials  = squeeze(mean(rasters(:, evkwdw, :),[2])); % evoked firing rate
evkbsl_trials = evoke_trials - bslmean; % evk - bsl score using overall baseline 

resp_trial_col = cellfun(@(idx)squeeze(evoke_trials(:,idx)),img_idxarr,'Un',0); 

% cell array of trial scores each img, shape imageN-by-1, each cell has chanN-by-trialN mat.
% compute mean, std, sem over the trial dimension. 
resp_mean_col = cellfun(@(resp)mean(resp,2),resp_trial_col,'Unif',0); 
resp_std_col = cellfun(@(resp)std(resp,0,2),resp_trial_col,'Unif',0); 
resp_sem_col = cellfun(@(resp)sem(resp,2),resp_trial_col,'Unif',0); 
%cellfun(@(idx)squeeze(mean(evoke_trials(:,idx),2)),img_idxarr,'Un',0);
resp_meanMat = cat(2,resp_mean_col{:}); % chan N-by-image N 
resp_stdMat = cat(2,resp_std_col{:}); % chan N-by-image N 
resp_semMat = cat(2,resp_sem_col{:}); % chan N-by-image N 

SS(iTr).resp.trial_col = resp_trial_col;
SS(iTr).resp.meanMat = resp_meanMat;
SS(iTr).resp.stdMat = resp_stdMat;
SS(iTr).resp.semMat = resp_semMat;
SS(iTr).resp.bslmean = bslmean;
SS(iTr).resp.bslsem = bslsem;
SS(iTr).resp.meanvec_pref = resp_meanMat(pref_chan_id, :);
SS(iTr).resp.stdvec_pref = resp_stdMat(pref_chan_id, :);
SS(iTr).resp.semvec_pref = resp_semMat(pref_chan_id, :);
%% psth of pref chan , TODO



%% format into a tensor , TODO

end
end 