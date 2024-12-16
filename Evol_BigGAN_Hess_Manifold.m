
!ExpRecordBackup.bat
%%
% saveroot = "E:\OneDrive - Washington University in St. Louis\Evol_Cosine";
saveroot = "E:\OneDrive - Harvard University\BigGAN_Hessian";
Set_Path;
%%
Animal = "Caos";
ExpRecord = loadBackupExpRecord(Animal);
currows = find(contains(ExpRecord.expControlFN,["241209","241210"]));  % ,"241202","241204"
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(currows(1:end), Animal, false);
bhvfns = ExpRecord.expControlFN(currows);
%%
Animal = "Diablito"; Set_Path;
currows = find(contains(ExpRecord.expControlFN,["241211", "241212"])); % "241203"
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(currows(1:end), Animal, false);
bhvfns = ExpRecord.expControlFN(currows);
%%
Animal = "Caos";
ExpRecord = loadBackupExpRecord(Animal);
currows = find(contains(ExpRecord.expControlFN,["241202","241204","241209","241210"]));  % ,
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(currows(1:end), Animal, false);
bhvfns = ExpRecord.expControlFN(currows);
%%
Animal = "Diablito"; 
ExpRecord = loadBackupExpRecord(Animal);
currows = find(contains(ExpRecord.expControlFN,["241203", "241211", "241212"]) & ...
               contains(ExpRecord.Exp_collection, "BigGAN_Hessian") ); % 
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(currows, Animal, false);
bhvfns = ExpRecord.expControlFN(currows);
%%  Final loading script for Hessian experiments
ExpRecord = [loadBackupExpRecord("Caos");...
             loadBackupExpRecord("Diablito")];
currows = find(contains(ExpRecord.expControlFN,["241202","241204","241209","241210","241203", "241211", "241212"]) & ...
               contains(ExpRecord.Exp_collection, "BigGAN_Hessian") ); % 
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments_tabrows(ExpRecord(currows, :), false);
bhvfns = ExpRecord.expControlFN(currows);

%%
hess_mask = contains(ExpRecord.Exp_collection(currows), "BigGAN_Hessian");
disp(ExpRecord(currows(hess_mask),:))
HEStats = Hess_BigGAN_collect_Stats_fun(meta_new(hess_mask), rasters_new(hess_mask), Trials_new(hess_mask));

%%  Final loading script for BigGAN-FC6 experiments
ExpRecord = [loadBackupExpRecord("Caos");...
             loadBackupExpRecord("Diablito")];
currows = find(contains(ExpRecord.expControlFN,["241202","241204","241209","241210","241203", "241211", "241212"]) & ...
               contains(ExpRecord.Exp_collection, "BigGAN_FC6") & ...
                ~isnan(ExpRecord.Expi)); % 
[meta_exp, rasters_exp, ~, Trials_exp] = loadExperiments_tabrows(ExpRecord(currows, :), false);

%% batch collect stats
BFEStats = Evol_BigGAN_FC6_Collect_Stats_fun(meta_exp(1:end), rasters_exp(1:end), Trials_exp(1:end));
%% batch animation
Evol_BigGAN_FC6_Animation_fun(BFEStats)
%% save stats for all experiments
matdir = "E:\OneDrive - Harvard University\Mat_Statistics";
savefast(fullfile(matdir, "CD_BigGAN_Hessian_Manifold_Stats.mat"), "HEStats");
savefast(fullfile(matdir, "CD_BigGAN_FC6_Evolution_Stats.mat"), "BFEStats");
%%
hess_exp_record = ExpRecord(currows(hess_mask), :);
disp(hess_exp_record)



%%
saveroot = "E:\OneDrive - Harvard University\BigGAN_Hessian";
for Expi = 1:size(hess_exp_record,1)
% pref_chan = 71;
ephys_name = HEStats(Expi).meta.ephysFN;
pref_chan = hess_exp_record.pref_chan(Expi);
fprintf("ephys %s Exp %d Pref Chan: %d\n", ephys_name, Expi, pref_chan);
figdir = fullfile(saveroot, ephys_name);
mkdir(figdir);
%%
stimroot = HEStats(Expi).meta.stimuli;
% Load noise images into cell array
noise_imgs = arrayfun(@(imgnm)imread(fullfile(stimroot, imgnm+".jpg")), ...
            HEStats(Expi).noise.imgnm_arr, 'UniformOutput', false);
class_imgs = arrayfun(@(imgnm)imread(fullfile(stimroot, imgnm+".jpg")), ...
            HEStats(Expi).class.imgnm_arr, 'UniformOutput', false);
figure(2);clf;set(2,'pos',[ 81         100       1720         750])
subplot(1,2,1);
if numel(class_imgs) > 0
    montage(class_imgs','Size',size(HEStats(Expi).class.imgnm_arr),...
                'ThumbnailSize',[],'BorderSize',1); % Montage fills columns first, then rows
    title('Class Space Images');
end
subplot(1,2,2); 
if numel(noise_imgs) > 0
    montage(noise_imgs','Size',size(HEStats(Expi).noise.imgnm_arr),...
                'ThumbnailSize',[],'BorderSize',1); % Montage fills columns first, then rows
    title('Noise Space Images');
end
sgtitle(compose("Image Montages for %s (Pref chan %d)", ephys_name, pref_chan));
saveallform(figdir, compose("hess_class_noise_images"),2,["png","pdf"]);
%%
pref_chan_id = find(HEStats(Expi).units.spikeID == pref_chan);
for chan_id = pref_chan_id'
    unit_str = HEStats(Expi).units.unit_name_arr{chan_id};
    if isempty(HEStats(Expi).class.resp_mat)
        fprintf("No class response matrix for channel %s (idx: %d)\n", unit_str, chan_id);
        class_resp_mat = [];
        class_eig_id_arr = [];
        class_dist_arr = [];
    else
        class_resp_mat = squeeze(HEStats(Expi).class.resp_mat(chan_id,:,:));
        class_eig_id_arr = HEStats(Expi).class.eig_id_arr;
        class_dist_arr = HEStats(Expi).class.dist_arr;
    end
    if isempty(HEStats(Expi).noise.resp_mat)
        fprintf("No noise response matrix for channel %s (idx: %d)\n", unit_str, chan_id);
        noise_resp_mat = [];
        noise_eig_id_arr = [];
        noise_dist_arr = [];
    else
        noise_resp_mat = squeeze(HEStats(Expi).noise.resp_mat(chan_id,:,:));
        noise_eig_id_arr = HEStats(Expi).noise.eig_id_arr;
        noise_dist_arr = HEStats(Expi).noise.dist_arr;
    end
    % plot the heatmap of the response matrix
    figure(1);clf;
    set(1, 'Position', [100, 100, 1300, 600]);
    % Calculate shared color limits from both matrices
    CLIM = prctile([class_resp_mat(:); noise_resp_mat(:)],[2,98],'all')';
    subplot(1,2,1);
    imagesc(class_resp_mat);
    clim(CLIM);
    axis image
    colorbar;
    yticks(1:length(class_eig_id_arr));
    yticklabels(HEStats(Expi).class.eig_id_arr);
    ylabel('Eigenvector ID');
    xticks(1:length(class_dist_arr));
    xticklabels(class_dist_arr);
    xlabel('Distance');
    title(sprintf("Class Response Matrix for Channel %s (idx: %d)", unit_str, chan_id));
    subplot(1,2,2);
    imagesc(noise_resp_mat);
    clim(CLIM);
    axis image
    colorbar;
    yticks(1:length(noise_eig_id_arr));
    yticklabels(noise_eig_id_arr);
    ylabel('Eigenvector ID');
    xticks(1:length(noise_dist_arr));
    xticklabels(noise_dist_arr);
    xlabel('Distance');
    title(sprintf("Noise Response Matrix for Channel %s (idx: %d)", unit_str, chan_id));
    sgtitle(sprintf("Ephys: %s, Pref Chan: %d Exp: %d", ephys_name, pref_chan, Expi));
    saveallform(figdir, compose("hess_class_noise_resp_mat_pref_chan%s", unit_str),1,["png","pdf"]);


    % CLIM = prctile([class_resp_mat(:); noise_resp_mat(:)],[2,98],'all')';
    frame_noise_imgs_col = score_frame_image_arr(noise_imgs, noise_resp_mat, CLIM);
    frame_class_imgs_col = score_frame_image_arr(class_imgs, class_resp_mat, CLIM);
    figure(3);clf;set(3,'pos',[ 81         100       1720         750])
    subplot(1,2,1);
    if numel(class_imgs) > 0
        imagesc(class_resp_mat);
        colorbar;
        montage(frame_class_imgs_col', 'Size', size(class_imgs), 'ThumbnailSize', [], 'BorderSize', 1);
        title('Class Response');
    end
    subplot(1,2,2);
    if numel(noise_imgs) > 0
        imagesc(noise_resp_mat);
        colorbar;
        montage(frame_noise_imgs_col', 'Size', size(noise_imgs), 'ThumbnailSize', [], 'BorderSize', 1);
        title('Noise Response');
    end
    sgtitle(compose("Ephys: %s, Pref Chan: %d Exp: %d", ephys_name, pref_chan, Expi));
    saveallform(figdir, compose("hess_class_noise_resp_mat_pref_chan%s_coloredframe", unit_str),3,["png","pdf"]);
end
end
%%