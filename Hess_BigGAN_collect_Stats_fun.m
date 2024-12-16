function HessBGStats = Hess_BigGAN_collect_Stats_fun(meta_new, rasters_new, Trials_new)

% %% Hessian BigGAN Collect Statistics
% 
% %%
% Animal="Beto";Set_Path;
% ftridx = find(contains(ExpRecord.Exp_collection,"BigGAN_Hessian") & contains(ExpRecord.expControlFN,"select") );
% [meta_new,rasters_new,~,Trials_new] = loadExperiments(ftridx, Animal);
% %%
MatStatsDir = "E:\OneDrive - Harvard University\Mat_Statistics";
figroot = "E:\OneDrive - Harvard University\HessBigGANTune";
HessBGStats = repmat(struct(),1,numel(meta_new));
% %%
% D = torchImDist();
D = [];
%%
for Triali = 1:numel(meta_new)
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};

if contains(meta.ephysFN,"Alfa"),Animal = "Alfa";
elseif contains(meta.ephysFN,"Beto"),Animal = "Beto";
elseif contains(meta.ephysFN,"Caos"),Animal = "Caos";
elseif contains(meta.ephysFN,"Diablito"),Animal = "Diablito";
end
if false
    exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
    Expi = ExpRecord.Expi(exp_rowi);
    fprintf("Processing  Exp %d:\n",Expi)
    fprintf([ExpRecord.comments{exp_rowi},'\n'])
    figdir = fullfile("E:\OneDrive - Harvard University\HessBigGANTune",compose("%s_Exp%02d",Animal,Expi));
    if meta.ephysFN == "Alfa-10082020-004", figdir = figdir + "_post"; end
end
figdir = fullfile(figroot, meta.ephysFN);%compose("%s_Exp%02d",Animal,Expi));
% make sure figdir exists
if ~exist(figdir, 'dir'), mkdir(figdir); end
Expi = nan;

prefchan = Trials.TrialRecord.User.prefChan;
if numel(prefchan) > 1
    % check all prefchan are the same
    if ~all(prefchan == prefchan(1))
        disp(prefchan)
        error("Multiple prefchan found, please check the data %s", meta.ephysFN);
    end
    prefchan = prefchan(1);
end
prefchan_ids = find(meta.spikeID == prefchan);
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% Record info in Stats 
HessBGStats(Triali).Animal = Animal;
% HessBGStats(Triali).Expi = Expi;
HessBGStats(Triali).imageName = Trials.imageName;
HessBGStats(Triali).meta = meta;
HessBGStats(Triali).figdir=figdir;
HessBGStats(Triali).units.pref_chan = prefchan;
HessBGStats(Triali).units.unit_name_arr = unit_name_arr;
HessBGStats(Triali).units.unit_num_arr = unit_num_arr;
HessBGStats(Triali).units.activ_msk = activ_msk;
HessBGStats(Triali).units.spikeID = meta.spikeID;
HessBGStats(Triali).units.pref_chan_id = prefchan_ids;
%%
imgname_uniq = unique( Trials.imageName );
[noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq);
% Extract parameters of images in noise space from img name
namepart_uniq = regexp(imgname_uniq, noise_pattern, 'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
eig_id_arr_nos = unique(cellfun(@(U)str2num(U.eig_id),namepart_uniq));
dist_arr_nos = unique(cellfun(@(U)str2num(U.dist),namepart_uniq));
% collect image name and trial indices in cell
imgnm_arr_nos = strings(length(eig_id_arr_nos), length(dist_arr_nos));
idx_arr_nos = cell(length(eig_id_arr_nos), length(dist_arr_nos));
for i = 1:length(eig_id_arr_nos)
    eig_id = eig_id_arr_nos(i);
    for j = 1:length(dist_arr_nos)
    dist = dist_arr_nos(j);
    imgnm = compose(noise_imgnm,eig_id,dist); % e.g. "noise_eig%d_lin%.1f"
    idx_arr_nos{i,j} = find(contains(Trials.imageName, imgnm));
    imgnm_arr_nos(i,j) = imgnm;
    if numel(idx_arr_nos{i,j}) == 0 && dist== 0 
    % Recently, we deleted images like "noise_eig2_lin0.0" so should use "noise_eig1_lin0.0" instead
    imgnm_arr_nos(i,j) = compose(noise_imgnm,eig_id_arr_nos(1),0);
    end
    end
end
img_noise = arrayfun(@(imgnm)imread(fullfile(meta.stimuli, imgnm+".jpg")),imgnm_arr_nos,"Un",false);
% Extract parameters of images in class space from img name
namepart_uniq = regexp(imgname_uniq,class_pattern,'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
eig_id_arr_cls = unique(cellfun(@(U)str2num(U.eig_id),namepart_uniq));
dist_arr_cls = unique(cellfun(@(U)str2num(U.dist),namepart_uniq));
% collect image name and trial indices in cell
imgnm_arr_cls = strings(length(eig_id_arr_cls), length(dist_arr_cls));
idx_arr_cls = cell(length(eig_id_arr_cls), length(dist_arr_cls));
for i = 1:length(eig_id_arr_cls)
    eig_id = eig_id_arr_cls(i);
    for j = 1:length(dist_arr_cls)
    dist = dist_arr_cls(j);
    imgnm = compose(class_imgnm,eig_id,dist); % "class_eig%d_lin%.1f"
    idx_arr_cls{i,j} = find(contains(Trials.imageName, imgnm));
    imgnm_arr_cls(i,j) = imgnm;
    if numel(idx_arr_cls{i,j}) == 0 && dist== 0 
    % Recently, we deleted images like "noise_eig2_lin0.0" so should use "noise_eig1_lin0.0" instead
    imgnm_arr_cls(i,j) = compose(class_imgnm,eig_id_arr_cls(1),0);
    end
    end
end
img_class = arrayfun(@(imgnm)imread(fullfile(meta.stimuli, imgnm+".jpg")),imgnm_arr_cls,"Un",false);

nrow_cls = length(eig_id_arr_cls);
ncol_cls = length(dist_arr_cls);
nrow_nos = length(eig_id_arr_nos);
ncol_nos = length(dist_arr_nos);
% image name info recorded here
HessBGStats(Triali).noise.RE = noise_pattern; % An RE pattern that could be used to search for noise space images
HessBGStats(Triali).noise.imgnm_tmpl = noise_imgnm;
HessBGStats(Triali).noise.eig_id_arr = eig_id_arr_nos;
HessBGStats(Triali).noise.dist_arr = dist_arr_nos;
HessBGStats(Triali).noise.imgnm_arr = imgnm_arr_nos;
HessBGStats(Triali).noise.idx_arr = idx_arr_nos;

HessBGStats(Triali).class.RE = class_pattern; % An RE pattern that could be used to search for class space images
HessBGStats(Triali).class.imgnm_tmpl = class_imgnm;
HessBGStats(Triali).class.eig_id_arr = eig_id_arr_cls; % eig_id and dist are annotation for the matrix 
HessBGStats(Triali).class.dist_arr = dist_arr_cls;
HessBGStats(Triali).class.imgnm_arr = imgnm_arr_cls; % imgnm_arr 
HessBGStats(Triali).class.idx_arr = idx_arr_cls; % idx_arr are the most useful. 
fprintf("Image name infor finished\n")
%% Reference images 
noise_parts = regexp(imgname_uniq, noise_pattern, 'names');
class_parts = regexp(imgname_uniq,class_pattern,'names');
ref_nmidx = find(cellfun(@isempty,noise_parts) .* cellfun(@isempty,class_parts));
ref_imgnms = imgname_uniq(ref_nmidx);
HessBGStats(Triali).ref.imgnm_arr = ref_imgnms;
HessBGStats(Triali).ref.idx_arr = arrayfun(@(nm) find(contains(Trials.imageName, nm)), ref_imgnms, 'Uni',0);
fprintf("Reference Image finished\n")

%% Compute image distances per row or for all
if isempty(D)
    fprintf("No distance matrix found, skip computing image distances\n")
else
    imdist_class = zeros(size(img_class));
    cent_i = (ncol_cls + 1) / 2;
    for rowi = 1:size(img_class,1)
    imrow = cat(4,img_class{rowi,:});
    distrow = D.distance(imrow, imrow(:,:,:,cent_i));
    imdist_class(rowi, :) = distrow;
    end
    imall = cat(4,img_class{:});
    distmat_class = D.distmat_B(imall);
    figure(6);imagesc(distmat_class);axis image;colorbar()
    title(compose("%s Exp %d class\nLPIPS squeezenet", Animal, Expi))
    saveas(6,fullfile(figdir, "class_img_dissim.jpg")); 

    imdist_noise = zeros(size(img_noise));
    cent_i = (ncol_nos + 1) / 2;
    for rowi = 1:size(img_noise,1)
    imrow = cat(4,img_noise{rowi,:});
    distrow = D.distance(imrow, imrow(:,:,:,cent_i));
    imdist_noise(rowi, :) = distrow;
    end
    imall = cat(4,img_noise{:,:});
    distmat_noise = D.distmat_B(imall);
    figure(9);imagesc(distmat_noise);axis image;colorbar()
    title(compose("%s Exp %d noise\nLPIPS squeezenet", Animal, Expi))
    saveas(9,fullfile(figdir, "noise_img_dissim.jpg")); 

    HessBGStats(Triali).class.imdist_row = imdist_class;
    HessBGStats(Triali).class.distmat = distmat_class;
    HessBGStats(Triali).noise.imdist_row = imdist_noise;
    HessBGStats(Triali).noise.distmat = distmat_noise;
    fprintf("Image Dissimilarity finished\n")
end
%% Collect response and compute statistics 
for iCh = prefchan_ids %1:size(rasters,1)
    % collect response and PSTH in cell or array
    resp_mat_nois = cellfun(@(idx)mean(rasters(iCh,51:200,idx),'all'), idx_arr_nos); % mean scaler RSP
    resp_mat_cls = cellfun(@(idx)mean(rasters(iCh,51:200,idx),'all'), idx_arr_cls);
    psth_mean_nois = cellfun(@(idx)mean(rasters(iCh,:,idx),3), idx_arr_nos, 'uni',0); % mean PSTH
    psth_mean_cls = cellfun(@(idx)mean(rasters(iCh,:,idx),3), idx_arr_cls, 'uni',0);
    psth_col_nois = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_nos, 'uni',0); % All trials
    psth_col_clas = cellfun(@(idx)rasters(iCh,:,idx), idx_arr_cls, 'uni',0);
end
    % Compute Statistics for psth. 
    stats_nois = calc_tune_stats(psth_col_nois);
    stats_row_nois = arrayfun(@(i)calc_tune_stats(psth_col_nois(i,:)), 1:nrow_nos);
    stats_col_nois = arrayfun(@(i)calc_tune_stats(psth_col_nois(:,i)), 1:ncol_nos);
    % arrayfun(@(S)S.F_P,stats_row_nois)'
    stats_clas = calc_tune_stats(psth_col_clas);
    stats_row_clas = arrayfun(@(i)calc_tune_stats(psth_col_clas(i,:)), 1:nrow_cls);
    stats_col_clas = arrayfun(@(i)calc_tune_stats(psth_col_clas(:,i)), 1:ncol_cls);
    % arrayfun(@(S)S.F_P,stats_row_clas)'
    HessBGStats(Triali).class.stats_row = stats_row_clas;
    HessBGStats(Triali).class.stats_col = stats_col_clas;
    HessBGStats(Triali).noise.stats_row = stats_row_nois;
    HessBGStats(Triali).noise.stats_col = stats_col_nois;
    % population average response matrix for each image 
    HessBGStats(Triali).class.resp_mat = cell2mat(shiftdim(cellfun(@(idx)mean(rasters(:,51:200,idx),[2,3]), idx_arr_cls, 'Un', 0),-1));
    HessBGStats(Triali).noise.resp_mat = cell2mat(shiftdim(cellfun(@(idx)mean(rasters(:,51:200,idx),[2,3]), idx_arr_nos, 'Un', 0),-1));
    % population single trial response for each image, it's a cell array, as the trial number is different for each image
    HessBGStats(Triali).class.resp_sgtr_col = cellfun(@(idx)squeeze(mean(rasters(:,51:200,idx),[2])), idx_arr_cls, 'Un', 0);
    HessBGStats(Triali).noise.resp_sgtr_col = cellfun(@(idx)squeeze(mean(rasters(:,51:200,idx),[2])), idx_arr_nos, 'Un', 0);
    % population average PSTH for each image 
    HessBGStats(Triali).class.psth_col = cell2mat(shiftdim(cellfun(@(idx)mean(rasters(:,:,idx),[3]), idx_arr_cls, 'Un', 0),-2));
    HessBGStats(Triali).noise.psth_col = cell2mat(shiftdim(cellfun(@(idx)mean(rasters(:,:,idx),[3]), idx_arr_nos, 'Un', 0),-2));
end
%%
savefast(fullfile(MatStatsDir, Animal+"HessBigGANStats.mat"),"HessBGStats")
end


%%
function [noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq)
% see if there is `noise` or `class` in the imgname_uniq, if not, skip that part 
is_noise_exist = contains(imgname_uniq, "noise");
is_class_exist = contains(imgname_uniq, "class");
if sum(is_noise_exist) == 0
    fprintf("No noise space images found\n")
    noise_pattern = "";
    noise_imgnm = {};
else
% noise part
namepart_uniq = regexp(imgname_uniq,"noise_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)",'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
if ~isempty(namepart_uniq)
    decim = cellfun(@(nm)split(nm.dist,"."), namepart_uniq, 'uni', false);
    float_L = unique(cellfun(@(dcm)length(dcm{2}),decim));
    if float_L == 1
        fprintf("Old Image Naming Convention in Noise Space\n")
        noise_pattern = "noise_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
        noise_imgnm = "noise_eig%d_lin%.1f";
    elseif float_L == 2
        fprintf("Latest Line Search Image Naming Convention in Noise Space\n")
        noise_pattern = "noise_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
        noise_imgnm = "noise_eig%d_lin%.2f";
    end
else
    fprintf("New Image Naming Convention in Noise Space\n")
    namepart_uniq = regexp(imgname_uniq,"noise_eig(?<eig_id>\d*)_exp(?<expon>[.\d]*)_lin(?<dist>[-.\d]*)",'names');
    namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
    if ~isempty(namepart_uniq)
        expon_str = string(unique(cellfun(@(U)U.expon,namepart_uniq,"Uni",false)));
        noise_pattern = "noise_eig(?<eig_id>\d*)_exp(?<expon>[.\d]*)_lin(?<dist>[-.\d]*)";
        noise_imgnm = "noise_eig%d_exp"+expon_str+"_lin%.1f";
    else
        error("Name pattern not recognized.")
    end
end
end

if sum(is_class_exist) == 0
    fprintf("No class space images found\n")
    class_pattern = "";
    class_imgnm = {};
else
% class part
namepart_uniq = regexp(imgname_uniq,"class_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)",'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
if ~isempty(namepart_uniq)
    decim = cellfun(@(nm)split(nm.dist,"."), namepart_uniq, 'uni', false);
    float_L = unique(cellfun(@(dcm)length(dcm{2}),decim));
    if float_L == 1
        fprintf("Old Image Naming Convention in Class Space\n")
        class_pattern = "class_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
        class_imgnm = "class_eig%d_lin%.1f";
    elseif float_L == 2
        fprintf("Latest Line Search Image Naming Convention in Class Space\n")
        class_pattern = "class_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
        class_imgnm = "class_eig%d_lin%.2f";
    end
else
    fprintf("New Image Naming Convention in Class Space\n")
    namepart_uniq = regexp(imgname_uniq,"class_eig(?<eig_id>\d*)_exp(?<expon>[.\d]*)_lin(?<dist>[-.\d]*)",'names');
    namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
    if ~isempty(namepart_uniq)
        expon_str = string(unique(cellfun(@(U)U.expon,namepart_uniq,"Uni",false)));
        class_pattern = "class_eig(?<eig_id>\d*)_exp(?<expon>[.\d]*)_lin(?<dist>[-.\d]*)";
        class_imgnm = "class_eig%d_exp"+expon_str+"_lin%.1f";
    else
        error("Name pattern not recognized.")
    end
end
end
end