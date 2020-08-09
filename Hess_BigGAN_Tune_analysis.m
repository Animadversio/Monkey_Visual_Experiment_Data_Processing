%% Adapted from Hess_Tune_Analysis
Animal = "Both";Set_Path;
%"200803","200804","200805"
expftr = contains(ExpRecord.expControlFN,["200806"]); %& contains(ExpRecord.Exp_collection,"BigGAN_Hessian");% & contains(ExpRecord.Exp_collection,"BigGAN");
fllist = find(expftr);no_return=false;
[meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(fllist(1:end),Animal,no_return);

figdir = "E:\OneDrive - Washington University in St. Louis\HessBigGANTune\Alfa_Exp02";
%%
Triali = 3;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
% lfps = lfps_new{Triali};
Trials = Trials_new{Triali};

if contains(meta.ephysFN,"Alfa"),Animal = "Alfa";elseif contains(meta.ephysFN,"Beto"),Animal = "Beto";end
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
imgname_uniq = unique(Trials.imageName); 
%% Load the images in the class space and noise space
% noise space
[noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq);

namepart_uniq = regexp(imgname_uniq, noise_pattern, 'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
eig_id_arr_nos = unique(cellfun(@(U)str2num(U.eig_id),namepart_uniq));
dist_arr_nos = unique(cellfun(@(U)str2num(U.dist),namepart_uniq));

imgnm_arr_nos = strings(length(eig_id_arr_nos), length(dist_arr_nos));
idx_arr_nos = cell(length(eig_id_arr_nos), length(dist_arr_nos));
for i = 1:length(eig_id_arr_nos)
    eig_id = eig_id_arr_nos(i);
    for j = 1:length(dist_arr_nos)
    dist = dist_arr_nos(j);
    imgnm = compose(noise_imgnm,eig_id,dist); % "noise_eig%d_lin%.1f"
    idx_arr_nos{i,j} = find(contains(Trials.imageName, imgnm));
    imgnm_arr_nos(i,j) = imgnm;
    end
end
img_noise = arrayfun(@(imgnm)imread(fullfile(meta.stimuli, imgnm+".jpg")),imgnm_arr_nos,"Un",false);
% class space
namepart_uniq = regexp(imgname_uniq,class_pattern,'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
eig_id_arr_cls = unique(cellfun(@(U)str2num(U.eig_id),namepart_uniq));
dist_arr_cls = unique(cellfun(@(U)str2num(U.dist),namepart_uniq));

imgnm_arr_cls = strings(length(eig_id_arr_cls), length(dist_arr_cls));
idx_arr_cls = cell(length(eig_id_arr_cls), length(dist_arr_cls));
for i = 1:length(eig_id_arr_cls)
    eig_id = eig_id_arr_cls(i);
    for j = 1:length(dist_arr_cls)
    dist = dist_arr_cls(j);
    imgnm = compose(class_imgnm,eig_id,dist); % "class_eig%d_lin%.1f"
    idx_arr_cls{i,j} = find(contains(Trials.imageName, imgnm));
    imgnm_arr_cls(i,j) = imgnm;
    end
end
img_class = arrayfun(@(imgnm)imread(fullfile(meta.stimuli, imgnm+".jpg")),imgnm_arr_cls,"Un",false);
%% Plot the image tile in the experiment
cls_tile = imtile(img_class','GridSize',size(img_class));
imwrite(cls_tile, fullfile(figdir, compose("Class_Images_Tile.jpg")))
nos_tile = imtile(img_noise','GridSize',size(img_noise));
imwrite(nos_tile, fullfile(figdir, compose("Noise_Images_Tile.jpg")))
%% Plot the neural response to image tiles 
prefchan_ids = find(meta.spikeID==Trials.TrialRecord.User.prefChan);
for iCh = 1:size(rasters,1) % prefchan_ids
resp_mat_nois = cellfun(@(idx)mean(rasters(iCh,51:200,idx),'all'), idx_arr_nos);
resp_mat_cls = cellfun(@(idx)mean(rasters(iCh,51:200,idx),'all'), idx_arr_cls);
% Score Heatmap in Noise space
figure(12);set(12,'position',[315   506   560   444])
imagesc(resp_mat_nois)
colorbar();axis image
ylabel("eigen id")
xlabel("Linear Distance")
xticks(1:length(dist_arr_nos));xticklabels(dist_arr_nos)
yticks(1:length(eig_id_arr_nos));yticklabels(eig_id_arr_nos)
title(compose("%s %s Pilot Exp\nBigGAN Noise Space Tuning Along Hessian EigenVectors",Animal, unit_name_arr(iCh)))
saveas(12, fullfile(figdir, compose("Noise_TuneMap_%s.jpg",unit_name_arr(iCh))))
% Score Heatmap in Class space
figure(13);set(13,'position',[315   159   560   608])
imagesc(resp_mat_cls)
colorbar();axis image
ylabel("eigen id")
xlabel("Linear Distance")
xticks(1:length(dist_arr_cls));xticklabels(dist_arr_cls)
yticks(1:length(eig_id_arr_cls));yticklabels(eig_id_arr_cls)
title(compose("%s %s Pilot Exp\nBigGAN Class Space Tuning Along Hessian EigenVectors",Animal, unit_name_arr(iCh)))
saveas(13, fullfile(figdir, compose("Class_TuneMap_%s.jpg",unit_name_arr(iCh))))
% score framed image tile in class and noise space
frame_img_class = score_frame_image_arr(img_class, resp_mat_cls, prctile(resp_mat_cls,[0,100],'all')', parula, 10);
frame_img_noise = score_frame_image_arr(img_noise, resp_mat_nois, prctile(resp_mat_nois,[0,100],'all')', parula, 10);
cls_scr_tile = imtile(frame_img_class','GridSize',[length(eig_id_arr_cls),11]);
imwrite(cls_scr_tile, fullfile(figdir, compose("Class_TuneTile_%s.jpg",unit_name_arr(iCh))))
nos_scr_tile = imtile(frame_img_noise','GridSize',[length(eig_id_arr_nos),11]);
imwrite(nos_scr_tile, fullfile(figdir, compose("Noise_TuneTile_%s.jpg",unit_name_arr(iCh))))
end
%%
% frame_img_noise = score_frame_image_arr(img_noise, resp_mat_nois, prctile(resp_mat_nois,[0,100],'all')', parula, 10);
% figure;
% montage(frame_img_noise','Size',[length(eig_id_arr_nos),11])
% 
% frame_img_class = score_frame_image_arr(img_class, resp_mat_cls, prctile(resp_mat_cls,[0,100],'all')', parula, 10);
% figure;
% montage(frame_img_class','Size',[length(eig_id_arr_cls),11])

% figure(3);
% montage(img_class','Size',[14,11])
% title(compose("%s Pilot Exp\nBigGAN Class Space Interpolation Along Hessian EigenVectors",Animal))
% saveas(3, fullfile(figdir, compose("Class_Images.jpg")))

% figure(4);
% montage(img_noise','Size',[10,11])
% title(compose("%s Pilot Exp\nBigGAN Noise Space Interpolation Along Hessian EigenVectors",Animal))
% saveas(4, fullfile(figdir, compose("Noise_Images.jpg")))
%%
[noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq);
namepart_uniq = regexp(imgname_uniq,noise_pattern,'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));


function [noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq)
% noise part
namepart_uniq = regexp(imgname_uniq,"noise_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)",'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
if ~isempty(namepart_uniq)
    fprintf("Old Image Naming Convention in Noise Space\n")
    noise_pattern = "noise_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
    noise_imgnm = "noise_eig%d_lin%.1f";
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
% class part
namepart_uniq = regexp(imgname_uniq,"class_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)",'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
if ~isempty(namepart_uniq)
    fprintf("Old Image Naming Convention in Class Space\n")
    class_pattern = "class_eig(?<eig_id>\d*)_lin(?<dist>[-.\d]*)";
    class_imgnm = "class_eig%d_lin%.1f";
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