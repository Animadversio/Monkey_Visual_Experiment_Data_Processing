function [idx_arr_nos, imgnm_arr_nos, idx_arr_cls, imgnm_arr_cls] = parse_image_idx_arr_hess(imageName)
% Adapt for code in Hess_BigGAN_Tune_analysis.m 
% Parse the `Trials.imageName` into 2 matrices of image name and 2 cell array
%  of indices of which trial show this image. 
imgname_uniq = unique(imageName);
[noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq);
% Extract parameters of images in noise space from img name
namepart_uniq = regexp(imgname_uniq, noise_pattern, 'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
eig_id_arr_nos = unique(cellfun(@(U)str2num(U.eig_id),namepart_uniq));
dist_arr_nos = unique(cellfun(@(U)str2num(U.dist),namepart_uniq));
% collect image name and trial indices in cell
imgnm_arr_nos = strings(length(eig_id_arr_nos), length(dist_arr_nos));
idx_arr_nos = cell(length(eig_id_arr_nos), length(dist_arr_nos));
for i = 1:numel(eig_id_arr_nos)
    eig_id = eig_id_arr_nos(i);
    for j = 1:numel(dist_arr_nos)
    dist = dist_arr_nos(j);
    imgnm = compose(noise_imgnm,eig_id,dist); % e.g. "noise_eig%d_lin%.1f"
    idx_arr_nos{i,j} = find(contains(imageName, imgnm));
    imgnm_arr_nos(i,j) = imgnm;
    if numel(idx_arr_nos{i,j}) == 0 && dist== 0 
    % Recently, we deleted images like "noise_eig2_lin0.0" so should use "noise_eig1_lin0.0" instead
    fprintf("Image %s doesn't exist, use the same image %s to substitute.\n",imgnm,compose(noise_imgnm,eig_id_arr_nos(1),0))
    imgnm = compose(noise_imgnm,eig_id_arr_nos(1),0);
    imgnm_arr_nos(i,j) = imgnm;
    idx_arr_nos{i,j} = find(contains(imageName, imgnm));
    end
    end
end
% Extract parameters of images in class space from img name
namepart_uniq = regexp(imgname_uniq,class_pattern,'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
eig_id_arr_cls = unique(cellfun(@(U)str2num(U.eig_id),namepart_uniq));
dist_arr_cls = unique(cellfun(@(U)str2num(U.dist),namepart_uniq));
% collect image name and trial indices in cell
imgnm_arr_cls = strings(length(eig_id_arr_cls), length(dist_arr_cls));
idx_arr_cls = cell(length(eig_id_arr_cls), length(dist_arr_cls));
for i = 1:numel(eig_id_arr_cls)
    eig_id = eig_id_arr_cls(i);
    for j = 1:numel(dist_arr_cls)
    dist = dist_arr_cls(j);
    imgnm = compose(class_imgnm,eig_id,dist); % "class_eig%d_lin%.1f"
    idx_arr_cls{i,j} = find(contains(imageName, imgnm));
    imgnm_arr_cls(i,j) = imgnm;
    if numel(idx_arr_cls{i,j}) == 0 && dist== 0 
    % Recently, we deleted images like "noise_eig2_lin0.0" so should use "noise_eig1_lin0.0" instead
    fprintf("Image %s doesn't exist, use the same image %s to substitute.\n",imgnm,compose(class_imgnm,eig_id_arr_cls(1),0))
    imgnm = compose(class_imgnm,eig_id_arr_cls(1),0);
    imgnm_arr_cls(i,j) = imgnm;
    idx_arr_cls{i,j} = find(contains(imageName, imgnm));
    end
    end
end
%  Plot the image tile in the experiment
% img_noise = arrayfun(@(imgnm)imread(fullfile(meta.stimuli, imgnm+".jpg")),imgnm_arr_nos,"Un",false);
% img_class = arrayfun(@(imgnm)imread(fullfile(meta.stimuli, imgnm+".jpg")),imgnm_arr_cls,"Un",false);
% cls_tile = imtile(img_class','GridSize',size(img_class));
% imwrite(cls_tile, fullfile(figdir, compose("Class_Images_Tile.jpg")))
% nos_tile = imtile(img_noise','GridSize',size(img_noise));
% imwrite(nos_tile, fullfile(figdir, compose("Noise_Images_Tile.jpg")))

nrow_cls = length(eig_id_arr_cls);
ncol_cls = length(dist_arr_cls);
nrow_nos = length(eig_id_arr_nos);
ncol_nos = length(dist_arr_nos);

end

function [noise_pattern, noise_imgnm, class_pattern, class_imgnm] = parse_naming_convention(imgname_uniq)
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