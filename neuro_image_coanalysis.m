%% Hard Image neural co-factorization

%%
image_dir = 'D:\Monkey_Data\2019-06-Evolutions\beto-190610a\BACKUP_3';
code_dir = 'D:\Monkey_Data\2019-06-Evolutions\beto-190610a';
% imglist = ls([image_dir,'\block*.jpg']);
% codelist = ls([code_dir,'\block*.mat']);
imglist = dir(fullfile(image_dir, 'block*.jpg'));
img_fn = {imglist.name};
img_fullfn = fullfile(image_dir, {imglist.name});
codelist = dir(fullfile(code_dir, 'block*.mat'));
code_fn = {codelist.name};
code_fullfn = fullfile(code_dir, {codelist.name});
%%
score_arr = []; % Note 
id_arr = {};
code_arr = [];
for i =1:numel(code_fullfn) 
   data = load(code_fullfn{i}); 
   block_i = i - 1;
   while ~contains(code_fullfn{i}, sprintf("block%03d", block_i))
       block_i = block_i + 1;
   end
   % Note the scores are scores in last generation and codes and ids are  
   % for next generation so the code and score do not match to each other
   img_ids = cellfun(@(id) strcat(sprintf("block%03d_gen_", block_i), id), data.ids); 
   code_arr = cat(1, code_arr, data.codes);
   id_arr = cat(2, id_arr, img_ids); % note there 
   if i~=1
       score_arr = cat(1, score_arr, data.scores);
   end
   if i == numel(code_fullfn)
       score_arr = cat(1, score_arr, nan([length(data.ids), 1]));
   end
end
clear img_ids block_i
%% Get Trial and 
if isempty(rasters)
    load('D:\Monkey_Data\monkey64chan-10062019-001_formatted.mat')
end
rasters = permute(rasters,[3,2,1]);
lfps = permute(lfps,[3,2,1]);
%% Sort the trial image information
unsort_evov_idx = find(strncmp('block', Trials.imageName, 5));
%%
[Sort_Name, sort_idx] = sort(Trials.imageName);
gen_num_i = zeros([length(Sort_Name),1]);
block_i = 0;
natural_stim_i = zeros([length(Sort_Name),1]);
stim_i = 1;
for i = 1:length(Sort_Name)
    if ~contains(Sort_Name{i}, "block")
        while ~contains(Sort_Name{i}, sprintf("[%02d]", stim_i))
            stim_i = stim_i + 1;
        end
        natural_stim_i(i) = stim_i;
        gen_num_i(i) = -1;
    else
        while ~contains(Sort_Name{i}, sprintf("block%03d", block_i))
            block_i = block_i + 1;
        end
        gen_num_i(i) = block_i;
    end
end
%% Find the code correspond to each trial. 
evolv_mask = (gen_num_i > -1);
evolv_img_names = Sort_Name(evolv_mask);
showed_code_arr = zeros(sum(evolv_mask), size(code_arr, 2));
for i=1:numel(evolv_img_names)
    row_i = find(strcmp(id_arr, evolv_img_names{i}));
    showed_code_arr(i, :) = code_arr(row_i, :);
end
clear row_i evolv_mask
%%
save("showed_code_arr.mat", 'showed_code_arr', 'evolv_img_names');
%%
channel_j = 4;
% cluster_input = rasters(sort_idx(gen_num_i~=-1), :, channel_j);
cluster_input = rasters(sort_idx(gen_num_i~=-1), :, channel_j);
Z = linkage(cluster_input, 'average', 'correlation');%''correlation
figure
[~, T , perm_indx] = dendrogram(Z,0);
figure("Position",[0,0,1500,500])
imagesc(cluster_input(perm_indx, :)',[0,500])
colorbar()
title(sprintf("Evolving Image PSTH Sorted by Hierachical clustering Channel %d",channel_j))
ylabel("Time")
xlabel("Sorted Image id")
%%
figure("Position",[0,0,1500,500])
imagesc(showed_code_arr(:, :)',[-10,10])
colorbar()
%%
figure("Position",[0,0,1500,500])
imagesc(squeeze(rasters(unsort_evov_idx, :, channel_j))')
colorbar()

%%
figure("Position",[0,0,1500,500])
imagesc(cluster_input(:, :)',[0,500])
%%
n_compon = 30 ;
% [img_mix,psth_basis,D] = nnmf(NMF_input, n_compon, 'algorithm','mult') ;
[psth_code_basis,img_mix,D] = nnmf([cluster_input, showed_code_arr]', n_compon);
%%
evolv_mask = (gen_num_i > -1);
showed_image_arr = zeros(sum(evolv_mask), 256*256*3,'single');
for i=1:numel(evolv_img_names)
    row_i = find(strncmp(evolv_img_names{i}, img_fn, length(evolv_img_names{i})));
    img = imread(img_fullfn{i});
    showed_image_arr(i, :) = img(:);
end
%%
%%
n_compon = 10;
[psth_img_basis,img_mix,D] = nnmf([cluster_input, showed_image_arr]', n_compon);
%%
figure;
imagesc(psth_img_basis(1:200,:))
%%
figure;
for i =1:n_compon
    subplot(1,n_compon,i)
    imshow(reshape(uint8(psth_img_basis(201:end, i)*mean(img_mix(i, :))), size(img)))
end
%%
nmf_img = reshape(uint8(psth_img_basis(201:end,1)), size(img));
%%
figure;
for i =1:n_compon
    subplot(3,n_compon,i)
    nmf_img = psth_img_basis(201:end, i)*0.7*max(img_mix(i, :));
    imshow(reshape(uint8(nmf_img), size(img)))
    subplot(3,n_compon,i+n_compon)
    plot(psth_img_basis(1:200, i))
    subplot(3,n_compon,i+2*n_compon)
    plot(img_mix(i, :)')
end
%% Sanity check how well it reconstruct image 
reconstr_psth = psth_img_basis(1:200, :) * img_mix(:, perm_indx);
figure;
imagesc(reconstr_psth)
%%
%%
tic
n_compon = 30;
% [psth_img_basis,cof_img_mix,D] = nnmf([cluster_input, showed_image_arr/256^2]', n_compon);
[psth_img_basis,cof_img_mix,iter] = nmf([cluster_input, showed_image_arr/256]', n_compon,...
    'verbose',0,'method','hals');
toc
save('image_psth_nmf30_norm.mat', 'psth_img_basis', 'cof_img_mix', 'iter');
%%
figure;
n_compon = 8;
for sub_i =1:n_compon
    i = sub_i + 22;
    subplot(3,n_compon,sub_i)
    nmf_img = psth_img_basis(201:end, i)*0.7*max(cof_img_mix(i, :))*256;
    imshow(reshape(uint8(nmf_img), size(img)))
    subplot(3,n_compon,sub_i+n_compon)
    plot(psth_img_basis(1:200, i))
    subplot(3,n_compon,sub_i+2*n_compon)
    plot(cof_img_mix(i, :)')
end
toc
%%
reconstr_psth2 = psth_img_basis(1:200, :) * cof_img_mix(:, perm_indx);
figure;
imagesc(reconstr_psth2)
%%
tic
n_compon = 20;
% [psth_img_basis,cof_img_mix,D] = nnmf([cluster_input, showed_image_arr/256^2]', n_compon);
[psth_img_basis,cof_img_mix,iter] = nmf([cluster_input, showed_image_arr/32]', n_compon,...
    'verbose',0,'method','hals');
toc
save('image_psth_nmf20_norm32.mat', 'psth_img_basis', 'cof_img_mix', 'iter');
%%
figure;
n_compon = 10;
for sub_i =1:n_compon
    i = sub_i+10;
    subplot(3,n_compon,sub_i)
    nmf_img = psth_img_basis(201:end, i)*0.7*max(cof_img_mix(i, :))*32;
    imshow(reshape(uint8(nmf_img), size(img)))
    subplot(3,n_compon,sub_i+n_compon)
    plot(psth_img_basis(1:200, i))
    subplot(3,n_compon,sub_i+2*n_compon)
    plot(cof_img_mix(i, :)')
end
toc
%%
reconstr_psth2 = psth_img_basis(1:200, :) * cof_img_mix(:, perm_indx);
figure;
imagesc(reconstr_psth2)
title("Reconstruction PSTH from 20 component NMF, Image normalization 32, HALS")
%%
tic
n_compon = 15;
% [psth_img_basis,cof_img_mix,D] = nnmf([cluster_input, showed_image_arr/256^2]', n_compon);
[psth_img_basis,cof_img_mix,iter] = nmf([cluster_input, showed_image_arr/16]', n_compon,...
    'verbose',0,'method','hals');
toc
save('image_psth_nmf15_norm16.mat', 'psth_img_basis', 'cof_img_mix', 'iter');
%%
figure('Position',[0         500        1916         450]);
pos = get(gca, 'Position');
pos(1) = 0.055;
pos(3) = 0.9;
set(gca, 'Position', pos)
n_compon = 15;
for sub_i =1:n_compon
    %i = sub_i+10;
    i = sub_i ;
    subplot(3,n_compon,sub_i)
    nmf_img = psth_img_basis(201:end, i)*0.7*max(cof_img_mix(i, :))*16;
    imshow(reshape(uint8(nmf_img), size(img)))
    subplot(3,n_compon,sub_i+n_compon)
    plot(psth_img_basis(1:200, i))
    subplot(3,n_compon,sub_i+2*n_compon)
    plot(cof_img_mix(i, :)')
    xlim([0,length(perm_indx)])
end
toc
%%
reconstr_psth2 = psth_img_basis(1:200, :) * cof_img_mix(:, perm_indx);
figure('Position',[0         500        1915         650]);
imagesc(reconstr_psth2)
title("Reconstruction PSTH from 15 component NMF, Image normalization 16, HALS")

%%
tic
n_compon = 15;
% [psth_img_basis,cof_img_mix,D] = nnmf([cluster_input, showed_image_arr/256^2]', n_compon);
[psth_img_basis,cof_img_mix,iter] = nmf([cluster_input, showed_image_arr/64]', n_compon,...
    'verbose',0,'method','hals');
toc
save('image_psth_nmf15_norm64.mat', 'psth_img_basis', 'cof_img_mix', 'iter');
%
figure('Position',[0         500        1916         450]);
pos = get(gca, 'Position');
pos(1) = 0.055;
pos(3) = 0.9;
set(gca, 'Position', pos)
n_compon = 15;
for sub_i =1:n_compon
    %i = sub_i+10;
    i = sub_i ;
    subplot(3,n_compon,sub_i)
    nmf_img = psth_img_basis(201:end, i)*0.7*max(cof_img_mix(i, :))*64;
    imshow(reshape(uint8(nmf_img), size(img)))
    subplot(3,n_compon,sub_i+n_compon)
    plot(psth_img_basis(1:200, i))
    subplot(3,n_compon,sub_i+2*n_compon)
    plot(cof_img_mix(i, :)')
    xlim([0,length(perm_indx)])
end
suptitle("Image and response Motif from 15 component NMF, Image normalization 64, HALS")
toc
%
reconstr_psth2 = psth_img_basis(1:200, :) * cof_img_mix(:, perm_indx);
figure('Position',[0         500        1915         650]);
imagesc(reconstr_psth2)
title("Reconstruction PSTH from 15 component NMF, Image normalization 64, HALS")

%%
tic
n_compon = 15;
% [psth_img_basis,cof_img_mix,D] = nnmf([cluster_input, showed_image_arr/256^2]', n_compon);
[psth_img_basis,cof_img_mix,iter] = nmf([cluster_input, showed_image_arr/32]', n_compon,...
    'verbose',0,'method','hals');
toc
save('image_psth_nmf15_norm32.mat', 'psth_img_basis', 'cof_img_mix', 'iter');
%
figure('Position',[0         500        1916         450]);
pos = get(gca, 'Position');
pos(1) = 0.055;
pos(3) = 0.9;
set(gca, 'Position', pos)
n_compon = 15;
for sub_i =1:n_compon
    %i = sub_i+10;
    i = sub_i ;
    subplot(3,n_compon,sub_i)
    nmf_img = psth_img_basis(201:end, i)*0.7*max(cof_img_mix(i, :))*32;
    imshow(reshape(uint8(nmf_img), size(img)))
    subplot(3,n_compon,sub_i+n_compon)
    plot(psth_img_basis(1:200, i))
    subplot(3,n_compon,sub_i+2*n_compon)
    plot(cof_img_mix(i, :)')
    xlim([0,length(perm_indx)])
end
suptitle("Image and response Motif from 15 component NMF, Image normalization 32, HALS")
toc
%
reconstr_psth2 = psth_img_basis(1:200, :) * cof_img_mix(:, perm_indx);
figure('Position',[0         500        1915         650]);
imagesc(reconstr_psth2)
title("Reconstruction PSTH from 15 component NMF, Image normalization 32, HALS")