%% pasu_images_montage
% demo code for montage and visualizign the pasupathy images (Included in PC_space_Pasu_tuning_analysis.m)
% part of this code goes to score_frame_image_arr.m

%'pasu_%02d_ori_%02d_wg_f'
img_dir = 'S:\Stimuli\2019-Manifold\pasupathy-wg-f-4-ori';
img_list = {};
for j = 1:4
    for i = 1:51
        cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
        img_idx = find(contains(Trials.imageName, cur_fn));
        if length(img_idx) == 0
            img_list{i,j} = [];
        else
            img_list{i,j} = imread(fullfile(img_dir, [cur_fn, '.jpg']));
        end
    end
end
figure;
montage(img_list, 'Size', [4, 51]);
%% Montage with scors decorating the montage
hl_mat = nanmean(score_mat,3);
figure;
subplot(2,1,1)
h = imagesc(hl_mat');
axis image;
xticks(1:51)
yticks(0:90:270)
cmap = colormap();
[Cmin, Cmax] = caxis();
mod_img_list = cell(51,4);
for j = 1:4
    for i = 1:51
        cur_fn = sprintf('pasu_%02d_ori_%02d_wg_f', i, 2*j-1);
        img_idx = find(contains(Trials.imageName, cur_fn));
        if length(img_idx) == 0
            mod_img_list{i,j} = [];
        else
            scale_val = (hl_mat(i,j) - Cmin) / (Cmax - Cmin);
            c = interp1(cmap, scale_val * (size(cmap, 1) - 1) + 1);
            LineWidth = 50;
            pad_img = padarray(img_list{i,j}, [2*LineWidth, 2*LineWidth], 0);
            tmp_img = insertShape(pad_img, ...
                'Rectangle', [LineWidth,LineWidth,size(img_list{i,j})+2*LineWidth], ...
                'LineWidth', 2 * LineWidth, 'Color', 256 * c);
            mod_img_list{i,j} = tmp_img;
        end
    end
end
subplot(2,1,2)
montage(mod_img_list, 'Size', [4, 51]);