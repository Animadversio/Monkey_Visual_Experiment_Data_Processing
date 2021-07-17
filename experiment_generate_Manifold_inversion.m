G = FC6Generator();
% Input the experimental backup folder containing the mat codes files.
global newimg_dir
newimg_dir = "N:\Stimuli\2021-ProjectPFC\Manifold_Demo\Banana_manifold";
mkdir(newimg_dir)
fprintf("Save new images to folder %s\n", newimg_dir)
%% Spherical interpolation
orig_img = imread("N:\Stimuli\2021-ProjectPFC\2021-Selectivity\bigSet\fruit (8).jpg");
target_img = imresize(orig_img, [256,256]);
[code_invert, img_fit, loss] = GAN_invert_fun(G, target_img, 200);
code_invert = double(code_invert);
orthonormalize2vect(code_invert, randn(2, 4096))
%%
[img_list] = explore_from_code(G,code_invert, rand_vec2, "RND12", ["RND1","RND2"]);
%%

function vects = orthonormalize2vect(centvect, vects)
% tmp = orthonormalize2vect([1,zeros(1,4095)],randn(5,4096));
% tmp = orthonormalize2vect(randn(5,4096));
if ~isempty(centvect) && nargin == 2
vects = vects - (vects * centvect') * centvect / vecnorm(centvect,2,2).^2;
elseif ~isempty(centvect) && nargin==1
vects=centvect;
end
for i = 2:size(vects,1)
    vects(i:end,:) = vects(i:end,:) - (vects(i:end, :) * vects(i-1, :)') * vects(i-1, :) / vecnorm(vects(i-1, :),2,2).^2;
end
vects = vects ./ vecnorm(vects,2,2);
end
% rand_vec2 = randn(2, 4096);
% rand_vec2 = rand_vec2 - (rand_vec2 * code_invert') * code_invert / vecnorm(code_invert,2,2).^2;
% rand_vec2(2,:) = rand_vec2(2,:) - (rand_vec2(2,:) * rand_vec2(1,:)') * rand_vec2(1,:) / vecnorm(rand_vec2(1,:),2,2).^2;
% rand_vec2 = rand_vec2 ./ vecnorm(rand_vec2,2,2); %np.sqrt((rand_vec2**2).sum(axis=1))[:, np.newaxis]

function [img_list] = explore_from_code(G, code_invert, tang_vecs, space_str, axes_str)
global newimg_dir
ANGLE_SPAN = 180;
imgN_per_arc = 11; 
assert(length(code_invert)==4096)
code_invert = reshape(code_invert,1,4096);
if size(tang_vecs,1)==4096
    tang_vecs = tang_vecs';
end
assert(size(tang_vecs,2)==4096)
Theta_step = ANGLE_SPAN / (imgN_per_arc-1);
Phi_step = ANGLE_SPAN / (imgN_per_arc-1);
sphere_norm = vecnorm(code_invert,2,2);
fprintf("norm of the inverted code %.2f\n", sphere_norm)
fprintf("Set sphere norm to the last generations norm!\n")
basis1 = code_invert / sphere_norm;
basis = [basis1; tang_vecs];
fprintf("Generating images on %s sphere (radius = %.1f)\n", space_str, sphere_norm)
img_list = {};
for j = -5:5
    for k = -5:5
        theta = Theta_step * j / 180 * pi;
        phi = Phi_step * k / 180 * pi;
        code_vec = [cos(theta)*cos(phi),...
                    sin(theta)*cos(phi),...
                    sin(phi)] * basis;
        code_vec = code_vec / norm(code_vec) * sphere_norm;
        img = G.visualize(code_vec);
        img_list{end+1} = img;
        imwrite(img, fullfile(newimg_dir, compose("norm_%d_%s_%d_%s_%d.jpg", ...
            sphere_norm, axes_str(1), Theta_step * j, axes_str(2), Phi_step* k)));
    end
end
mtg = imtile(img_list,'GridSize',[11,11], 'BorderSize',4,'ThumbnailSize',[256,256]);
imwrite(mtg,fullfile(newimg_dir,compose("%s_Montage.png", space_str)));
save(fullfile(newimg_dir, compose("%s_tan_vect_data.mat", space_str)), ...
    "code_invert", "tang_vecs", "sphere_norm", "Theta_step", "Phi_step");
end
