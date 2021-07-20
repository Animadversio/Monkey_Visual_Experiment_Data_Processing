G = FC6Generator();
[eigvals,eigvect] = loadH();
% Input the experimental backup folder containing the mat codes files.
global newimg_dir
%%

%% Spherical interpolation
manifold_root = "C:\Users\ponce\Documents\monkeylogic2\selectivityBasic\manifold";
imgname = "Apple";
orig_img = imread("C:\Users\ponce\Documents\monkeylogic2\selectivityBasic\catsA\fruit (1).jpeg");
target_img = imresize(orig_img, [256,256]);
newimg_dir = fullfile(manifold_root, compose("%s_manifold",imgname));
mkdir(newimg_dir);fprintf("Save new images to folder %s\n", newimg_dir);
[code_invert, img_fit, loss] = GAN_invert_fun(G, target_img, 300);
figure;imshow(imtile({target_img,img_fit}))
code_invert = double(code_invert);
%%
rand_tanvec2 = orthonormalize2vect(code_invert, randn(2, 4096));
[img_list] = explore_from_code(G,code_invert, rand_tanvec2, "RND12", ["RND1","RND2"],148,9);
%%
EIG_tanvec2 = orthonormalize2vect(code_invert, eigvect(:, 1:2)');
[img_list] = explore_from_code(G,code_invert, EIG_tanvec2, "eig12", ["eig1","eig2"],30,7);
%%
EIG_tanvec2 = orthonormalize2vect(code_invert, eigvect(:, [8,10])');
[img_list] = explore_from_code(G,code_invert, EIG_tanvec2, "eig810", ["eig8","eig10"],40,7,false);

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

function [img_list] = explore_from_code(G, code_invert, tang_vecs, space_str, axes_str, ...
    ANGLE_SPAN, imgN_per_arc, DOSAVE)
% code_invert: center code, that you invert.
% tang_vecs: tangent vectors to explore along. 
% space_str: name of the space, for file labelling. e.g. "eig12"
% axes_str: string array containing the names of the axes. ["eig1","eig2"]
% ANGLE_SPAN: angle span of exploration, two sided. e.g. 148
% imgN_per_arc: final grid will be imgN_per_arc-by-imgN_per_arc, step size
%               will be ANGLE_SPAN/(imgN_per_arc - 1). e.g. 7
% DOSAVE: if you are trying out parameters set DOSAVE to false.
global newimg_dir
if nargin<=5, ANGLE_SPAN = 180; end
if nargin<=6, imgN_per_arc = 11; end
if nargin<=7, DOSAVE = true; end
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
for j = -(imgN_per_arc-1)/2:(imgN_per_arc-1)/2
    for k = -(imgN_per_arc-1)/2:(imgN_per_arc-1)/2
        theta = Theta_step * j / 180 * pi;
        phi = Phi_step * k / 180 * pi;
        code_vec = [cos(theta)*cos(phi),...
                    sin(theta)*cos(phi),...
                    sin(phi)] * basis;
        code_vec = code_vec / norm(code_vec) * sphere_norm;
        img = G.visualize(code_vec);
        img_list{end+1} = img;
        if DOSAVE
        imwrite(img, fullfile(newimg_dir, compose("norm_%d_%s_%d_%s_%d.jpg", ...
            int32(sphere_norm), axes_str(1), int32(Theta_step * j), axes_str(2), int32(Phi_step* k))));
        end
    end
end
mtg = imtile(img_list,'GridSize',[imgN_per_arc,imgN_per_arc], 'BorderSize',4,'ThumbnailSize',[256,256]);
figure;imshow(mtg)
if DOSAVE
imwrite(mtg,fullfile(newimg_dir, compose("norm_%d_Span%d_%dperArc_%s_Montage.png", int32(sphere_norm), ANGLE_SPAN, imgN_per_arc, space_str)));
save(fullfile(newimg_dir, compose("%s_tan_vect_data.mat", space_str)), ...
    "code_invert", "tang_vecs", "sphere_norm", "Theta_step", "Phi_step");
end
end

function [eigvals,eigvect] = loadH()
py.importlib.import_module("numpy");
data = py.numpy.load("N:\Data-Computational\Evolution_Avg_Hess.npz");
eigvect = data.get("eigvect_avg").double; % each column is an eigen vector. 
eigvals = data.get("eigv_avg").double;
eigvals = eigvals(end:-1:1);
eigvect = eigvect(:,end:-1:1);
end