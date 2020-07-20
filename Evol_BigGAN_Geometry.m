%% BigGAN Evolution Geometry
pyenv("Version","C:\Users\ponce\.conda\envs\caffe36\python.exe")

%%
G = torchBigGAN("biggan-deep-256");
Embed_mat = G.get_embedding();
[BGcodes_all, ~, BGcode_geni] = load_codes_all("N:\Stimuli\2020-BigGAN\2020-07-20-Beto-01\2020-07-20-12-40-51",1);
[data,MLConfig,TrialRecord,~] = mlread("N:\Data-Behavior (BHV2)\200720_Beto_generate_BigGAN(1).bhv2");
space = TrialRecord.User.space_cfg{1};
if strcmp(space, 'BigGAN_class')
noise_vec = TrialRecord.User.space_cfg{2};
end
truncnorm = truncate(makedist("Normal"),-2,2);
%%
net = alexnet;
ImgNetLabel = net.Layers(end).Classes;
clear net
%%
noise_vec = truncnorm.random([1,128]);
%%
G = G.select_space(space, noise_vec);
imgs = G.visualize(BGcodes_all(BGcode_geni==max(BGcode_geni),:));figure;montage(imgs)
%%
traj_distmat = pdist2(BGcodes_all, Embed_mat',"correlation"); % "euclidean"
[sort_dist, sort_idx] = sort(traj_distmat, 2, "Ascend");
NNClass = ImgNetLabel(sort_idx(:,1:5));
%%
for geni=1:max(BGcode_geni)
for rowi = find(BGcode_geni==geni)
fprintf("\n%02d ", 1 + rowi - min(find(BGcode_geni==geni)))
for i = 1:5
    fprintf("%s %.1f\t\t",NNClass(rowi,i),sort_dist(rowi,i))
end
fprintf("\n")
end
imgs = G.visualize(BGcodes_all(BGcode_geni==geni,:));
figure(8);montage(imgs)
pause;
end

