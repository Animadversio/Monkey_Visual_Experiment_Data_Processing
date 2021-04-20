%% Evol_Cosine_summary
saveroot = "O:\Evol_Cosine";
mkdir(saveroot)
refcoldir = "N:\Stimuli\2020-CosineEvol\RefCollection";
targmap = get_refimg_map(refcoldir);
%% Summary the evolution trajectory
subfdrs = string(ls(saveroot+"\*Alfa*"));
% figure;montage(fullfile(saveroot,deblank(subfdrs(6:end)),"online_scoretraj.png"),...
%         'Size', [5,8], 'BorderSize', 4,'ThumbnailSize',[656   875])
mtg = imtile(fullfile(saveroot,deblank(subfdrs(6:end)),"online_scoretraj.png"),'GridSize', [4,9],'BorderSize', 4);
imwrite(mtg,fullfile(saveroot,"CosineSummary","CosineTraj_summary.png"))
%%
function refnmMap = get_refimg_map(stim_dir,parent)
if nargin==1, parent=false; end
if parent,
parent_dir = fileparts(stim_dir);
else
parent_dir = stim_dir;
end
refimgs = dir(parent_dir);
refimgs = refimgs(~arrayfun(@(R)R.isdir,refimgs));
refnmMap = containers.Map();
for i = 1:numel(refimgs)
    [~,fn,ext] = fileparts(refimgs(i).name);
    refnmMap(string(fn)) = string(fullfile(refimgs(i).folder,refimgs(i).name));
end
end