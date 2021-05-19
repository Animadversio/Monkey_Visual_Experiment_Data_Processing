function refnmMap = get_refimg_map(stim_dir,parent)
% given a folder `stim_dir` return the map from image names to the full paths
% handy function to load the reference images from their name without suffix 
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