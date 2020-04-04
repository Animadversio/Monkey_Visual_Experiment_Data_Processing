%% Getting Stats
function sphere_norm = infer_norm_from_imgname(imgnames, pattern)
    if nargin == 1
        pattern = "PC2"; % by default it's the PC2 PC3 space
    end
    tmp = cellfun(@(c) regexp(c,strcat("norm_(?<norm>\d*)_", pattern),'names'), imgnames, 'UniformOutput', false);      
    extnorms = cellfun(@(c) str2num(c.norm), tmp(~cellfun('isempty',tmp)));
    sphere_norm = mode(extnorms);
end