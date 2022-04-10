function [chanX_broke, vec_broke] = prepare_broken_target(chanX_targ, vec_targ)
% Add nan to gaps in the chanX array such that plotting will not add
% straigt lines across the missing gaps.
% Example 
%   [chanX_broke, vec_broke] = prepare_broken_target([1,2,5,6,8], [10,5,56,4,5]);
% 
%   assert(isequaln(chanX_broke, [1     2   nan     5     6   nan     8]));
%   assert(isequaln(vec_broke, [10     5   nan    56     4   nan     5]));

delta = diff(chanX_targ);
if all(delta == 1)
chanX_broke = chanX_targ;
vec_broke = vec_targ;
else
breakidx = find(delta ~= 1);
chanX_broke = [];
vec_broke = [];
csr = 1;
for i = 1:numel(breakidx)
    idx = breakidx(i);
    chanX_broke = [chanX_broke, reshape(chanX_targ(csr:idx),1,[]), nan];
    vec_broke = [vec_broke, reshape(vec_targ(csr:idx),1,[]), nan];
    csr = breakidx(i) + 1;
end
chanX_broke = [chanX_broke, reshape(chanX_targ(csr:end),1,[])];
vec_broke = [vec_broke, reshape(vec_targ(csr:end),1,[])];
end
end