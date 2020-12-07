function [codes_all, img_ids, code_geni] = load_codes_all(stim_path, threadi, loadblocks)
% stim_path = "N:\Stimuli\2019-12-Evolutions\2020-03-10-Alfa-01\2020-03-10-13-50-57";
% threadi = 1;
% Note threadi starts from 1, so 2nd thread is threadi=2
if nargin == 1, threadi = 1; loadblocks=[]; end
if nargin == 2, loadblocks=[]; end
data_fn  = ls(fullfile(stim_path, sprintf("*_thread%03d_code.mat", threadi - 1)));
data_fn = sort(string(data_fn)); 
if contains(data_fn{1}, "block000"), offset = -1; elseif contains(data_fn{1}, "block001"), offset = 0; else error; end
% note older version starts from block000
codes_all = [];
code_geni = [];
img_ids = {};
if isempty(loadblocks)
blocks = 1:length(data_fn);
elseif strcmp(loadblocks,"last")
blocks = length(data_fn);
else
blocks = loadblocks;
end
for block_k = blocks
D = load(fullfile(stim_path, data_fn{block_k}),"codes", "ids");
assert(contains(data_fn{block_k}, sprintf("block%03d_thread%03d_code.mat", block_k + offset, threadi - 1)),...
    sprintf("Missing code mat file %s", sprintf("block%03d_thread%03d_code.mat", block_k, threadi - 1)))
codes_all = [codes_all; D.codes];
code_geni = [code_geni, block_k * ones(1, size(D.codes, 1))];
img_ids = [img_ids, D.ids];
end
img_ids = string(img_ids);
end