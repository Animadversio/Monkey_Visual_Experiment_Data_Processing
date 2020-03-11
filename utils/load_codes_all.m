function [codes_all, img_ids, code_geni] = load_codes_all(stim_path, threadi)
% stim_path = "N:\Stimuli\2019-12-Evolutions\2020-03-10-Alfa-01\2020-03-10-13-50-57";
% threadi = 1;
if nargin == 1
    threadi = 1;
end
data_fn  = ls(fullfile(stim_path, sprintf("*_thread%03d_code.mat", threadi - 1)));
data_fn = sort(string(data_fn)); 
codes_all = [];
code_geni = [];
img_ids = {};
for block_k = 1:length(data_fn)
D = load(fullfile(stim_path, data_fn{block_k}),"codes", "ids");
assert(contains(data_fn{block_k}, sprintf("block%03d_thread%03d_code.mat", block_k, threadi - 1)),...
    sprintf("Missing code mat file %s", sprintf("block%03d_thread%03d_code.mat", block_k, threadi - 1)))
codes_all = [codes_all; D.codes];
code_geni = [code_geni, block_k * ones(1, size(D.codes, 1))];
img_ids = [img_ids, D.ids];
end
img_ids = string(img_ids);
end