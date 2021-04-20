function [codes_all, img_ids, code_geni] = load_codes_all(stim_path, threadi, loadblocks)
% Parameters: 
% stim_path: backuped stimuli path to the evolution experiment 
%            e.g. "N:\Stimuli\2019-12-Evolutions\2020-03-10-Alfa-01\2020-03-10-13-50-57";
% threadi: thread number to load; default as 1. 
%            Note threadi starts from 1, so 2nd thread is threadi=2
% loadblocks: default to be empty, means every block. 
%            "last" means the final block 
%            it can be a list of block id or a single number
%            For numbering convention, it usually starts from block001,
%            unless in some case there is a block000 then starts from 0,
%            offset=-1
%            There is test to make sure the file name you are loading is 
%            "block%03d_thread%03d_code.mat"%( block_k + offset, threadi - 1)
if nargin == 1, threadi = 1; loadblocks=[]; end
if nargin == 2, loadblocks=[]; end
data_fn  = ls(fullfile(stim_path, sprintf("*_thread%03d_code.mat", threadi - 1)));
data_fn = sort(string(data_fn)); 
if contains(data_fn{1}, "block000"), offset = -1; fprintf("Saved codes starts from 000. offset -1\n");
elseif contains(data_fn{1}, "block001"), offset = 0; else error; end
% note older version starts from block000
codes_all = [];
code_geni = [];
img_ids = {};
if isempty(loadblocks)
blocks = 1:length(data_fn);
elseif strcmp(loadblocks,"last")
blocks = length(data_fn);
else
blocks = loadblocks; % blocks is basically the list of blocks to load.
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