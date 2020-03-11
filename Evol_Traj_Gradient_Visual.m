clearvars -except meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug ExpRecord
%% Visualizing Code Evolution Traj


% Loading the Exp data

% Sort image name and get scores in the given window

% Find and load all codes 
[codes_all, img_ids, code_geni] = load_codes_all(stim_path, threadi);
% load Generator
G = FC6Generator("matlabGANfc6");
% For each generation in the experiment 

% get all codes in this gen, basis 

% do local PCA on codes? Or norm ? 

% Plot the image at the given location 

% Visualize the gradient / the next sample 
