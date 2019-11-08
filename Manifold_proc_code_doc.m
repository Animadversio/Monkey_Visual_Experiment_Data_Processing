%% Analysis Pipeline Documentation

%% Pre-Processing 
% Directly Process and Load the raw files and give back the structured data
Project_Manifold_Beto_loadRaw;
% Arrange the structures loaded, into a larger structure saved in mat. 
Append_Exp_to_Collection; 
%% Basic functions for processing
Set_Exp_Specs; % Set the sphere_norm and prefered channel for each experiment in order. May be done in automatic fashion
generate_unit_labels; % Generate labels from meta.spikeID