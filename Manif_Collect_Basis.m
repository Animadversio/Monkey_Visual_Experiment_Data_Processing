%% Collect Statisitcs of the basis. 
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa"; 
load(fullfile(mat_dir, Animal+'_Evol_stats.mat')) 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat')) 
% load(fullfile(mat_dir, Animal+'_ManifPopDynamics.mat'),'ManifDyn')
%%
ManifBasisStats = repmat(struct(),numel(Stats),1);
for Expi = 18:numel(Stats)
spaceN = numel(Stats(Expi).manif.psth);
basis_path = fullfile(EStats(Expi).meta.stimuli,"PC_imgs","PC_vector_data.npz");
f = py.numpy.load(basis_path);
PC_Vec = f.get('PC_vecs').double;
sphere_norm = f.get('sphere_norm').double;
rnd_Vec = f.get("rand_vec2").double;
f.close();
% Get the mat containing all the codes of the last generation. To see
% whether we should inverse PC1 
matfns = string(ls(fullfile(EStats(Expi).meta.stimuli,"*.mat")));
code_tmp = load(fullfile(EStats(Expi).meta.stimuli,matfns(end)));
proj_coord = mean(code_tmp.codes,1) * PC_Vec';
basis = cell(1,spaceN);
basis{1} = [sign(proj_coord(1)),1,1]' .* PC_Vec(1:3,:);
if proj_coord(1) < 0
    % Note the final PC may need to reverse! not always the same dir! 
    fprintf("Exp %d The evolution direction is inverse to the PC1 direction of PCA. Inverse PC1 as basis\n", Expi)
end
if spaceN == 3
    basis{2} = [basis{1}(1,:); PC_Vec(49:50,:)];
    basis{3} = [basis{1}(1,:); rnd_Vec];
end
clear code_tmp matfns
ManifBasisStats(Expi).units.pref_chan = EStats(Expi).units.pref_chan;
ManifBasisStats(Expi).basis_path = basis_path;
ManifBasisStats(Expi).sphere_norm = sphere_norm;
ManifBasisStats(Expi).basis = basis; % cell array the same length of space
ManifBasisStats(Expi).subsp_n = spaceN;
end
%%
savefast(fullfile(mat_dir, Animal+'_ManifBasis.mat'), 'ManifBasisStats')
%%
% if Animal=="Alfa",EStats(19).meta.stimuli = "N:\\Stimuli\\2019-Manifold\\alfa-191210a\\backup_12_10_2019_13_07_57";end
% save(fullfile(mat_dir, Animal+"_Evol_stats.mat"), 'EStats')