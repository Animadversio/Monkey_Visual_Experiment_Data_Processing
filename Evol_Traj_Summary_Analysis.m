%% Beto Evolution Trajectory Synthesis
all_genes=[];
for i = 1:numel(StatsB)
    all_genes = [all_genes;StatsB{i}.genes.all];
end
%%
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(all_genes, 'Centered', false, 'NumComponents', 50);
save('D:\\PC_all_evol.mat', 'SCORE', 'COEFF', 'EXPLAINED')
codes = COEFF;
save('D:\\PC_codes_all_evol.mat', 'codes')