%% Analyze the evolution of the norm of the code
Norm_arr = {};
%%
codes_all = StatsB{1}.genes.all;
generations = StatsB{1}.genes.gen;
gen_num = max(generations) - min(generations) + 1;
%%
norms = sqrt(sum(codes_all.^2, 2));
mean_code_norm = zeros(1, gen_num);
norm_avg = zeros(1, gen_num);
for geni = min(generations):max(generations)
    norm_gen_avg = norm(mean(codes_all(generations==geni, :), 1));
    if min(generations) == 0
        mean_code_norm(geni+1) = norm_gen_avg;
        norm_avg(geni+1) = mean(norms(generations==geni, :));
    end
end
%%
code_in_gen_PC1_std = zeros(1, gen_num); 
code_in_gen_PC1_range = zeros(1, gen_num);
code_in_gen_PC2_std = zeros(1, gen_num); 
code_in_gen_PC2_range = zeros(1, gen_num);
for geni = min(generations):max(generations)
    [COEFF, SCORE] = pca(codes_all(generations==geni, :));
    if min(generations) == 0
        code_in_gen_PC1_std(geni+1) = std(SCORE(:,1));
        code_in_gen_PC1_range(geni+1) = range(SCORE(:,1));
        code_in_gen_PC2_std(geni+1) = std(SCORE(:,2));
        code_in_gen_PC2_range(geni+1) = range(SCORE(:,2));
    end
end
%%
figure;hold on 
scatter(generations, norms)
plot(min(generations):max(generations), norm_avg)
plot(min(generations):max(generations), mean_code_norm)
plot(min(generations):max(generations), code_in_gen_PC1_std)
plot(min(generations):max(generations), code_in_gen_PC1_range)
plot(min(generations):max(generations), code_in_gen_PC2_std)
plot(min(generations):max(generations), code_in_gen_PC2_range)
legend({'Norm of Code Sample', 'Mean of Norm', 'Norm of Mean Code', 'STD of PC1 proj',... 
    'Range of PC1 proj', 'STD of PC2 proj', 'Range of PC2 proj'}, 'FontSize', 10)
ylabel("Distance in Code Space")
xlabel("Generations")
hold off

%%
[COEFF, SCORE] = pca(codes_all(generations==geni, :));
% norm_gen_std = 1;
%%
figure;hold on 
plot(min(generations):max(generations), StatsB{1}.meanResp_syn)
% plot(min(generations):max(generations), StatsB{1}.maxResp_syn)
ylabel("Firing Rate")
xlabel("Generations")
legend({'Mean Response', 'Baseline'}, 'FontSize', 10)
hold off
%%
