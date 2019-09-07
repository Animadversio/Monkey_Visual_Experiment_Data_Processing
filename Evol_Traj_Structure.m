exp_num = 4;
codes_all = StatsB{exp_num}.genes.all; % size(codes_all) = [img_num, code_dim]
generations = StatsB{exp_num}.genes.gen;
%%
figure(1)
imagesc(codes_all')
xlabel("image by generation")
ylabel("Code")
title(sprintf("Evolution %d: code overview", exp_num))
%% Norm vs Generation
code_norm = sum(codes_all.^2, 2);
%
figure(2)
scatter(generations, code_norm)
title("Code Norm grows with generation", 'Fontsize', 16)
xlabel("Generation #", 'Fontsize', 14)
ylabel("Norm of code", 'Fontsize', 14)
%%
tic
[PC_axiss,code_PC,PC_var,~,var_explained] = pca(codes_all,'NumComponents',10);
toc
%% See the cumulative variance curve. 
figure
plot(cumsum(var_explained))
%% 
figure(3)
for i=1:10
   subplot(10,1,i)
   scatter(generations, code_PC(:, i))
end
suptitle("PC Component Weights across Generation")%, 'Fontsize', 16)
%%
figure(4);hold on 
for i=1:10
   scatter(generations, code_PC(:, i),'filled')
end
xlabel("Generation #", 'Fontsize', 14)
ylabel("Weight of the PC vector", 'Fontsize', 14)
hold off
title("PC Component Weights across Generation", 'Fontsize', 16)
%%
figure(5);hold on 
for i=1:10
   plot(code_PC(:, i))
end
hold off
title("PC Component Weights across Generation", 'Fontsize', 16)