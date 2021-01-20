% This code target for 
Data = load("N:\Code\Generator_DB_Windows\init_population\init_code.mat");
initcodes = Data.codes;
%%
ctrlcodes_col = [];
for triali=1:300
optim = CMAES_simple(4096,[],struct());
% initcodes = codes_all(1:30,:); %randn(30,4096);%
codes = initcodes;
% ctrl_codes_all = [];
for geni = 1:60
% ctrl_codes_all = [ctrl_codes_all; codes];
scores = randn(1,size(codes,1));
newcodes = optim.doScoring(codes,scores,1,struct());
codes = newcodes;
end
ctrlcodes_col = [ctrlcodes_col; mean(newcodes,1)];
end
save(fullfile(summarydir,"ctrl_code_col.mat"),'ctrlcodes_col')
%%
ctrlcodes_shf = cell2mat(arrayfun(@(idx)ctrlcodes_col(idx,randperm(4096)),1:300,'uni',0)');
ctrl_proj_cc = ctrlcodes_col * eigvect;
ctrl_perm_proj_cc = ctrlcodes_shf * eigvect;
ctrl_proj_cc = ctrl_proj_cc./norm_axis(ctrl_proj_cc,2).*300;
ctrl_perm_proj_cc = ctrl_perm_proj_cc./norm_axis(ctrl_perm_proj_cc,2).*300;
%%
rksm_col_ctrl = arrayfun(@(id)ranksum(abs(ctrl_proj_cc(:,id)),abs(perm_proj_cc(:,id))),1:4096,'uni',1);
krsm_col_ctrl = arrayfun(@(id)kstest2(abs(ctrl_proj_cc(:,id)),abs(perm_proj_cc(:,id)),'alpha',0.01),1:4096,'uni',1);
%%
[neurctrl_rksm_P, ~, neurctrl_rksm_Stat] = arrayfun(@(id)ranksum(abs(code_proj_cc(:,id)),abs(ctrl_proj_cc(:,id))),1:4096,'uni',1);
[~,neurctrl_krsm_P,neurctrl_krsm_Stat] = arrayfun(@(id)kstest2(abs(code_proj_cc(:,id)),abs(ctrl_proj_cc(:,id)),'alpha',0.01),1:4096,'uni',1);

%%
figure(9);clf;hold on
eigid = 5;
histogram(code_proj_cc(:,eigid),40,'FaceAlpha',0.3,'EdgeColor','none') % ,-20:2:30
histogram(ctrl_proj_cc(:,eigid),40,'FaceAlpha',0.5,'EdgeColor','none') % ,-20:2:30
histogram(perm_proj_cc(:,eigid),40,'FaceAlpha',0.5,'EdgeColor','none') % ,-20:2:30
% histogram(ctrl_perm_proj_cc(:,1),40,'FaceAlpha',0.5,'EdgeColor','none')
legend(["Neural Evol", "Noise Evol", "Neural Shuffled", "Noise Shuffled"])
title(compose("Hessian Eigenvector %d projection",eigid));xlabel("Projection Coefficient on Eigvect")
%%
ctrlcodes = (ctrl_codes_all(end-40+1:end, :)*eigvect);%,1);
cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((ctrlcodes>cutoff | ctrlcodes<-cutoff) & eigid>=eigthr-200+1 & eigid<=eigthr,2)),...
            prctile(abs(ctrlcodes(:)),[80,90,95,99,99.5,99.9,99.95,99.99],'all')')./min(eigthr,200),...[0,0.2,1,2,4,6,8,10,13,15,18]
        [25,50,100,200,400,800,1000,2000,4096]','Uni',0))

cell2mat(arrayfun(@(eigthr)...
            arrayfun(@(cutoff)mean(sum((ctrlcodes>cutoff | ctrlcodes<-cutoff) & eigid<=eigthr,2)),...
            prctile(abs(ctrlcodes(:)),[80,90,95,99,99.5,99.9,99.95,99.99],'all')')./eigthr,...[0,0.2,1,2,4,6,8,10,13,15,18]
        [25,50,100,200,400,600,800,1000,2000,4096]','Uni',0))

    
%%
Data = load("N:\Code\Generator_DB_Windows\init_population\init_code.mat");
initcodes = Data.codes;
%%
drft_codes_col = cell(numel(evo_meta),1);
for rowi = 1:numel(evo_meta)
for triali = 1:5
tic;
target_gene = evol_codes(rowi,randperm(4096));
target_norm = norm(target_gene);
genes_mean = mean(initcodes,1);
C =  CMAES_simple(4096,[],struct());

mean_activation = nan(1) ;
iGen = 1 ; 
genes = initcodes;
current_norm = mean(norm_axis(genes,2));
while current_norm < target_norm
    %         fprintf('%d ',iGen) ; if ~mod(iGen,30), fprintf('\n'); end
    % score is a function of the distance from current genes to target gene
    act_unit_distance = pdist2(target_gene,genes) ; 
    act_unit = 1./(act_unit_distance)*100 ; 
    mean_activation(iGen) = mean(act_unit) ;
    [genes,~] = C.doScoring(genes,act_unit,true);
    current_norm = mean( sqrt(sum( genes.^2 , 2 ) ) ) ;  
    mean_norm(iGen) = current_norm; 
    genes_mean = cat(1,genes_mean,mean(genes,1));
    iGen = iGen + 1 ; 
end % of iGen
drft_codes_col{rowi} = [drft_codes_col{rowi};genes_mean(end,:)];
toc
end
end
%%
save(fullfile(mat_dir,"evol_ctrl_codes.mat"),'drft_codes_col','evo_meta')
%%
drft_codes_gen_col = cell(numel(evo_meta),1);
for rowi = 1:numel(evo_meta)
for triali = 1:5
tic;
target_gene = evol_codes(rowi,randperm(4096));
genes_mean = mean(initcodes,1);
C =  CMAES_simple(4096,[],struct());
mean_activation = nan(1) ;
genes = initcodes;
for iGen = 1:evo_meta(rowi).nGen
    %         fprintf('%d ',iGen) ; if ~mod(iGen,30), fprintf('\n'); end
    % score is a function of the distance from current genes to target gene
    act_unit_distance = pdist2(target_gene,genes) ; 
    act_unit = 1./(act_unit_distance)*100 ; 
    mean_activation(iGen) = mean(act_unit) ;
    [genes,~] = C.doScoring(genes,act_unit,true);
    genes_mean = cat(1,genes_mean,mean(genes,1));
end % of iGen
drft_codes_gen_col{rowi} = [drft_codes_gen_col{rowi};genes_mean(end,:)];
toc
end
end
save(fullfile(mat_dir,"evol_ctrl_gen_codes.mat"),'drft_codes_gen_col','evo_meta')