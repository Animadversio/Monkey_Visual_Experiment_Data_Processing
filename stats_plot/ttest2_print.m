function [tval,pval,sumstr,mean_arr,sem_arr] = ttest2_print(group1,group2,g1label,g2label)
if nargin==2, g1label="Group 1";g2label="Group 2";end
g1_mean = mean(group1);
g2_mean = mean(group2);
mean_arr = [g1_mean,g2_mean];
sem_arr = [sem(group1),sem(group2)];
[~,P,CI,TST] = ttest2(group1,group2);
tval = TST.tstat;
pval = P;
sumstr = sprintf("%s (%.1f) - %s (%.1f): t=%.3f(df=%d), P=%.1e, CI=[%.1f,%.1f]\n",g1label,g1_mean,g2label,g2_mean,tval,TST.df,P,CI(1),CI(2));
fprintf(sumstr)
end