function [tval,pval,sumstr,mean_arr,sem_arr] = ttest_print(group1,g1label)
% Signature:
%   ttest_print(group1,g1label)
%    
% Examples:
%    ttest_print( val_arr );
if nargin==1, g1label="Group 1"; end
g1_mean = mean(group1);
mean_arr = [g1_mean];
sem_arr = [sem(group1)];
[~,P,CI,TST] = ttest(group1);
tval = TST.tstat; pval = P;
sumstr = sprintf("%s (%.1f) - 0: t=%.3f(df=%d), P=%.1e, CI=[%.1f,%.1f]\n",g1label,g1_mean,tval,TST.df,pval,CI(1),CI(2));
fprintf(sumstr)
end