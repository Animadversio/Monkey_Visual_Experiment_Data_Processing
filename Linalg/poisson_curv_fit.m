%% Poisson Likelihood Curve Fit a demo to fit discrete output.


real_score = predStats(11).score_Evol(:,end); % real value firing rate
intscore = round(real_score * 0.15); % get the discrete spike number here.
lfit_score = predStats(11).E2M.lfitscore(:,end); % the input x to be transformed
%%
BReLUfun = @(x,param) max(0, param(1).* (x - param(2))) + param(3);
[param,fval] = fminsearchbnd(@(param) sum( - intscore.*log(BReLUfun(lfit_score,param)) + BReLUfun(lfit_score,param)),...
    double([1,min(lfit_score),min(intscore)+0.1]),[0,min(lfit_score),1E-5],[10000,max(lfit_score),max(intscore)*2])
disp(corr(intscore, BReLUfun(lfit_score, param)))
%% Iterative Gradient Descent
param0 = [1,min(lfit_score),min(intscore)+0.1];
lr = 0.01;
param = param0;
for i = 1:5000
activ_msk = lfit_score > param(2);
partigrad = [activ_msk.*lfit_score, - lfit_score * param(1), ones(length(lfit_score),1)];
grad = mean((intscore ./ BReLUfun(lfit_score, param) - 1) .* partigrad,1);
param = param + lr * grad;
end
disp(param) % parameter
disp(sum( intscore.*log(BReLUfun(lfit_score,param)) - BReLUfun(lfit_score,param)))% likelihood function
disp(corr(intscore, BReLUfun(lfit_score, param)))% correlation