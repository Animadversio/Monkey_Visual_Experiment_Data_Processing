% Test the properties of CCA and OLS and Cov as weights.
X = rand(1000,100);
Ampl = 1+abs(randn(1,size(X,2)));%0.1+rand(1,size(X,2));
X = Ampl.*X;
% w = randn(size(X,2),1);
w = rand(size(X,2),1);
X = X - mean(X,1);
err = randn(size(X,1),1)*0.1;
y = X * w + err;
[CCAb,~,r,~,~,stats] = canoncorr(X,y);
[OLSb,OLSbintv,r] = regress(y,X);
CovMat = cov([y,X]);
cvyX = CovMat(2:end, 1);
cvyX = cvyX * std(y)/std(X*cvyX);
CorrMat = corrcoef([y,X]);
crryX = CorrMat(2:end, 1);
crryX = crryX * std(y)/std(X*crryX);
fprintf("\n")
fprintf("Xy Correlation vs True W corr: %.3f\n",corr(crryX, w))
fprintf("Xy Covarince vs True W corr: %.3f\n",corr(cvyX, w))
fprintf("OLS est vs True W corr: %.3f\n",corr(OLSb, w))
fprintf("CCA est vs True W corr: %.3f\n",corr(CCAb, w))
fprintf("CCA est vs OLS est corr: %.3f\n",corr(OLSb, CCAb))
fprintf("Reconstr with CovXy vs y corr: %.3f\n",corr(X*cvyX, y))
fprintf("Reconstr with CorrXy vs y corr: %.3f\n",corr(X*crryX, y))
fprintf("Reconstr with OLS vs y corr: %.3f\n",corr(X*OLSb, y))
fprintf("Reconstr with CCA vs y corr: %.3f\n",corr(X*CCAb, y))
fprintf("Reconstr with True w vs y corr: %.3f\n",corr(X*w, y))
%%
X = rand(100,10180);
w = rand(size(X,2),1);
X = X - mean(X,1);
err = randn(size(X,1),1)*0.1;
y = X * w + err;
[CCAb,~,r,~,~,stats] = canoncorr(X,y);
ridgeB = ridge(y,X,0.01);
OLSb = ridgeB;
CovMat = cov([y,X]);
cvyX = CovMat(2:end, 1);
CorrMat = corrcoef([y,X]);
crryX = CorrMat(2:end, 1);
fprintf("\nProblem samp N=%d; feat p=%d\n",size(X,2),size(X,1))
fprintf("Xy Correlation vs True W corr: %.3f\n",corr(crryX, w))
fprintf("Xy Covarince vs True W corr: %.3f\n",corr(cvyX, w))
fprintf("Xy Covarince vs Xy Correlation corr: %.3f\n",corr(cvyX, crryX))
fprintf("Ridge est vs True W corr: %.3f\n",corr(OLSb, w))
fprintf("CCA est vs True W corr: %.3f\n",corr(CCAb, w))
fprintf("CCA est vs Ridge est corr: %.3f\n",corr(OLSb, CCAb))
fprintf("Reconstr with CovXy vs y corr: %.3f\n",corr(X*cvyX, y))
fprintf("Reconstr with CorrXy vs y corr: %.3f\n",corr(X*crryX, y))
fprintf("Reconstr with Ridge vs y corr: %.3f\n",corr(X*OLSb, y))
fprintf("Reconstr with CCA vs y corr: %.3f\n",corr(X*CCAb, y))
fprintf("Reconstr with True w vs y corr: %.3f\n",corr(X*w, y))
%%
thresh = 3 ;
yXtval = r2tval(crryX,size(X,1));
sum(r2tval(crryX,size(X,1))>thresh);
feat_msk = abs(yXtval)>thresh;
pred = X(:,feat_msk)*cvyX(feat_msk,:);
fprintf("Reconstr with feature selected correlation vs y corr: %.3f\n",corr(pred, y))
function t = r2tval(r,N)
t = sqrt(N-2)*r./sqrt(1-r.^2);
end