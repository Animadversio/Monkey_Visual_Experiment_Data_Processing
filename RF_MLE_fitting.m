%% RF MLE fitting
global img_maskstack gridX gridY
ntick = 201;
visualField = [-10 10]; 
coli = linspace(visualField(1),visualField(2),ntick);
rowi = linspace(visualField(1),visualField(2),ntick);
[gridX,gridY]  = meshgrid(coli,rowi);
%% infer the RF parameters through MLE fitting....
A=1; rfX=0; rfY=0; sigma=1;
RFmask = exp(- ( (gridX - rfX).^2 + (gridY - rfY).^2 ) / sigma^2);
convRF = A * squeeze(sum(img_maskstack .* RFmask,[1,2]));%einsum(img_maskstack, RFmask, 'ijk,ij->k');
%%
% param = fit(x,y,fitType); % creates the fit to the data in x and y with the model specified by fitType.
% loss = @(A, rfX, rfY, sigma) sum((score_mat(:,end) - RFpred(A, rfX, rfY, sigma)).^2);
target = score_mat(:,76);
loss = @(x) sum((target - RFpred(x(1), x(2), x(3), x(4))).^2);
[x,fval,exitflag,output] = fminsearch(loss, [1,0,0,1])
%%
predRsp = RFpred(x(1), x(2), x(3), x(4));
%%
function convRF = RFpred(A, rfX, rfY, sigma)
global img_maskstack gridX gridY
RFmask = exp(- ( (gridX - rfX).^2 + (gridY - rfY).^2 ) / sigma^2);
convRF = A * squeeze(sum(img_maskstack .* RFmask,[1,2]));%einsum(img_maskstack, RFmask, 'ijk,ij->k');
end
