figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning\Schematics";
%%
[XX, YY]= meshgrid(-50:50,-50:50);
ZZ = exp(-0.005*(XX.^2 + YY.^2-XX.*YY));
figure(1);clf;hold on 
p = surf(XX,YY,ZZ,'FaceAlpha',0.4,'EdgeAlpha',0.1);
msk = randi(numel(XX),100,1);
Xsamp = XX(msk);
Ysamp = YY(msk);
datasamp = ZZ(msk) .* (1+0.2*randn(numel(msk),1));
scatter3(XX(msk),YY(msk),datasamp,'red','filled')
zlabel("Activation")
saveas(1,figdir+"\TuningPeak.png")
%%
XYdistmat = squareform(pdist([Xsamp,Ysamp]));
[cc_vec, p_vec] = corrcoef_vec(datasamp, XYdistmat);
[min_cc,min_idx] = min(cc_vec);

figure(2); 
subplot("Position",[0.05,0.05,0.87,0.9])
imagesc(XYdistmat); axis image; title("Stimuli Dissimilarity")
vline([min_idx-0.5,min_idx+0.5],{'r-.','r-.'})
subplot("Position",[0.93,0.05,0.05,0.9])
imagesc(datasamp); yticklabels([]);xticklabels([]); title("Response")
saveas(2,figdir+"\Distmat_Rsp.png")
%%
%%
figure(3);clf;hold on 
scatter(XYdistmat(:,min_idx), datasamp)
plotGPRfitting(XYdistmat(:,min_idx), datasamp)
xlabel("Distance to reference")
ylabel("Response")
title("Radial Tuning Curve")
saveas(3,figdir+"\Distmat_Rsp_Corr.png")
%%
function plotGPRfitting(X, Y, vargin)
if nargin == 2, vargin={};end
gprMdl = fitrgp(X, Y, 'SigmaLowerBound', 1E-3); % This lower bound is to fix a problem! 
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit, vargin{:})
end