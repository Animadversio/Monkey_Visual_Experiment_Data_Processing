%% FC6GAN spectrum if linearized
%  This code demonstrate there is cut off effect to spectrum as we go
%  through the pretrained GAN. Thus components in many directions are
%  suppressed by the linear mapping layers thus not affect the subsequent
%  layers or image. 
%  As a result the GANs will have quite a few directions with small Hessian
%  eigenvalue. 
%% Extract the layers from GAN 
Wdefc7 = G.LinNet.Layers(2).Weights; % out_feat by in_feat
Wdefc6 = G.LinNet.Layers(4).Weights;
Wdefc5 = G.LinNet.Layers(6).Weights;
Bdefc7 = G.LinNet.Layers(2).Bias; % out_feat by in_feat
Bdefc6 = G.LinNet.Layers(4).Bias;
Bdefc5 = G.LinNet.Layers(6).Bias;
%% SVD to the weight matrices
tic
[U7,S7,D7] = svd(Wdefc7); S7=diag(S7);
[U67,S67,D67] = svd(Wdefc6 * Wdefc7); S67=diag(S67);
[U567,S567,D567] = svd(Wdefc5 * Wdefc6 * Wdefc7); S567=diag(S567); % the first 3 layers have no negative spectrum 
toc
%% Plot the Spectrum 
savedir = "E:\OneDrive - Washington University in St. Louis\HessTune\Spectrum";
figure(7);clf
subplot(121);hold on 
plot(S7,'LineWidth',1.5);
plot(S67,'LineWidth',1.5);
plot(S567,'LineWidth',1.5);
title("raw singular value")
legend(["defc7", "defc6,7", "defc5,6,7"])
subplot(122);hold on 
plot(S7 ./ prctile(S7,95),'LineWidth',1.5); 
plot(S67 ./ prctile(S67,95),'LineWidth',1.5);
plot(S567 ./ prctile(S567,95),'LineWidth',1.5);
title("singular value normalized to 95% percentile")
legend(["defc7", "defc6,7", "defc5,6,7"])
suptitle("SVD analysis for first 3 fc Layer of FC6GAN")
savefig(7,fullfile(savedir,"SV_fc_layers.fig"))
saveas(7,fullfile(savedir,"SV_fc_layers.jpg"))
%%
savedir = "E:\OneDrive - Washington University in St. Louis\HessTune\Spectrum";
figure(8);clf;set(8,'position',[1000         145         975         833])
corrplot(log10([S7,S67,S567]),'varNames',["fc7","fc67","fc567"]);
annotation(figure(8),'textbox', [0.103 0.965 0.785 0.04],...
    'String',{'Correlation of Singlular Value of first 3 fc Layer of FC6GAN','Acting on Right Eigenvectors of each Operator'},...
    'FontSize',18, 'FitBoxToText','off', 'EdgeColor','none');
savefig(8,fullfile(savedir,"SV_fc_layers_corr.fig"))
saveas(8,fullfile(savedir,"SV_fc_layers_corr.jpg"))
%%
S7a = norm_axis(Wdefc7 * D7,1)';
S67a = norm_axis(Wdefc6 * Wdefc7 * D7,1)';
S567a = norm_axis(Wdefc5 * Wdefc6 * Wdefc7 * D7,1)';
% S6a = norm_axis(Wdefc6 * D7,1)'; % this doesn't work
% S5a = norm_axis(Wdefc5 * D7,1)'; % this doesn't work, since it's a diff space
figure(5);clf
corrplot(log10([S7a,S67a,S567a]),'varNames',["fc7","fc67","fc567"]);%,S6a,S5a,"fc6","fc5"
annotation(figure(5),'textbox', [0.103 0.965 0.785 0.04],...
    'String',{'Correlation of Amplification Factor of first 3 fc Layer of FC6GAN','Acting on Right Eigenvectors of defc7'},...
    'FontSize',18, 'FitBoxToText','off', 'EdgeColor','none');
savefig(5,fullfile(savedir,"SVact_fc_layers_corr.fig"))
saveas(5,fullfile(savedir,"SVact_fc_layers_corr.jpg"))
%% Load the average hessian Matrix from python calculation. 
py.importlib.import_module("numpy");
out_dir = "E:\OneDrive - Washington University in St. Louis\HessTune\NullSpace";
data = py.numpy.load(fullfile(out_dir, "Evolution_Avg_Hess.npz"));
H_evo_avg = data.get("H_avg").double;
eigval_evo_avg = data.get("eigv_avg").double;
eigvec_evo_avg = data.get("eigvect_avg").double;
%%
savedir = "E:\OneDrive - Washington University in St. Louis\HessTune\Spectrum";
%% Confirm that the Hessian eigenvalues are correlated with the action of first 3 layers on it 
linact_eigvec = Wdefc5 * Wdefc6 * Wdefc7 * eigvec_evo_avg;
linact_ampl = norm_axis(linact_eigvec,1);
figure(3);
subplot(121)
scatter(eigval_evo_avg,linact_ampl)
ylabel("Amplification of Linear Operator")
xlabel("Eigenvalue to the Image Difference Hessian")
subplot(122)
scatter(log10(eigval_evo_avg),log10(linact_ampl))
ylabel("Amplification of Linear Operator (log)")
xlabel("Eigenvalue to the Image Difference Hessian (log)")
suptitle(["Hessian EigVal to Image Difference Metric Correlates with Linear Amplification of First 3 Layers on Eigenvectors"])
saveas(3, fullfile(savedir, "spect_linact_align.jpg"))
%% nonlinear action with relu
nlnact_eigvec = max(0, Wdefc5 * max(0, Wdefc6 * max(0, Wdefc7 * eigvec_evo_avg)));% + Bdefc7) + Bdefc6) + Bdefc5);
nlnact_ampl = norm_axis(nlnact_eigvec,1);
figure(6);
subplot(121)
scatter(eigval_evo_avg,nlnact_ampl)
ylabel("Amplification of NonLinear Operator")
xlabel("Eigenvalue to the Image Difference Hessian")
subplot(122)
scatter(log10(eigval_evo_avg),log10(nlnact_ampl))
ylabel("Amplification of NonLinear Operator (log)")
xlabel("Eigenvalue to the Image Difference Hessian (log)")
suptitle(["Hessian EigVal to Image Difference Metric Correlates with Amplification of First 3 Layers on Eigenvectors","with Relu NonLinear"])
saveas(6, fullfile(savedir, "spect_reluact_align.jpg"))
%% nonlinear action with relu and biases (need to scale input to suit bias)
scale = 300;
nlnact_eigvec = max(0, Wdefc5 * max(0, Wdefc6 * max(0, Wdefc7 * scale * eigvec_evo_avg + Bdefc7) + Bdefc6) + Bdefc5);
nlnact_ampl = norm_axis(nlnact_eigvec - Bdefc5, 1) / scale;
figure(4);
subplot(121)
scatter(eigval_evo_avg,nlnact_ampl)
ylabel("Amplification of NonLinear Operator")
xlabel("Eigenvalue to the Image Difference Hessian")
subplot(122)
scatter(log10(eigval_evo_avg),log10(nlnact_ampl))
ylabel("Amplification of NonLinear Operator (log)")
xlabel("Eigenvalue to the Image Difference Hessian (log)")
suptitle(["Hessian EigVal to Image Difference Metric Correlates with NonLinear Amplification of First 3 Layers on Eigenvectors",compose("with Relu and Bias scale %d",scale)])
saveas(4, fullfile(savedir, "spect_reluBact300_align.jpg"))