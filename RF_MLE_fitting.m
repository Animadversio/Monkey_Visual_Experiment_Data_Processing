% function RF_View_fun()
Animal = "Alfa";
MatStats_path = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Manif_RFstats.mat", Animal)), 'RFStats')
% load(fullfile(MatStats_path, compose("%s_Manif_RFMaps.mat", Animal)), 'MaskStats')
% who('-file', fullfile(MatStats_path, compose("%s_Manif_RFMaps.mat", Animal)));
% end
%% Comparison of fitting and inference results for multiple methods
Expi = 2;
npos = size(RFStats(Expi).stim.uniqpos,1);
stim_info = [RFStats(Expi).stim.uniqpos, ones(npos,1) * RFStats(Expi).stim.imgsize_deg];
scores_mat = reshape(RFStats(Expi).psth.score_mean, size(RFStats(Expi).psth.score_mean,1), [])';
xpos = RFStats(Expi).stim.xpos; 
ypos = RFStats(Expi).stim.ypos;
uniqpos = RFStats(Expi).stim.uniqpos;
for iCh = 36:size(scores_mat,2)
% Log likelihood with regression to infer Amp and baseline (may converge faster)
[~, maxIdx] = max(scores_mat(:,iCh));
peakX = stim_info(maxIdx,1); peakY = stim_info(maxIdx,2);
initWid = (xpos(2)-xpos(1))*1.5;%(max(xpos) - min(xpos)) / 5;
% un bounded search can result in negative values
% [x, minL] = fminsearch(@(x) - logLikeReg(stim_info(:,1:2), stim_info(:,3), scores_mat(:,iCh), [x(1),x(2)], x(3)),...
%     [peakX, peakY, initWid],optimset("MaxFunEvals",5000));
LB=[min(xpos)-2,min(ypos)-2,0.01];
UB=[max(xpos)+2,max(ypos)+2, 6];
[x, minL] = fminsearchbnd(@(x) - logLikeReg(stim_info(:,1:2), stim_info(:,3), scores_mat(:,iCh), [x(1),x(2)], x(3)),...
    [peakX, peakY, initWid],LB,UB,optimset("MaxFunEvals",5000));
[x, predRsp] = RegressRsp(stim_info(:,1:2), stim_info(:,3), scores_mat(:,iCh), [x(1),x(2)], x(3));
disp(minL);disp(x)
% Use the MLE inferenced parameter to get the pdf
[XX,YY]=meshgrid(-10:0.1:10,-10:0.1:10);
estimRF = mvnpdf([XX(:),YY(:)],x(1:2),eye(2)*x(3)) * x(4) + x(5);
% Fit the response to a 2d Gaussian function
param = fmgaussfit(uniqpos(:,1),uniqpos(:,2),scores_mat(:,iCh));
% fitresult: 1 amplitude, 2 angle, 3 stdx, 4 stdy, 5 x0, 6 y0, 7 baseline(z0)
gaussfitRF = gaussian2D(param,XX,YY);
gaussfitRsp = gaussian2D(param,reshape(uniqpos(:,1), length(ypos),length(xpos)),...
                reshape(uniqpos(:,2), length(ypos),length(xpos)));
%% Plotting 
figure(12); % show inference result of RF
subplot(1,3,1)
imagesc(-10:0.1:10,-10:0.1:10,reshape(estimRF,size(XX)));axis image;set(gca,"Ydir",'normal')
title(compose("MLE RF center [%.1f,%.1f] std %.1f",x(1),x(2),x(3)))
subplot(1,3,2)
interpScore = griddata(uniqpos(:,1),uniqpos(:,2),scores_mat(:,iCh),XX,YY);
imagesc(-10:0.1:10,-10:0.1:10,interpScore);axis image;set(gca,"Ydir",'normal')
title("Interpolation of Response")
subplot(1,3,3);
imagesc(-10:0.1:10,-10:0.1:10,gaussfitRF);axis image;set(gca,"Ydir",'normal')
title("Gauss fit of RF")

figure(6);clf % show fitting result of response
subplot(1,3,1)
imagesc(xpos, ypos, reshape(scores_mat(:,iCh), length(ypos),length(xpos)))
CLIM = caxis();colorbar();axis image;set(gca,"Ydir",'normal')
title(compose("Mean Rsp, F:%.1f(%.E) T:%.1f(%.E)",RFStats(Expi).stats.F(iCh),RFStats(Expi).stats.anovaP(iCh),...
    RFStats(Expi).stats.T(iCh),RFStats(Expi).stats.ttestP(iCh)));
subplot(1,3,2)
imagesc(xpos, ypos, reshape(predRsp, length(ypos),length(xpos)))
caxis(CLIM);colorbar();axis image;set(gca,"Ydir",'normal')
title(compose("Fitted Rsp LLh %.1f nsd %.1f",-minL,x(end)))
subplot(1,3,3);
imagesc(xpos, ypos, gaussfitRsp);
caxis(CLIM);axis image;colorbar();set(gca,"Ydir",'normal')
title("Gauss fit of 2d data")
suptitle(compose("%s RFMap Exp %d, kernel fitting Ch %s", Animal, Expi, RFStats(Expi).unit.unit_name_arr(iCh)))
pause;
end
%%
function L = logLike(stimpos, stimsize, rsp, rfpos, rfstd, A, bsl, noisesd)
% Assume Gaussian kernel with square shape image
LB = ((stimpos - rfpos) - stimsize/2) / rfstd;
UB = ((stimpos - rfpos) + stimsize/2) / rfstd;
prob1d = abs(normcdf(LB) - normcdf(UB)); % compute the probability mass along 2 axis 
predrsp = prod(prob1d, 2) * A + bsl; % compute the 2d probability mass 
LLH = log(normpdf(rsp, predrsp, noisesd));
L = sum(LLH);
% sum((rsp - predrsp).^2 / noisesd.^2)
end

function L = logLikeReg(stimpos, stimsize, rsp, rfpos, rfstd)
% Assume Gaussian kernel with square shape image
LB = ((stimpos - rfpos) - stimsize/2) / rfstd;
UB = ((stimpos - rfpos) + stimsize/2) / rfstd;
prob1d = abs(normcdf(LB) - normcdf(UB));
prob = prod(prob1d, 2); % compute the 2d probability mass 
[coef,~,residue] = regress(rsp, [ones(length(prob),1),prob]);
A = coef(2) ;bsl = coef(1); % Inference baseline and amplitude
if A < 0, bsl = mean(rsp);residue = rsp - bsl; end % Inference baseline and amplitude
noisesd = std(residue);
% LLH = log(normpdf(residue / noisesd));
LLH = log(normpdf(residue, 0, noisesd));
L = sum(LLH);
% sum((rsp - predrsp).^2 / noisesd.^2)
end

function [param, predRsp] = RegressRsp(stimpos, stimsize, rsp, rfpos, rfstd)
LB = ((stimpos - rfpos) - stimsize/2) / rfstd;
UB = ((stimpos - rfpos) + stimsize/2) / rfstd;
prob1d = abs(normcdf(LB) - normcdf(UB));
prob = prod(prob1d, 2);
[coef,~,residue] = regress(rsp, [ones(length(prob),1),prob]);
A = coef(2); bsl = coef(1);
if A < 0, A=0; bsl = mean(rsp);residue = rsp - bsl; end
noisesd = std(residue); 
param = [rfpos, rfstd, A, bsl, noisesd];
predRsp = prob * A + bsl;
end

function z = gaussian2D(par,xx,yy)
% compute 2D gaussian
z = par(7) + ...
    par(1)*exp(-(((xx-par(5)).*cosd(par(2))+(yy-par(6)).*sind(par(2)))./par(3)).^2-...
    ((-(xx-par(5)).*sind(par(2))+(yy-par(6)).*cosd(par(2)))./par(4)).^2);
end