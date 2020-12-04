%% Smoothness of Manifold Tuning
%  This script is written to measure the smoothness of a 1d or 2d tuning
%  function. 
%  Note the sampling and neuronal noise will affect this roughness, thus
%  during testing, we need to generate hypothesis taking this sampling noise 
%  into account, based on this roughness measurement. 

Animal = "Beto";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(mat_dir,"gab_imdist.mat"),'gab_imdist')
load(fullfile(mat_dir,"pasu_imdist.mat"),'pasu_imdist')
load(fullfile(mat_dir,Animal+"_Manif_ImDist.mat"),"ManifImDistStat")
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
%% Load the map 
% tunemap_dirichilet_energy
Expi=10;si=1; ui=1;
unitlab = Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui));
actmap = cellfun(@(P)mean(P(ui,51:200,:),'all'),Stats(Expi).manif.psth{si});
act_col = cellfun(@(P)squeeze(mean(P(ui,51:200,:),2)),Stats(Expi).manif.psth{si},'uni',0);
cntmap = cellfun(@(P)size(P,3),Stats(Expi).manif.psth{si});
%%
kerd2 = [[0,-1,0];[-1,4,-1];[0,-1,0]]/4;
laplsMap = conv2(actmap, kerd2);
[actGx,actGy] = gradient(actmap);
gradAmpMap = actGx.^2 + actGy.^2;
dirEng = sum(gradAmpMap,'all');
%%
[PHI, THETA] = meshgrid(-90:18:90,-90:18:90); % PHI is the angle towards PC3, THETA is the angle towards PC2 
%%
figure(3);
subplot(131)
imagesc(actmap)
axis image;colorbar
title("Neural Activation")
subplot(132)
imagesc(sqrt(gradAmpMap))
axis image;colorbar
title("Gradient Amplitude")
subplot(133)
imagesc(laplsMap)
axis image;colorbar
title("Laplacian Filtered Map")
%%
actmap_shfl = reshape(actmap(randperm(numel(actmap))),size(actmap));
%%
actrng = range(actmap,'all');
[D2E,TVE,D1E,D2map,GPmap] = map_energy(actmap);
[D2E_shfl,TVE_shfl,D1E_shfl,D2map_shfl,GPmap_shfl] = map_energy(actmap_shfl);
%%
btrp_Scol = [];
for rep = 1:1000
actmap_btrp = cellfun(@(A)mean(A(datasample(1:numel(A),numel(A)))),act_col);
[D2E_btrp,TVE_btrp,D1E_btrp,D2map_btrp,GPmap_btrp] = map_energy(actmap_btrp);
btrp_Scol = [btrp_Scol;[D2E_btrp,TVE_btrp,D1E_btrp]];
end
shfl_Scol = [];
for rep = 1:1000
actmap_shfl = reshape(actmap(randperm(numel(actmap))),size(actmap));
[D2E_shfl,TVE_shfl,D1E_shfl,D2map_shfl,GPmap_shfl] = map_energy(actmap_shfl);
shfl_Scol = [shfl_Scol;[D2E_shfl,TVE_shfl,D1E_shfl]];
end
%%
h=figure(4);clf;
subplot(131);hold on
histogram(btrp_Scol(:,1),'Norm','pdf');
histogram(shfl_Scol(:,1),'Norm','pdf');
title("Laplacian Filtered Energy")
subplot(132);hold on
histogram(btrp_Scol(:,2),'Norm','pdf');
histogram(shfl_Scol(:,2),'Norm','pdf');
title("Total Variation Energy")
subplot(133);hold on
histogram(btrp_Scol(:,3),'Norm','pdf');
histogram(shfl_Scol(:,3),'Norm','pdf');
title("Dirichlet Energy")
legend(["Trial Resampled", "Matrix Shuffled"])
title(h,"Smooth Energy Compared")
%%
%%
[PHI, THETA] = meshgrid(-90:18:90,-90:18:90); 
detJac = abs(cosd(PHI));
Ginv11 = 1./cosd(PHI).^2;Ginv11(:,1)=0;Ginv11(:,end)=0;
Ginv22 = ones(size(PHI));
% Sample points for the gradients are far away.
G_W_map = cosd([mean(PHI(:,[1,2]),2),PHI(:,2:end-1),mean(PHI(:,[end-1,end]),2)]);
actmap_sp = actmap;
actmap_sp(:,1)=mean(actmap(:,1));
actmap_sp(:,end)=mean(actmap(:,end));
kerd2 = [[0,-1,0];[-1,4,-1];[0,-1,0]]/4;
laplsMap_sp = del2(actmap_sp);%conv2(actmap_sp, kerd2);
[actGx_sp,actGy_sp] = gradient(actmap_sp);
gradNormMap_sp = actGx_sp.^2 .* Ginv22 + actGy_sp.^2 .*Ginv11;
dirEng_sp = sum(gradNormMap_sp.*G_w_map,'all');
%%
figure(2);
subplot(141)
imagesc(actmap)
axis image;colorbar
title("Neural Activation")
subplot(142)
imagesc(gradNormMap_sp)
axis image;colorbar
title("Gradient Power")
subplot(143)
imagesc(gradNormMap_sp.*G_w_map)
axis image;colorbar
title("Gradient Power * Grid Weights")
subplot(144)
imagesc(laplsMap)
axis image;colorbar
title("Laplacian Filtered Map")
%%
figdir = "O:\Manif_MapSmooth";
spacelabels = ["PC23","PC4950","RND12"];
Expi=10;si=1; ui=1;
SmthStatTab = repmat(struct(),1,0);cnt=1;
h=figure(5);h2=figure(6);h3=figure(7);h4=figure(8);
set(5,'pos',[275         423        1616         437]);
set(7,'pos',[275         423        1616         437]);
set(6,'pos',[386         135        1703         766]);
set(8,'pos',[386         135        1703         766]);
for Expi = 1:numel(Stats)
for si = 1:numel(Stats(Expi).manif.psth)
splab = spacelabels(si);
for ui = 1:numel(Stats(Expi).units.pref_chan_id)
unitlab = Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui));
SmthStatTab(cnt).Expi = Expi;
SmthStatTab(cnt).chan = Stats(Expi).units.spikeID(Stats(Expi).units.pref_chan_id(ui));
SmthStatTab(cnt).unitlab = unitlab;
SmthStatTab(cnt).unitnum = ui;
SmthStatTab(cnt).spaceId = si; 
% Get the activities 
actmap = cellfun(@(P)mean(P(ui,51:200,:),'all'),Stats(Expi).manif.psth{si});
act_col = cellfun(@(P)squeeze(mean(P(ui,51:200,:),2)),Stats(Expi).manif.psth{si},'uni',0);
cntmap = cellfun(@(P)size(P,3),Stats(Expi).manif.psth{si});
% Spherical Version
[D2E_sp, TVE_sp, D1E_sp, D2map_sp, GPmap_sp, detJac_map] = sph_map_energy(actmap);
btrp_Scol_sp = [];
for rep = 1:1000
actmap_btrp = cellfun(@(A)mean(A(datasample(1:numel(A),numel(A)))),act_col); % Resample trial responses
[D2E_sp_btrp, TVE_sp_btrp, D1E_sp_btrp, D2map_sp_btrp, GPmap_sp_btrp, detJac_map] = sph_map_energy(actmap_btrp);
btrp_Scol_sp = [btrp_Scol_sp;[D2E_sp_btrp, TVE_sp_btrp, D1E_sp_btrp]];
end
shfl_Scol_sp = [];
for rep = 1:1000
actmap_shfl = reshape(actmap(randperm(numel(actmap))),size(actmap));
[D2E_sp_shfl, TVE_sp_shfl, D1E_sp_shfl, D2map_sp_shfl, GPmap_sp_shfl, detJac_map] = sph_map_energy(actmap_shfl);
shfl_Scol_sp = [shfl_Scol_sp;[D2E_sp_shfl, TVE_sp_shfl, D1E_sp_shfl]];
end
[~,P_D2E_sp,~,tSt_D2E_sp] = ttest2(btrp_Scol_sp(:,1), shfl_Scol_sp(:,1));
Dpr_D2E_sp = computeCohen_d(btrp_Scol_sp(:,1), shfl_Scol_sp(:,1));
[~,P_TVE_sp,~,tSt_TVE_sp] = ttest2(btrp_Scol_sp(:,2), shfl_Scol_sp(:,2));
Dpr_TVE_sp = computeCohen_d(btrp_Scol_sp(:,2), shfl_Scol_sp(:,2));
[~,P_D1E_sp,~,tSt_D1E_sp] = ttest2(btrp_Scol_sp(:,3), shfl_Scol_sp(:,3));
Dpr_D1E_sp = computeCohen_d(btrp_Scol_sp(:,3), shfl_Scol_sp(:,3));
SmthStatTab(cnt).D2E_sp = D2E_sp; 
SmthStatTab(cnt).TVE_sp = TVE_sp; 
SmthStatTab(cnt).D1E_sp = D1E_sp; 
SmthStatTab(cnt).D2map_sp = D2map_sp; 
SmthStatTab(cnt).GPmap_sp = GPmap_sp; 
SmthStatTab(cnt).D2E_sp_btrp_mean = mean(btrp_Scol_sp(:,1));
SmthStatTab(cnt).D2E_sp_shfl_mean = mean(shfl_Scol_sp(:,1));
SmthStatTab(cnt).D2E_sp_btrp_std = std(btrp_Scol_sp(:,1));
SmthStatTab(cnt).D2E_sp_shfl_std = std(shfl_Scol_sp(:,1));
SmthStatTab(cnt).TVE_sp_btrp_mean = mean(btrp_Scol_sp(:,2));
SmthStatTab(cnt).TVE_sp_shfl_mean = mean(shfl_Scol_sp(:,2));
SmthStatTab(cnt).TVE_sp_btrp_std = std(btrp_Scol_sp(:,2));
SmthStatTab(cnt).TVE_sp_shfl_std = std(shfl_Scol_sp(:,2));
SmthStatTab(cnt).D1E_sp_btrp_mean = mean(btrp_Scol_sp(:,3));
SmthStatTab(cnt).D1E_sp_shfl_mean = mean(shfl_Scol_sp(:,3));
SmthStatTab(cnt).D1E_sp_btrp_std = std(btrp_Scol_sp(:,3));
SmthStatTab(cnt).D1E_sp_shfl_std = std(shfl_Scol_sp(:,3));
SmthStatTab(cnt).D2E_sp_t = tSt_D2E_sp.tstat;
SmthStatTab(cnt).D2E_sp_P = P_D2E_sp;
SmthStatTab(cnt).D2E_sp_Dpr = Dpr_D2E_sp;
SmthStatTab(cnt).TVE_sp_t = tSt_TVE_sp.tstat;
SmthStatTab(cnt).TVE_sp_P = P_TVE_sp;
SmthStatTab(cnt).TVE_sp_Dpr = Dpr_TVE_sp;
SmthStatTab(cnt).D1E_sp_t = tSt_D1E_sp.tstat;
SmthStatTab(cnt).D1E_sp_P = P_D1E_sp;
SmthStatTab(cnt).D1E_sp_Dpr = Dpr_D1E_sp;
% Euclidean Version
[D2E, TVE, D1E, D2map, GPmap] = map_energy(actmap);
btrp_Scol = [];
for rep = 1:1000
actmap_btrp = cellfun(@(A)mean(A(datasample(1:numel(A),numel(A)))),act_col); % Resample trial responses
[D2E_btrp, TVE_btrp, D1E_btrp, D2map_btrp, GPmap_btrp] = map_energy(actmap_btrp);
btrp_Scol = [btrp_Scol;[D2E_btrp, TVE_btrp, D1E_btrp]];
end
shfl_Scol = [];
for rep = 1:1000
actmap_shfl = reshape(actmap(randperm(numel(actmap))),size(actmap));
[D2E_shfl, TVE_shfl, D1E_shfl, D2map_shfl, GPmap_shfl] = map_energy(actmap_shfl);
shfl_Scol = [shfl_Scol;[D2E_shfl, TVE_shfl, D1E_shfl]];
end
[~,P_D2E,~,tSt_D2E] = ttest2(btrp_Scol(:,1), shfl_Scol(:,1));
Dpr_D2E = computeCohen_d(btrp_Scol(:,1), shfl_Scol(:,1));
[~,P_TVE,~,tSt_TVE] = ttest2(btrp_Scol(:,2), shfl_Scol(:,2));
Dpr_TVE = computeCohen_d(btrp_Scol(:,2), shfl_Scol(:,2));
[~,P_D1E,~,tSt_D1E] = ttest2(btrp_Scol(:,3), shfl_Scol(:,3));
Dpr_D1E = computeCohen_d(btrp_Scol(:,3), shfl_Scol(:,3));
SmthStatTab(cnt).D2E = D2E; 
SmthStatTab(cnt).TVE = TVE; 
SmthStatTab(cnt).D1E = D1E; 
SmthStatTab(cnt).D2map = D2map; 
SmthStatTab(cnt).GPmap = GPmap; 
SmthStatTab(cnt).D2E_btrp_mean = mean(btrp_Scol(:,1));
SmthStatTab(cnt).D2E_shfl_mean = mean(shfl_Scol(:,1));
SmthStatTab(cnt).D2E_btrp_std = std(btrp_Scol(:,1));
SmthStatTab(cnt).D2E_shfl_std = std(shfl_Scol(:,1));
SmthStatTab(cnt).TVE_btrp_mean = mean(btrp_Scol(:,2));
SmthStatTab(cnt).TVE_shfl_mean = mean(shfl_Scol(:,2));
SmthStatTab(cnt).TVE_btrp_std = std(btrp_Scol(:,2));
SmthStatTab(cnt).TVE_shfl_std = std(shfl_Scol(:,2));
SmthStatTab(cnt).D1E_btrp_mean = mean(btrp_Scol(:,3));
SmthStatTab(cnt).D1E_shfl_mean = mean(shfl_Scol(:,3));
SmthStatTab(cnt).D1E_btrp_std = std(btrp_Scol(:,3));
SmthStatTab(cnt).D1E_shfl_std = std(shfl_Scol(:,3));
SmthStatTab(cnt).D2E_t = tSt_D2E.tstat;
SmthStatTab(cnt).D2E_P = P_D2E;
SmthStatTab(cnt).D2E_Dpr = Dpr_D2E;
SmthStatTab(cnt).TVE_t = tSt_TVE.tstat;
SmthStatTab(cnt).TVE_P = P_TVE;
SmthStatTab(cnt).TVE_Dpr = Dpr_TVE;
SmthStatTab(cnt).D1E_t = tSt_D1E.tstat;
SmthStatTab(cnt).D1E_P = P_D1E;
SmthStatTab(cnt).D1E_Dpr = Dpr_D1E;
cnt = cnt+1;
%%
set(0,'CurrentFigure',h);clf;
T = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile(T,1);hold on
histogram(btrp_Scol_sp(:,1),'Norm','pdf');
histogram(shfl_Scol_sp(:,1),'Norm','pdf');
vline(D2E_sp)
title(compose("Laplacian Filtered Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D2E_sp.tstat,P_D2E_sp,Dpr_D2E_sp))
nexttile(T,2);hold on
histogram(btrp_Scol_sp(:,2),'Norm','pdf');
histogram(shfl_Scol_sp(:,2),'Norm','pdf');
vline(TVE_sp)
title(compose("Total Variation Energy\n t=%.2f(%.1e) d'=%.2f",tSt_TVE_sp.tstat,P_TVE_sp,Dpr_TVE_sp))
nexttile(T,3);hold on
histogram(btrp_Scol_sp(:,3),'Norm','pdf');
histogram(shfl_Scol_sp(:,3),'Norm','pdf');
vline(D1E_sp)
title(compose("Dirichlet Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D1E_sp.tstat,P_D1E_sp,Dpr_D1E_sp))
legend(["Trial Resampled", "Matrix Shuffled"])
title(T,["Map Smoothness Energy (Sphere) Compared",compose("Exp %d Space%d(%s) Unit %s",Expi,si,splab,unitlab)])
saveas(h, fullfile(figdir, compose("Smooth_DE_sph_hist_cmp_%s_Exp%d_%s_%s.png",Animal,Expi,splab,unitlab)))
%%
set(0,'CurrentFigure',h2);clf;
T = tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
nexttile(T,1);
imagesc(actmap)
axis image;colorbar
title("Neural Activation")
ylabel("Trial Resampling")
xlabel("PC3 Phi")
nexttile(T,2);
imagesc(GPmap_sp_btrp)
axis image;colorbar
title("Gradient Power")
xlabel("PC3 Phi")
nexttile(T,3);
imagesc(GPmap_sp_btrp.*detJac_map)
axis image;colorbar
title(["Gradient Power * Grid Weights",compose("Dirichlet Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D1E_sp.tstat,P_D1E_sp,Dpr_D1E_sp)])
xlabel("PC3 Phi")
nexttile(T,4);
imagesc(abs(D2map_sp_btrp))
axis image;colorbar
title(["Laplacian Filtered Amplitude",compose("Laplacian Filtered Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D2E_sp.tstat,P_D2E_sp,Dpr_D2E_sp)])
xlabel("PC3 Phi")

nexttile(T,5);
imagesc(actmap_shfl)
axis image;colorbar
title("Neural Activation")
ylabel("Shuffled Control")
xlabel("PC3 Phi")
nexttile(T,6);
imagesc(GPmap_sp_shfl)
axis image;colorbar
title("Gradient Power")
xlabel("PC3 Phi")
nexttile(T,7);
imagesc(GPmap_sp_shfl.*detJac_map)
axis image;colorbar
title("Gradient Power * Grid Weights")
xlabel("PC3 Phi")
nexttile(T,8);
imagesc(abs(D2map_sp_shfl))
axis image;colorbar
title("Laplacian Filtered Amplitude")
xlabel("PC3 Phi")
title(T,compose("Filtered map Compared (Sphere)\nExp %d Space%d(%s) Unit %s",Expi,si,splab,unitlab))
saveas(h2, fullfile(figdir, compose("Smooth_Map_sph_cmp_%s_Exp%d_%s_%s.png",Animal,Expi,unitlab,splab)))


set(0,'CurrentFigure',h3);clf;
T = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile(T,1);hold on
histogram(btrp_Scol(:,1),'Norm','pdf');
histogram(shfl_Scol(:,1),'Norm','pdf');
vline(D2E)
title(compose("Laplacian Filtered Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D2E.tstat,P_D2E,Dpr_D2E))
nexttile(T,2);hold on
histogram(btrp_Scol(:,2),'Norm','pdf');
histogram(shfl_Scol(:,2),'Norm','pdf');
vline(TVE)
title(compose("Total Variation Energy\n t=%.2f(%.1e) d'=%.2f",tSt_TVE.tstat,P_TVE,Dpr_TVE))
nexttile(T,3);hold on
histogram(btrp_Scol(:,3),'Norm','pdf');
histogram(shfl_Scol(:,3),'Norm','pdf');
vline(D1E)
title(compose("Dirichlet Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D1E.tstat,P_D1E,Dpr_D1E))
legend(["Trial Resampled", "Matrix Shuffled"])
title(T,["Map Smoothness Energy Compared",compose("Exp %d Space%d(%s) Unit %s",Expi,si,splab,unitlab)])
saveas(h3, fullfile(figdir, compose("Smooth_DE_hist_cmp_%s_Exp%d_%s_%s.png",Animal,Expi,splab,unitlab)))
%%
set(0,'CurrentFigure',h4);clf;
T = tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
nexttile(T,1);
imagesc(actmap)
axis image;colorbar
title("Neural Activation")
ylabel("Trial Resampling")
xlabel("PC3 Phi")
nexttile(T,2);
imagesc(GPmap_btrp)
axis image;colorbar
title("Gradient Power")
xlabel("PC3 Phi")
nexttile(T,3);
imagesc(GPmap_btrp.*detJac_map)
axis image;colorbar
title(["Gradient Power * Grid Weights",compose("Dirichlet Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D1E.tstat,P_D1E,Dpr_D1E)])
xlabel("PC3 Phi")
nexttile(T,4);
imagesc(abs(D2map_btrp))
axis image;colorbar
title(["Laplacian Filtered Amplitude",compose("Laplacian Filtered Energy\n t=%.2f(%.1e) d'=%.2f",tSt_D2E.tstat,P_D2E,Dpr_D2E)])
xlabel("PC3 Phi")

nexttile(T,5);
imagesc(actmap_shfl)
axis image;colorbar
title("Neural Activation")
ylabel("Shuffled Control")
xlabel("PC3 Phi")
nexttile(T,6);
imagesc(GPmap_shfl)
axis image;colorbar
title("Gradient Power")
xlabel("PC3 Phi")
nexttile(T,7);
imagesc(GPmap_shfl.*detJac_map)
axis image;colorbar
title("Gradient Power * Grid Weights")
xlabel("PC3 Phi")
nexttile(T,8);
imagesc(abs(D2map_shfl))
axis image;colorbar
title("Laplacian Filtered Amplitude")
xlabel("PC3 Phi")
title(T,compose("Exp %d Space%d(%s) Unit %s",Expi,si,splab,unitlab))
saveas(h4, fullfile(figdir, compose("Smooth_Map_cmp_%s_Exp%d_%s_%s.png",Animal,Expi,unitlab,splab)))

end
end
end
%%
writetable(struct2table(SmthStatTab),fullfile(figdir,"Beto_SmoothStatTab.csv"))
%%
SmthStatTable = struct2table(SmthStatTab);
%%
D1E_Dpr_V1 = SmthStatTable.D1E_Dpr((SmthStatTable.chan<=48)&(SmthStatTable.chan>=33));
D1E_Dpr_V4 = SmthStatTable.D1E_Dpr((SmthStatTable.chan>48)&(SmthStatTable.chan>=33));
D1E_Dpr_IT = SmthStatTable.D1E_Dpr((SmthStatTable.chan<=48)&(SmthStatTable.chan<33));
%% Visualize effect of bootstrapping on the activation maps. 
figure(7);
s = RandStream('mrg32k3a');
for i=1:500
actmap_btrp = cellfun(@(A)mean(A(datasample(s,1:numel(A),numel(A)))),act_col);
imagesc(actmap_btrp);
axis image
pause(0.05)
end
%%
function [laplsEng,TVEng,dirEng,laplsMap,gradPowMap] = map_energy(actmap)
kerd2 = [[0,-1,0];[-1,4,-1];[0,-1,0]]/4;
% laplsMap = conv2(actmap, kerd2);
laplsMap = conv2(padarray(actmap,[1,1],'replicate'), kerd2,'valid');
[actGx,actGy] = gradient(actmap);
gradPowMap = actGx.^2 + actGy.^2;
dirEng = sum(gradPowMap,'all');
TVEng = sum(sqrt(gradPowMap),'all');
laplsEng = sum(abs(laplsMap),'all');
end
function [laplsEng_sp, TVEng_sp, dirEng_sp, laplsMap_sp, gradNormMap_sp, detJac_map] = sph_map_energy(actmap)
[PHI, THETA] = meshgrid(-90:18:90,-90:18:90); 
kerd2 = [[0,-1,0];[-1,4,-1];[0,-1,0]]/4;
% detJac = abs(cosd(PHI));
Ginv11 = 1./cosd(PHI).^2;Ginv11(:,1)=0;Ginv11(:,end)=0; % These components explode, set to 0. 
Ginv22 = ones(size(PHI));
% Sample points for the gradients are far away.
detJac_map = cosd([mean(PHI(:,[1,2]),2),PHI(:,2:end-1),mean(PHI(:,[end-1,end]),2)]);
actmap_sp = actmap;
actmap_sp(:,1)=mean(actmap(:,1));
actmap_sp(:,end)=mean(actmap(:,end));
laplsMap_sp = conv2(padarray(actmap_sp,[1,1],'replicate'), kerd2,'valid');%del2(actmap_sp);%
[actGx_sp,actGy_sp] = gradient(actmap_sp);
gradNormMap_sp = actGx_sp.^2 .* Ginv22 + actGy_sp.^2 .*Ginv11;
dirEng_sp = sum(gradNormMap_sp.*detJac_map,'all');
TVEng_sp = sum(sqrt(gradNormMap_sp).*detJac_map,'all');
laplsEng_sp = sum(abs(laplsMap_sp).*detJac_map,'all');
end