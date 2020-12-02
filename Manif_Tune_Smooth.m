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
btrp_Scol_sp = [];
for rep = 1:1000
actmap_btrp = cellfun(@(A)mean(A(datasample(1:numel(A),numel(A)))),act_col);
[D2E_sp_btrp, TVE_sp_btrp, D1E_sp_btrp, D2map_sp_btrp, GPmap_sp_btrp, detJac_map] = sph_map_energy(actmap_btrp);
btrp_Scol_sp = [btrp_Scol_sp;[D2E_sp_btrp, TVE_sp_btrp, D1E_sp_btrp]];
end
shfl_Scol_sp = [];
for rep = 1:1000
actmap_shfl = reshape(actmap(randperm(numel(actmap))),size(actmap));
[D2E_sp_shfl, TVE_sp_shfl, D1E_sp_shfl, D2map_sp_shfl, GPmap_sp_shfl, detJac_map] = sph_map_energy(actmap_shfl);
shfl_Scol_sp = [shfl_Scol_sp;[D2E_sp_shfl, TVE_sp_shfl, D1E_sp_shfl]];
end
%%
h=figure(5);clf;
T = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile(T,1);hold on
histogram(btrp_Scol_sp(:,1),'Norm','pdf');
histogram(shfl_Scol_sp(:,1),'Norm','pdf');
[H,P,~,tSt] = ttest2(btrp_Scol_sp(:,1), shfl_Scol_sp(:,1));
Dpr = computeCohen_d(btrp_Scol_sp(:,1), shfl_Scol_sp(:,1));
title(compose("Laplacian Filtered Energy\n t=%.2f(%.1e) d'=%.2f",tSt.tstat,P,Dpr))
nexttile(T,2);hold on
histogram(btrp_Scol_sp(:,2),'Norm','pdf');
histogram(shfl_Scol_sp(:,2),'Norm','pdf');
[H,P,~,tSt] = ttest2(btrp_Scol_sp(:,2), shfl_Scol_sp(:,2));
Dpr = computeCohen_d(btrp_Scol_sp(:,2), shfl_Scol_sp(:,2));
title(compose("Total Variation Energy\n t=%.2f(%.1e) d'=%.2f",tSt.tstat,P,Dpr))
nexttile(T,3);hold on
histogram(btrp_Scol_sp(:,3),'Norm','pdf');
histogram(shfl_Scol_sp(:,3),'Norm','pdf');
[H,P,~,tSt] = ttest2(btrp_Scol_sp(:,3), shfl_Scol_sp(:,3));
Dpr = computeCohen_d(btrp_Scol_sp(:,3), shfl_Scol_sp(:,3));
title(compose("Dirichlet Energy\n t=%.2f(%.1e) d'=%.2f",tSt.tstat,P,Dpr))
legend(["Trial Resampled", "Matrix Shuffled"])
title(T,"Map Smoothness Energy Compared")
%%
figure(6);
T = tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
nexttile(T,1);
imagesc(actmap_sp)
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
title("Gradient Power * Grid Weights")
xlabel("PC3 Phi")
nexttile(T,4);
imagesc(abs(D2map_sp_btrp))
axis image;colorbar
title("Laplacian Filtered Amplitude")
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
title(T,compose("Exp %d Space%d Unit %s",Expi,si,Stats(Expi).units.unit_name_arr(Stats(Expi).units.pref_chan_id(ui))))

%%
function [laplsEng,TVEng,dirEng,laplsMap,gradPowMap] = map_energy(actmap)
kerd2 = [[0,-1,0];[-1,4,-1];[0,-1,0]]/4;
laplsMap = conv2(actmap, kerd2);
[actGx,actGy] = gradient(actmap);
gradPowMap = actGx.^2 + actGy.^2;
dirEng = sum(gradPowMap,'all');
TVEng = sum(sqrt(gradPowMap),'all');
laplsEng = sum(abs(laplsMap),'all');
end
function [laplsEng_sp, TVEng_sp, dirEng_sp, laplsMap_sp, gradNormMap_sp, detJac_map] = sph_map_energy(actmap)
[PHI, THETA] = meshgrid(-90:18:90,-90:18:90); 
% detJac = abs(cosd(PHI));
Ginv11 = 1./cosd(PHI).^2;Ginv11(:,1)=0;Ginv11(:,end)=0; % These components explode, set to 0. 
Ginv22 = ones(size(PHI));
% Sample points for the gradients are far away.
detJac_map = cosd([mean(PHI(:,[1,2]),2),PHI(:,2:end-1),mean(PHI(:,[end-1,end]),2)]);
actmap_sp = actmap;
actmap_sp(:,1)=mean(actmap(:,1));
actmap_sp(:,end)=mean(actmap(:,end));
kerd2 = [[0,-1,0];[-1,4,-1];[0,-1,0]]/4;
laplsMap_sp = conv2(padarray(actmap_sp,[1,1],'replicate'), kerd2,'valid');%del2(actmap_sp);%
[actGx_sp,actGy_sp] = gradient(actmap_sp);
gradNormMap_sp = actGx_sp.^2 .* Ginv22 + actGy_sp.^2 .*Ginv11;
dirEng_sp = sum(gradNormMap_sp.*detJac_map,'all');
TVEng_sp = sum(sqrt(gradNormMap_sp).*detJac_map,'all');
laplsEng_sp = sum(abs(laplsMap_sp).*detJac_map,'all');
end