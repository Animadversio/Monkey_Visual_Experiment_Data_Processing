%% PC structure of Evol Trajectories (Manifold paper supp figure 5)
%%
Animal = "Both"; Set_Path;
mat_dir = "O:\Mat_Statistics";
pyenv("Version","C:\ProgramData\Anaconda3\envs\tf\python.exe")
py.importlib.import_module("numpy");
%%
Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
% load(fullfile(mat_dir, Animal+"_ManifMapVarStats.mat"),'MapVarStats')
%%
ETrajStats = repmat(struct(),1,numel(EStats));
for Expi = 1:numel(EStats)
    stim_path = EStats(Expi).meta.stimuli;
    stim_path = strrep(stim_path,'N:','S:');
    [codes_all, img_ids, generations] = load_codes_all(stim_path, 1);
    data = py.numpy.load(fullfile(stim_path, "PC_imgs", "PC_vector_data.npz"));
    PC_vecs = data.get('PC_vecs').single';
    PCprojected = codes_all*PC_vecs;
    
%     [PCaxes, PCprojected, latent] = pca(codes_all,'NumCom',50);
    ETrajStats(Expi).generations = generations;
    ETrajStats(Expi).PCaxes = PC_vecs;
    ETrajStats(Expi).PCprojected = PCprojected;
    
    initgen_PCprojs_m = mean(PCprojected(generations==min(generations),:),1);
    initgen_PCprojs_s = sem(PCprojected(generations==min(generations),:),1);
    lastgen_PCprojs_m = mean(PCprojected(generations==max(generations),:),1);
    lastgen_PCprojs_s = sem(PCprojected(generations==max(generations),:),1);
    ETrajStats(Expi).initgen_PCprojs_m = initgen_PCprojs_m;
    ETrajStats(Expi).initgen_PCprojs_s = initgen_PCprojs_s;
    ETrajStats(Expi).lastgen_PCprojs_m = lastgen_PCprojs_m;
    ETrajStats(Expi).lastgen_PCprojs_s = lastgen_PCprojs_s;
%     colorseq = brewermap(5,'Spectral');
%     figure(1);clf;set(1,'pos',[500,200,900,450])
%     subplot(121);hold on
%     for i = 1:3
%     plot(generations,PCprojected(:,i),'o','color',[colorseq(i,:),0.3])
%     end
%     subplot(122);
%     plot(abs(mean(PCprojected(generations==max(generations),:),1)))
%     sgtitle(compose("%s Exp%02d",Animal,Expi))
%     pause
end
%%
save(fullfile(mat_dir,Animal+"_EvolTrajPCcoefs.mat"),'ETrajStats')

%% Find the generation with the highest score / activation and project it on the Manifold plane
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_EvolTrajPCcoefs.mat"),'ETrajStats');
load(fullfile(mat_dir,Animal+"_Evol_stats.mat"),'EStats');
for Expi = 1:numel(EStats)
wdw = 51:200;
scoretraj = cellfun(@(P)mean(P(:,wdw,:),'all'),EStats(Expi).evol.psth);
scoretraj_smth = movmean(scoretraj(1:end-1), 3);
generations = ETrajStats(Expi).generations;
PCprojected = ETrajStats(Expi).PCprojected;
[maxScore,maxGeni] = max(scoretraj_smth);
bestgen_PCprojs_m = mean(PCprojected(generations==maxGeni,:),1);
bestgen_PCprojs_s = sem(PCprojected(generations==maxGeni,:),1);
ETrajStats(Expi).bestgen_PCprojs_m = bestgen_PCprojs_m;
ETrajStats(Expi).bestgen_PCprojs_s = bestgen_PCprojs_s;
end
save(fullfile(mat_dir,Animal+"_EvolTrajPCcoefs.mat"),'ETrajStats');
end

%% Plot the PC23 plane end points on the coordinate plane or in 3d space 
PC123projs_mat = cell2mat(arrayfun(@(E)E.lastgen_PCprojs_m(1:3),ETrajStats','uni',0));
[THE,PHI,Rvec] = cart2sph(-PC123projs_mat(:,1),PC123projs_mat(:,2),PC123projs_mat(:,3));
%%
figure;
scatter3(-PC123projs_mat(:,1),PC123projs_mat(:,2),PC123projs_mat(:,3))
axis equal
%%
figure;
scatter(THE,PHI);% very clustered
xlim([-pi/2,pi/2]);ylim([-pi/2,pi/2]);axis equal;
%% get the basis from python saving
% py.importlib.import_module("torch");
% data = py.numpy.load(fullfile(stim_path,"PC_imgs","PC_vector_data.npz"));
% PC_vecs = data.get('PC_vecs').double;


%% Reload the data and do statistics 
ETrajStats = [];
for Animal = ["Alfa","Beto"]
D = load(fullfile(mat_dir,Animal+"_EvolTrajPCcoefs.mat"),'ETrajStats');
ETrajStats = [ETrajStats,D.ETrajStats];
end
%%
PC123projs_last = cell2mat(arrayfun(@(E)E.lastgen_PCprojs_m(1:3),ETrajStats','uni',0));
[THE_L,PHI_L,Rvec_L] = cart2sph(-PC123projs_last(:,1),PC123projs_last(:,2),PC123projs_last(:,3));
PC123projs_best = cell2mat(arrayfun(@(E)E.bestgen_PCprojs_m(1:3),ETrajStats','uni',0));
[THE_B,PHI_B,Rvec_B] = cart2sph(-PC123projs_best(:,1),PC123projs_best(:,2),PC123projs_best(:,3));
%%
EStats_all = struct();
for Animal = ["Alfa","Beto"]
load(fullfile(mat_dir,Animal+"_Evol_stats.mat"),'EStats')
EStats_all.(Animal) = EStats;
end
poptabdir = "O:\Manif_Fitting\popstats";
poptab = [readtable(fullfile(poptabdir,"Alfa_Exp_all_KentStat_bsl_pole.csv"));...
          readtable(fullfile(poptabdir,"Beto_Exp_all_KentStat_bsl_pole.csv"))];
% Create the masks
drivermsk = zeros(size(poptab,1),1,'logical'); % Masks of real driver units instead of using the first. 
for i = 1:numel(drivermsk)
    driver_unit = EStats_all.(poptab.Animal{i})(poptab.Expi(i)).evol.unit_in_pref_chan;
    drivermsk(i) = (poptab.unitnum(i) == driver_unit) & (poptab.chan(i) == poptab.prefchan(i));
end
driver_fittab = poptab(drivermsk, :);
driverPC23_fittab = driver_fittab(driver_fittab.space==1,:);
%%
%% Visualize the PC trajectory
figdir = "O:\Manuscript_Manifold\FigureS2\TrajPCStruct";
%% Visualize individual Experiment's PC coefficient. 
Expi = 5;
colorseq = brewermap(9,'Spectral');%'Spectral');
[PC_mean,PC_std,blockarr] = block_summarize(ETrajStats(Expi).PCprojected,ETrajStats(Expi).generations);
figure(3);clf;set(gcf,'pos',[500,200,500,350]);hold on
for i = 1:8
shadedErrorBar(blockarr,PC_mean(:,i),PC_std(:,i),'lineprop',{'color',[colorseq(i,:),0.6],'DisplayName',"PC"+num2str(i)})
% plot(blockarr,PC_mean(:,i),'color',[colorseq(i,:),0.6],'DisplayName',"PC"+num2str(i))
% scatter(ETrajStats(Expi).generations,ETrajStats(Expi).PCprojected(:,i),'o','markeredgecolor',[colorseq(i,:)],'markeredgealpha',0.6,...
%     'markerfacecolor',[colorseq(i,:)],'markerfacealpha',0.6,'DisplayName',"PC"+num2str(i))
end
legend('location','best')
ylabel("Proj Coefficient");xlabel("Generation")
title(compose("%s Exp%02d latent code PC Trajectory","Alfa",Expi))
saveallform(figdir,compose("%s_Exp%02d_PCtraj","Alfa",Expi))
%% Summary of all evolution trajectory in PC spaces.
figure(4);T=tiledlayout(3,1,'pad','compact','tilesp','compact');
for Expi = 1:numel(ETrajStats)
[PC_mean,PC_std,blockarr] = block_summarize(ETrajStats(Expi).PCprojected,ETrajStats(Expi).generations);
nexttile(T,1);hold on;plot(blockarr,PC_mean(:,1),'color',[0,0,0,0.25],'LineWidth',1)
%shadedErrorBar(blockarr,PC_mean(:,1),PC_std(:,1),'lineprop',{'color',[0,0,0,0.3]})
nexttile(T,2);hold on;plot(blockarr,PC_mean(:,2),'color',[0,0,0,0.25],'LineWidth',1)
%shadedErrorBar(blockarr,PC_mean(:,2),PC_std(:,2),'lineprop',{'color',[0,0,0,0.3]})
nexttile(T,3);hold on;plot(blockarr,PC_mean(:,3),'color',[0,0,0,0.25],'LineWidth',1)
%shadedErrorBar(blockarr,PC_mean(:,3),PC_std(:,3),'lineprop',{'color',[0,0,0,0.3]})
end
nexttile(T,1);xlim([0,80]);ylabel("PC1")
nexttile(T,2);xlim([0,80]);ylabel("PC2")
nexttile(T,3);xlim([0,80]);ylabel("PC3");xlabel("Generation")
title(T,"All Evoltion Trajectory in PC space")
saveallform(figdir,compose("AllEvol_PC123coef"))


%%
figure;
subplot(121);scatter(PHI_B', driverPC23_fittab.phi);xlim([-pi/2,pi/2]);ylim([-pi/2,pi/2]);axis equal square
subplot(122);scatter(THE_B', driverPC23_fittab.theta);xlim([-pi/2,pi/2]);ylim([-pi/2,pi/2]);axis equal square
%%
[cval,pval] = corr(PHI_B(msk), driverPC23_fittab.phi(msk),'Rows','complete')
[cval,pval] = corr(THE_B(msk), driverPC23_fittab.theta(msk),'Rows','complete')
%% Figure Compare the final and best generation latent code projection and the peak location in Manifold exp.
msk = driverPC23_fittab.R2>0.5;
Amsk = driverPC23_fittab.Animal=="Alfa";
Bmsk = driverPC23_fittab.Animal=="Beto";
Cord = colororder;
figure(5);set(5,'pos',[1000         378         630         415])
T=tiledlayout(1,2,'pad','compact','tilesp','compact');
nexttile(T,1);
scatter(THE_B(msk), PHI_B(msk),'DisplayName',"Best Gen");hold on
scatter(THE_L(msk), PHI_L(msk),'DisplayName',"Last Gen");
xlabel("PC2 (Theta)");ylabel("PC3 (Phi)");title("Best Evol Codes Projection")
xlim([-pi/2,pi/2]);ylim([-pi/2,pi/2]);axis square
xticks(-1.5:0.5:1.5);yticks(-1.5:0.5:1.5);
plot([0,0],[-pi/2,pi/2],'-.k','HandleVisibility','off');plot([-pi/2,pi/2],[0,0],'-.k','HandleVisibility','off')
legend('location','best')
nexttile(T,2);hold on
% scatter(driverPC23_fittab.theta(msk), driverPC23_fittab.phi(msk));
scatter(driverPC23_fittab.theta(msk&Amsk), driverPC23_fittab.phi(msk&Amsk),'Marker','o','MarkerEdgeColor',Cord(1,:),'DisplayName',"Monkey A");
scatter(driverPC23_fittab.theta(msk&Bmsk), driverPC23_fittab.phi(msk&Bmsk),'Marker','*','MarkerEdgeColor',Cord(1,:),'DisplayName',"Monkey B");
xlabel("PC2 (Theta)");ylabel("PC3 (Phi)");title("Manif Peak (Kent fit)")
xlim([-pi/2,pi/2]);ylim([-pi/2,pi/2]);axis square
xticks(-1.5:0.5:1.5);yticks(-1.5:0.5:1.5);
plot([0,0],[-pi/2,pi/2],'-.k','HandleVisibility','off');plot([-pi/2,pi/2],[0,0],'-.k','HandleVisibility','off')
legend('location','best')
[cval_PHI_B,pval_PHI_B] = corr(PHI_B(msk), driverPC23_fittab.phi(msk),'Rows','complete');
[cval_THE_B,pval_THE_B] = corr(THE_B(msk), driverPC23_fittab.theta(msk),'Rows','complete');
[tval,pval,sumstr_Phi] = ttest_print(driverPC23_fittab.phi(msk),'Phi');
[tval,pval,sumstr_The] = ttest_print(driverPC23_fittab.theta(msk),'Theta');
title(T,compose("Correlation of Coordinate of Best Codes in Evol and Manifold Peak\nTheta %.3f(P=%.1e) Phi %.3f(P=%.1e)\n%s\n%s",...
        cval_THE_B,pval_THE_B,cval_PHI_B,pval_PHI_B,sumstr_The,sumstr_Phi))
saveallform(figdir,compose("Evol-Manif_BestImg_coord_cmp"))
%%
fprintf("Last generation Theta(PC2) deviates by %.3f+-%.3f rad (%.1f deg), Phi(PC3) deviated by %.3f+-%.3f rad (%.1f deg)\n",...
    mean(abs(THE_L)),sem(abs(THE_L)),rad2deg(mean(abs(THE_L))),mean(abs(PHI_L)),sem(abs(PHI_L)),rad2deg(mean(abs(PHI_L))))
%%
ttest_print(PHI_B,'Phi Best');
ttest_print(THE_B,'Theta Best');
[tval,pval,sumstr_Phi] = ttest_print(driverPC23_fittab.phi(msk),'Phi fit');
[tval,pval,sumstr_The] = ttest_print(driverPC23_fittab.theta(msk),'Theta fit');
function [PC_mean,PC_std,blockarr] = block_summarize(PCcoefs,gens)
assert(numel(gens) == size(PCcoefs,1))
PC_mean = [];
PC_std = [];
blockarr = min(gens):max(gens);
for blocki = blockarr
PC_mean(blocki,:) = mean(PCcoefs(gens==blocki,:),1);
PC_std(blocki,:) = std(PCcoefs(gens==blocki,:),1,1);
end
end