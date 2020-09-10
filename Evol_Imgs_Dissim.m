%% 
Animal="Beto"; Set_Path;
ftrrows = find(contains(ExpRecord.expControlFN,"generate_") & ExpRecord.Exp_collection=="Manifold");
%%
D = torchImDist();
%%
stimuli_path = ExpRecord.stimuli{22}; %ExpRecord.stimuli{241}; % Manif Expi 11
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1); % each row is an evolved code
%%
avg_codes = arrayfun(@(geni)mean(codes_all(code_geni==geni,:)), 1:max(code_geni),'Uni',0);
avg_codes = cell2mat(avg_codes');
avg_imgs = G.visualize(avg_codes);
distmat = D.distmat(avg_imgs);
%%
figure;imagesc(distmat)
%%
all_imgs = G.visualize(codes_all);
%%
tic
distmat_all = D.distmat_B(all_imgs);
toc % 2577 sec for the accelerated distmat.
%%
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'))
%%
Evol_distmat = repmat(struct(),1,numel(EStats));
%%
for Expi = 1:numel(EStats)
    tic;
    fprintf("Expi %d\t",Expi)
    stimuli_path = EStats(Expi).meta.stimuli; %ExpRecord.stimuli{241}; % Manif Expi 11
    [codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1); % each row is an evolved code
    fprintf("Finish loading codes %.1fs\t",toc)
    avg_codes = arrayfun(@(geni)mean(codes_all(code_geni==geni,:)), 1:max(code_geni),'Uni',0);
    avg_codes = cell2mat(avg_codes');
    fprintf("(%d gens)\t",size(avg_codes,1))
    avg_imgs = G.visualize(avg_codes);
    fprintf("Finish visualize %.1fs\t",toc)
    D = D.select_metric("squeeze");
    distmat_s = D.distmat_B(avg_imgs);
    fprintf("Finish distmat %.1fs\t",toc)
    D = D.select_metric("alex");
    distmat_a = D.distmat_B(avg_imgs);
    fprintf("Finish distmat %.1fs\t",toc)
%     D = D.select_metric("vgg");
%     distmat_v = D.distmat_B(avg_imgs);
%     fprintf("Finish distmat %.1f\n",toc)
    Evol_distmat(Expi).avg.squ = distmat_s;
    Evol_distmat(Expi).avg.alex = distmat_a;
%     Evol_distmat(Expi).avg.vgg = distmat_v;
    L2dist = squareform(pdist(reshape(avg_imgs,256*256*3,[])'));
    Evol_distmat(Expi).avg.L2 = L2dist;
    fprintf("\n")
end
%%
save(fullfile(mat_dir, Animal+"_Evol_ImDist.mat"), 'Evol_distmat')