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
%%
figdir = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning";
Animal = "Alfa";
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'))
load(fullfile(mat_dir, Animal+"_Evol_ImDist.mat"), 'Evol_distmat')
load(fullfile(mat_dir, Animal+"_Manif_ImDist.mat"), "ManifImDistStat")
load(fullfile(mat_dir, "gab_imdist.mat"),'gab_imdist')
load(fullfile(mat_dir, "pasu_imdist.mat"),'pasu_imdist')
Manif_distmat = ManifImDistStat;
%%
pasu_nm_grid = cellfun(@(idx)unique(Stats(1).imageName(idx)),Stats(1).ref.pasu_idx_grid,'Uni',0);
pasu_nm_grid = reshape(pasu_nm_grid',[],1); % reshape into one row. But transpose to make similar object adjacent
pasu_val_msk = ~cellfun(@isempty,pasu_nm_grid);
%% Plot that tuning averaging over each block
for Expi = 15:numel(EStats)
    ui =1;si=1;
    score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),EStats(Expi).evol.psth)';
    evkrt_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all')-mean(psth(ui,1:50,:),'all'),EStats(Expi).evol.psth)';
    [sortScore,sortId]=sort(score_vec,'Descend');
    [maxScore,maxId]=max(score_vec);
    if size(Evol_distmat(Expi).avg.squ,1)==numel(score_vec)+1
    Evol_distmat(Expi).avg.squ = Evol_distmat(Expi).avg.squ(1:end-1,1:end-1);
    Evol_distmat(Expi).avg.L2 = Evol_distmat(Expi).avg.L2(1:end-1,1:end-1);
    Evol_distmat(Expi).avg.alex = Evol_distmat(Expi).avg.alex(1:end-1,1:end-1);
    end
    figure(1);
    T=tiledlayout(1,4,'TileSpacing','compact');
    title(T,compose("%s Exp %d prefchan %d",Animal,Expi,EStats(Expi).units.pref_chan),'FontSize',16)
    % Evolution Data
    nexttile(1)%subplot(141);
    scatter(Evol_distmat(Expi).avg.squ(:,maxId),score_vec)
    xlabel("Perceptual Similarity(SqueezeNet)");ylabel("Neural Activation")
    titstr1 = { compose("Evol: pear %.3f spear %.3f",corr(Evol_distmat(Expi).avg.squ(:,maxId),score_vec),corr(Evol_distmat(Expi).avg.squ(:,maxId),score_vec,'Type','Spearman'))};
    % nexttile(2)
    % scatter(ManifImDistStat(Expi).SSIM(:,maxId),score_vec)
    % xlabel("SSIM")
    nexttile(2)
    scatter(Evol_distmat(Expi).avg.alex(:,maxId),score_vec)
    titstr2 = { compose("Evol: pear %.3f spear %.3f",corr(Evol_distmat(Expi).avg.alex(:,maxId),score_vec),corr(Evol_distmat(Expi).avg.alex(:,maxId),score_vec,'Type','Spearman'))};
    xlabel("Perceptual Similarity(alexnet)")
    nexttile(3)
    scatter(Evol_distmat(Expi).avg.L2(:,maxId),score_vec)
    titstr3 = { compose("Evol: pear %.3f spear %.3f",corr(Evol_distmat(Expi).avg.L2(:,maxId),score_vec),corr(Evol_distmat(Expi).avg.L2(:,maxId),score_vec,'Type','Spearman'))};
    xlabel("L2")
    % Manifold Data
    score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
    evkrt_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all')-mean(psth(ui,1:50,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
    [sortScore,sortId]=sort(score_vec,'Descend');
    [maxScore,maxId]=max(score_vec);
    nexttile(1);hold on%subplot(141);
    scatter(Manif_distmat(Expi).squ(:,maxId), score_vec)
    titstr1{end+1} = compose("Manif: pear %.3f spear %.3f",corr(Manif_distmat(Expi).squ(:,maxId),score_vec),corr(Manif_distmat(Expi).squ(:,maxId),score_vec,'Type','Spearman'));
    % Pasupathy Data
    if Stats(Expi).ref.didPasu
    pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
    pasu_vec(~pasu_val_msk) = [];
    [pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
    [pasu_maxScore,pasu_maxId]=max(pasu_vec);
    nexttile(1);hold on %subplot(141);
    scatter(pasu_imdist.squ(:,pasu_maxId),pasu_vec)
    titstr1{end+1} = compose("Pasu: pear %.3f spear %.3f",corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Rows','complete'),corr(pasu_imdist.squ(:,pasu_maxId),pasu_vec,'Type','Spearman','Rows','complete'));
    end
    if Stats(Expi).ref.didGabor
    gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
    [gab_sortScore,sortId]=sort(gab_vec,'Descend');
    [gab_maxScore,gab_maxId]=max(gab_vec);
    nexttile(1);hold on
    scatter(gab_imdist.squ(:,gab_maxId),gab_vec)%,'g'
    titstr1{end+1} = compose("Gabor: pear %.3f spear %.3f",corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Rows','complete'),corr(gab_imdist.squ(:,gab_maxId),gab_vec,'Type','Spearman','Rows','complete'));
    end
    nexttile(1);title(titstr1);legend(["Evol","Manif","Pasupathy","Gabor"])
    nexttile(2);title(titstr2)
    nexttile(3);title(titstr3)
    
    pause
end