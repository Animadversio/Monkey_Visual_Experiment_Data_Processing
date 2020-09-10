%% Hess_Frame_Evol_Analysis_batch
Animal="Beto"; Set_Path;
% evol_ids = find(contains(ExpRecord.expControlFN,"generate_integrated") & ExpRecord.Exp_collection=="ReducDimen_Evol") ;
evol_ids = find(contains(ExpRecord.expControlFN,"generate_integrated") & ExpRecord.Exp_collection=="Seq_Evol") ;
[meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(evol_ids,"Beto",true,false);
%%
global final_gen_xexp G
G = FC6Generator();
%% load the Hessian Matrices
py.importlib.import_module("numpy")
data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\HessTune\NullSpace\Evolution_Avg_Hess.npz");
% data = py.numpy.load("E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy\Nullspace\Pasu_Space_Avg_Hess.npz");
eigvect = data.get("eigvect_avg").double; % each column is an eigen vector. 
eigvals = data.get("eigv_avg").double;
eigvals = eigvals(end:-1:1);
eigvect = eigvect(:,end:-1:1);
%% load the evolution codes
for evoli = evol_ids(31:end)'
disp(ExpRecord.comments(evoli))
stimuli_path = ExpRecord.stimuli{evoli}; %ExpRecord.stimuli{241}; % Manif Expi 11
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1); % each row is an evolved code
% [codes_all_rd, img_ids_rd, code_geni_rd] = load_codes_all(stimuli_path, 2); % each row is an evolved code
%
Hproj_coef = codes_all * eigvect;
Hproj_coef_rd = codes_all_rd * eigvect;
%
geni = max(code_geni);
[~,P,~,STATS] = ttest2(abs(mean(Hproj_coef(code_geni==geni,1:100),1)), abs(mean(Hproj_coef(code_geni==geni,300:400),1)));
fprintf("1-100 vs 300-400, t=%.3f (p=%.3f)\n",STATS.tstat,P)
[~,P,~,STATS] = ttest2(abs(mean(Hproj_coef(code_geni==geni,1:100),1)), abs(mean(Hproj_coef(code_geni==geni,3996:4096),1)));
fprintf("1-100 vs 3996-4096, t=%.3f (p=%.3f)\n",STATS.tstat,P)
[~,P,~,STATS] = ttest2(abs(mean(Hproj_coef(code_geni==geni,300:400),1)), abs(mean(Hproj_coef(code_geni==geni,3996:4096),1)));
fprintf("300-400 vs 3996-4096, t=%.3f (p=%.3f)\n",STATS.tstat,P)
end
%% Collect evolution experiment on both monkey 
Animal="Both"; Set_Path;
final_gen_xexp = []; 
perm_gen_xexp = []; 
evol_ids = find(contains(ExpRecord.expControlFN,"generate_inte") | contains(ExpRecord.expControlFN,"generate_para")) ; % & ExpRecord.Exp_collection=="Seq_Evol"
%%
for evoli = 1:numel(evol_ids)
fprintf("%d/%d\n",evoli,numel(evol_ids))
stimuli_path = ExpRecord.stimuli{evol_ids(evoli)}; %ExpRecord.stimuli{241}; % Manif Expi 11
try
if ~exist(stimuli_path,'dir'), continue;end
[codes_all, img_ids, code_geni] = load_codes_all(stimuli_path, 1, 'last'); % each row is an evolved code
mean_gene = mean(codes_all,1);
final_gen_xexp = [final_gen_xexp; mean_gene];
for shufi = 1:1
perm_gen_xexp = [perm_gen_xexp; mean_gene(randperm(4096))];
end
catch ME
    disp(ExpRecord(evol_ids(evoli),:))
    disp(ME)
    continue
end
end
%%
save(fullfile(summarydir,"evollastgencode_col.mat"),'final_gen_xexp')
%%
load(fullfile(summarydir,"ctrl_code_col.mat"),'ctrlcodes_col')
%%
ctrlcodes_shf = cell2mat(arrayfun(@(idx)ctrlcodes_col(idx,randperm(4096)),1:300,'uni',0)');
ctrl_proj_cc = ctrlcodes_col * eigvect;
ctrl_perm_proj_cc = ctrlcodes_shf * eigvect;
% ctrl_proj_cc = ctrl_proj_cc./norm_axis(ctrl_proj_cc,2).*300;
% ctrl_perm_proj_cc = ctrl_perm_proj_cc./norm_axis(ctrl_perm_proj_cc,2).*300;
%%
init_proj_cc = initcodes * eigvect;
code_proj_cc = final_gen_xexp * eigvect;
perm_proj_cc = perm_gen_xexp * eigvect;
% code_proj_cc = code_proj_cc./norm_axis(code_proj_cc,2).*300;
% perm_proj_cc = perm_proj_cc./norm_axis(perm_proj_cc,2).*300;
%%
figure(2);clf;hold on;set(2,'pos',[568         236        1325         666])
scatter([1:4096],mean(code_proj_cc,1)-mean(init_proj_cc,1),12,'k')
scatter([1:4096],mean(ctrl_proj_cc,1)-mean(init_proj_cc,1),5,'blue')
legend(["Evolved","NoiseDriven"],'FontSize',14);xlim([-50,4150])
ylabel("Projection Coefficient Residue (-Init code)",'FontSize',15);xlabel("Eigen Idx",'FontSize',15)
title("Scatter of Residue Mean Projection Coefficient (~ Init Code) on Each Axes",'FontSize',18)
saveas(2,fullfile(summarydir,"Proj_cc_Residue_cmp_scatter.png"))
savefig(2,fullfile(summarydir,"Proj_cc_Residue_cmp_scatter.fig"))
%%
figure(2);clf;hold on;set(2,'pos',[568         236        1325         666])
errorbar([1:4096],mean(code_proj_cc,1)-mean(init_proj_cc,1),sem(code_proj_cc,1),'.','Color','k','CapSize',1)
errorbar([1:4096],mean(ctrl_proj_cc,1)-mean(init_proj_cc,1),sem(ctrl_proj_cc,1),'.','Color','blue','CapSize',1)
legend(["Evolved","NoiseDriven"],'FontSize',14);xlim([-50,4150])
ylabel("Projection Coefficient Residue (-Init code)",'FontSize',15);xlabel("Eigen Idx",'FontSize',15)
title("Scatter of Residue Mean Projection Coefficient (~ Init Code) on Each Axes",'FontSize',18)
saveas(2,fullfile(summarydir,"Proj_cc_Residue_cmp_errbar.png"))
savefig(2,fullfile(summarydir,"Proj_cc_Residue_cmp_errbar.fig"))
%%
figure(1);clf;hold on
scatter([1:4096],mean(code_proj_cc,1),12,'k')
scatter([1:4096],mean(perm_proj_cc,1),5,'g')
scatter([1:4096],mean(ctrl_proj_cc,1),5,'blue')
scatter([1:4096],mean(init_proj_cc,1),5,'r')
legend(["Evolved","Shuffle","NoiseDriven","Initial Gen"],'FontSize',14);xlim([-50,4150])
ylabel("Projection Coefficient",'FontSize',15);xlabel("Eigen Idx",'FontSize',15)
title("Scatter of Mean Projection Coefficient on Each Axes",'FontSize',18)
saveas(1,fullfile(summarydir,"Proj_cc_cmp_scatter.png"))
savefig(1,fullfile(summarydir,"Proj_cc_cmp_scatter.fig"))
%%
figure(4);clf;hold on
errorbar([1:4096],mean(code_proj_cc,1),sem(code_proj_cc,1),'.','Color','k','CapSize',1)
errorbar([1:4096],mean(perm_proj_cc,1),sem(perm_proj_cc,1),'.','Color','g','CapSize',1)
errorbar([1:4096],mean(ctrl_proj_cc,1),sem(ctrl_proj_cc,1),'.','Color','blue','CapSize',1)
% errorbar([1:4096],mean(init_proj_cc,1),sem(init_proj_cc,1),'.','Color','r','CapSize',1)
legend(["Evolved","Shuffle","NoiseDriven","Initial Gen"],'FontSize',14);xlim([-50,4150])
ylabel("Projection Coefficient",'FontSize',15);xlabel("Eigen Idx",'FontSize',15)
title("Scatter of Mean Projection Coefficient on Each Axes",'FontSize',18)
% saveas(4,fullfile(summarydir,"Proj_cc_cmp_errbar.png"))
% savefig(4,fullfile(summarydir,"Proj_cc_cmp_errbar.fig"))
%%
[~, ttst_P, ~, ttst_Stat] = arrayfun(@(id)ttest2(abs(code_proj_cc(:,id)),abs(perm_proj_cc(:,id))),1:4096,'uni',0);
ttst_Stat = cell2mat(ttst_Stat); ttst_P = cell2mat(ttst_P); 
[rksm_P, ~, rksm_Stat] = arrayfun(@(id)ranksum(abs(code_proj_cc(:,id)),abs(perm_proj_cc(:,id))),1:4096,'uni',1);
[~, krsm_P, krsm_Stat] = arrayfun(@(id)kstest2(abs(code_proj_cc(:,id)),abs(perm_proj_cc(:,id)),'alpha',0.01),1:4096,'uni',1);
%%
[~, ctrl_ttst_P, ~, ctrl_ttst_Stat] = arrayfun(@(id)ttest2(abs(code_proj_cc(:,id)),abs(ctrl_proj_cc(:,id))),1:4096,'uni',0);
ctrl_ttst_Stat = cell2mat(ctrl_ttst_Stat); ctrl_ttst_P = cell2mat(ctrl_ttst_P); 

[~, ctrl_ttst_abs_P, ~, ctrl_ttst_abs_Stat] = arrayfun(@(id)ttest2(code_proj_cc(:,id),ctrl_proj_cc(:,id)),1:4096,'uni',0);
ctrl_ttst_abs_Stat = cell2mat(ctrl_ttst_abs_Stat); ctrl_ttst_abs_P = cell2mat(ctrl_ttst_abs_P); 

[ctrl_rksm_P, ~, ctrl_rksm_Stat] = arrayfun(@(id)ranksum(abs(code_proj_cc(:,id)),abs(ctrl_proj_cc(:,id))),1:4096,'uni',1);
[~, ctrl_krsm_P, ctrl_krsm_Stat] = arrayfun(@(id)kstest2(abs(code_proj_cc(:,id)),abs(ctrl_proj_cc(:,id)),'alpha',0.01),1:4096,'uni',1);
%%
[ctrl_rksm_P, ~, ctrl_rksm_Stat] = arrayfun(@(id)ranksum(abs(code_proj_cc(:,id)),abs(ctrl_proj_cc(:,id))),1:4096,'uni',1);
[~, ctrl_krsm_P, ctrl_krsm_Stat] = arrayfun(@(id)kstest2(abs(code_proj_cc(:,id)),abs(ctrl_proj_cc(:,id)),'alpha',0.01),1:4096,'uni',1);
%%
[ctrl_rksm_P, ~, ctrl_rksm_Stat] = arrayfun(@(id)ranksum(code_proj_cc(:,id),ctrl_proj_cc(:,id)),1:4096,'uni',1);
[~, ctrl_krsm_P, ctrl_krsm_Stat] = arrayfun(@(id)kstest2(code_proj_cc(:,id),ctrl_proj_cc(:,id),'alpha',0.01),1:4096,'uni',1);
%%
% [ctrl_sgrk_P, ~, ctrl_sgrk_Stat] = arrayfun(@(id)signrank(code_proj_cc(:,id),ctrl_proj_cc(:,id)),1:4096,'uni',0);
% ctrl_sgrk_Stat = cell2mat(ctrl_sgrk_Stat); ctrl_sgrk_P = cell2mat(ctrl_sgrk_P); 
% [ctrl_sgrk_abs_P, ~, ctrl_sgrk_abs_Stat] = arrayfun(@(id)signrank(abs(code_proj_cc(:,id)),abs(ctrl_proj_cc(:,id)),'tail','right'),1:4096,'uni',0);
% ctrl_sgrk_abs_Stat = cell2mat(ctrl_sgrk_abs_Stat); ctrl_sgrk_abs_P = cell2mat(ctrl_sgrk_abs_P); 
%%
corr_P = mafdr(ctrl_krsm_P,'BHFDR',true); % nothing fall out of it
% should we test the absolute value or the raw value. 
%%
candid_id = find(corr_P<0.05);
%%
% candid_id = find(ctrl_rksm_P<0.05 & ctrl_rksm_Stat<0.05 & rksm_P<0.05 & krsm_P<0.05);

for eigid = candid_id
figure(9);clf;hold on
histogram(code_proj_cc(:,eigid),40,'FaceAlpha',0.3,'EdgeColor','none') % ,-20:2:30
histogram(ctrl_proj_cc(:,eigid),40,'FaceAlpha',0.5,'EdgeColor','none') % ,-20:2:30
histogram(perm_proj_cc(:,eigid),40,'FaceAlpha',0.5,'EdgeColor','none') % ,-20:2:30
% histogram(ctrl_perm_proj_cc(:,1),40,'FaceAlpha',0.5,'EdgeColor','none')
legend(["Neural Evol", "Noise Evol", "Neural Shuffled", "Noise Shuffled"])
title(compose("Hessian Eigenvector %d projection\nKS test ~shuffle p=%.1e ~noise p=%.1e\nranksum ~shuffle p=%.1e ~noise p=%.1e",eigid,...
    krsm_P(eigid),neurctrl_krsm_P(eigid),rksm_P(eigid),neurctrl_rksm_P(eigid)));xlabel("Projection Coefficient on Eigvect")
saveas(9,fullfile(summarydir,compose("proj_cc_cmp_eig%d.jpg",eigid)))
pause;
end
%%
for eigid = candid_id
imgs = G.visualize(eigvect(:,eigid)*linspace(-30,30,7));
figure(11);
imshow(imtile(imgs,'GridSize',[1,7]))
title(compose("Image Change Along Eigen Vector %d Projection CC 5-95 percentile [%.1f,%.1f]",eigid,prctile(code_proj_cc(:,eigid),[5,95])))
xlabel("Deviation from center [-30,30]")
saveas(11,fullfile(summarydir,compose("G_img_axis_eig%d.jpg",eigid)))
pause;
end
%%
proj_cc_msk = (abs(code_proj_cc) > prctile(abs(code_proj_cc),99.5,2));
figure(3);stem(sum(proj_cc_msk,1),'marker','none')
xlim([-25,490]);xlabel("Eig Id");
ylabel({'Number of Being Top0.5% of the Coefficent'});
title({'Frequency of Extreme Coefficients Occuring on Each Eigen Axes'});
saveas(3,fullfile(summarydir,"ExtrmValDist.png"))
savefig(3,fullfile(summarydir,"ExtrmValDist.fig"))
%%
code_proj_cc_orig = final_gen_xexp * eigvect;
[sorted_proj_cc, sorted_idx] = sort(code_proj_cc_orig,1,'Descend'); % code_proj_cc_orig
[sorted_proj_cc_abs, sorted_idx_abs] = sort(abs(code_proj_cc_orig),1,'Descend'); % code_proj_cc_orig
%%
code_proj_cc_orig = final_gen_xexp * eigvect;
[sorted_proj_cc, sorted_idx] = sort(code_proj_cc,1,'Descend'); % code_proj_cc_orig
[sorted_proj_cc_abs, sorted_idx_abs] = sort(abs(code_proj_cc),1,'Descend'); % code_proj_cc_orig
%%
eigid=12;
top_cc = median(sorted_proj_cc(1:8,eigid));
bot_cc = median(sorted_proj_cc(end-7:end,eigid));
med_cc = median(sorted_proj_cc_abs(end-7:end,eigid));
img_axis = G.visualize(eigvect(:,eigid)*[top_cc,bot_cc,med_cc]*5);
imgs_top = G.visualize(final_gen_xexp(sorted_idx(1:8,eigid),:));
imgs_bot = G.visualize(final_gen_xexp(sorted_idx(end-7:end,eigid),:));
imgs_med = G.visualize(final_gen_xexp(sorted_idx_abs(end-7:end,eigid),:));
imgs_top_m = G.visualize(mean(final_gen_xexp(sorted_idx(1:8,eigid),:),1));
imgs_bot_m = G.visualize(mean(final_gen_xexp(sorted_idx(end-7:end,eigid),:),1));
imgs_med_m = G.visualize(mean(final_gen_xexp(sorted_idx_abs(end-7:end,eigid),:),1));
figure(5);
montage(cat(4,img_axis(:,:,:,1),imgs_top,imgs_top_m,...
              img_axis(:,:,:,2),imgs_bot,imgs_bot_m,...
              img_axis(:,:,:,3),imgs_med,imgs_med_m),'size',[3,10])
title(compose("Maximal, Minimal, and Min Abs Projected value on Eigen Axis %d Max %.1f Min %.1f Mid %.1f", eigid, top_cc, bot_cc, med_cc))
%%
figure();
montage(cat(4,imgs_top_m,imgs_bot_m,imgs_med_m),'size',[3,1])

%%
% code_proj_cc_orig = final_gen_xexp * eigvect;
%%
eigid=30;
imgs_top = G.visualize(final_gen_xexp(sorted_idx_abs(1:8,eigid),:));
imgs_bot = G.visualize(final_gen_xexp(sorted_idx_abs(end-7:end,eigid),:));
figure(6);
montage(cat(4,imgs_top,imgs_bot),'size',[2,8])
title(compose("Maximal and Minimal Projected ABOLUTE value on Eigen Axis %d Max %.1f Min %.1f", eigid, median(sorted_proj_cc_abs(1:8,eigid)), median(sorted_proj_cc_abs(end-7:end,eigid))))
%%
code_proj_cc_orig
%%
%%
[red_cc,umap,clusterId,extras] = run_umap(double(code_proj_cc_orig(:,1:800).*eigvals(1:800)),'metric','correlation','n_neighbors',10,'n_epochs',1000);
%%
figure(10);
% scatter(red_cc(:,1),red_cc(:,2))
S = scatter(red_cc(:,1),red_cc(:,2));
S.ButtonDownFcn = @VisCode;
xlabel("U1");ylabel("U2");
title("Evolved code Dimension Reduced to 2d")
%%
figure(19);
S = scatter(code_proj_cc_orig(:,1),code_proj_cc_orig(:,5));
S.ButtonDownFcn = @VisCode;
xlabel("Eigvec1");ylabel("Eigvec5");
title("Evolved code scattered on Eig Spaces")
% S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Subj", ExpTabPool.Animal(tunedrows));
% S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp", ExpTabPool.Expi(tunedrows));
% S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Area", ExpTabPool.area(tunedrows));
% S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Chan", ExpTabPool.chan_str(tunedrows));
% S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Driver", ExpTabPool.prefchan(tunedrows));
% S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Pos", compose("(%.1f %.1f)",ExpTabPool.imgpos(tunedrows,:)));
% S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Size", ExpTabPool.imgsize(tunedrows));
% S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("F", ExpTabPool.F(tunedrows));
% S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("T", ExpTabPool.T(tunedrows));
% S.DataTipTemplate.DataTipRows(1:2)=[];
function [coordinateSelected, rowi] = VisCode(hObj, event)
global final_gen_xexp G
x = hObj.XData; 
y = hObj.YData; 
pt = event.IntersectionPoint(1:2);       % The (x0,y0) coordinate you just selected
coordinates = [x(:),y(:)];     % matrix of your input coordinates
dist = pdist2(pt,coordinates); 
[~, expi] = min(dist);
% rowi = tunedrows(rel_rowi);
% Expi = ExpTabPool.Expi(rowi); iCh = ExpTabPool.iCh(rowi); Subj = ExpTabPool.Animal(rowi);
% titlestr = compose("row %d %s Exp %d prefCh %d Ch %s\n pos %s size %d F=%.1f T=%.1f",rowi,...
%     ExpTabPool.Animal(rowi),ExpTabPool.Expi(rowi),ExpTabPool.prefchan(rowi),ExpTabPool.chan_str(rowi),compose("(%.1f,%.1f)",ExpTabPool.imgpos(rowi,:)),ExpTabPool.imgsize(rowi),ExpTabPool.F(rowi),ExpTabPool.T(rowi));
% fprintf("Row %d in table, %s\n",rowi,titlestr)
figure(13);
imshow(G.visualize(final_gen_xexp(expi,:)))
end