Animal="Alfa"; Set_Path;
mat_dir = "O:\Mat_Statistics";
tabdir = "O:\Manif_Fitting\popstats";
% load(fullfile(Matdir, "Beto_ManifPopDynamics.mat"))
load(fullfile(mat_dir,Animal+"_ManifMapVarStats.mat"),'MapVarStats')
load(fullfile(mat_dir,Animal+"_Manif_stats.mat"),'Stats')
%%
addpath D:\Github\Fit_Spherical_Tuning
addpath e:\Github_Projects\Fit_Spherical_Tuning
pool = parpool(4);
T00 = tic;
parfor Expi = 1:numel(Stats)
% Basic information to collect into the Stats
nCh = numel(MapVarStats(Expi).units.spikeID);
Animal_tab = array2table(repmat(Animal,nCh,1),'VariableNames',{'Animal'});
Expi_tab = array2table(repmat(Expi,nCh,1),'VariableNames',{'Expi'});
unitstr_tab = array2table(MapVarStats(Expi).units.unit_name_arr,'VariableNames',{'unitstr'});
unit_tab = array2table(MapVarStats(Expi).units.unit_num_arr,'VariableNames',{'unitnum'});
chan_tab = array2table(MapVarStats(Expi).units.spikeID,'VariableNames',{'chan'});
prefchan_tab = array2table(repmat(MapVarStats(Expi).units.pref_chan,nCh,1),'VariableNames',{'prefchan'});
tic
for si = 1:MapVarStats(Expi).manif.subsp_n
space_tab = array2table(repmat(si,nCh,1),'VariableNames',{'space'});
actmap_col = arrayfun(@(iCh)cellfun(@(A)...
    squeeze(A(iCh,:)),MapVarStats(Expi).manif.act_col{si},'uni',0),...
    1:nCh,'uni',0);
actmap_mean = arrayfun(@(iCh)cellfun(@(A)...
    mean(A(iCh,:),2),MapVarStats(Expi).manif.act_col{si},'uni',1),...
    1:nCh,'uni',0);
toc
KentStats_exp = cellfun(@(actmap)fit_Kent_stats(actmap),actmap_mean);
anovaStats_exp = cellfun(@anova_cells,actmap_col);
keyboard
toc
anovaStats_tab = struct2table(anovaStats_exp);
KentStats_tab = struct2table(KentStats_exp);
StatsTab = [Animal_tab, Expi_tab, unitstr_tab, unit_tab, chan_tab, prefchan_tab, space_tab, anovaStats_tab, smthStats_tab];
writetable(StatsTab,fullfile(tabdir,compose("%s_Exp%d_sp%d.csv",Animal,Expi,si)));
% StatsTab_sum = [StatsTab_sum;StatsTab];
% disp(size(StatsTab_sum))
end
end
toc(T00);

%%
Kent_stats = repmat(struct(), 1, length(Stats));
for i = 1:length(Stats)
% there can be multiple units, or multiple subspace, so we have this cell
% array
stats_grid = cell(numel(Stats(i).units.pref_chan_id), Stats(i).manif.subsp_n);
scr_stats_grid = cell(numel(Stats(i).units.pref_chan_id), Stats(i).manif.subsp_n);
subsp_str = ['PC23',"PC4950",'RND12'];
subsp_axis = ["PC2","PC3";"PC49","PC50";"RND1","RND2"];
for subsp_i = 1:Stats(i).manif.subsp_n
imgnm_grid = cellfun(@(idx) unique([Stats(i).imageName(idx)]), Stats(i).manif.idx_grid{subsp_i});
psths = Stats(i).manif.psth{subsp_i};
fprintf(subsp_str{subsp_i}+"\n")
for uniti = 1:1:numel(Stats(i).units.pref_chan_id)
channel_j = Stats(i).units.pref_chan_id(uniti);
chan_label_str = sprintf("%s Exp%d Channel %s", Stats(i).Animal, Stats(i).Expi, ...
            Stats(i).units.unit_name_arr{channel_j});

mean_score_map = cellfun(@(psth) mean(psth(uniti, 51:200, :),[2,3])-mean(psth(uniti, 1:40, :),[2,3]), psths);
[Parameter, gof] = fit_Kent(mean_score_map);
scr_stats_grid{uniti, subsp_i} = gof;
% V = coeffvalues(Parameter);
% CI = confint(Parameter);
% param_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
%                 V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
mean_act_map = cellfun(@(psth) mean(psth(uniti, 51:200, :),[2,3]), psths);
[Parameter, gof] = fit_Kent(mean_act_map);
stats_grid{uniti, subsp_i} = gof;
V = coeffvalues(Parameter);
CI = confint(Parameter);
param_str = sprintf("theta=%.2f (%.2f, %.2f)  phi=%.2f (%.2f, %.2f)\n psi=%.2f (%.2f, %.2f)  A=%.2f (%.2f, %.2f)\n kappa=%.2f (%.2f, %.2f)  beta=%.2f (%.2f, %.2f) ", ...
                V(1), CI(1,1), CI(2,1), V(2), CI(1,2), CI(2,2), V(3), CI(1,3), CI(2,3), V(6), CI(1,6), CI(2,6), V(4), CI(1,4), CI(2,4), V(5), CI(1,5), CI(2,5));
fprintf(chan_label_str+"\n "+param_str+"\n")
end
end
Kent_stats(i).act_fit = stats_grid;
Kent_stats(i).scr_fit = scr_stats_grid;
end
toc

%%

function S = fit_Kent_stats(funcval)
[Parameter, gof] = fit_Kent(funcval);
S.sse = gof.sse;
S.R2 = gof.rsquare;
S.dfe = gof.dfe;
S.adjR2 = gof.adjrsquare;
S.rmse = gof.rmse;
for i = 1:numel(gof.coefname)
S.(gof.coefname(i)) = gof.coef(i);
S.(gof.coefname(i)+"_CI") = gof.confint(:,i);
end
end



