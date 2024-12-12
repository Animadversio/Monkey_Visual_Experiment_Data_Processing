%% Relative increase 
%  statistics and plots to test the relative activation increase for each thread. 
%  and also the difference in initial generation. 
figdir = "E:\OneDrive - Harvard University\Manuscript_BigGAN\Figures\Evol_activation_stats";
%%
Animal = "Both"; Set_Path;
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
%%
wdw = [51:200]; bslwdw = [1:50];
GANstr = ["fc6","BigGAN"];
STS_arr = [];
for iTr = 1:numel(BFEStats)
S = BFEStats(iTr);
% ui = BFEStats(iTr).evol.unit_in_pref_chan(1);
if isempty(S.evol), continue; end
if ~(contains(S.evol.space_names{1},"fc6")) && (contains(S.evol.space_names{2},"BigGAN"))
    % assert(contains(S.evol.space_names{1},"fc6"))
    % assert(contains(S.evol.space_names{2},"BigGAN"))
    fprintf('Exp%d %s %s %s\n', iTr, S.meta.ephysFN, ...
        S.evol.space_names{1}, S.evol.space_names{2})
    continue; 
end
ui=1;
activ_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:), [2])),S.evol.psth,"uni",0); % cell arr, vector of score per gen. 
score_col = cellfun(@(psth)squeeze(mean(psth(ui,wdw,:), [2]))-mean(psth(ui,bslwdw,:),[2,3]), S.evol.psth,"uni",0);
block_col = cellfun(@(idx) S.evol.block_arr(idx), S.evol.idx_seq,"uni",0); % cell arr, vector of block id per gen. 
STS = struct();
STS.Expi = iTr;
STS.Animal = S.Animal;
STS.prefchan = S.evol.pref_chan(1);
expdate = parseExpDate(S.meta.ephysFN, S.meta.expControlFN);
if (S.Animal == "Beto") && expdate > datetime(2021,08,01)
    array_layout = "Beto_new";
else
    array_layout = S.Animal;
end
visarea = area_map(S.evol.pref_chan(1), array_layout);
STS.area = visarea;
STS.ephysFN = S.meta.ephysFN;
for GANi = 1:2
% Drop the last block as we assume it's non-complete
% score_traces{iTr,GANi} = cat(1,score_col{GANi, 1:end-1}); 
% block_traces{iTr,GANi} = cat(1,block_col{GANi, 1:end-1});
% zscore_traces{iTr,GANi} = zscores_tsr(pref_chan_id, row_gen & thread_msks{threadi});
end_acts = cat(1,score_col{GANi, end-2:end-1});
init_acts = cat(1,score_col{GANi, 1:2});
score_m = cellfun(@mean, score_col(GANi, 1:end-1));
[~,maxN] = max(score_m);
if maxN == 1, maxN=maxN+1; end
max_acts = cat(1,score_col{GANi, maxN-1:maxN});
SM_init = mean_CI(init_acts, GANstr(GANi)+"_init", 0.05);
SM_end  = mean_CI(end_acts,  GANstr(GANi)+"_end", 0.05);
SM_max  = mean_CI(max_acts,  GANstr(GANi)+"_max", 0.05);
STS_end = ttest2struct(end_acts,init_acts,GANstr(GANi)+"_endinit",0.05);
STS_max = ttest2struct(max_acts,init_acts,GANstr(GANi)+"_maxinit",0.05);
STS = catstruct(STS,SM_init,SM_end,SM_max,STS_end,STS_max);
end
STS_arr = [STS_arr, STS];
end
acttab = struct2table(STS_arr);
%%
writetable(acttab,fullfile(figdir,"act_increase_summary.csv"))
%%
Evol_act_stats = STS_arr;
save(fullfile(figdir, "Evol_act_increase_stats.mat"),"Evol_act_stats")
%%
Alfamsk = acttab.Animal =="Alfa";
Betomsk = acttab.Animal =="Beto";
V4msk = acttab.area =="V4";
ITmsk = acttab.area =="IT";
%%
paired_stripe_plot({acttab.fc6_endinit_m,acttab.BigGAN_endinit_m},["fc6","BigGAN"],{},[])
%%
paired_stripe_plot({acttab.fc6_endinit_m,acttab.BigGAN_endinit_m},["fc6","BigGAN"],{Alfamsk,Betomsk},["Alfa","Beto"])
%%
paired_stripe_plot({acttab.fc6_endinit_m,acttab.BigGAN_endinit_m},["fc6","BigGAN"],{V4msk},["V4"])
%%
paired_stripe_plot({acttab.fc6_endinit_m,acttab.BigGAN_endinit_m},["fc6","BigGAN"],{ITmsk},["IT"])
%%
paired_stripe_plot({acttab.fc6_init_m,acttab.BigGAN_init_m},["fc6","BigGAN"],{Alfamsk,Betomsk},["Alfa","Beto"])
%%
paired_stripe_error_plot({acttab.fc6_init_m,acttab.BigGAN_init_m},...
    {acttab.fc6_init_CI, acttab.BigGAN_init_CI}, ["fc6","BigGAN"], {}, [])
%%
paired_stripe_error_plot({acttab.fc6_endinit_m,acttab.BigGAN_endinit_m},...
    {acttab.fc6_endinit_CI, acttab.BigGAN_endinit_CI}, ["fc6","BigGAN"], {}, [])
%%
figure(2);set(2,'pos',[ 680   631   374   347])
errorbar(acttab.fc6_endinit_m(Expi2cmp(:,1)), acttab.BigGAN_endinit_m(Expi2cmp(:,2)),...
        fina_Z_sem(Expi2cmp(:,1)), fina_Z_sem(Expi2cmp(:,1)),...
        fina_Z_sem(Expi2cmp(:,2)), fina_Z_sem(Expi2cmp(:,2)),'o')
xlabel("1 deg Evol");ylabel("3 deg Evol");
title("Final Activation z score Comparison across Image size")
add_diagonal();xlim([0,0.85]);ylim([0,0.85])
%%
h = xscatter_error_plot(acttab, "fc6_endinit_m", "BigGAN_endinit_m", "fc6_endinit_CI", "BigGAN_endinit_CI", ...
    {}, ["all"], "", false);
axis equal
add_diagonal(gca,'HandleVisibility','off');
%xlim([0,0.85]);ylim([0,0.85])
%%
h = xscatter_error_plot(acttab, "fc6_endinit_m", "BigGAN_endinit_m", "fc6_endinit_CI", "BigGAN_endinit_CI", ...
    {V4msk, ITmsk}, ["V4","IT"], "", false,"",'CapSize',1);
add_diagonal(gca,'HandleVisibility','off','Color','k');
hline(0,'k:')
vline(0,'k:')
axis equal
%%
alpha = 0.1;
ax = gca();
for i = 1:2
erb = ax.Children(i);
set([erb.Bar, erb.Line], 'ColorType', 'truecoloralpha', 'ColorData', [erb.Line.ColorData(1:3); 255*alpha])
set(erb.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [erb.Cap.EdgeColorData(1:3); 255*alpha])
end
%% 
function SM = mean_CI(x,prefix,alpha)
mean_x = mean(x);               % convert cell array to matrix, mean by row
std_x = std(x);  %mean_x/size(x,2);                     % std dev of mean(y) is mean(y)/nObs
SEM = (std_x / sqrt(length(x)));
ts = tinv([alpha/2, 1-alpha/2], length(x)-1);      % T-Score
CI95 = mean_x + ts.*(std_x/sqrt(length(x)));
SM = struct();
SM.(prefix+"_m") = mean_x;
SM.(prefix+"_sem") = SEM;
SM.(prefix+"_CI") = CI95;
end

function STS = ttest2struct(arr1,arr2,prefix,alpha)
STS = struct();
[~,p,CI,stats] = ttest2(arr1,arr2,'alpha',0.05);
STS.(prefix+"_m") = mean(arr1) - mean(arr2);
STS.(prefix+"_CI") = reshape(CI,1,[]);
STS.(prefix+"_tval") = stats.tstat;
STS.(prefix+"_t_p") = p;
STS.(prefix+"_t_df") = stats.df;
end
function [c] = appendStruct(a,b)
%APPENDSTRUCT Appends two structures ignoring duplicates
%   Developed to append two structs while handling cases of non-unique
%   fieldnames.  The default keeps the last occurance of the duplicates in
%   the appended structure.
ab = [struct2cell(a); struct2cell(b)];
abNames = [fieldnames(a); fieldnames(b)];
[~,iab] = unique(abNames,'last');
c = cell2struct(ab(iab),abNames(iab));
end