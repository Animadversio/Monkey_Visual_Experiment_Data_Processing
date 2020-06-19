%% This code collect the SNR and related Stats for manifold and evolution
%  and then do correlation among them to find the relationship
%%  Load data
Animal = "Alfa";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
load(fullfile(MatStats_path, compose("%s_Manif_stats.mat", Animal)), 'Stats')
Kent_path = "E:\OneDrive - Washington University in St. Louis\Manif_SUHash\summary";
load(fullfile(Kent_path, Animal+"_Manif_Kent_Fit.mat"),"Kent_stats")
%% Prepare mask for experiments
prefchan_arr = arrayfun(@(E) E.units.pref_chan,EStats);
prefname_arr = arrayfun(@(E) E.units.unit_name_arr(E.units.pref_chan_id),EStats);
V1msk = prefchan_arr <=48 & prefchan_arr>=33;
V4msk = prefchan_arr <=64 & prefchan_arr>=49;
ITmsk = prefchan_arr <=32 & prefchan_arr>=1;
msk.V1 = V1msk;msk.V4 = V4msk; msk.IT = ITmsk; msk.all = ones(1,length(Stats),'logical');
ColorArr = [V1msk',V4msk',ITmsk'];

%% Evolution Ref Signal to Noise Ratio
%  Visualize the scatter of mean and std of firing rate towards reference
%  image.
for Expi = 22%45:length(Stats)
refimgmean = cellfun(@(psth)mean(psth(1,window,:),'all'),EStats(Expi).ref.psth_arr);
refimgbsl = cellfun(@(psth)mean(psth(1,1:50,:),'all'),EStats(Expi).ref.psth_arr);
refimgstd = cellfun(@(psth)std(mean(psth(1,window,:),2)),EStats(Expi).ref.psth_arr);
figure;
scatter(refimgmean,refimgstd)
ylabel("std");xlabel("mean")
title(compose("Exp %d pref chan %s", Expi, prefname_arr(Expi)));
pause
end
%%
ChanQualStats = repmat(struct(),length(Stats),1);
%% Compute the Channel Quality statistics (T and F for tunings)
nanstruct = struct('F',nan,'F_P',nan,'F_bsl',nan,'F_P_bsl',nan,'T',nan,'t_P',nan);
for Expi = 1:length(Stats)
fprintf("Processing %s Evol Exp %d\n",Animal,Expi)
ChanQualStats(Expi).evol.ref = calc_tune_stats(EStats(Expi).ref.psth_arr);
ChanQualStats(Expi).evol.gen = calc_tune_stats(EStats(Expi).evol.psth);
ui = EStats(Expi).units.unit_num_arr(EStats(Expi).units.pref_chan_id);
getunit = @(pstharr)cellfun(@(psth)psth(ui,:,:),pstharr,'UniformOutput', false);
for si = 1:Stats(Expi).manif.subsp_n
ChanQualStats(Expi).manif.gen(si) = calc_tune_stats(getunit(Stats(Expi).manif.psth{si}));
end
if Stats(Expi).ref.didGabor
gabStats = calc_tune_stats(getunit(Stats(Expi).ref.gab_psths));
ChanQualStats(Expi).manif.gabor = gabStats;
else
ChanQualStats(Expi).manif.gabor = nanstruct;    
end
if Stats(Expi).ref.didPasu
pasuStats = calc_tune_stats(getunit(Stats(Expi).ref.pasu_psths));
ChanQualStats(Expi).manif.pasu = pasuStats;
else
ChanQualStats(Expi).manif.pasu = nanstruct;
end
end
%% Get the kappa values for evolving channel
drivechan_tune_kappa = zeros(length(Stats),1);
for Expi = 1:length(Stats)
ui = EStats(Expi).units.unit_num_arr(EStats(Expi).units.pref_chan_id);
si = 1; % number of subspace, in case there are multiple subspace there.
kappa = Kent_stats(Expi).act_fit{ui,si}.coef(4);
drivechan_tune_kappa(Expi) = kappa;
ChanQualStats(Expi).manif.kappa = kappa;
end
%% Evolution Statistics Extraction
EvolSuccStats = repmat(struct(),1,length(EStats));
for Expi = 1:length(Stats)
blockn = EStats(Expi).evol.block_n;
block_mean_score = cellfun(@(psth)mean(psth(1,51:200,:),'all'),EStats(Expi).evol.psth(1:blockn-1));
[~,peakBlock]=max(block_mean_score);
if peakBlock == blockn-1, peakBlock = blockn-2; end % avoid the last generation.
endspsths = cell2mat(reshape(EStats(Expi).evol.psth(peakBlock:peakBlock+1),1,1,[]));
initpsths = cell2mat(reshape(EStats(Expi).evol.psth(2:3),1,1,[]));
window = [51:200];
initacts = squeeze(mean(initpsths(1,window,:),2));
endsacts = squeeze(mean(endspsths(1,window,:),2));
[H,P,CI,STATS] = ttest2(endsacts, initacts);
EvolSuccStats(Expi).tstat = STATS.tstat;
EvolSuccStats(Expi).DAOA = (mean(endsacts) - mean(initacts)) / mean(initacts);
ChanQualStats(Expi).evol.tstat = STATS.tstat;
ChanQualStats(Expi).evol.DAOA = (mean(endsacts) - mean(initacts)) / mean(initacts);
end
DAOA_arr = arrayfun(@(E)E.DAOA, EvolSuccStats)';
tstat_arr = arrayfun(@(E)E.tstat, EvolSuccStats)';
%% Save or load the statistics
save(fullfile(MatStats_path,Animal+"_ChanQualStat.mat"),'ChanQualStats')
% load(fullfile(MatStats_path,Animal+"_ChanQualStat.mat"),'ChanQualStats')
%% Put stats into tabel to do corrplot
QualTab = [arrayfun(@(C)C.evol.ref.F,ChanQualStats),... 
           arrayfun(@(C)C.evol.ref.T,ChanQualStats), ...
           arrayfun(@(C)C.manif.gabor.F,ChanQualStats), ...
           arrayfun(@(C)C.manif.gabor.T,ChanQualStats), ...
           arrayfun(@(C)C.manif.pasu.F,ChanQualStats), ...
           arrayfun(@(C)C.manif.pasu.T,ChanQualStats)];
SU_index = arrayfun(@(S)S.meta.SUidx,Stats)'; % EStats is the same number. 
DAOA_arr = arrayfun(@(C)C.evol.DAOA, ChanQualStats);
tstat_arr = arrayfun(@(C)C.evol.tstat, ChanQualStats);
drivechan_tune_kappa = arrayfun(@(C)C.manif.kappa,ChanQualStats);
%% struct to collect stats of different area
rho_mats = struct();PVal_mats = struct();
%%
figure(3);clf;set(3,'position',[405          42        1155         954]);
ax = subtightplot(1,1,1); % make the image tighter
[rho_mat, PVal_mat,~] = corrplot(ax,[drivechan_tune_kappa, DAOA_arr, tstat_arr, SU_index, QualTab, ],'testR','on','alpha',0.005,...
    'varNames',["kappa","DAOA","evolT", "SUidx","eRefF","eRefT","gabF","gabT","pasuF","pasuT"]);
axis equal
annotation(figure(3),'textbox',[0.150 0.982 0.700 0.031],...
    'String',compose("%s all experiment correlation of statistics of driving channel",Animal),...
    'FontSize',18,'FitBoxToText','off', 'EdgeColor','none');
rho_mats.all = rho_mat; PVal_mats.all = PVal_mat;
savedir = "E:\OneDrive - Washington University in St. Louis\Evol_Succ_Manif_Kappa";
saveas(3,fullfile(savedir, Animal+"_Kappa_ChanQual_corrplot.jpg"))
savefig(3,fullfile(savedir, Animal+"_Kappa_ChanQual_corrplot.fig"))

%% The same thing for experiment driven by 3 different areas
corrdata = [drivechan_tune_kappa, DAOA_arr, tstat_arr, SU_index, QualTab, ];
for area = ["V1", "V4", "IT"]
figure(5);clf;set(5,'position',[405          42        1155         954]);
ax = subtightplot(1,1,1);
[rho_mat, PVal_mat,~] = corrplot(ax,corrdata(msk.(area),:),'testR','on','alpha',0.01,...
    'varNames',["kappa","DAOA","evolT", "SUidx","eRefF","eRefT","gabF","gabT","pasuF","pasuT"]);
axis equal
annotation(figure(5),'textbox',[0.150 0.982 0.700 0.031],...
    'String',compose("%s %s experiment correlation of statistics of driving channel",Animal,area),...
    'FontSize',18,'FitBoxToText','off', 'EdgeColor','none');
% rho_mats.(area) = rho_mat; PVal_mats.(area) = PVal_mat;
saveas(5,fullfile(savedir, compose("%s_Kappa_ChanQual_corrplot_%s.jpg",Animal,area)))
savefig(5,fullfile(savedir, compose("%s_Kappa_ChanQual_corrplot_%s.fig",Animal,area)))
% pause
end
%%
varNames = ["kappa","DAOA","evolT","SUidx","evolRefF","evolRefT","manifGabF","manifGabT","manifPasuF","manifPasuT"];
save(fullfile(savedir, Animal+"_ChanQualcorrMatrix.mat"),'rho_mats', 'PVal_mats', 'varNames', 'QualTab', 'drivechan_tune_kappa', 'DAOA_arr', 'tstat_arr', 'SU_index')
%% What should be the criterion to exclude a channel.
msk.good =~((arrayfun(@(C)C.manif.pasu.t_P,ChanQualStats) > 0.01) & ...
            (arrayfun(@(C)C.manif.gabor.t_P,ChanQualStats) > 0.01) & ...
            (arrayfun(@(C)C.evol.ref.t_P,ChanQualStats) > 0.01));%find()
msk.nooutlier =~(DAOA_arr>12);
varNames = ["kappa","DAOA","evflT","maniF","maniT","evolF","evolT","SUidx","eRefF","eRefT","gabF","gabT","pasuF","pasuT"];
%% Correlate Everything
rho_mats = struct();PVal_mats = struct();
corrdata = [drivechan_tune_kappa, DAOA_arr, tstat_arr, ...
            arrayfun(@(C)C.manif.gen.F,ChanQualStats), arrayfun(@(C)C.manif.gen.T,ChanQualStats),...
            arrayfun(@(C)C.evol.gen.F,ChanQualStats), arrayfun(@(C)C.evol.gen.T,ChanQualStats),...
            SU_index, QualTab, ];
for area = "nooutlier"%["all","good","V1","V4","IT"]
figure(6);clf;set(6,'position',[405          42        1155         954]);
ax = subtightplot(1,1,1);
[rho_mat, PVal_mat,~] = corrplot(ax,corrdata(msk.(area),:),'testR','on','alpha',0.01,...
    'varNames',varNames);
axis equal
annotation(figure(6),'textbox',[0.150 0.982 0.700 0.031],...
    'String',compose("%s %s experiment correlation of statistics of driving channel",Animal,area),..."Good Signal"
    'FontSize',18,'FitBoxToText','off', 'EdgeColor','none');
rho_mats.(area) = rho_mat; PVal_mats.(area) = PVal_mat;
saveas(6,fullfile(savedir, compose("%s_All_corrplot_%s.jpg",Animal,area)))
savefig(6,fullfile(savedir, compose("%s_All_corrplot_%s.fig",Animal,area)))
end
%% Correlate the T and F statistics with tuning kappa
figure;
subtightplot(1,2,1,0.05,0.08,0.05)
S1 = scatter(arrayfun(@(C)C.evol.ref.T,ChanQualStats), drivechan_tune_kappa, 25, ColorArr);
S1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp",1:length(Stats));
S1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("prefchan",prefname_arr);
ylabel("Manif Tuning Kappa");xlabel("t stat for Evol")
title(compose("Corr Coef %.3f",corr(tstat_arr,drivechan_tune_kappa)))
subtightplot(1,2,2,0.05,0.10,0.05)
S2 = scatter(arrayfun(@(C)C.evol.ref.F,ChanQualStats),  drivechan_tune_kappa, 25, ColorArr);
S2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp",1:length(Stats));
S2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("prefchan",prefname_arr);
msk = DAOA_arr < 2 & drivechan_tune_kappa <2 ;
xlabel("Delta Activation / Activation_0")
title(compose("Corr Coef %.3f (< 2 subregion Corr Coef %.3f)",...
    corr(DAOA_arr,drivechan_tune_kappa),corr(DAOA_arr(msk),drivechan_tune_kappa(msk))))
arrayfun(@(C)C.evol.ref.F,ChanQualStats)
%%

%%
% function reportStats = calc_Stats(psth_cells)
% % PSTH Chacteristics
% reportStats = struct();
% groupsize = cellfun(@(psth) size(psth,3), psth_cells);
% indices = reshape(1:numel(psth_cells),size(psth_cells));
% idx_vect = arrayfun(@(L, idx) idx*ones(L,1), groupsize, indices, 'UniformOutput', false);
% 
% act_wdw = 51:200; bsl_wdw=1:50;
% score_vect = cellfun(@(psth)squeeze(mean(psth(1,act_wdw,:),2)),psth_cells,'UniformOutput',false);
% basel_vect = cellfun(@(psth)squeeze(mean(psth(1,bsl_wdw,:),2)),psth_cells,'UniformOutput',false);
% idx_vect = cell2mat(reshape(idx_vect,[],1));
% score_vect = cell2mat(reshape(score_vect,[],1));
% basel_vect = cell2mat(reshape(basel_vect,[],1));
% [P,ANOVATAB,STATS] = anova1(score_vect,idx_vect,'off');
% reportStats.F = ANOVATAB{2,5};
% reportStats.F_P = P;
% [P,ANOVATAB,STATS] = anova1(basel_vect,idx_vect,'off');
% reportStats.F_bsl = ANOVATAB{2,5};
% reportStats.F_P_bsl = P;
% [H,P,CI,STATS] = ttest2(score_vect,basel_vect);
% reportStats.T = STATS.tstat;
% reportStats.t_P = P;
% end