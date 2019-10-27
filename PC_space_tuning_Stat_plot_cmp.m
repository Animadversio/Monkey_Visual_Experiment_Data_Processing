%% Line Plot for the statistics F and kappa 
%  Plot the channels' statistics over different experiment and subspace! 
Set_Exp_Specs
global sphere_norm Trials channel rasters ang_step Reps meta
ExpNum = 10;
%% 
Fstat_arr = [];
kstat_arr = [];
RootDir = "C:\\Users\\ponce\\OneDrive\\Desktop\\OneDrive_Binxu\\OneDrive\\PC_space_tuning";
figure(13);clf;set(13, 'position', [1          41        2560         963]);
for Expi=1:10%1:6
%     rasters = storedStruct.rasters{Expi};
    Trials = storedStruct.Trials{Expi};
    meta = storedStruct.meta{Expi};
    padded_name_arr = generate_unit_labels(meta.spikeID, [], '%03d');
    pref_chan = pref_chan_arr(Expi);
    sphere_norm = norm_arr(Expi);
    ang_step = 18;
    savepath = sprintf("C:\\Users\\ponce\\OneDrive\\Desktop\\OneDrive_Binxu\\OneDrive\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
    load(fullfile(savepath,"Basic_Stats.mat")) % Stat_summary
    load(fullfile(savepath,"KentFit_Stats.mat")) % Param_summary
    
    pref_chan_msk = contains(padded_name_arr, num2str(pref_chan, '%03d'));
    pref_chan_idx = find(pref_chan_msk);
    ax1 = subplot(2, ExpNum, Expi);hold on 
    Fstat_arr = [Fstat_arr; plotStat(ax1, Stat_summary, pref_chan_idx, 'anova_F')];
    ylim([0,25])
    title(['Exp',num2str(Expi),'  pref channel',num2str(pref_chan)])
    if Expi == 1, ylabel('anova F'), end
    ax2 = subplot(2, ExpNum, 1*ExpNum + Expi);hold on 
    kstat_arr = [kstat_arr; plotStat(ax2, Param_summary, pref_chan_idx, 'kappa')];
    ylim([0,16])
    if Expi == 1, ylabel('Kent Fit kappa'), end
%     
%     Full_Data_table = cell2table([cellfun(@(c) {c.t_CI}, Stat_summary),...
%                               cellfun(@(c) {c.t_p}, Stat_summary), ...
%                               cellfun(@(c) {c.anova_F}, Stat_summary), ...
%                               cellfun(@(c) {c.anova_p}, Stat_summary)], ...
%                               'VariableNames',{'PC23_CI' 'PC4950_CI' 'RND12_CI' 'PC23_p' 'PC4950_p' 'RND12_p' 'PC23_ANOVA_F' 'PC4950_ANOVA_F' 'RND12_ANOVA_F' 'PC23_ANOVA_p' 'PC4950_ANOVA_p' 'RND12_ANOVA_p'}, ...
%                               'RowNames', unit_name_arr); % Can be appended to more data! 
%     writetable(Full_Data_table, fullfile(savepath,'StatsTable.csv'),'QuoteStrings',true,'WriteRowNames',true);

    
%     imagesc(cell2mat(cellfun(@(c) {c.anova_F}, Stat_summary)))
%     title(sprintf("Exp%d F Stats, ANOVA1", Expi))
%     colorbar()
%     yticks(0.5:1:length(unit_name_arr)-0.5);yticklabels(unit_name_arr)
%     xticks(1:3)
%     saveas(8, fullfile(savepath, sprintf("Exp%d_ANOVA_stat.bmp", Expi)))
end
set(13, 'position', [1          41        2560         963])
saveas(13, fullfile(RootDir, sprintf("TuningStatCmp_EolveUnit")))
% Doing ANOVA on the tuning statistics values 
anova1(kstat_arr, {'PC23', 'PC4950', 'RND12'})
ylabel("\kappa value Kent fitting")
title({"Comparison of tuning Statistics","between Subspaces"})
set(gcf,'position',[1328         130         373         548])
saveas(gcf, fullfile(RootDir, sprintf("kappa_value_cmp_EvolveUnit.png")))
anova1(Fstat_arr, {'PC23', 'PC4950', 'RND12'})
ylabel("F value ANOVA")
title({"Comparison of tuning Statistics","between Subspaces"})
set(gcf,'position',[1328         130         373         548])
saveas(gcf, fullfile(RootDir, sprintf("F_value_cmp_EvolveUnit.png")))

%%
function stat_arr = plotStat(ax, stat_array, chan_idx, stat_name)
    stat_arr = [];
    for idx = chan_idx
        plot_stat = [];
        for j =1:3
            plot_stat(j) = getfield(stat_array{idx,j}, stat_name);
        end
        plot(ax, plot_stat)
        stat_arr = [stat_arr; plot_stat ];
    end
    xlim([0.5,3.5])
    xticks(1:3)
    xticklabels({'PC23','PC4950',"RND12"})
end