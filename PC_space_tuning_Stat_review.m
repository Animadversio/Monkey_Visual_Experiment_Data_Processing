%% PC space tuning Stats Review
% read the statstics in matlab / export it in readable format! 

Set_Exp_Specs
global sphere_norm Trials channel rasters ang_step Reps meta
%%
for Expi=10%1:6
%     rasters = storedStruct.rasters{Expi};
    Trials = storedStruct.Trials{Expi};
    meta = storedStruct.meta{Expi};
    unit_name_arr = generate_unit_labels(meta.spikeID);
    [activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
    pref_chan = pref_chan_arr(Expi);
    sphere_norm = norm_arr(Expi);%328; % 269 Day3 % 326 Daye % 328 Day 1 [328, 326, 269, 329, 401]
    ang_step = 18;
    savepath = sprintf("C:\\Users\\ponce\\OneDrive\\Desktop\\OneDrive_Binxu\\OneDrive\\PC_space_tuning\\Exp%d_chan%02d", Expi, pref_chan);
    load(fullfile(savepath,"Basic_Stats.mat"))
    Full_Data_table = cell2table([cellfun(@(c) {c.t_CI}, Stat_summary),...
                              cellfun(@(c) {c.t_p}, Stat_summary), ...
                              cellfun(@(c) {c.anova_F}, Stat_summary), ...
                              cellfun(@(c) {c.anova_p}, Stat_summary)], ...
                              'VariableNames',{'PC23_CI' 'PC4950_CI' 'RND12_CI' 'PC23_p' 'PC4950_p' 'RND12_p' 'PC23_ANOVA_F' 'PC4950_ANOVA_F' 'RND12_ANOVA_F' 'PC23_ANOVA_p' 'PC4950_ANOVA_p' 'RND12_ANOVA_p'}, ...
                              'RowNames', unit_name_arr); % Can be appended to more data! 
    writetable(Full_Data_table, fullfile(savepath,'StatsTable.csv'),'QuoteStrings',true,'WriteRowNames',true);
    
    
    figure(8);clf;set(8, 'position', [977    42   392   954]);
    imagesc(cell2mat(cellfun(@(c) {c.anova_F}, Stat_summary)))
    title(sprintf("Exp%d F Stats, ANOVA1", Expi))
    colorbar()
    yticks(0.5:1:length(unit_name_arr)-0.5);yticklabels(unit_name_arr)
    xticks(1:3)
    saveas(8, fullfile(savepath, sprintf("Exp%d_ANOVA_stat.bmp", Expi)))
end
%%
Data_table = cell2table(cellfun(@(c) {mean(c.anova_F)}, Stat_summary), ...
    'VariableNames',{'PC23' 'PC4950' 'RND12'}, 'RowNames', unit_name_arr);
%%
Full_Data_table = cell2table([cellfun(@(c) {c.t_CI}, Stat_summary),...
                              cellfun(@(c) {c.t_p}, Stat_summary), ...
                              cellfun(@(c) {c.anova_F}, Stat_summary), ...
                              cellfun(@(c) {c.anova_p}, Stat_summary)], ...
                              'VariableNames',{'PC23_CI' 'PC4950_CI' 'RND12_CI' 'PC23_p' 'PC4950_p' 'RND12_p' 'PC23_ANOVA_F' 'PC4950_ANOVA_F' 'RND12_ANOVA_F' 'PC23_ANOVA_p' 'PC4950_ANOVA_p' 'RND12_ANOVA_p'}, ...
                              'RowNames', unit_name_arr);
writetable(Full_Data_table, fullfile(savepath,'StatsTable.csv'),'QuoteStrings',true,'WriteRowNames',true);
%%
function unit_name_arr = generate_unit_labels()
% Generate the unit labels 17B from the spikeID variable
global meta
Unit_id = meta.spikeID;
unit_name_arr = {}; % name tag for each unit 
for i = 1:length(Unit_id)
    cur_chan = Unit_id(i);
    if sum(Unit_id == cur_chan) == 1
        unit_name_arr{i} = num2str(cur_chan);
    else
        cur_chan = Unit_id(i);
        rel_idx = find(find(Unit_id == cur_chan) == i);
        unit_name_arr{i} = [num2str(cur_chan), char(64+rel_idx)];
    end
end
end

