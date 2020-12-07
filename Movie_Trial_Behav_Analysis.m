%% 
% This script dedicates to quantify and plot the behavior part of Movie
% viewing exp. Majorly see the fixation breaking decay through time, and
% the fixation breaking events for different movies.
Animal = "Alfa";Set_Path;
setMatlabTitle("Movie Behavior Quantification")
%%
ftr = find(contains(ExpRecord.expControlFN,"Alfa") & contains(ExpRecord.Exp_collection,"Movie"));
ExpRecord(ftr,:)
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftr(6:7), Animal);
%%

%%
% blcidx = arrayfun(@(iblc)find(errornums == 0 & block_arr == iblc & valid_mask), 1:max(block_arr), 'uni', 0);
%%
figdirs = ["E:\OneDrive - Washington University in St. Louis\MovieDynamics\2020-10-27-Alfa-Chan09-1\behav",...
           "E:\OneDrive - Washington University in St. Louis\MovieDynamics\2020-10-29-Alfa-Chan09-1\behav"];
for Triali = 1:2
Trials = Trials_new{Triali};
meta = meta_new{Triali};
figdir = figdirs(Triali);

mkdir(figdir)
errornums = arrayfun(@(B)B.TrialError, Trials.B);
valid_mask = arrayfun(@(B)numel(B.TaskObject.Attribute)>1, Trials.B); % Trial with fix and a movie in it. 
block_arr = arrayfun(@(B)B.Block, Trials.B);
% 
corr_trial_num = arrayfun(@(iblc)sum(errornums == 0 & valid_mask & block_arr == iblc), 1:max(block_arr)); % correct Trial
brkfix_trial_num = arrayfun(@(iblc)sum(errornums == 3 & valid_mask & block_arr == iblc), 1:max(block_arr)); % break fixation
none_trial_num = arrayfun(@(iblc)sum(errornums == 4 & valid_mask & block_arr == iblc), 1:max(block_arr)); % no fixation, lazy or sleepy
tot_trial_num = arrayfun(@(iblc)sum(valid_mask & block_arr == iblc), 1:max(block_arr));
%
trialmovs = arrayfun(@(B)string(B.TaskObject.Attribute{1}{2}),Trials.B);
uniqmovs = unique(trialmovs);
uniqmovs = uniqmovs(uniqmovs~="0"); % for the trials that only contains a fixation point this will extrat 0 instead of the movie name.
stim_msks = arrayfun(@(mv)trialmovs==mv, uniqmovs,'Uni',0);
brkfix_4_stim = cellfun(@(msk)sum(msk & errornums==3),stim_msks);
brkfix_4stim_block = arrayfun(@(iblc) ...
                        cellfun(@(msk)sum(msk & errornums==3 & block_arr == iblc), stim_msks), ...
                            [1:max(block_arr)]','Uni',0);
brkfix_4stim_block = cell2mat(brkfix_4stim_block);
%%
mvnm_short = cellfun(@(mv)split(mv,"\"),uniqmovs,'Uni',0);
mvnm_short = string(cellfun(@(parts)strrep(parts{end}(1:end-20),'_',' '),mvnm_short,'Uni',0));
%%
h = figure;
bar([brkfix_trial_num;none_trial_num;corr_trial_num]')
legend(["Break fixation","No Fixation","Correct"]);
ylabel("Trial #");xlabel("Block #");box off
title(compose("%s %s\nBehavior Summary", meta.ephysFN, meta.expControlFN),'Interp','none')
savefig(h, fullfile(figdir, "behavior_correct_summary.fig"))
saveas(h, fullfile(figdir, "behavior_correct_summary.png"))

%%
h2 = figure;
bar(brkfix_4stim_block','stacked')
legend([compose("%d",1:6)]);
ylabel("Trial #");xlabel("Stimuli #");box off;xticklabels(mvnm_short);xtickangle(30)
title(compose("%s %s\n Break Fixation Trial # for Each Movie", meta.ephysFN, meta.expControlFN),'Interp','none')
savefig(h2, fullfile(figdir, "movie_indiv_summary.fig"))
saveas(h2, fullfile(figdir, "movie_indiv_summary.png"))

end

